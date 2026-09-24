"""pyGCC backend.

Implementation notes, all verified against pyGCC 1.5.3:

* pyGCC exposes no arbitrary-reaction API. ``calcRxnlogK`` only evaluates a
  species' database-defined formation reaction, so this backend supplies
  per-species formation energies via ``species_eos.supcrtaq`` and the reaction
  layer above assembles ``sum(nu_i * dGf_i)`` itself.
* ``supcrtaq`` returns an ``ndarray`` of shape (1,) in calories per mole, not a
  float.
* Water is *not* in ``dbaccessdic`` -- HKF parameterises solutes, not the
  solvent -- so H2O is a required special case handled through ``iapws95``.
* At exactly 100 C, P = 1 bar lies below the saturation pressure, water would be
  steam, and ``supcrtaq`` returns NaN. Only pyGCC's ``P='T'`` saturation
  sentinel is reliable on that boundary; a literal 1.013 also fails.
* A single evaluation costs roughly 28 ms, and passing a vector of temperatures
  is *slower* than looping. Hence scalar calls behind a cache, and precomputed
  grids for anything interactive.
"""

from __future__ import annotations

import difflib
import os
import warnings
from typing import Any

import numpy as np

from ..exceptions import MissingDataError, OutOfRangeError, SpeciesNotFoundError
from ..units import CAL_PER_MOL, KJ_PER_MOL_STR, Quantity
from .base import SpeciesRecord, ThermoBackend

#: Above this temperature, pressure must follow the saturation curve.
_SATURATION_THRESHOLD_C = 99.0

#: The IAPWS-95 water EOS fails at exactly 0 C.
_MIN_TEMPERATURE_C = 0.01
_MAX_TEMPERATURE_C = 100.0

#: Aqueous HKF entries carry 13 parameters; gases and minerals carry a
#: different, shorter set and must be routed through ``heatcap`` instead.
_HKF_ENTRY_LENGTH = 13


#: Names for liquid water, which never appears in the species database.
_WATER_NAMES = {"H2O", "H2O(l)", "water", "WATER"}

#: Formation energy of liquid water at 25 C, 1 bar, on the SUPCRT92 convention
#: (Helgeson & Kirkham), in cal/mol.
#:
#: pyGCC's IAPWS-95 routine returns -56677.9 cal/mol instead. The 9 cal/mol
#: difference is a reference-state convention, not an error: -237.14 kJ/mol is
#: the modern CODATA value while -237.18 is the older SUPCRT one.
#:
#: It matters because every other number in this library descends from the
#: SUPCRT lineage. Back-calculating the water energy implied by log K values in
#: thermo.com.dat -- independently from the OH-, Fe+++ and CO2(aq) reactions --
#: gives -56687.7 cal/mol every time, agreeing with the SUPCRT convention to
#: within the rounding of the tabulated log K. Using the IAPWS value instead
#: leaves a systematic 0.041 kJ/mol error per mole of water in every reaction.
#:
#: So water is shifted by a constant onto the SUPCRT reference. The constant
#: preserves the IAPWS-95 temperature dependence exactly; only the datum moves.
SUPCRT_WATER_GIBBS_CAL_25C = -56687.0

#: Helgeson ion-size parameters (Angstrom) for the extended Debye-Huckel term.
#: Only ions this library routinely encounters are tabulated; anything else
#: falls back to _DEFAULT_ION_SIZE, which is the usual practice.
ION_SIZE_ANGSTROM: dict[str, float] = {
    "H+": 9.0,
    "OH-": 3.5,
    "Na+": 4.0,
    "K+": 3.0,
    "Ca++": 6.0,
    "Mg++": 8.0,
    "Cl-": 3.0,
    "SO4--": 4.0,
    "HCO3-": 4.0,
    "CO3--": 4.5,
    "HS-": 3.5,
    "S--": 5.0,
    "NO3-": 3.0,
    "NO2-": 3.0,
    "NH4+": 2.5,
    "Fe++": 6.0,
    "Fe+++": 9.0,
    "Mn++": 6.0,
    "HPO4--": 4.0,
    "H2PO4-": 4.0,
    "Acetate": 4.5,
    "Formate(aq)": 3.5,
    "Lactate(aq)": 4.5,
}
_DEFAULT_ION_SIZE = 4.5


#: The HKF database this library pins. pyGCC's own default is speq21.dat;
#: speq23 is a strict superset of it -- 1594 -> 1597 species, nothing removed,
#: and not one formation energy of the species we expose moves by so much as
#: 1e-9 kJ/mol (checked across all 85). It adds epsomite, hexahydrite and
#: kieserite. Pinning it by name rather than relying on pyGCC's default also
#: means the provenance record can say which file was read.
DEFAULT_HKF_DATABASE = "speq23.dat"

#: Extra HKF sources consulted after the primary one, in order, for species it
#: does not carry. supcrtbl.dat is SUPCRTBL (Zimmer et al. 2016), which revised
#: SUPCRT92's mineral end-members against Holland & Powell (2011) and added
#: arsenic minerals. It holds only 444 species and lacks 1458 that speq23 has,
#: so it can only ever be a supplement -- never a base.
SUPPLEMENTARY_HKF_DATABASES = {"supcrtbl.dat": "HP11"}


def bundled_database(filename: str) -> str:
    """Absolute path to a database file shipped inside pyGCC."""
    import pygcc

    return os.path.join(os.path.dirname(pygcc.__file__), "default_db", filename)


class PygccBackend(ThermoBackend):
    """Formation energies and solvent properties from pyGCC."""

    name = "pygcc"

    def __init__(
        self,
        database: str | None = None,
        dielectric_method: str = "JN91",
        mineral_database: str | None = None,
        use_minerals: bool = True,
        water_convention: str = "supcrt",
        allow_extrapolation: bool = False,
    ):
        self._allow_extrapolation = allow_extrapolation
        if water_convention not in ("supcrt", "iapws"):
            raise ValueError(
                "water_convention must be 'supcrt' (consistent with the rest of "
                "the database) or 'iapws' (pyGCC's raw IAPWS-95 value)"
            )
        self._water_convention = water_convention
        self._water_offset_cal: float | None = None
        self._database = (
            database if database is not None else bundled_database(DEFAULT_HKF_DATABASE)
        )
        self._dielectric_method = dielectric_method
        self._mineral_database = mineral_database
        self._use_minerals = use_minerals
        self._db: Any = None
        self._species: dict[str, Any] | None = None
        self._minerals: dict[str, Any] | None = None
        self._supplementary: dict[str, Any] | None = None
        self._gibbs_cache: dict[tuple[str, float, Any], float] = {}
        self._solvent_cache: dict[tuple[float, Any], dict] = {}

    # --- lazy database loading -------------------------------------------------

    @property
    def species_dict(self) -> dict[str, Any]:
        if self._species is None:
            from pygcc.pygcc_utils import db_reader

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._db = db_reader(dbaccess=self._database)
            self._species = self._db.dbaccessdic
        return self._species

    @property
    def supplementary_species(self) -> dict:
        """Species from the supplementary HKF files, keyed name -> (entry, method).

        Consulted **last**, after the primary database and the GWB minerals, so
        adding a file here can only ever add species -- it never changes a
        number that already resolved. That matters because SUPCRTBL is not a
        newer edition of the same data: it revised SUPCRT92's mineral
        end-members against Holland & Powell (2011), and the two disagree by
        real amounts (hematite 1.7, goethite 2.4, magnetite 3.4 kJ/mol).
        Preferring it silently would move published numbers.

        Its entries also need a different equation of state. A SUPCRT92 mineral
        carries Maier-Kelley coefficients in calories; a SUPCRTBL one carries
        Holland & Powell parameters in kJ, in a longer record. pyGCC can
        evaluate both, but only if told which -- hence the method tag travelling
        with each entry.
        """
        if self._supplementary is None:
            from pygcc.pygcc_utils import db_reader

            self._supplementary = {}
            for filename, method in SUPPLEMENTARY_HKF_DATABASES.items():
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        entries = db_reader(dbaccess=bundled_database(filename)).dbaccessdic
                except Exception:
                    continue
                for name, entry in entries.items():
                    self._supplementary.setdefault(name, (entry, method, filename))
        return self._supplementary

    @property
    def minerals(self) -> dict:
        """Mineral entries from the GWB database, loaded on first use.

        Supplies phases the direct-access database lacks -- notably every
        manganese oxide -- by deriving their formation energies from tabulated
        log K values and the basis species. See :mod:`.gwb`.
        """
        if not self._use_minerals:
            return {}
        if self._minerals is None:
            from .gwb import default_gwb_path, load_gwb_minerals

            path = self._mineral_database or default_gwb_path()
            try:
                self._minerals = load_gwb_minerals(path)
            except OSError:
                self._minerals = {}
        return self._minerals

    @property
    def version(self) -> str:
        import importlib.metadata as md

        try:
            return md.version("pygcc")
        except md.PackageNotFoundError:  # pragma: no cover
            return "unknown"

    @property
    def database_name(self) -> str:
        return os.path.basename(self._database)

    # --- species resolution ----------------------------------------------------

    def resolve(self, name: str) -> SpeciesRecord:
        """Where this species' formation energy comes from.

        Follows the same dispatch order as :meth:`_compute_gibbs_cal` --
        water, HKF, GWB log K, hand-entered -- so anything ``delta_Gf`` can
        evaluate, this can describe. Checking only the HKF dictionary, as this
        once did, made it raise for Goethite and ``As(OH)3(aq)`` while those
        species worked perfectly well in reactions, and then suggest the name
        just passed, which reads as a bug rather than a scope limit.
        """
        if name in _WATER_NAMES:
            return SpeciesRecord(name=name, backend_name="H2O", source="IAPWS-95 via pygcc.iapws95")

        if name in self.species_dict:
            return SpeciesRecord(name=name, backend_name=name, source=f"pyGCC {self.database_name}")

        if name in self.minerals:
            return SpeciesRecord(
                name=name,
                backend_name=name,
                source=f"{self.mineral_database_name} (log K route)",
            )

        if name in self.supplementary_species:
            _, method, filename = self.supplementary_species[name]
            return SpeciesRecord(name=name, backend_name=name, source=f"{filename} ({method})")

        from ..supplemental import supplemental_species

        entry = supplemental_species().get(name)
        if entry is not None:
            return SpeciesRecord(
                name=name,
                backend_name=name,
                source=f"hand-entered supplemental table ({entry.provenance_summary})",
                verified=entry.verified,
            )

        raise SpeciesNotFoundError(name, suggestions=self.suggest(name))

    @property
    def mineral_database_name(self) -> str:
        from .gwb import default_gwb_path

        return os.path.basename(self._mineral_database or default_gwb_path())

    def suggest(self, name: str, n: int = 5) -> list[str]:
        """Closest species names, for error messages."""
        from ..supplemental import supplemental_species

        candidates = (
            list(self.species_dict)
            + list(self.minerals)
            + list(self.supplementary_species)
            + list(supplemental_species())
        )
        return difflib.get_close_matches(name, candidates, n=n, cutoff=0.6)

    def available_species(self) -> list[str]:
        from ..supplemental import supplemental_species

        return sorted(
            set(self.species_dict)
            | set(self.minerals)
            | set(self.supplementary_species)
            | set(supplemental_species())
        )

    def search(self, pattern: str) -> list[str]:
        """Case-insensitive substring search over species names."""
        needle = pattern.lower()
        return sorted(s for s in self.species_dict if needle in s.lower())

    # --- temperature and pressure policy ---------------------------------------

    def _check_temperature(self, temperature_c: float) -> float:
        t = float(temperature_c)
        if not (_MIN_TEMPERATURE_C <= t <= _MAX_TEMPERATURE_C):
            raise OutOfRangeError(
                f"temperature {t} C is outside the supported range "
                f"{_MIN_TEMPERATURE_C}-{_MAX_TEMPERATURE_C} C. "
                "The IAPWS-95 water equation of state fails at exactly 0 C."
            )
        return t

    def _pressure_for(self, temperature_c: float, pressure_bar) -> Any:
        """Resolve the pressure argument passed to pyGCC.

        Defaults to 1 bar, switching to pyGCC's saturation sentinel at and above
        99 C where 1 bar would put water in the vapour field and yield NaN.
        Below that threshold the two differ by well under 1 cal/mol.
        """
        if pressure_bar is not None:
            return pressure_bar
        return "T" if temperature_c >= _SATURATION_THRESHOLD_C else 1.0

    # --- the single numerical primitive ----------------------------------------

    def delta_Gf(self, name: str, temperature_c: float, pressure_bar=None) -> Quantity:
        """Standard-state formation free energy of ``name`` at T, P, in kJ/mol."""
        t = self._check_temperature(temperature_c)
        p = self._pressure_for(t, pressure_bar)
        key = (name, round(t, 6), p)
        if key not in self._gibbs_cache:
            self._gibbs_cache[key] = self._compute_gibbs_cal(name, t, p)
        value_cal = self._gibbs_cache[key]
        return (Quantity(value_cal, CAL_PER_MOL)).to(KJ_PER_MOL_STR)

    def _compute_gibbs_cal(self, name: str, temperature_c: float, pressure) -> float:
        if name in _WATER_NAMES:
            return self._water_gibbs_cal(temperature_c, pressure)

        species = self.species_dict
        if name not in species:
            if name in self.minerals:
                return self._mineral_gibbs_cal(name, temperature_c, pressure)
            if name in self.supplementary_species:
                return self._supplementary_gibbs_cal(name, temperature_c, pressure)
            supplemented = self._supplemental_gibbs_cal(name, temperature_c)
            if supplemented is not None:
                return supplemented
            raise SpeciesNotFoundError(name, suggestions=self.suggest(name))

        entry = species[name]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if len(entry) == _HKF_ENTRY_LENGTH:
                from pygcc.species_eos import supcrtaq

                raw = supcrtaq(temperature_c, pressure, entry)
            else:
                # Gases and minerals are not HKF solutes; they carry
                # Maier-Kelley heat-capacity coefficients instead and go
                # through pyGCC's heatcap class.
                from pygcc import heatcap

                raw = heatcap(
                    T=temperature_c,
                    P=pressure,
                    Species_ppt=entry,
                    Species=name,
                    method="SUPCRT",
                ).dG
        value = float(np.asarray(raw, dtype=float).ravel()[0])
        if not np.isfinite(value):
            raise MissingDataError(
                f"pyGCC returned a non-finite free energy for {name!r} at "
                f"{temperature_c} C, P={pressure!r}. This usually means the "
                "state point falls outside the equation of state's valid region."
            )
        return value

    def _mineral_gibbs_cal(self, name: str, temperature_c: float, pressure) -> float:
        """Formation energy of a mineral, derived from its log K.

        The basis species are looked up through the normal path, so the result
        sits on the same standard-state convention as everything else.
        """
        from .gwb import mineral_gibbs_cal

        def basis(species_name: str) -> float:
            return self._compute_gibbs_cal(species_name, temperature_c, pressure)

        return mineral_gibbs_cal(self.minerals[name], temperature_c, basis)

    def _supplementary_gibbs_cal(self, name: str, temperature_c: float, pressure) -> float:
        """Formation energy from a supplementary HKF file, in cal/mol.

        These carry a different equation of state from the primary database --
        Holland & Powell rather than Maier-Kelley -- so the method has to be
        passed explicitly; pyGCC defaults to SUPCRT and would misread the
        record entirely.

        The unit needs no handling, which is worth stating because it looks as
        though it should. A SUPCRTBL record stores kJ where a SUPCRT92 one
        stores calories -- Pyrite is -160.16 against -38293.0 -- but
        ``heatcap`` converts internally and returns calories either way.
        Converting again here lands you a factor of 4.184 out, on a number
        that still looks entirely plausible.
        """
        from pygcc import heatcap

        entry, method, _ = self.supplementary_species[name]
        raw = heatcap(
            T=temperature_c,
            P=pressure,
            Species_ppt=entry,
            Species=name,
            method=method,
        ).dG
        return float(raw if not hasattr(raw, "__len__") else raw[0])

    def _supplemental_gibbs_cal(self, name: str, temperature_c: float):
        """Hand-entered value for a species no database carries, or None.

        Unverified entries warn on every use; see :mod:`..supplemental`.
        """
        from ..supplemental import supplemental_species

        entry = supplemental_species().get(name)
        if entry is None:
            return None
        entry.warn_if_unverified()
        value_kj = entry.gibbs_kJ_mol(temperature_c, allow_extrapolation=self._allow_extrapolation)
        return value_kj * 1000.0 / 4.184  # kJ/mol -> cal/mol

    def _water_gibbs_cal(self, temperature_c: float, pressure) -> float:
        value, density = self._raw_water(temperature_c, pressure)
        if density < 500.0:
            raise OutOfRangeError(
                f"water is not liquid at {temperature_c} C, P={pressure!r} "
                f"(density {density:.1f} kg/m3)"
            )
        return value + self._water_offset()

    def _raw_water(self, temperature_c: float, pressure):
        """IAPWS-95 Gibbs energy and density, before any datum shift."""
        from pygcc import iapws95

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            water = iapws95(T=temperature_c, P=pressure)
        return (
            float(np.asarray(water.G, dtype=float).ravel()[0]),
            float(np.asarray(water.rho, dtype=float).ravel()[0]),
        )

    def _water_offset(self) -> float:
        """Constant shift putting water on the same datum as everything else.

        Computed once from the 25 C reference point, so the IAPWS-95
        temperature dependence is preserved and only the datum moves. See
        :data:`SUPCRT_WATER_GIBBS_CAL_25C`.
        """
        if self._water_convention == "iapws":
            return 0.0
        if self._water_offset_cal is None:
            raw, _ = self._raw_water(25.0, 1.0)
            self._water_offset_cal = SUPCRT_WATER_GIBBS_CAL_25C - raw
        return self._water_offset_cal

    # --- solvent properties ----------------------------------------------------

    def solvent_properties(self, temperature_c: float, pressure_bar=None) -> dict:
        """Debye-Huckel A, B, the b-dot parameter, dielectric constant, density."""
        t = self._check_temperature(temperature_c)
        p = self._pressure_for(t, pressure_bar)
        key = (round(t, 6), p)
        if key in self._solvent_cache:
            return self._solvent_cache[key]

        from pygcc import water_dielec

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = water_dielec(T=t, P=p, Dielec_method=self._dielectric_method)

        def scalar(attr):
            return float(np.asarray(getattr(w, attr), dtype=float).ravel()[0])

        props = {
            "A": scalar("Ah"),
            "B": scalar("Bh"),
            "bdot": scalar("bdot"),
            "dielectric_constant": scalar("E"),
            "density_kg_m3": scalar("rho"),
            "temperature_c": t,
            "pressure": p,
        }
        self._solvent_cache[key] = props
        return props

    def activity_coefficient(
        self,
        charge: int,
        ionic_strength: float,
        temperature_c: float,
        pressure_bar=None,
        ion_size: float | None = None,
        species_name: str | None = None,
    ) -> float:
        """Extended Debye-Huckel (B-dot) activity coefficient.

        .. math::
            \\log_{10}\\gamma = \\frac{-A z^2 \\sqrt{I}}{1 + B \\mathring{a} \\sqrt{I}}
            + \\dot{b} I

        Neutral species are assigned gamma = 1. That is the usual Helgeson
        convention and avoids inventing Setchenow coefficients we do not have;
        it does mean dissolved-gas activities are treated as ideal.
        """
        if charge == 0:
            return 1.0
        if ionic_strength <= 0:
            return 1.0

        props = self.solvent_properties(temperature_c, pressure_bar)
        if ion_size is None:
            ion_size = ION_SIZE_ANGSTROM.get(species_name or "", _DEFAULT_ION_SIZE)

        sqrt_i = float(np.sqrt(ionic_strength))
        log_gamma = (
            -props["A"] * charge**2 * sqrt_i / (1.0 + props["B"] * ion_size * sqrt_i)
            + props["bdot"] * ionic_strength
        )
        return float(10.0**log_gamma)
