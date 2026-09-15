"""Hand-entered formation energies for species pyGCC does not carry.

Glucose, pyruvate and hydroxylamine are absent from every database pyGCC
ships, and hydroxylamine in particular is hard to do without -- it is the
central nitrification intermediate.

The values here are of a different kind from everything else in this library.
They are literature figures typed in by hand rather than computed from an
equation of state, and they have not been cross-checked against a second
source the way the mineral data was. So each carries a ``verified`` flag, and
an unverified one warns every time it is used. A figure built on one says so.

That is deliberate. For a teaching tool, a number of unknown provenance that
looks like all the others is worse than no number at all.

Temperature: these are single-temperature values with no heat-capacity data.
They are honoured at their stated temperature and refused elsewhere, unless an
enthalpy is given, in which case a van 't Hoff correction is applied and a
warning issued.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import yaml

from .exceptions import OutOfRangeError
from .units import celsius_to_kelvin

_DATA_FILE = "supplemental_gibbs.yaml"

#: How far from its stated temperature a value may be used without correction.
_TOLERANCE_C = 0.5


@dataclass(frozen=True)
class SupplementalSpecies:
    """One hand-entered formation energy, with its provenance."""

    backend: str
    formula: str
    delta_Gf_kJ_mol: float
    temperature_c: float = 25.0
    delta_Hf_kJ_mol: float | None = None
    display: str = ""
    phase: str = "aq"
    smiles: str | None = None
    aliases: tuple = ()
    verified: bool = False
    provenance: str = ""

    def gibbs_kJ_mol(self, temperature_c: float, allow_extrapolation: bool = False) -> float:
        """Formation energy at ``temperature_c``.

        Returns the tabulated value at its own temperature. Elsewhere it
        refuses, unless an enthalpy is available and extrapolation is allowed,
        in which case a van 't Hoff correction is applied.
        """
        if abs(temperature_c - self.temperature_c) <= _TOLERANCE_C:
            return self.delta_Gf_kJ_mol

        if not allow_extrapolation or self.delta_Hf_kJ_mol is None:
            detail = (
                "no enthalpy is tabulated for it, so there is nothing to extrapolate with"
                if self.delta_Hf_kJ_mol is None
                else "pass allow_extrapolation=True to apply a van 't Hoff correction"
            )
            raise OutOfRangeError(
                f"{self.backend} is a supplemental value tabulated at "
                f"{self.temperature_c:g} C only, and you asked for "
                f"{temperature_c:g} C; {detail}"
            )

        # van 't Hoff, assuming dH is constant over the interval.
        t0 = celsius_to_kelvin(self.temperature_c)
        t1 = celsius_to_kelvin(temperature_c)
        entropy_term = (self.delta_Hf_kJ_mol - self.delta_Gf_kJ_mol) / t0
        value = self.delta_Hf_kJ_mol - t1 * entropy_term
        warnings.warn(
            f"{self.backend} was extrapolated from {self.temperature_c:g} C to "
            f"{temperature_c:g} C by van 't Hoff, assuming a constant enthalpy. "
            "Treat the result as indicative only.",
            UserWarning,
            stacklevel=3,
        )
        return value

    def warn_if_unverified(self) -> None:
        if self.verified:
            return
        warnings.warn(
            f"{self.backend} uses a hand-entered formation energy that has not "
            f"been traced to a primary source ({self.provenance.strip()}). "
            "Check it before relying on the result.",
            UserWarning,
            stacklevel=4,
        )


@lru_cache(maxsize=1)
def supplemental_species() -> dict:
    """The supplemental table, keyed by backend name."""
    text = resources.files("microbial_thermo.data").joinpath(_DATA_FILE).read_text()
    document = yaml.safe_load(text) or {}
    entries = {}
    for item in document.get("species", []):
        entry = SupplementalSpecies(
            backend=item["backend"],
            formula=item["formula"],
            delta_Gf_kJ_mol=float(item["delta_Gf_kJ_mol"]),
            temperature_c=float(item.get("temperature_c", 25.0)),
            delta_Hf_kJ_mol=(float(item["delta_Hf_kJ_mol"]) if "delta_Hf_kJ_mol" in item else None),
            display=item.get("display", ""),
            phase=item.get("phase", "aq"),
            smiles=item.get("smiles"),
            aliases=tuple(str(a) for a in item.get("aliases", [])),
            verified=bool(item.get("verified", False)),
            provenance=item.get("provenance", ""),
        )
        entries[entry.backend] = entry
    return entries


def unverified_names() -> list[str]:
    """Names whose values have not been traced to a primary source."""
    return sorted(name for name, e in supplemental_species().items() if not e.verified)


def uses_unverified(names) -> list[str]:
    """Which of ``names`` are backed by an unverified supplemental value."""
    table = supplemental_species()
    return sorted({n for n in names if n in table and not table[n].verified})
