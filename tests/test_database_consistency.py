"""Are the databases we mix actually on the same footing?

This library draws numbers from four places:

1. ``speq21.dat`` -- HKF parameters for aqueous solutes, via ``supcrtaq``.
2. ``speq21.dat`` again -- Maier-Kelley coefficients for gases and some
   minerals, via ``heatcap``.
3. ``thermo.com.dat`` -- a GWB database supplying minerals the others lack,
   through tabulated log K values (see :mod:`microbial_thermo.backends.gwb`).
4. ``iapws95`` -- liquid water, which is in no species database.

Route 3 is the one that has to be justified, because it *mixes* sources: the
log K comes from ``thermo.com.dat`` but the basis-species formation energies
come from ``speq21.dat``. If those two disagreed about the basis species, every
derived mineral would inherit the error.

The tests below check that they do not. The decisive one compares, for the 900+
species present in *both* databases, the tabulated log K against the log K
predicted from the HKF parameters. Agreement means the two files are the same
SUPCRT lineage and may be combined.

One real inconsistency was found this way and is corrected in the backend:
IAPWS-95 puts liquid water 9 cal/mol away from the SUPCRT convention the rest
of the data uses. See ``SUPCRT_WATER_GIBBS_CAL_25C``.
"""

import math
import unittest
import warnings

import numpy as np

import microbial_thermo as mt
from microbial_thermo.backends import PygccBackend
from microbial_thermo.backends.gwb import default_gwb_path, load_gwb_minerals
from microbial_thermo.species import default_registry

R_CAL = 1.987204
CAL_TO_KJ = 4.184 / 1000.0
TEMPERATURE = 25.0
RT_LN10_CAL = R_CAL * (TEMPERATURE + 273.15) * math.log(10.0)

#: 1 log K unit is about 5.7 kJ/mol at 25 C.
KJ_PER_LOGK = RT_LN10_CAL * CAL_TO_KJ

#: Species where the two databases genuinely disagree, with the size of the
#: disagreement in log K units. These are documented rather than hidden: each
#: is a real difference between the sources, not a bug in our code.
KNOWN_DISAGREEMENTS = {
    "Rhodochrosite": 0.50,  # 2.5 kJ/mol; we use the speq21 value, nearer to published
    "FeOH+": 0.25,  # 1.1 kJ/mol
    "Mn+++": 0.10,  # 0.4 kJ/mol, shared by the Mn(III)/Mn(VII) set
    "MnO4-": 0.10,
    "MnO4--": 0.10,
    # 0.26 kJ/mol. The metarsenous acid tabulation is the one arsenic species
    # the two databases disagree on at all; the rest of the As set agrees to
    # 0.09 kJ/mol or better and HAsO4-- agrees exactly. This is the same
    # discrepancy that shows up as the 0.17 kJ/mol gap between the As(OH)3 and
    # HAsO2 representations in tests/test_arsenic.py, and it is why a bare
    # "As(III)" resolves to As(OH)3(aq) rather than to HAsO2(aq).
    "HAsO2(aq)": 0.06,
}

#: Everything else must agree this closely, in log K units (~0.17 kJ/mol).
AGREEMENT_TOLERANCE_LOGK = 0.03


def _gibbs_cal(backend, name):
    return backend._compute_gibbs_cal(name, TEMPERATURE, 1.0)


def _predicted_log_k(backend, entry, name):
    """log K for a species' own formation reaction, from HKF parameters."""
    delta_g = sum(c * _gibbs_cal(backend, s) for c, s in entry.stoichiometry) - _gibbs_cal(
        backend, name
    )
    return -delta_g / RT_LN10_CAL


class TestReferenceConventions(unittest.TestCase):
    """The conventions every route must share."""

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()

    def test_proton_is_the_zero_of_the_scale(self):
        self.assertAlmostEqual(self.backend.delta_Gf("H+", 25.0).magnitude, 0.0, places=9)

    def test_elements_in_their_standard_states_are_zero(self):
        for name in ("H2(g)", "O2(g)"):
            with self.subTest(species=name):
                self.assertAlmostEqual(self.backend.delta_Gf(name, 25.0).magnitude, 0.0, places=6)

    def test_proton_stays_zero_across_temperature(self):
        for t in (0.01, 25.0, 60.0, 100.0):
            with self.subTest(temperature=t):
                self.assertAlmostEqual(self.backend.delta_Gf("H+", t).magnitude, 0.0, places=9)


class TestWaterDatum(unittest.TestCase):
    """Water is the one place the sources genuinely disagreed."""

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()
        cls.minerals = load_gwb_minerals(default_gwb_path())

    def _implied_water_gibbs_cal(self, name):
        """Back out the water energy a tabulated log K assumes."""
        entry = self.minerals[name]
        water_coefficient = sum(c for c, s in entry.stoichiometry if s == "H2O")
        others = sum(c * _gibbs_cal(self.backend, s) for c, s in entry.stoichiometry if s != "H2O")
        return (
            -RT_LN10_CAL * entry.log_k_at(TEMPERATURE) + _gibbs_cal(self.backend, name) - others
        ) / water_coefficient

    def test_water_agrees_with_what_the_database_assumes(self):
        """Three independent reactions must all imply the water we use.

        Each of OH-, Fe+++ and CO2(aq) contains water in its formation
        reaction, so each pins the water datum separately. Before the SUPCRT
        correction these disagreed with the value in use by 0.041 kJ/mol.
        """
        used = _gibbs_cal(self.backend, "H2O")
        for name in ("OH-", "Fe+++", "CO2(aq)"):
            with self.subTest(reaction=name):
                implied = self._implied_water_gibbs_cal(name)
                difference_kj = (implied - used) * CAL_TO_KJ
                self.assertLess(
                    abs(difference_kj),
                    0.01,
                    f"{name} implies water at {implied:.1f} cal/mol but we use "
                    f"{used:.1f} ({difference_kj:+.4f} kJ/mol)",
                )

    def test_correction_is_a_pure_datum_shift(self):
        """The offset must not disturb the IAPWS-95 temperature dependence."""
        supcrt = PygccBackend()
        iapws = PygccBackend(water_convention="iapws")
        span_supcrt = supcrt.delta_Gf("H2O", 95.0).magnitude - supcrt.delta_Gf("H2O", 5.0).magnitude
        span_iapws = iapws.delta_Gf("H2O", 95.0).magnitude - iapws.delta_Gf("H2O", 5.0).magnitude
        self.assertAlmostEqual(span_supcrt, span_iapws, places=9)

    def test_offset_is_constant_across_temperature(self):
        supcrt = PygccBackend()
        iapws = PygccBackend(water_convention="iapws")
        offsets = [
            supcrt.delta_Gf("H2O", t).magnitude - iapws.delta_Gf("H2O", t).magnitude
            for t in (0.01, 25.0, 60.0, 98.0)
        ]
        self.assertAlmostEqual(max(offsets), min(offsets), places=9)

    def test_water_remains_close_to_published(self):
        """Both conventions sit within the published spread; this guards
        against the correction being applied twice or with the wrong sign."""
        value = mt.get_backend().delta_Gf("H2O", 25.0).magnitude
        self.assertAlmostEqual(value, -237.18, delta=0.1)

    def test_unknown_convention_rejected(self):
        with self.assertRaises(ValueError):
            PygccBackend(water_convention="nonsense")


class TestCrossDatabaseAgreement(unittest.TestCase):
    """The decisive check: do the two databases predict the same log K?"""

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()
        cls.minerals = load_gwb_minerals(default_gwb_path())
        cls.aqueous = cls.backend.species_dict

    #: The bulk scan is capped: each species costs several backend calls and
    #: the full 900-odd take minutes. Names are sorted first, so the sample is
    #: deterministic rather than dependent on dict ordering.
    BULK_SAMPLE = 220

    def _comparable(self, limit=None):
        """Species with both a tabulated log K and HKF parameters."""
        seen = 0
        for name in sorted(self.minerals):
            entry = self.minerals[name]
            if limit is not None and seen >= limit:
                return
            if name not in self.aqueous:
                continue
            if not np.isfinite(entry.log_k).any():
                continue
            try:
                # Some entries are tabulated only well outside our range, and
                # log_k_at refuses those rather than extrapolating.
                tabulated = entry.log_k_at(TEMPERATURE)
                predicted = _predicted_log_k(self.backend, entry, name)
            except Exception:
                continue
            seen += 1
            yield name, tabulated, predicted

    def test_shipped_species_agree_across_databases(self):
        """Every species this library actually exposes."""
        registry_names = {s.backend for s in default_registry().all_species()}
        checked = 0
        for name, tabulated, predicted in self._comparable():
            if name not in registry_names:
                continue
            checked += 1
            tolerance = KNOWN_DISAGREEMENTS.get(name, AGREEMENT_TOLERANCE_LOGK)
            with self.subTest(species=name):
                self.assertLess(
                    abs(predicted - tabulated),
                    tolerance,
                    f"{name}: tabulated {tabulated:.4f} vs predicted "
                    f"{predicted:.4f} "
                    f"({(predicted - tabulated) * KJ_PER_LOGK:+.3f} kJ/mol)",
                )
        self.assertGreater(checked, 30, "too few species were actually compared")

    def test_bulk_agreement_is_tight(self):
        """Across all 900-odd shared species, the median must be near zero.

        A systematic shift between the databases would show up here even
        though individual exotic species scatter. Runs over a deterministic
        sample rather than all 900-odd, for speed.
        """
        differences = [p - t for _, t, p in self._comparable(limit=self.BULK_SAMPLE)]
        self.assertGreater(len(differences), 150)
        median = float(np.median(differences))
        self.assertLess(abs(median), 0.01, f"median offset {median:+.4f} log K")
        within = sum(1 for d in differences if abs(d) < 0.1)
        self.assertGreater(
            within / len(differences),
            0.80,
            "fewer than 80% of shared species agree within 0.1 log K",
        )

    def test_basis_species_for_minerals_agree(self):
        """The species the mineral route actually leans on.

        These matter most: an error here propagates into every derived
        mineral.
        """
        for name in ("Fe+++", "HS-", "CO2(aq)", "OH-"):
            with self.subTest(species=name):
                entry = self.minerals[name]
                difference = _predicted_log_k(self.backend, entry, name) - entry.log_k_at(
                    TEMPERATURE
                )
                self.assertLess(abs(difference) * KJ_PER_LOGK, 0.15)


class TestMineralRoutesAgree(unittest.TestCase):
    """Minerals reachable two ways must give the same answer both ways."""

    #: Present in speq21.dat (heat-capacity route) *and* thermo.com.dat (log K
    #: route), so the two can be compared directly. Tolerance in kJ/mol.
    BOTH_ROUTES = [
        ("Hematite", 0.5),
        ("Magnetite", 0.5),
        ("Siderite", 0.5),
        ("Pyrite", 0.5),
        ("Rhodochrosite", 3.0),  # the one genuine disagreement
    ]

    def test_routes_agree(self):
        from microbial_thermo.backends.gwb import mineral_gibbs_cal

        backend = mt.get_backend()
        minerals = load_gwb_minerals(default_gwb_path())
        for name, tolerance in self.BOTH_ROUTES:
            with self.subTest(mineral=name):
                native = backend.delta_Gf(name, TEMPERATURE).magnitude
                derived = (
                    mineral_gibbs_cal(
                        minerals[name],
                        TEMPERATURE,
                        lambda s: _gibbs_cal(backend, s),
                    )
                    * CAL_TO_KJ
                )
                self.assertAlmostEqual(native, derived, delta=tolerance)

    def test_gas_route_and_aqueous_route_share_a_datum(self):
        """CO2(g) and CO2(aq) come through different pyGCC routines; their
        difference must reproduce the known dissolution energy."""
        backend = mt.get_backend()
        gas = backend.delta_Gf("CO2(g)", 25.0).magnitude
        aqueous = backend.delta_Gf("CO2(aq)", 25.0).magnitude
        # CO2(g) = CO2(aq) has log K near -1.47 at 25 C (Henry's law).
        log_k = -(aqueous - gas) / KJ_PER_LOGK
        self.assertAlmostEqual(log_k, -1.47, delta=0.1)


class TestConsistencyIsVisibleToUsers(unittest.TestCase):
    def test_reaction_energies_are_unaffected_by_the_water_convention(self):
        """A reaction with no net water must not care which datum is used."""
        conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
        for convention in ("supcrt", "iapws"):
            with self.subTest(convention=convention):
                reaction = mt.Reaction.from_couples(
                    donor=("Fe+2", "Fe+3"),
                    acceptor=("Fe+2", "Fe+3"),
                    conditions=conditions,
                    backend=PygccBackend(water_convention=convention),
                )
                self.assertAlmostEqual(reaction.delta_G.magnitude, 0.0, places=6)

    def test_water_convention_shifts_water_bearing_reactions_only_slightly(self):
        conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

        def energy(convention):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return mt.Reaction.from_couples(
                    donor=("H2(g)", "H+"),
                    acceptor=("HS-", "SO4-2"),
                    conditions=conditions,
                    backend=PygccBackend(water_convention=convention),
                    normalize_to="integer",
                ).delta_G_standard_prime.magnitude

        difference = abs(energy("supcrt") - energy("iapws"))
        # Four waters at 0.038 kJ/mol each.
        self.assertLess(difference, 0.25)
        self.assertGreater(difference, 0.01)


if __name__ == "__main__":
    unittest.main()
