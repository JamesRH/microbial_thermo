"""Acid-base speciation.

p*K*a values are computed from the backend, not tabulated, so these tests
double as a check that the backend's formation energies are internally
consistent: getting acetate's p*K*a right to two decimals requires both
acetate and acetic acid to be right.
"""

import unittest
import warnings

import microbial_thermo as mt
from microbial_thermo.exceptions import MicrobialThermoError
from microbial_thermo.speciation import (
    default_families,
    dominant,
    family_by_name,
    family_for,
    fraction_of,
    fractions,
    near_pKa,
    pKa,
    pKa_ladder,
    speciation_table,
)


class TestPkaValues(unittest.TestCase):
    """Computed p*K*a against published values at 25 C."""

    PUBLISHED = {
        "sulfide": [7.00],
        "carbonate": [6.35, 10.33],
        "ammonia": [9.25],
        "phosphate": [2.15, 7.20, 12.35],
        "acetate": [4.76],
        "lactate": [3.86],
        "formate": [3.75],
        "propanoate": [4.87],
        "butanoate": [4.82],
        "sulfite": [7.20],
        "nitrite": [3.25],
        "sulfate": [1.99],
    }

    def test_every_family_matches_published_values(self):
        for name, expected in self.PUBLISHED.items():
            with self.subTest(family=name):
                computed = pKa_ladder(name, temperature_c=25.0)
                self.assertEqual(len(computed), len(expected))
                for got, want in zip(computed, expected, strict=True):
                    self.assertAlmostEqual(got, want, delta=0.1)

    def test_all_families_load(self):
        self.assertEqual(len(default_families()), len(self.PUBLISHED))

    def test_carbonate_first_step_is_balanced_with_water(self):
        """CO2(aq) + H2O = HCO3- + H+, not a bare proton loss.

        Ignoring the water puts pKa1 tens of units out, so this is the case
        that would silently break a naive implementation.
        """
        value = pKa("CO2(aq)", "HCO3-", temperature_c=25.0)
        self.assertAlmostEqual(value, 6.35, delta=0.1)

    def test_pka_shifts_with_temperature(self):
        cold = pKa("H2S(aq)", "HS-", temperature_c=5.0)
        hot = pKa("H2S(aq)", "HS-", temperature_c=95.0)
        self.assertNotAlmostEqual(cold, hot, places=1)

    def test_mismatched_pair_raises(self):
        with self.assertRaises(MicrobialThermoError):
            pKa("H2S(aq)", "SO4-2")


class TestFractions(unittest.TestCase):
    def test_fractions_sum_to_one(self):
        for name in ("sulfide", "carbonate", "phosphate", "ammonia"):
            for ph in (1.0, 4.0, 7.0, 10.0, 13.0):
                with self.subTest(family=name, pH=ph):
                    total = sum(fractions(name, ph).values())
                    self.assertAlmostEqual(total, 1.0, places=9)

    def test_sulfide_is_evenly_split_at_its_pka(self):
        """The motivating case: at pH 7 sulfide is close to 50/50, so naming
        a single form is very nearly a coin flip."""
        distribution = fractions("sulfide", 6.99)
        for value in distribution.values():
            self.assertAlmostEqual(value, 0.5, delta=0.01)

    def test_extreme_ph_does_not_overflow(self):
        """Relative abundances are built in log space; a triprotic acid at
        pH 0 would overflow otherwise."""
        distribution = fractions("phosphate", 0.0)
        self.assertAlmostEqual(sum(distribution.values()), 1.0, places=9)
        self.assertTrue(all(0.0 <= v <= 1.0 for v in distribution.values()))

    def test_dominant_form_tracks_ph(self):
        expected = {
            5.0: "H2S(aq)",
            6.0: "H2S(aq)",
            8.0: "HS-",
            9.0: "HS-",
        }
        for ph, backend_name in expected.items():
            with self.subTest(pH=ph):
                member, _ = dominant("sulfide", ph)
                self.assertEqual(member.backend, backend_name)

    def test_bicarbonate_dominates_at_circumneutral_ph(self):
        member, value = dominant("carbonate", 7.5)
        self.assertEqual(member.backend, "HCO3-")
        self.assertGreater(value, 0.9)

    def test_fraction_of_species_outside_any_family_is_one(self):
        self.assertEqual(fraction_of("O2(aq)", 7.0), 1.0)
        self.assertEqual(fraction_of("Fe+2", 7.0), 1.0)

    def test_speciation_table_columns(self):
        frame = speciation_table("carbonate", 7.5)
        self.assertEqual(len(frame), 3)
        self.assertAlmostEqual(frame["fraction"].sum(), 1.0, places=9)


class TestFamilyLookup(unittest.TestCase):
    def test_lookup_by_member_alias_and_name(self):
        self.assertEqual(family_by_name("sulfide").name, "sulfide")
        self.assertEqual(family_by_name("DIC").name, "carbonate")
        self.assertEqual(family_by_name("H2S").name, "sulfide")
        self.assertEqual(family_for("HCO3-").name, "carbonate")

    def test_unknown_family_returns_none(self):
        self.assertIsNone(family_by_name("unobtainium"))

    def test_near_pka_detection(self):
        self.assertTrue(near_pKa("sulfide", 7.0))
        self.assertFalse(near_pKa("sulfide", 10.0))


class TestReactionIntegration(unittest.TestCase):
    """A measured total becomes the activity of the written species."""

    @staticmethod
    def _activity_of(backend_name, conditions):
        from microbial_thermo.reaction import _activity

        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"), acceptor=("HS-", "SO4-2"), conditions=conditions
        )
        species = next(s for s in reaction.coefficients if s.backend == backend_name)
        return _activity(species, conditions, reaction.backend)

    def _conditions(self, ph):
        return mt.Conditions(
            temperature_c=25.0,
            pH=ph,
            activity_model="ideal",
            total_concentrations={"sulfide": 1e-6},
            partial_pressures={"H2(g)": 1e-5},
            concentrations={"SO4-2": 2.8e-2},
        )

    def test_total_is_split_by_ph(self):
        low = self._activity_of("HS-", self._conditions(6.0))
        high = self._activity_of("HS-", self._conditions(8.0))
        self.assertLess(low, high)
        self.assertAlmostEqual(low, 1e-6 * fraction_of("HS-", 6.0), places=12)
        self.assertAlmostEqual(high, 1e-6 * fraction_of("HS-", 8.0), places=12)

    def test_total_never_exceeds_the_measured_amount(self):
        for ph in (4.0, 7.0, 10.0):
            with self.subTest(pH=ph):
                self.assertLessEqual(self._activity_of("HS-", self._conditions(ph)), 1e-6)

    def test_explicit_concentration_takes_precedence(self):
        conditions = mt.Conditions(
            temperature_c=25.0,
            pH=8.0,
            activity_model="ideal",
            total_concentrations={"sulfide": 1e-6},
            concentrations={"HS-": 5e-9},
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertAlmostEqual(self._activity_of("HS-", conditions), 5e-9, places=12)

    def test_warns_when_naming_one_form_near_a_pka(self):
        conditions = mt.Conditions(
            temperature_c=25.0,
            pH=7.0,
            activity_model="ideal",
            concentrations={"HS-": 1e-6},
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._activity_of("HS-", conditions)
        messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        self.assertTrue(any("pKa" in m for m in messages), messages)

    def test_no_warning_far_from_any_pka(self):
        conditions = mt.Conditions(
            temperature_c=25.0,
            pH=11.0,
            activity_model="ideal",
            concentrations={"HS-": 1e-6},
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._activity_of("HS-", conditions)
        messages = [str(w.message) for w in caught if "pKa" in str(w.message)]
        self.assertEqual(messages, [])

    def test_free_energy_responds_to_speciation(self):
        def energy(ph):
            return mt.Reaction.from_couples(
                donor=("H2(g)", "H+"),
                acceptor=("HS-", "SO4-2"),
                conditions=self._conditions(ph),
            ).delta_G.magnitude

        self.assertNotAlmostEqual(energy(6.0), energy(8.0), places=1)

    def test_cross_check_still_holds_with_speciation(self):
        mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self._conditions(7.0),
        ).verify_consistency()


if __name__ == "__main__":
    unittest.main()
