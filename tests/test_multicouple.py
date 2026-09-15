"""Couples whose sides name more than one species.

A single-species couple cannot express an *incomplete* oxidation: syntrophic
propionate oxidation yields acetate **and** bicarbonate together. The ratio
between those products is not fixed by conservation -- propionate could in
principle go to one and a half acetate, or to three bicarbonate -- so the
caller states it and conservation then fixes the scale between the two sides.
"""

import unittest
import warnings
from fractions import Fraction

import microbial_thermo as mt
from microbial_thermo.balance import (
    balance_half_reaction,
    normalize_side,
    side_oxidation_state,
    verify_conservation,
)
from microbial_thermo.exceptions import BalancingError
from microbial_thermo.library import reaction
from microbial_thermo.species import resolve


def coefficient_map(coefficients):
    return {s.backend: v for s, v in coefficients.items()}


class TestNormalizeSide(unittest.TestCase):
    def test_a_bare_name_becomes_one_mole(self):
        side = normalize_side("Acetate")
        self.assertEqual(len(side), 1)
        self.assertEqual(side[0][1], Fraction(1))

    def test_a_species_object_is_accepted(self):
        side = normalize_side(resolve("acetate"))
        self.assertEqual(side[0][0].backend, "Acetate")

    def test_a_list_means_equal_proportions(self):
        side = normalize_side(["Acetate", "HCO3-"])
        self.assertEqual([c for _, c in side], [Fraction(1), Fraction(1)])

    def test_explicit_proportions_by_pairs(self):
        side = normalize_side([("Acetate", 2)])
        self.assertEqual(side[0][1], Fraction(2))

    def test_explicit_proportions_by_mapping(self):
        side = normalize_side({"Acetate": 2, "HCO3-": 1})
        self.assertEqual(dict((s.backend, c) for s, c in side)["Acetate"], Fraction(2))

    def test_empty_side_rejected(self):
        with self.assertRaises(BalancingError):
            normalize_side([])

    def test_non_positive_coefficient_rejected(self):
        with self.assertRaises(BalancingError):
            normalize_side([("Acetate", 0)])
        with self.assertRaises(BalancingError):
            normalize_side([("Acetate", -1)])


class TestMultiProductBalancing(unittest.TestCase):
    def test_propionate_to_acetate_and_bicarbonate(self):
        """The textbook syntrophic half reaction."""
        half = balance_half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        self.assertEqual(half.n_electrons, 6)
        self.assertEqual(
            coefficient_map(half.coefficients),
            {
                "Acetate": -1,
                "HCO3-": -1,
                "H+": -7,
                "e-": -6,
                "Propanoate(aq)": 1,
                "H2O": 3,
            },
        )
        verify_conservation(half.coefficients)

    def test_butyrate_to_two_acetate(self):
        half = balance_half_reaction("Butanoate(aq)", [("Acetate", 2)])
        self.assertEqual(half.n_electrons, 4)
        verify_conservation(half.coefficients)

    def test_side_is_flagged_as_not_simple(self):
        half = balance_half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        self.assertFalse(half.is_simple)
        self.assertTrue(balance_half_reaction("HS-", "SO4-2").is_simple)

    def test_singular_accessors_refuse_a_multi_product_side(self):
        half = balance_half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        with self.assertRaises(BalancingError):
            _ = half.oxidized
        # The reduced side is still single, so that one works.
        self.assertEqual(half.reduced.backend, "Propanoate(aq)")

    def test_oxidation_state_is_the_atom_weighted_mean(self):
        """Two acetate carbons at 0 and one bicarbonate carbon at +4 average
        to +4/3, not to +2."""
        half = balance_half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        oxidized, reduced = half.oxidation_states()
        self.assertEqual(oxidized, Fraction(4, 3))
        self.assertEqual(reduced, Fraction(-2, 3))

    def test_weighted_mean_helper_directly(self):
        side = normalize_side(["Acetate", "HCO3-"])
        self.assertEqual(side_oxidation_state(side, "C"), Fraction(4, 3))

    def test_unbalanced_elements_are_caught(self):
        """Adding nitrogen to one side only must fail conservation."""
        with self.assertRaises(BalancingError):
            balance_half_reaction("Propanoate(aq)", ["Acetate", "NO3-"])

    def test_proportions_change_the_answer(self):
        """A different product split is a different reaction, and both
        balance -- which is exactly why the caller must supply the ratio."""
        one = balance_half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        other = balance_half_reaction("Propanoate(aq)", [("HCO3-", 3)])
        self.assertNotEqual(one.n_electrons, other.n_electrons)
        verify_conservation(other.coefficients)


class TestSimpleCouplesUnaffected(unittest.TestCase):
    """Regression: generalising the sides must not disturb the common case."""

    CASES = [
        ("HS-", "SO4-2", 8),
        ("H2O", "O2(aq)", 2),
        ("H2(g)", "H+", 2),
        ("NH4+", "NO3-", 8),
        ("Fe+2", "Fe+3", 1),
        ("NO2-", "NO3-", 2),
    ]

    def test_reference_half_reactions_unchanged(self):
        for reduced, oxidized, electrons in self.CASES:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                half = balance_half_reaction(reduced, oxidized)
                self.assertEqual(half.n_electrons, electrons)
                self.assertTrue(half.is_simple)
                self.assertEqual(half.oxidized.backend, resolve(oxidized).backend)


class TestSyntrophyEnergetics(unittest.TestCase):
    """Published dG0' values, and the reason syntrophy exists at all."""

    @classmethod
    def setUpClass(cls):
        cls.standard = mt.Conditions(temperature_c=25.0, pH=7.0)

    def _energy(self, name, conditions):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return reaction(
                name, conditions, normalize_to="integer"
            ).delta_G_standard_prime.magnitude

    def test_propionate_matches_published(self):
        """Propionate + 3 H2O -> acetate + HCO3- + H+ + 3 H2, near +76 kJ/mol."""
        value = self._energy("syntrophic_propionate_oxidation", self.standard)
        self.assertAlmostEqual(value, 76.1, delta=6.0)

    def test_butyrate_matches_published(self):
        """Butyrate + 2 H2O -> 2 acetate + H+ + 2 H2, near +48 kJ/mol."""
        value = self._energy("syntrophic_butyrate_oxidation", self.standard)
        self.assertAlmostEqual(value, 48.3, delta=6.0)

    def test_endergonic_at_standard_state(self):
        for name in (
            "syntrophic_propionate_oxidation",
            "syntrophic_butyrate_oxidation",
            "syntrophic_ethanol_oxidation",
        ):
            with self.subTest(metabolism=name):
                self.assertGreater(self._energy(name, self.standard), 0.0)

    def test_low_hydrogen_makes_syntrophy_pay(self):
        """The whole point: the reaction only proceeds once a partner has
        drawn the hydrogen down."""

        def energy(p_h2):
            conditions = mt.Conditions(
                temperature_c=25.0,
                pH=7.0,
                activity_model="ideal",
                partial_pressures={"H2(g)": p_h2},
                concentrations={"Propanoate(aq)": 1e-4, "Acetate": 1e-4, "HCO3-": 2e-3},
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return reaction("syntrophic_propionate_oxidation", conditions).delta_G.magnitude

        self.assertGreater(energy(1.0), 0.0)  # 1 bar: hopeless
        self.assertLess(energy(1e-6), 0.0)  # 1 ppm: pays
        self.assertLess(energy(1e-6), energy(1e-3))

    def test_a_syntrophic_window_exists(self):
        """Propionate oxidation needs hydrogen low; hydrogenotrophic
        methanogenesis needs it high. Both are exergonic only inside a narrow
        overlap, which is why the partners must stay coupled."""

        def energies(p_h2):
            conditions = mt.Conditions(
                temperature_c=25.0,
                pH=7.0,
                activity_model="ideal",
                partial_pressures={"H2(g)": p_h2},
                concentrations={
                    "Propanoate(aq)": 1e-4,
                    "Acetate": 1e-4,
                    "HCO3-": 2e-3,
                    "CO2(aq)": 1e-3,
                    "Methane(aq)": 1e-4,
                },
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                producer = reaction("syntrophic_propionate_oxidation", conditions).delta_G.magnitude
                consumer = reaction("hydrogenotrophic_methanogenesis", conditions).delta_G.magnitude
            return producer, consumer

        inside = [p for p in (1e-6, 1e-5, 1e-4) if all(v < 0 for v in energies(p))]
        self.assertTrue(inside, "no hydrogen pressure suited both partners")
        # Outside the window one partner or the other fails.
        self.assertGreater(energies(1.0)[0], 0.0)
        self.assertGreater(energies(1e-12)[1], 0.0)


class TestLibraryIntegration(unittest.TestCase):
    def test_syntrophy_entries_are_catalogued(self):
        from microbial_thermo.library import default_library

        names = default_library().names()
        for expected in (
            "syntrophic_propionate_oxidation",
            "syntrophic_butyrate_oxidation",
            "syntrophic_ethanol_oxidation",
        ):
            self.assertIn(expected, names)

    def test_show_work_handles_a_multi_product_couple(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            text = (
                reaction(
                    "syntrophic_propionate_oxidation",
                    mt.Conditions(temperature_c=25.0, pH=7.0),
                )
                .show_work()
                .to_text()
            )
        self.assertIn("Acetate", text)
        self.assertIn("HCO3-", text)

    def test_figure_renders_a_multi_product_couple(self):
        import tempfile
        from pathlib import Path

        import matplotlib

        matplotlib.use("Agg")
        from microbial_thermo.figures import plot_half_reactions

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            built = reaction(
                "syntrophic_propionate_oxidation",
                mt.Conditions(temperature_c=25.0, pH=7.0),
            )
            with tempfile.TemporaryDirectory() as directory:
                base = Path(directory) / "syntrophy"
                plot_half_reactions(built, save=base)
                self.assertTrue(base.with_suffix(".svg").exists())


if __name__ == "__main__":
    unittest.main()
