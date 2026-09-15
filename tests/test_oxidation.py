"""Oxidation states and NOSC, against published reference values.

The NOSC fixtures are the worked examples from LaRowe & Van Cappellen (2011).
The per-atom fixtures exercise the conservation invariant that RDKit's own
experimental CalcOxidationNumbers violates.
"""

import unittest
from fractions import Fraction

from microbial_thermo.exceptions import MicrobialThermoError
from microbial_thermo.oxidation import (
    atom_oxidation_states,
    format_oxidation_state,
    mean_oxidation_state,
    nosc,
)


class TestMeanOxidationState(unittest.TestCase):
    CASES = [
        ("S", "SO4-2", 6),
        ("S", "SO3-2", 4),
        ("S", "S2O3-2", 2),
        ("S", "HS-", -2),
        ("S", "H2S", -2),
        ("S", "S", 0),
        ("N", "NO3-", 5),
        ("N", "NO2-", 3),
        ("N", "N2", 0),
        ("N", "NH4+", -3),
        ("C", "CO2", 4),
        ("C", "CH4", -4),
        ("C", "C2H3O2-", 0),
        ("C", "CH2O", 0),
        ("Fe", "Fe+2", 2),
        ("Fe", "Fe+3", 3),
        ("Mn", "MnO2", 4),
        ("O", "O2", 0),
        ("O", "H2O", -2),
        ("H", "H2", 0),
        ("H", "H+", 1),
    ]

    def test_reference_values(self):
        for element, formula, expected in self.CASES:
            with self.subTest(formula=formula):
                self.assertEqual(mean_oxidation_state(element, formula), expected)

    def test_fractional_state_is_exact(self):
        """Magnetite's iron is +8/3, not 2.6666667."""
        self.assertEqual(mean_oxidation_state("Fe", "Fe3O4"), Fraction(8, 3))

    def test_missing_element_raises(self):
        with self.assertRaises(MicrobialThermoError):
            mean_oxidation_state("N", "CO2")

    def test_formatting(self):
        self.assertEqual(format_oxidation_state(Fraction(6)), "+6")
        self.assertEqual(format_oxidation_state(Fraction(-2)), "-2")
        self.assertEqual(format_oxidation_state(Fraction(0)), "0")
        self.assertEqual(format_oxidation_state(Fraction(8, 3)), "+8/3")


class TestNOSC(unittest.TestCase):
    CASES = [
        ("CO2", 4),
        ("CH4", -4),
        ("C2H3O2-", 0),
        ("C6H12O6", 0),
        ("CH2O", 0),
        ("CH4O", -2),
        ("CHO2-", 2),
    ]

    def test_reference_values(self):
        for formula, expected in self.CASES:
            with self.subTest(formula=formula):
                self.assertEqual(nosc(formula), expected)

    def test_agrees_with_mean_carbon_state(self):
        """NOSC is by construction the mean oxidation state of carbon."""
        for formula, _ in self.CASES:
            with self.subTest(formula=formula):
                self.assertEqual(nosc(formula), mean_oxidation_state("C", formula))

    def test_carbon_free_raises(self):
        with self.assertRaises(MicrobialThermoError):
            nosc("SO4-2")


class TestAtomOxidationStates(unittest.TestCase):
    def test_ethanol_conserves_charge(self):
        """The case where RDKit's own experimental routine gives +6, not 0."""
        states = atom_oxidation_states("CCO")
        self.assertEqual(sum(o for _, o in states), 0)
        carbons = [o for s, o in states if s == "C"]
        self.assertEqual(sorted(carbons), [-3, -1])

    def test_acetate_distinguishes_its_two_carbons(self):
        states = atom_oxidation_states("CC(=O)[O-]")
        self.assertEqual(sum(o for _, o in states), -1)
        self.assertEqual(sorted(o for s, o in states if s == "C"), [-3, 3])

    def test_mean_of_per_atom_matches_nosc(self):
        """Per-atom states averaged over carbon must reproduce NOSC."""
        for smiles, formula in [
            ("CCO", "C2H6O"),
            ("CC(=O)[O-]", "C2H3O2-"),
            ("C", "CH4"),
            ("O=C=O", "CO2"),
            ("CO", "CH4O"),
        ]:
            with self.subTest(smiles=smiles):
                carbons = [o for s, o in atom_oxidation_states(smiles) if s == "C"]
                mean = Fraction(sum(carbons), len(carbons))
                self.assertEqual(mean, nosc(formula))


if __name__ == "__main__":
    unittest.main()
