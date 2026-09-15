"""Formula parsing, including both charge notations pyGCC and users each use."""

import unittest

from microbial_thermo.exceptions import MicrobialThermoError
from microbial_thermo.formula import (
    format_charge,
    normalise_charge_notation,
    parse_formula,
    strip_charge,
)


class TestChargeNotation(unittest.TestCase):
    def test_doubled_signs_normalise(self):
        self.assertEqual(normalise_charge_notation("SO4--"), "SO4-2")
        self.assertEqual(normalise_charge_notation("Fe+++"), "Fe+3")

    def test_single_sign_unchanged(self):
        self.assertEqual(normalise_charge_notation("NH4+"), "NH4+")

    def test_both_notations_agree(self):
        """The legacy database spelling and the modern one must parse alike."""
        for legacy, modern in [("SO4--", "SO4-2"), ("Fe++", "Fe+2"), ("CO3--", "CO3-2")]:
            with self.subTest(legacy=legacy):
                self.assertEqual(parse_formula(legacy).elements, parse_formula(modern).elements)
                self.assertEqual(parse_formula(legacy).charge, parse_formula(modern).charge)


class TestParseFormula(unittest.TestCase):
    def test_neutral(self):
        parsed = parse_formula("C6H12O6")
        self.assertEqual(parsed.elements, {"C": 6, "H": 12, "O": 6})
        self.assertEqual(parsed.charge, 0)

    def test_anion(self):
        parsed = parse_formula("SO4-2")
        self.assertEqual(parsed.elements, {"O": 4, "S": 1})
        self.assertEqual(parsed.charge, -2)

    def test_electron(self):
        parsed = parse_formula("e-")
        self.assertTrue(parsed.is_electron)
        self.assertEqual(parsed.charge, -1)
        self.assertEqual(parsed.elements, {})

    def test_empty_raises(self):
        with self.assertRaises(MicrobialThermoError):
            parse_formula("   ")

    def test_nonsense_raises(self):
        with self.assertRaises(MicrobialThermoError):
            parse_formula("Methane(aq)")  # a database name, not a formula

    def test_strip_charge(self):
        self.assertEqual(strip_charge("SO4--"), "SO4")
        self.assertEqual(strip_charge("H2O"), "H2O")

    def test_format_charge(self):
        self.assertEqual(format_charge(-2), "2-")
        self.assertEqual(format_charge(1), "+")
        self.assertEqual(format_charge(0), "")


if __name__ == "__main__":
    unittest.main()
