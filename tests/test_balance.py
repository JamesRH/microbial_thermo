"""Half-reaction and full-reaction balancing against hand-balanced references."""

import unittest
from fractions import Fraction

from microbial_thermo.balance import (
    ELECTRON,
    balance_equation,
    balance_half_reaction,
    combine_half_reactions,
    verify_conservation,
)
from microbial_thermo.exceptions import BalancingError


def coefficient_map(coefficients):
    return {s.backend: v for s, v in coefficients.items()}


class TestHalfReactionBalancing(unittest.TestCase):
    """Reference half reactions, balanced by hand, written as reductions."""

    CASES = [
        # (reduced, oxidized, n electrons, {backend name: signed coefficient})
        ("HS-", "SO4-2", 8, {"SO4--": -1, "H+": -9, "e-": -8, "HS-": 1, "H2O": 4}),
        ("H2O", "O2(aq)", 2, {"O2(aq)": Fraction(-1, 2), "H+": -2, "e-": -2, "H2O": 1}),
        ("H2(g)", "H+", 2, {"H+": -2, "e-": -2, "H2(g)": 1}),
        ("NH4+", "NO3-", 8, {"NO3-": -1, "H+": -10, "e-": -8, "NH4+": 1, "H2O": 3}),
        ("N2(aq)", "NO3-", 10, {"NO3-": -2, "H+": -12, "e-": -10, "N2(aq)": 1, "H2O": 6}),
        ("Fe+2", "Fe+3", 1, {"Fe+++": -1, "e-": -1, "Fe++": 1}),
        ("NO2-", "NO3-", 2, {"NO3-": -1, "H+": -2, "e-": -2, "NO2-": 1, "H2O": 1}),
    ]

    def test_reference_half_reactions(self):
        for reduced, oxidized, n, expected in self.CASES:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                half = balance_half_reaction(reduced, oxidized)
                self.assertEqual(half.n_electrons, n)
                self.assertEqual(coefficient_map(half.coefficients), expected)

    def test_all_conserve_atoms_and_charge(self):
        for reduced, oxidized, _, _ in self.CASES:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                half = balance_half_reaction(reduced, oxidized)
                verify_conservation(half.coefficients)

    def test_acetate_half_reaction(self):
        """2 CO2 + 7 H+ + 8 e- -> acetate + 2 H2O, the acetogenesis couple."""
        half = balance_half_reaction("acetate", "CO2(aq)", key_element="C")
        self.assertEqual(half.n_electrons, 8)
        self.assertEqual(
            coefficient_map(half.coefficients),
            {"CO2(aq)": -2, "H+": -7, "e-": -8, "Acetate": 1, "H2O": 2},
        )

    def test_non_couple_raises(self):
        with self.assertRaises(BalancingError):
            balance_half_reaction("SO4-2", "NO3-")


class TestNormalization(unittest.TestCase):
    def test_electron_pair_default_is_exact(self):
        """Normalising to two electrons must give exact quarters, not floats."""
        half = balance_half_reaction("HS-", "SO4-2").normalized_to_electrons(2)
        self.assertEqual(half.n_electrons, 2)
        self.assertEqual(coefficient_map(half.coefficients)["SO4--"], Fraction(-1, 4))

    def test_single_electron(self):
        half = balance_half_reaction("HS-", "SO4-2").normalized_to_electrons(1)
        self.assertEqual(half.n_electrons, 1)

    def test_normalize_to_named_species(self):
        half = balance_half_reaction("HS-", "SO4-2").normalized_to("SO4-2")
        self.assertEqual(abs(coefficient_map(half.coefficients)["SO4--"]), 1)

    def test_normalize_to_absent_species_raises(self):
        with self.assertRaises(BalancingError):
            balance_half_reaction("Fe+2", "Fe+3").normalized_to("SO4-2")

    def test_reversed_flips_every_sign(self):
        half = balance_half_reaction("HS-", "SO4-2")
        flipped = half.reversed()
        self.assertEqual(flipped.n_electrons, -half.n_electrons)
        for species, value in half.coefficients.items():
            self.assertEqual(flipped.coefficients[species], -value)


class TestCombining(unittest.TestCase):
    def test_methane_oxidation(self):
        """CH4 + 2 O2 -> CO2 + 2 H2O."""
        donor = balance_half_reaction("methane", "CO2(aq)")
        acceptor = balance_half_reaction("H2O", "O2(aq)")
        combined = combine_half_reactions(donor, acceptor, n_electrons=8)
        self.assertEqual(
            coefficient_map(combined),
            {"Methane(aq)": -1, "O2(aq)": -2, "CO2(aq)": 1, "H2O": 2},
        )

    def test_electrons_cancel(self):
        donor = balance_half_reaction("H2(g)", "H+")
        acceptor = balance_half_reaction("HS-", "SO4-2")
        combined = combine_half_reactions(donor, acceptor, n_electrons=8)
        self.assertNotIn(ELECTRON, combined)
        verify_conservation(combined)


class TestEquationBalancing(unittest.TestCase):
    def test_balances_a_written_equation(self):
        coefficients = balance_equation("methane + O2(aq) -> CO2(aq) + H2O")
        self.assertEqual(
            coefficient_map(coefficients),
            {"Methane(aq)": -1, "O2(aq)": -2, "CO2(aq)": 1, "H2O": 2},
        )

    def test_missing_arrow_raises(self):
        with self.assertRaises(BalancingError):
            balance_equation("methane + O2(aq)")

    def test_accepts_unicode_arrow(self):
        coefficients = balance_equation("methane + O2(aq) → CO2(aq) + H2O")
        self.assertEqual(coefficient_map(coefficients)["Methane(aq)"], -1)


if __name__ == "__main__":
    unittest.main()
