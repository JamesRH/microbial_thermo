"""Two kinds of ambiguity the library refuses to resolve silently.

A name two species both claim, and an equation conservation cannot balance
uniquely. In both cases guessing would be easy and wrong, so the choice is
either declared in the data or asked of the caller.
"""

import unittest

import microbial_thermo as mt
from microbial_thermo.balance import balance_equation, verify_conservation
from microbial_thermo.exceptions import AmbiguousReactionError, BalancingError
from microbial_thermo.species import default_registry, resolve


class TestPhaseAmbiguity(unittest.TestCase):
    """A bare name shared by a dissolved and a gaseous entry."""

    @classmethod
    def setUpClass(cls):
        cls.registry = default_registry()

    def test_nothing_resolves_by_file_ordering(self):
        """Every contested name must have a declared winner. Adding a new
        aqueous/gas pair without saying which one a bare name means should
        fail here rather than depending on where it sits in the YAML."""
        undeclared = self.registry.undeclared_ambiguities
        self.assertEqual(undeclared, {}, f"undeclared ambiguities: {sorted(undeclared)}")

    def test_the_contested_names_are_the_ones_expected(self):
        self.assertEqual(sorted(self.registry.ambiguities()), ["ch4", "co2", "h2", "n2", "o2"])

    def test_a_bare_name_means_the_dissolved_form(self):
        """What a cell actually sees. The gas needs saying explicitly."""
        for name, expected in [
            ("H2", "H2(aq)"),
            ("O2", "O2(aq)"),
            ("N2", "N2(aq)"),
            ("CO2", "CO2(aq)"),
            ("CH4", "Methane(aq)"),
        ]:
            with self.subTest(name=name):
                self.assertEqual(resolve(name).backend, expected)
                self.assertEqual(resolve(name).phase, "aq")

    def test_phase_can_be_demanded(self):
        self.assertEqual(resolve("H2", phase="g").backend, "H2(g)")
        self.assertEqual(resolve("H2", phase="aq").backend, "H2(aq)")
        self.assertEqual(resolve("CO2", phase="g").backend, "CO2(g)")

    def test_an_impossible_phase_says_what_is_available(self):
        from microbial_thermo.exceptions import SpeciesNotFoundError

        with self.assertRaises(SpeciesNotFoundError) as caught:
            resolve("H2", phase="s")
        self.assertIn("aq, g", str(caught.exception))

    def test_the_gas_and_the_dissolved_form_differ_enough_to_matter(self):
        """91 mV apart at pH 7 -- worth getting right."""
        conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
        aqueous = mt.half_reaction("H2(aq)", "H+", conditions)
        gas = mt.half_reaction("H2(g)", "H+", conditions)
        gap = abs(
            aqueous.E_standard_prime.to("V").magnitude - gas.E_standard_prime.to("V").magnitude
        )
        self.assertGreater(gap, 0.05)


class TestUnderdeterminedEquations(unittest.TestCase):
    """Equations conservation alone cannot balance."""

    AMBIGUOUS = "HS- + O2 -> SO4-2 + Sulfur(s) + H2O + H+"

    def test_a_determined_equation_still_balances_unchanged(self):
        coefficients = balance_equation("methane + O2(aq) -> CO2(aq) + H2O")
        self.assertEqual(
            {s.backend: v for s, v in coefficients.items()},
            {"Methane(aq)": -1, "O2(aq)": -2, "CO2(aq)": 1, "H2O": 2},
        )

    def test_it_refuses_and_says_how_many_solutions(self):
        with self.assertRaises(AmbiguousReactionError) as caught:
            balance_equation(self.AMBIGUOUS)
        self.assertIn("2 independent balanced solutions", str(caught.exception))

    def test_the_error_carries_a_readable_basis(self):
        """Not raw sympy vectors -- coefficient sets you can look at."""
        with self.assertRaises(AmbiguousReactionError) as caught:
            balance_equation(self.AMBIGUOUS)
        basis = caught.exception.basis
        self.assertEqual(len(basis), 2)
        for option in basis:
            self.assertTrue(all(hasattr(s, "backend") for s in option))
            verify_conservation(option)

    def test_the_error_names_the_way_out(self):
        with self.assertRaises(AmbiguousReactionError) as caught:
            balance_equation(self.AMBIGUOUS)
        message = str(caught.exception)
        self.assertIn("fix=", message)
        self.assertIn("from_couples", message)

    def test_constraints_resolve_it(self):
        coefficients = balance_equation(self.AMBIGUOUS, fix={"HS-": 2, "Sulfur(s)": 1})
        verify_conservation(coefficients)
        by_name = {s.backend: v for s, v in coefficients.items()}
        self.assertEqual(by_name["HS-"], -2)  # reactant, so negative
        self.assertEqual(by_name["Sulfur(s)"], 1)

    def test_the_sign_follows_the_side_it_was_written_on(self):
        """Moles are given as positive numbers; the caller should not have to
        think about which side carries the minus."""
        coefficients = balance_equation(self.AMBIGUOUS, fix={"HS-": 2, "Sulfur(s)": 1})
        by_name = {s.backend: v for s, v in coefficients.items()}
        self.assertLess(by_name["HS-"], 0)  # written as a reactant
        self.assertGreater(by_name["Sulfur(s)"], 0)  # written as a product

    def test_too_few_constraints_says_how_many_are_needed(self):
        with self.assertRaises(AmbiguousReactionError) as caught:
            balance_equation(self.AMBIGUOUS, fix={"O2": 2})
        message = str(caught.exception)
        self.assertIn("degree(s) of freedom remain", message)
        self.assertIn("needs 2 constraints", message)

    def test_fixing_a_species_that_is_not_there_is_reported(self):
        with self.assertRaises(BalancingError) as caught:
            balance_equation(self.AMBIGUOUS, fix={"NO3-": 1, "HS-": 1})
        self.assertIn("does not appear", str(caught.exception))

    def test_an_inconsistent_constraint_is_reported(self):
        """Demanding a stoichiometry conservation forbids must fail, not
        silently return something else."""
        with self.assertRaises(BalancingError):
            balance_equation(self.AMBIGUOUS, fix={"O2": 2, "Sulfur(s)": 1})

    def test_different_constraints_give_different_chemistry(self):
        """Which is the whole point: the caller is choosing a reaction."""
        first = balance_equation(self.AMBIGUOUS, fix={"HS-": 2, "Sulfur(s)": 1})
        second = balance_equation(self.AMBIGUOUS, fix={"O2": 3, "Sulfur(s)": 1})
        self.assertNotEqual(
            {s.backend: v for s, v in first.items()},
            {s.backend: v for s, v in second.items()},
        )
        for option in (first, second):
            verify_conservation(option)


if __name__ == "__main__":
    unittest.main()
