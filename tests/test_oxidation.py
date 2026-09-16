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


class TestPerAtomDisplay(unittest.TestCase):
    """Per-atom states for the half-reaction figure.

    A mean hides real chemistry: acetate's two carbons are four units apart
    and average to zero, which describes neither of them.
    """

    def test_acetate_shows_both_carbons(self):
        from microbial_thermo.oxidation import format_per_atom_states

        self.assertEqual(format_per_atom_states("CC(=O)[O-]", "C"), "-3, +3")

    def test_repeats_are_collapsed_with_a_count(self):
        from microbial_thermo.oxidation import format_per_atom_states

        # butyrate: methyl, two methylenes, carboxyl
        self.assertEqual(format_per_atom_states("CCCC(=O)[O-]", "C"), "-3, -2(x2), +3")

    def test_a_single_carbon_reads_like_the_mean(self):
        from microbial_thermo.oxidation import format_per_atom_states

        self.assertEqual(format_per_atom_states("CO", "C"), "-2")

    def test_states_are_sorted(self):
        from microbial_thermo.oxidation import per_atom_states

        self.assertEqual(per_atom_states("CCC(=O)[O-]", "C"), [-3, -2, 3])

    def test_the_mean_of_the_per_atom_states_is_the_nosc(self):
        """The two views must agree: averaging the per-atom states has to give
        back what the formula-only route reports."""
        from fractions import Fraction

        from microbial_thermo.oxidation import nosc, per_atom_states
        from microbial_thermo.species import default_registry

        for species in default_registry().all_species():
            if not species.smiles or "C" not in species.parsed.elements:
                continue
            with self.subTest(species=species.backend):
                states = per_atom_states(species.smiles, "C")
                mean = Fraction(sum(states), len(states))
                self.assertEqual(mean, nosc(species.formula))

    def test_missing_element_raises(self):
        from microbial_thermo.oxidation import format_per_atom_states

        with self.assertRaises(MicrobialThermoError):
            format_per_atom_states("CO", "N")


class TestPerAtomOnTheFigure(unittest.TestCase):
    def setUp(self):
        import matplotlib

        matplotlib.use("Agg")
        import microbial_thermo as mt

        self.mt = mt
        self.reaction = mt.Reaction.from_couples(
            donor=("acetate", "CO2(aq)"),
            acceptor=("H2O", "O2(aq)"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
        )

    def _annotations(self, per_atom):
        from microbial_thermo.typeset import build_equation_tokens

        layout = build_equation_tokens(
            self.reaction.donor_half.half,
            "oxidation",
            annotate_element="C",
            per_atom=per_atom,
        )
        return {
            t.species.backend: t.oxidation_state
            for t in layout.species_tokens()
            if t.oxidation_state
        }

    def test_the_mean_is_still_the_default(self):
        self.assertEqual(self._annotations(False)["Acetate"], "0")

    def test_per_atom_replaces_it_when_asked(self):
        self.assertEqual(self._annotations(True)["Acetate"], "-3, +3")

    def test_species_without_a_smiles_keep_the_mean(self):
        """CO2 has no SMILES in the registry, so it cannot do better than
        the mean -- and must not break the figure."""
        self.assertEqual(self._annotations(True)["CO2(aq)"], "+4")

    def test_the_figure_renders_either_way(self):
        from microbial_thermo.figures import plot_half_reactions

        for per_atom in (False, True):
            with self.subTest(per_atom=per_atom):
                plot_half_reactions(self.reaction, per_atom=per_atom)


if __name__ == "__main__":
    unittest.main()
