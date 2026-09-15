"""Figure layout and rendering.

Layout is asserted on the token model rather than on pixels: the point of
building layout as a data structure was to make it testable without
rendering. A couple of smoke tests confirm the renderers actually produce
files.
"""

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import microbial_thermo as mt  # noqa: E402
from microbial_thermo.typeset import (  # noqa: E402
    ARROW,
    align_at_arrows,
    build_equation_tokens,
)


class TestEquationLayout(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
            n_electrons=8,
        )

    def test_layout_has_exactly_one_arrow(self):
        layout = build_equation_tokens(self.reaction.donor_half.half, "oxidation")
        arrows = [t for t in layout.tokens if t.kind == ARROW]
        self.assertEqual(len(arrows), 1)
        self.assertIs(arrows[0], layout.arrow)

    def test_oxidation_state_annotations_are_attached(self):
        half = self.reaction.acceptor_half.half
        layout = build_equation_tokens(half, "reduction", annotate_element="S")
        annotated = [t for t in layout.species_tokens() if t.oxidation_state]
        self.assertTrue(annotated)
        states = {t.oxidation_state for t in annotated}
        self.assertEqual(states, {"+6", "-2"})

    def test_only_species_with_the_element_are_annotated(self):
        half = self.reaction.acceptor_half.half
        layout = build_equation_tokens(half, "reduction", annotate_element="S")
        for token in layout.species_tokens():
            with self.subTest(species=token.species.backend):
                has_sulfur = "S" in token.species.parsed.elements
                self.assertEqual(bool(token.oxidation_state), has_sulfur)

    def test_arrows_align_across_stacked_equations(self):
        """The whole point of the layout model: arrows share an x position."""
        top = build_equation_tokens(self.reaction.donor_half.half, "oxidation")
        bottom = build_equation_tokens(self.reaction.acceptor_half.half, "reduction")
        for layout in (top, bottom):
            for index, token in enumerate(layout.tokens):
                token.width = 0.05 + 0.01 * index  # deterministic stand-in widths
        align_at_arrows([top, bottom])
        self.assertAlmostEqual(top.arrow.x, bottom.arrow.x, places=9)

    def test_alignment_leaves_no_equation_overflowing_left(self):
        top = build_equation_tokens(self.reaction.donor_half.half, "oxidation")
        bottom = build_equation_tokens(self.reaction.acceptor_half.half, "reduction")
        for layout in (top, bottom):
            for token in layout.tokens:
                token.width = 0.05
        align_at_arrows([top, bottom])
        for layout in (top, bottom):
            self.assertGreaterEqual(min(t.x for t in layout.tokens), -1e-9)

    def test_direction_reverses_the_equation(self):
        half = self.reaction.donor_half.half
        oxidative = build_equation_tokens(half, "oxidation")
        reductive = build_equation_tokens(half, "reduction")
        first_ox = oxidative.tokens[0]
        first_red = reductive.tokens[0]
        self.assertNotEqual(first_ox.text, first_red.text)


class TestRendering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
            normalize_to="integer",
        )

    def test_half_reaction_figure_writes_svg_and_png(self):
        from microbial_thermo.figures import plot_half_reactions

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "half"
            plot_half_reactions(self.reaction, save=base)
            self.assertTrue(base.with_suffix(".svg").exists())
            self.assertTrue(base.with_suffix(".png").exists())
            self.assertGreater(base.with_suffix(".svg").stat().st_size, 1000)

    def test_tower_figure_writes_svg_and_png(self):
        from microbial_thermo.figures import plot_redox_tower

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "tower"
            plot_redox_tower(self.reaction, save=base)
            self.assertTrue(base.with_suffix(".svg").exists())
            self.assertTrue(base.with_suffix(".png").exists())


class TestReferenceCouples(unittest.TestCase):
    def test_tower_is_ordered_and_spans_the_expected_range(self):
        from microbial_thermo.tower import reference_couples

        entries = reference_couples(mt.Conditions(temperature_c=25.0, pH=7.0))
        self.assertGreater(len(entries), 5)
        values = [e.potential_v for e in entries]
        self.assertEqual(values, sorted(values))
        # The classic tower runs from roughly -0.5 V to +0.9 V at pH 7.
        self.assertLess(values[0], -0.3)
        self.assertGreater(values[-1], 0.7)

    def test_oxygen_is_the_strongest_acceptor(self):
        from microbial_thermo.tower import reference_couples

        entries = reference_couples(mt.Conditions(temperature_c=25.0, pH=7.0))
        self.assertEqual(entries[-1].label, "O2/H2O")


class TestAtpScale(unittest.TestCase):
    def test_atp_equivalents_is_a_linear_rescaling(self):
        from microbial_thermo.figures.style import atp_equivalents
        from microbial_thermo.units import Quantity

        self.assertAlmostEqual(
            atp_equivalents(Quantity(-100.0, "kJ/mol"), Quantity(-50.0, "kJ/mol")),
            2.0,
        )

    def test_endergonic_reaction_gives_negative_atp(self):
        from microbial_thermo.figures.style import atp_equivalents
        from microbial_thermo.units import Quantity

        self.assertLess(atp_equivalents(Quantity(+30.0, "kJ/mol"), Quantity(-50.0, "kJ/mol")), 0)


if __name__ == "__main__":
    unittest.main()


class TestElectronPairNormalization(unittest.TestCase):
    """Figures normalise to an electron pair by default.

    The brief specifies the half-reaction diagram around a transferred
    electron pair, with everything else balanced to match, and the tower
    follows the same convention so its energy scale bar is meaningful.
    """

    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def _built_with_eight(self):
        return mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            n_electrons=8,
        )

    def test_renormalized_changes_stoichiometry_not_potentials(self):
        octet = self._built_with_eight()
        pair = octet.renormalized(2)
        self.assertEqual(pair.n_electrons, 2)
        self.assertAlmostEqual(
            pair.delta_E_standard_prime.to("V").magnitude,
            octet.delta_E_standard_prime.to("V").magnitude,
            places=9,
        )
        self.assertAlmostEqual(
            pair.delta_G_standard_prime.magnitude * 4.0,
            octet.delta_G_standard_prime.magnitude,
            places=6,
        )

    def test_half_reaction_figure_defaults_to_a_pair(self):
        from microbial_thermo.figures import plot_half_reactions

        figure = plot_half_reactions(self._built_with_eight())
        title = figure.axes[0].get_title()
        self.assertIn("2 e", title)

    def test_tower_defaults_to_a_pair(self):
        from microbial_thermo.figures import plot_redox_tower

        figure = plot_redox_tower(self._built_with_eight())
        self.assertIn("2 e", figure.axes[0].get_title())

    def test_passing_none_leaves_the_reaction_alone(self):
        from microbial_thermo.figures import plot_half_reactions

        figure = plot_half_reactions(self._built_with_eight(), n_electrons=None)
        self.assertIn("8 e", figure.axes[0].get_title())


class TestEnergyScaleBar(unittest.TestCase):
    """The tower's map-style scale bar.

    Valid only because n is pinned: with n fixed, a difference in potential
    maps to free energy by the constant nF.
    """

    def test_volts_per_kilojoule_at_two_electrons(self):
        from microbial_thermo.figures.tower import volts_per_kilojoule

        # 1 kJ/mol / (2 * 96.485 kJ/mol/V) = 5.182 mV
        self.assertAlmostEqual(volts_per_kilojoule(2), 0.005182, places=6)

    def test_energy_quantum_and_atp_marks(self):
        from microbial_thermo.figures.tower import volts_per_kilojoule

        per_kj = volts_per_kilojoule(2)
        self.assertAlmostEqual(20 * per_kj, 0.1036, places=4)
        self.assertAlmostEqual(50 * per_kj, 0.2591, places=4)

    def test_scale_halves_when_electrons_double(self):
        """Twice the electrons, half the potential span per kJ/mol."""
        from microbial_thermo.figures.tower import volts_per_kilojoule

        self.assertAlmostEqual(volts_per_kilojoule(4), volts_per_kilojoule(2) / 2.0, places=9)

    def test_scale_bar_is_consistent_with_the_reported_atp_yield(self):
        """Reading the donor-acceptor gap against the bar must agree with the
        ATP figure quoted in the energy panel."""
        from microbial_thermo.figures.style import atp_equivalents
        from microbial_thermo.figures.tower import volts_per_kilojoule

        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
        )
        gap_v = reaction.delta_E_standard_prime.to("V").magnitude
        from_bar = gap_v / volts_per_kilojoule(reaction.n_electrons) / 50.0
        from_panel = atp_equivalents(reaction.delta_G_standard_prime, mt.Quantity(-50.0, "kJ/mol"))
        self.assertAlmostEqual(from_bar, from_panel, places=6)
