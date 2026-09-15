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


class TestAffinityLadder(unittest.TestCase):
    """The ladder over the curated library."""

    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(
            temperature_c=12.0,
            pH=7.4,
            activity_model="ideal",
            total_concentrations={"sulfide": 1e-5, "DIC": 3e-3},
            concentrations={"SO4-2": 2.0e-2, "acetate": 1e-5, "O2(aq)": 1e-6},
            partial_pressures={"H2(g)": 5e-6, "CH4(g)": 1e-3},
        )

    def test_ladder_renders_and_saves(self):
        import warnings

        from microbial_thermo.figures import plot_affinity_ladder

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "ladder"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                plot_affinity_ladder(self.conditions, save=base)
            self.assertTrue(base.with_suffix(".svg").exists())
            self.assertTrue(base.with_suffix(".png").exists())

    def test_ladder_accepts_a_precomputed_table(self):
        import warnings

        from microbial_thermo.figures import plot_affinity_ladder
        from microbial_thermo.library import energy_table

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = energy_table(self.conditions)
            figure = plot_affinity_ladder(self.conditions, table=table)
        # One bar per evaluable metabolism. Count the bar container rather than
        # every patch, since the shaded energy-quantum band is a patch too.
        bars = figure.axes[0].containers[0]
        self.assertEqual(len(bars), len(table[table["problem"] == ""]))

    def test_ladder_has_a_true_atp_axis(self):
        """The x axis is kJ/mol here, so ATP really is a second axis --
        unlike the tower, where the axis is a potential."""
        import warnings

        from microbial_thermo.figures import plot_affinity_ladder

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            figure = plot_affinity_ladder(self.conditions)
        self.assertGreaterEqual(len(figure.axes), 2)
        twin = figure.axes[-1]
        self.assertIn("ATP", twin.get_xlabel())

    def test_empty_selection_raises(self):
        from microbial_thermo.figures import plot_affinity_ladder
        from microbial_thermo.library import energy_table

        empty = energy_table(self.conditions).iloc[0:0]
        with self.assertRaises(ValueError):
            plot_affinity_ladder(self.conditions, table=empty)


class TestScaleBarShowsBothConventions(unittest.TestCase):
    """The tower's scale bar carries per-reaction and per-electron readings.

    The potential axis is intensive -- every rung stays put under
    renormalisation -- but the reference quantities (one ATP, the energy
    quantum) are per-reaction, so their position on that axis depends on n.
    Showing both makes that visible rather than leaving it to a caption.
    """

    @staticmethod
    def _texts(reaction, **kwargs):
        import warnings as _w

        from microbial_thermo.figures import plot_redox_tower

        with _w.catch_warnings():
            _w.simplefilter("ignore")
            figure = plot_redox_tower(reaction, **kwargs)
        return [t.get_text() for ax in figure.axes for t in ax.texts]

    @classmethod
    def setUpClass(cls):
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
        )

    def test_per_electron_conversion_does_not_depend_on_n(self):
        """dG/n = F dE, so the per-electron reading is the same at any n."""
        from microbial_thermo.figures.tower import volts_per_kilojoule

        per_electron = volts_per_kilojoule(1)
        self.assertAlmostEqual(50 * per_electron, 0.5182, places=4)
        self.assertAlmostEqual(20 * per_electron, 0.2073, places=4)

    def test_both_labels_are_drawn(self):
        joined = "\n".join(self._texts(self.reaction))
        self.assertIn("per reaction", joined)
        self.assertIn("per electron", joined)

    def test_per_reaction_marks_move_with_n_but_per_electron_do_not(self):
        two = "\n".join(self._texts(self.reaction, n_electrons=2))
        six = "\n".join(self._texts(self.reaction, n_electrons=6))
        # Per-reaction ATP mark: 0.259 V at n=2, 0.086 V at n=6.
        self.assertIn("0.259", two)
        self.assertIn("0.086", six)
        # Per-electron reading is identical in both.
        for text in (two, six):
            self.assertIn("0.518", text)
            self.assertIn("0.207", text)

    def test_the_two_conventions_coincide_at_one_electron(self):
        from microbial_thermo.figures.tower import volts_per_kilojoule

        self.assertAlmostEqual(50 * volts_per_kilojoule(1), 50 * volts_per_kilojoule(1), places=9)
        one = "\n".join(self._texts(self.reaction, n_electrons=1))
        self.assertEqual(one.count("0.518"), 2)  # per reaction and per electron


if __name__ == "__main__":
    unittest.main()
