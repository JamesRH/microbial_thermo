"""Frost-Ebsworth diagrams, and the volt-equivalent convention behind them.

The convention for *negative* oxidation states was the thing that held this
up, so it is what most of these tests are about. The claim being checked is
that there is no special convention at all: the volt equivalent is the free
energy of formation from the element, per mole of the element, divided by the
Faraday, and a hydride falls out of that with no sign rule applied.

Nitrogen settles it, because its textbook Frost diagram is one of the most
reproduced in chemistry and ammonium sits at a negative state:

    NH4+   -3   -0.82 V
    N2      0    0
    HNO2   +3   +4.39 V
    NO3-   +5   +6.22 V

All four come out of this module to within a few millivolts, and finding the
nitrogen zero point wrong -- N2 dissolved rather than N2 gas, which shifts
every point by 0.094 V -- is how the reference question was found at all.
"""

import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from microbial_thermo.figures.basis import FARADAY_KJ, element_series
from microbial_thermo.figures.frost import frost_diagram, plot_frost
from microbial_thermo.figures.latimer import latimer_diagram

#: The published nitrogen Frost diagram at pH 0, in volts.
NITROGEN = {"NH4+": -0.82, "N2(g)": 0.0, "HNO2(aq)": 4.39, "NO3-": 6.22}

NITROGEN_SPECIES = ["N2(g)", "NO3-", "NH4+", "HNO2(aq)"]


class TestTheVoltEquivalentConvention(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.nitrogen = frost_diagram("N", species=NITROGEN_SPECIES, pH=0.0, activity=1.0)

    def test_every_nitrogen_point_matches_the_published_diagram(self):
        found = {p.backend: p.volt_equivalent for p in self.nitrogen.points}
        for name, expected in NITROGEN.items():
            with self.subTest(species=name):
                self.assertAlmostEqual(found[name], expected, delta=0.02)

    def test_ammonium_is_negative_with_no_sign_rule(self):
        """The case that blocked this module.

        A hydride has a *positive* electron count in its formation from the
        element -- N2 + 8H+ + 6e- -> 2NH4+ -- so its volt equivalent comes
        out negative on its own. Nothing anywhere treats z < 0 specially.
        """
        ammonium = [p for p in self.nitrogen.points if p.backend == "NH4+"][0]
        self.assertLess(ammonium.oxidation_state, 0)
        self.assertGreater(ammonium.entry.electrons, 0)
        self.assertAlmostEqual(ammonium.volt_equivalent, -0.82, delta=0.02)

    def test_sulfide_is_negative_too(self):
        # Published -0.288 V: S + 2H+ + 2e- -> H2S at E = +0.144 V, twice.
        sulfur = frost_diagram("S", pH=0.0, activity=1.0)
        sulfide = [p for p in sulfur.points if p.backend == "H2S(aq)"][0]
        self.assertAlmostEqual(sulfide.volt_equivalent, -0.288, delta=0.01)

    def test_the_element_is_the_zero(self):
        # Nitrogen is left out deliberately: its z = 0 species here is
        # *dissolved* N2, which really does sit 0.094 V above the gas that
        # sets the zero. That case is its own test below.
        for element in ("S", "As", "Se", "Cu"):
            with self.subTest(element=element):
                diagram = frost_diagram(element, pH=0.0, activity=1.0)
                zeros = [p for p in diagram.points if abs(p.oxidation_state) < 1e-9]
                self.assertTrue(zeros)
                self.assertAlmostEqual(zeros[0].volt_equivalent, 0.0, delta=0.03)

    def test_the_volt_equivalent_is_the_free_energy_over_the_faraday(self):
        for point in self.nitrogen.points:
            with self.subTest(species=point.backend):
                self.assertAlmostEqual(point.volt_equivalent, point.gibbs / FARADAY_KJ, places=9)


class TestTheReferenceMatters(unittest.TestCase):
    """On an Eh-pH diagram the element reference cancels. Here it does not."""

    def test_dissolved_nitrogen_as_the_zero_shifts_every_point(self):
        gas = frost_diagram("N", species=NITROGEN_SPECIES, pH=0.0, activity=1.0)
        dissolved = frost_diagram(
            "N",
            species=["N2(aq)", "NO3-", "NH4+", "HNO2(aq)"],
            reference="N2(aq)",
            pH=0.0,
            activity=1.0,
        )
        self.assertEqual(gas.reference, "N2(g)")
        self.assertEqual(dissolved.reference, "N2(aq)")

        # Compare the species the two have in common; the zero point itself
        # is a different species in each, so it is not one of them.
        by_name = {p.backend: p.volt_equivalent for p in gas.points}
        shifts = [
            by_name[p.backend] - p.volt_equivalent for p in dissolved.points if p.backend in by_name
        ]
        self.assertEqual(len(shifts), 3)
        # One common offset, and it is the dissolution energy of N2 halved.
        self.assertAlmostEqual(max(shifts), min(shifts), places=6)
        self.assertAlmostEqual(max(shifts), 0.094, delta=0.005)

    def test_slopes_survive_a_bad_reference_even_though_heights_do_not(self):
        gas = frost_diagram("N", species=NITROGEN_SPECIES, pH=0.0)
        dissolved = frost_diagram(
            "N", species=["N2(aq)", "NO3-", "NH4+", "HNO2(aq)"], reference="N2(aq)", pH=0.0
        )
        self.assertAlmostEqual(
            gas.slope("NO3-", "HNO2(aq)"), dissolved.slope("NO3-", "HNO2(aq)"), places=9
        )

    def test_a_non_elemental_reference_is_reported(self):
        # Chromium metal is in none of our databases, so the fallback picks an
        # arbitrary zero and the figure has to say so.
        diagram = frost_diagram("Cr", pH=0.0)
        self.assertFalse(diagram.reference_is_element)

    def test_an_elemental_reference_is_reported(self):
        self.assertTrue(frost_diagram("As", pH=0.0).reference_is_element)

    def test_dissolved_nitrogen_sits_above_the_gas_that_is_the_zero(self):
        """The default nitrogen diagram's z = 0 point is not at zero.

        ``N2(aq)`` is the species a biologist cares about and ``N2(g)`` is
        the standard state, and they differ by the dissolution energy. The
        figure states which zero it used rather than quietly moving it.
        """
        diagram = frost_diagram("N", pH=0.0, activity=1.0)
        self.assertEqual(diagram.reference, "N2(g)")
        dissolved = [p for p in diagram.points if p.backend == "N2(aq)"][0]
        self.assertAlmostEqual(dissolved.volt_equivalent, 0.094, delta=0.005)


class TestSlopesAreCouplePotentials(unittest.TestCase):
    def test_every_slope_matches_the_latimer_arrow(self):
        for element in ("As", "S", "Se", "Cu", "Mn", "Fe", "C"):
            frost = frost_diagram(element, pH=7.0)
            ladder = latimer_diagram(element, pH=7.0)
            for step in ladder.steps:
                with self.subTest(element=element, step=step.oxidized.backend):
                    self.assertAlmostEqual(
                        frost.slope(step.oxidized.backend, step.reduced.backend),
                        step.potential,
                        places=9,
                    )

    def test_the_same_state_twice_is_refused(self):
        diagram = frost_diagram("As", pH=0.0, predominant_only=False)
        with self.assertRaises(ValueError):
            diagram.slope("H3AsO4(aq)", "H2AsO4-")


class TestDisproportionation(unittest.TestCase):
    def test_copper_one_disproportionates(self):
        diagram = frost_diagram("Cu", pH=0.0, activity=1.0)
        found = {d.species: d for d in diagram.disproportionation()}
        self.assertIn("Cu+", found)
        self.assertEqual(set(found["Cu+"].into), {"Cu++", "Cu"})
        # 2Cu+ -> Cu2+ + Cu releases about 35 kJ, so half that per Cu+.
        self.assertAlmostEqual(found["Cu+"].delta_g, -17.2, delta=2.0)
        self.assertAlmostEqual(sum(found["Cu+"].fractions), 1.0, places=6)

    def test_it_is_exergonic_by_construction(self):
        for element in ("Cu", "Mn", "Fe", "S", "N"):
            for event in frost_diagram(element, pH=0.0).disproportionation():
                with self.subTest(element=element, species=event.species):
                    self.assertLess(event.delta_g, 0.0)

    def test_manganese_three_is_unstable_in_acid(self):
        # Textbook: Mn(III) disproportionates to Mn(II) and MnO2.
        diagram = frost_diagram("Mn", pH=0.0, activity=1.0)
        unstable = {d.species for d in diagram.disproportionation()}
        self.assertIn("Bixbyite", unstable)

    def test_nothing_on_the_hull_disproportionates(self):
        diagram = frost_diagram("Mn", pH=0.0)
        hull = {p.backend for p in diagram.stable}
        for event in diagram.disproportionation():
            self.assertNotIn(event.species, hull)

    def test_the_latimer_flag_implies_a_point_off_the_hull(self):
        """The two tests are the same test, applied locally and globally.

        A point above the chord of its immediate neighbours is necessarily
        above the hull, since the hull never rises above a chord. The
        converse does not hold, so this is an implication, not an equality.
        """
        for element in ("Cu", "Mn", "Fe", "S", "N", "C"):
            ladder = latimer_diagram(element, pH=0.0)
            hull = {p.backend for p in frost_diagram(element, pH=0.0).stable}
            for rung in ladder.disproportionating:
                with self.subTest(element=element, species=rung.backend):
                    self.assertNotIn(rung.backend, hull)


class TestTheHull(unittest.TestCase):
    def test_the_most_stable_point_is_the_lowest(self):
        diagram = frost_diagram("Fe", pH=0.0)
        self.assertEqual(diagram.most_stable.volt_equivalent, min(diagram.volt_equivalents))

    def test_the_hull_ends_at_the_extreme_states(self):
        diagram = frost_diagram("Mn", pH=0.0)
        self.assertEqual(diagram.stable[0].oxidation_state, min(diagram.states))
        self.assertEqual(diagram.stable[-1].oxidation_state, max(diagram.states))

    def test_the_hull_is_convex(self):
        diagram = frost_diagram("Mn", pH=0.0)
        hull = diagram.stable
        slopes = [
            (b.volt_equivalent - a.volt_equivalent) / (b.oxidation_state - a.oxidation_state)
            for a, b in zip(hull, hull[1:], strict=False)
        ]
        self.assertEqual(slopes, sorted(slopes))


class TestConditions(unittest.TestCase):
    def test_the_diagram_moves_with_pH(self):
        acid = frost_diagram("As", pH=0.0)
        neutral = frost_diagram("As", pH=7.0)
        by_state = {round(p.oxidation_state, 3): p.volt_equivalent for p in acid.points}
        for point in neutral.points:
            if abs(point.oxidation_state) < 1e-9:
                continue
            self.assertNotAlmostEqual(
                point.volt_equivalent, by_state[round(point.oxidation_state, 3)], places=3
            )

    def test_the_diagram_moves_with_temperature(self):
        cold = frost_diagram("As", pH=7.0, temperature_c=2.0)
        warm = frost_diagram("As", pH=7.0, temperature_c=60.0)
        self.assertNotAlmostEqual(
            cold.points[-1].volt_equivalent, warm.points[-1].volt_equivalent, places=3
        )

    def test_predominant_only_keeps_one_species_per_state(self):
        one = frost_diagram("As", pH=7.0, predominant_only=True)
        every = frost_diagram("As", pH=7.0, predominant_only=False)
        self.assertEqual(len(one.states), len(set(one.states)))
        self.assertGreater(len(every.points), len(one.points))

    def test_the_predominant_form_is_the_one_the_ladder_chose(self):
        frost = frost_diagram("As", pH=7.0)
        ladder = latimer_diagram("As", pH=7.0)
        self.assertEqual(sorted(p.backend for p in frost.points), sorted(ladder.species))

    def test_a_series_can_be_shared_between_the_two_diagrams(self):
        series = element_series("As", pH=7.0, activity=1.0)
        frost = frost_diagram("As", series=series, pH=7.0)
        ladder = latimer_diagram("As", series=series, pH=7.0)
        self.assertEqual(sorted(p.backend for p in frost.points), sorted(ladder.species))
        # At pH 7 arsenite, not the element, is the bottom of the diagram --
        # seven pH units are worth 1.24 V on a three-electron step.
        self.assertEqual(frost.most_stable.backend, "As(OH)3(aq)")


class TestEveryFormOfEveryState(unittest.TestCase):
    """``predominant_only=False`` stacks several species at one oxidation
    state, which is where the difference between *predominant* and *stable*
    stops being academic.

    Predominant is a comparison within one state -- which of the four
    arsenates, decided by pH. Stable is a comparison across states -- whether
    that arsenate survives at all, decided by the convex hull. An earlier
    version ran the hull over every point, which reported H3AsO4 as stable at
    pH 7 and drew a vertical hull segment down the As(V) column: an acid
    dissociation drawn as a redox step.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.full = frost_diagram("As", pH=7.0, activity=1e-6, predominant_only=False)
        cls.reduced = frost_diagram("As", pH=7.0, activity=1e-6)

    def test_one_predominant_form_per_state(self):
        states = [p.oxidation_state for p in self.full.predominant]
        self.assertEqual(len(states), len(set(states)))

    def test_the_predominant_form_is_the_lowest_of_its_state(self):
        for point in self.full.predominant:
            siblings = [
                other
                for other in self.full.points
                if abs(other.oxidation_state - point.oxidation_state) < 1e-9
            ]
            self.assertEqual(point.gibbs, min(s.gibbs for s in siblings))

    def test_it_agrees_with_the_diagram_built_predominant_only(self):
        self.assertEqual(
            [p.backend for p in self.full.predominant],
            [p.backend for p in self.reduced.points],
        )

    def test_the_hull_never_holds_two_forms_of_one_state(self):
        states = [p.oxidation_state for p in self.full.stable]
        self.assertEqual(len(states), len(set(states)))

    def test_a_minority_acid_form_is_not_called_stable(self):
        # H3AsO4 is the As(V) form at pH 0; at pH 7 it is a trace species that
        # sits above HAsO4(2-) at the same oxidation state.
        self.assertIn("H3AsO4(aq)", [p.backend for p in self.full.points])
        self.assertNotIn("H3AsO4(aq)", [p.backend for p in self.full.stable])

    def test_a_minority_acid_form_is_not_called_disproportionating(self):
        """It is not falling apart; it is just not the dominant acid form."""
        reported = {event.species for event in self.full.disproportionation()}
        minority = {p.backend for p in self.full.points} - {
            p.backend for p in self.full.predominant
        }
        self.assertTrue(minority)
        self.assertEqual(reported & minority, set())

    def test_the_verdicts_do_not_depend_on_how_the_diagram_was_built(self):
        for element in ("As", "S", "C", "Mn"):
            with self.subTest(element=element):
                full = frost_diagram(element, pH=7.0, predominant_only=False)
                reduced = frost_diagram(element, pH=7.0)
                self.assertEqual(
                    [p.backend for p in full.stable], [p.backend for p in reduced.stable]
                )
                self.assertEqual(
                    [str(e) for e in full.disproportionation()],
                    [str(e) for e in reduced.disproportionation()],
                )


class TestShowAndLabelOptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.diagram = frost_diagram("As", pH=7.0, activity=1e-6, predominant_only=False)

    def texts(self, **kwargs):
        import matplotlib.pyplot as plt

        figure, _ = plot_frost("As", diagram=self.diagram, **kwargs)
        try:
            return {
                text.get_text()
                for text in figure.axes[0].texts
                if text.get_text() and "most stable" not in text.get_text()
            }
        finally:
            plt.close(figure)

    def test_label_all_names_every_form(self):
        labelled = self.texts(label="all")
        self.assertIn("$H_3AsO_4$", labelled)
        self.assertIn("$HAsO_4^{2-}$", labelled)

    def test_label_predominant_names_one_per_state(self):
        labelled = self.texts(label="predominant")
        self.assertIn("$HAsO_4^{2-}$", labelled)
        self.assertNotIn("$H_3AsO_4$", labelled)

    def test_predominant_is_the_default(self):
        self.assertEqual(self.texts(), self.texts(label="predominant"))

    def test_prominent_is_accepted_as_a_spelling(self):
        self.assertEqual(self.texts(label="prominent"), self.texts(label="predominant"))

    def test_show_predominant_draws_fewer_markers(self):
        import matplotlib.pyplot as plt

        everything, _ = plot_frost("As", diagram=self.diagram, show="all")
        fewer, _ = plot_frost("As", diagram=self.diagram, show="predominant")
        try:
            self.assertGreater(len(everything.axes[0].lines), len(fewer.axes[0].lines))
        finally:
            plt.close("all")

    def test_show_predominant_matches_a_diagram_built_that_way(self):
        """The cheap path: narrow at draw time rather than recomputing."""
        import matplotlib.pyplot as plt

        narrowed, _ = plot_frost("As", diagram=self.diagram, show="predominant")
        built, _ = plot_frost("As", pH=7.0, activity=1e-6)
        try:
            self.assertEqual(
                {t.get_text() for t in narrowed.axes[0].texts},
                {t.get_text() for t in built.axes[0].texts},
            )
        finally:
            plt.close("all")

    def test_the_key_names_only_the_classes_present(self):
        """Arsenic has minority forms and nothing disproportionating.

        An earlier version always drew a red "disproportionates" swatch, so
        the arsenic figure carried a key to a marker that was not on it.
        """
        import matplotlib.pyplot as plt

        figure, _ = plot_frost("As", diagram=self.diagram, show="all")
        try:
            entries = [text.get_text() for text in figure.axes[0].get_legend().get_texts()]
            self.assertIn("other forms of the same state", entries)
            self.assertNotIn("disproportionates", entries)
        finally:
            plt.close(figure)

    def test_the_key_names_disproportionation_where_it_happens(self):
        import matplotlib.pyplot as plt

        crowded = frost_diagram("S", pH=7.0, activity=1e-5, predominant_only=False)
        figure, _ = plot_frost("S", diagram=crowded, show="all")
        try:
            entries = [text.get_text() for text in figure.axes[0].get_legend().get_texts()]
            self.assertIn("disproportionates", entries)
        finally:
            plt.close(figure)

    def test_an_unknown_selection_is_refused(self):
        import matplotlib.pyplot as plt

        with self.assertRaises(ValueError):
            plot_frost("As", diagram=self.diagram, show="dominant-ish")
        plt.close("all")

    def test_a_key_appears_only_when_there_is_something_to_explain(self):
        import matplotlib.pyplot as plt

        crowded, _ = plot_frost("As", diagram=self.diagram, show="all")
        plain, _ = plot_frost("As", pH=7.0, activity=1e-6)
        try:
            self.assertIsNotNone(crowded.axes[0].get_legend())
            self.assertIsNone(plain.axes[0].get_legend())
        finally:
            plt.close("all")


class TestFigure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_it_renders_and_saves(self):
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "n"
            figure, _ = plot_frost("N", pH=7.0, save=base, annotate_slopes=True)
            try:
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
            finally:
                plt.close(figure)

    def test_a_bad_reference_earns_a_warning_on_the_figure(self):
        import matplotlib.pyplot as plt

        figure, _ = plot_frost("Cr", pH=0.0)
        try:
            texts = " ".join(t.get_text() for t in figure.texts)
            self.assertIn("zero of the", texts)
        finally:
            plt.close(figure)

    def test_the_title_states_the_conditions_and_the_zero(self):
        import matplotlib.pyplot as plt

        figure, _ = plot_frost("As", pH=7.0)
        try:
            title = figure.axes[0].get_title()
            self.assertIn("pH 7", title)
            self.assertIn("zero at", title)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
