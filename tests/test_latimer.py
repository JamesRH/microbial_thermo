"""Latimer diagrams: the redox ladder as a chain of couple potentials.

The potentials here are computed from formation energies by a route that
never touches a table of standard potentials, so comparing them against
published ones is a real check on the arithmetic rather than a tautology.
Copper is the sharpest case: Cu(2+)/Cu(+) at +0.153 V and Cu(+)/Cu at
+0.521 V are textbook, they bracket in the awkward order, and the
consequence -- that Cu(+) disproportionates -- has to fall out on its own.

The other thing being tested is the rule that keeps pyrite off the iron
ladder. Its electron count against a sulfate basis is -12, which is a real
number about a real reaction and is not iron's oxidation state.
"""

import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from microbial_thermo.figures.basis import element_series, ladder_entries
from microbial_thermo.figures.latimer import (
    _coefficient,
    _state_label,
    latimer_diagram,
    plot_latimer,
)


class TestPotentialsAgainstPublishedValues(unittest.TestCase):
    """Standard-state conditions, where published numbers apply."""

    def step(self, element, oxidized, reduced, **kwargs):
        diagram = latimer_diagram(element, pH=0.0, activity=1.0, **kwargs)
        for step in diagram.steps:
            if step.oxidized.backend == oxidized and step.reduced.backend == reduced:
                return step
        raise AssertionError(
            f"no {oxidized} -> {reduced} step; got "
            + ", ".join(f"{s.oxidized.backend}->{s.reduced.backend}" for s in diagram.steps)
        )

    def test_copper_two_to_copper_one(self):
        # Published +0.153 V.
        self.assertAlmostEqual(self.step("Cu", "Cu++", "Cu+").potential, 0.153, delta=0.02)

    def test_copper_one_to_copper_metal(self):
        # Published +0.521 V.
        self.assertAlmostEqual(self.step("Cu", "Cu+", "Cu").potential, 0.521, delta=0.02)

    def test_sulfur_to_sulfide(self):
        # Published +0.144 V for S + 2H+ + 2e- -> H2S(aq).
        self.assertAlmostEqual(self.step("S", "Sulfur(s)", "H2S(aq)").potential, 0.144, delta=0.01)

    def test_arsenate_to_arsenite(self):
        # Published +0.560 V, quoted against HAsO2; As(OH)3 is the same species
        # plus a water, so the potential is the same.
        self.assertAlmostEqual(
            self.step("As", "H3AsO4(aq)", "As(OH)3(aq)").potential, 0.560, delta=0.02
        )

    def test_arsenite_to_arsenic(self):
        # Published +0.240 V.
        self.assertAlmostEqual(self.step("As", "As(OH)3(aq)", "As").potential, 0.240, delta=0.02)

    def test_manganese_two_to_the_metal(self):
        # Published -1.185 V. Elemental Mn carries a 2.5 kJ/mol database
        # offset, which is 13 mV over two electrons -- inside this tolerance
        # but worth knowing it is there.
        self.assertAlmostEqual(self.step("Mn", "Mn++", "Mn").potential, -1.185, delta=0.03)

    def test_nitrate_to_nitrite(self):
        # Published +0.835 V.
        self.assertAlmostEqual(self.step("N", "NO3-", "NO2-").potential, 0.835, delta=0.02)


class TestElectronCounts(unittest.TestCase):
    def test_a_step_carries_the_difference_in_oxidation_state(self):
        diagram = latimer_diagram("As", pH=0.0)
        for step in diagram.steps:
            self.assertAlmostEqual(
                step.n_electrons,
                step.oxidized.oxidation_state - step.reduced.oxidation_state,
            )

    def test_a_mixed_valence_oxide_gives_a_fractional_step(self):
        diagram = latimer_diagram("Fe", pH=0.0)
        magnetite = [s for s in diagram.steps if s.reduced.backend == "Magnetite"]
        self.assertTrue(magnetite)
        # Fe(III) to Fe(8/3) is a third of an electron per iron.
        self.assertAlmostEqual(magnetite[0].n_electrons, 1 / 3, places=6)

    def test_the_energy_matches_the_potential(self):
        diagram = latimer_diagram("S", pH=0.0)
        for step in diagram.steps:
            self.assertAlmostEqual(
                step.delta_g, -step.n_electrons * step.potential * 96.485, delta=0.02
            )


class TestWhichSpeciesSitsOnEachRung(unittest.TestCase):
    """The rung is the form that predominates at the working pH."""

    def test_arsenate_changes_form_with_pH(self):
        chosen = {ph: latimer_diagram("As", pH=ph).rungs[0].backend for ph in (0.0, 4.0, 7.0, 12.0)}
        self.assertEqual(chosen[0.0], "H3AsO4(aq)")
        self.assertEqual(chosen[4.0], "H2AsO4-")
        self.assertEqual(chosen[7.0], "HAsO4--")
        self.assertEqual(chosen[12.0], "AsO4---")

    def test_the_forms_not_chosen_are_kept(self):
        diagram = latimer_diagram("As", pH=7.0)
        rejected = {name for _, name in diagram.alternatives}
        self.assertIn("H3AsO4(aq)", rejected)
        self.assertNotIn("HAsO4--", rejected)

    def test_one_rung_per_oxidation_state(self):
        diagram = latimer_diagram("As", pH=7.0)
        self.assertEqual(len(diagram.states), len(set(diagram.states)))

    def test_potentials_move_with_pH(self):
        acid = latimer_diagram("As", pH=0.0).steps[-1].potential
        neutral = latimer_diagram("As", pH=7.0).steps[-1].potential
        # As(OH)3 + 3H+ + 3e- -> As consumes three protons per three
        # electrons, so the shift is a full -59 mV per pH unit.
        self.assertAlmostEqual(neutral - acid, -7 * 0.0592, delta=0.004)


class TestPyriteIsNotAnIronOxidationState(unittest.TestCase):
    def test_it_is_excluded_with_a_reason(self):
        diagram = latimer_diagram("Fe", pH=0.0)
        excluded = dict(diagram.excluded)
        self.assertIn("Pyrite", excluded)
        self.assertIn("oxidation state", excluded["Pyrite"])
        self.assertNotIn("Pyrite", diagram.species)

    def test_siderite_goes_too_although_its_carbon_is_fine(self):
        # Carbonate carbon is C(+IV) in both siderite and the bicarbonate
        # basis, so siderite's +2 would have been right. The rule refuses it
        # anyway, because telling it apart from pyrite needs an independent
        # oxidation-state assignment and the one available refuses sulfides.
        # It loses nothing: siderite is Fe(II) and Fe(2+) already holds that
        # rung. Asserted so the exclusion stays deliberate.
        series = element_series("Fe", activity=1.0)
        kept, dropped = ladder_entries(series)
        self.assertNotIn("Siderite", [e.backend for e in kept])
        self.assertIn("Siderite", dict(dropped))
        siderite = [e for e in series if e.backend == "Siderite"][0]
        self.assertAlmostEqual(siderite.oxidation_state, 2.0, places=6)

    def test_no_ladder_at_all_is_an_error_not_a_silent_empty_figure(self):
        with self.assertRaises(ValueError):
            latimer_diagram("Fe", species=["Pyrite"], pH=0.0)


class TestDisproportionation(unittest.TestCase):
    def test_copper_one_is_flagged(self):
        diagram = latimer_diagram("Cu", pH=0.0)
        self.assertIn("Cu+", [rung.backend for rung in diagram.disproportionating])

    def test_the_flag_is_the_two_arrow_comparison(self):
        diagram = latimer_diagram("Cu", pH=0.0)
        into, out = diagram.steps[0], diagram.steps[1]
        self.assertGreater(out.potential, into.potential)


class TestSkipStepPotentials(unittest.TestCase):
    def test_a_skip_is_the_electron_weighted_mean(self):
        diagram = latimer_diagram("As", pH=0.0)
        first, second = diagram.steps
        weighted = (first.n_electrons * first.potential + second.n_electrons * second.potential) / (
            first.n_electrons + second.n_electrons
        )
        self.assertAlmostEqual(diagram.potential("H3AsO4(aq)", "As"), weighted, places=9)

    def test_it_is_not_the_arithmetic_mean(self):
        diagram = latimer_diagram("As", pH=0.0)
        naive = sum(s.potential for s in diagram.steps) / len(diagram.steps)
        self.assertNotAlmostEqual(diagram.potential("H3AsO4(aq)", "As"), naive, places=3)

    def test_an_unknown_species_says_what_is_on_the_diagram(self):
        diagram = latimer_diagram("As", pH=0.0)
        with self.assertRaises(KeyError) as caught:
            diagram.potential("Hematite", "As")
        self.assertIn("As", str(caught.exception))


class TestHalfReactions(unittest.TestCase):
    def test_arsenate_reduction_is_balanced_as_written(self):
        diagram = latimer_diagram("As", pH=0.0)
        text = diagram.steps[0].half_reaction()
        self.assertEqual(text, "H3AsO4(aq) + 2 H+ + 2 e- -> As(OH)3(aq) + H2O")

    def test_fractional_coefficients_print_as_fractions(self):
        self.assertEqual(_coefficient(1 / 3), "1/3")
        self.assertEqual(_coefficient(8 / 3), "8/3")
        self.assertEqual(_coefficient(2.0), "2")
        self.assertEqual(_coefficient(1.0), "")

    def test_state_labels_are_roman(self):
        self.assertEqual(_state_label(5), "V")
        self.assertEqual(_state_label(-3), "-III")
        self.assertEqual(_state_label(0), "0")
        self.assertEqual(_state_label(8 / 3), "8/3")


class TestTemperature(unittest.TestCase):
    def test_the_ladder_recomputes_with_temperature(self):
        cold = latimer_diagram("As", pH=7.0, temperature_c=2.0).steps[0].potential
        warm = latimer_diagram("As", pH=7.0, temperature_c=60.0).steps[0].potential
        self.assertNotAlmostEqual(cold, warm, places=3)


class TestFigure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_it_renders_and_saves(self):
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "as"
            figure, diagram = plot_latimer("As", pH=7.0, save=base)
            try:
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
                self.assertIn("pH 7", figure.axes[0].get_title())
            finally:
                plt.close(figure)

    def test_it_takes_a_prebuilt_diagram(self):
        import matplotlib.pyplot as plt

        diagram = latimer_diagram("Cu", pH=0.0)
        figure, same = plot_latimer("Cu", diagram=diagram)
        try:
            self.assertIs(same, diagram)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
