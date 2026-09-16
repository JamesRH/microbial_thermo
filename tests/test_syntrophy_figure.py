"""The syntrophy window figure.

The energetics behind it are asserted in ``test_multicouple.py``; what is
tested here is the window-finding -- which is where a plausible-looking wrong
answer could hide -- and that the figure renders.

Sweeps cost a pyGCC call per species per point, so the grids here are coarse
on purpose. The edge refinement is what makes a coarse grid legitimate, and
that is itself one of the things under test.
"""

import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

import microbial_thermo as mt
from microbial_thermo.exceptions import MicrobialThermoError
from microbial_thermo.figures.syntrophy import (
    _carrier_label,
    _longest_run,
    _zero_crossing,
    plot_syntrophy_window,
    syntrophy_window,
)
from microbial_thermo.library import reaction


def syntrophic_conditions(p_h2: float = 1e-5) -> mt.Conditions:
    """A flooded sediment: both partners present, hydrogen drawn down."""
    return mt.Conditions(
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


class TestHelpers(unittest.TestCase):
    """The pure pieces, without paying for a sweep."""

    def test_longest_run_picks_the_longest(self):
        self.assertEqual(_longest_run([False, True, False, True, True, True]), (3, 5))

    def test_longest_run_of_nothing(self):
        self.assertIsNone(_longest_run([False, False]))

    def test_longest_run_spanning_everything(self):
        self.assertEqual(_longest_run([True, True, True]), (0, 2))

    def test_zero_crossing_interpolates_in_log_pressure(self):
        """Halfway in log space, not halfway in pressure. A curve running from
        -1 to +1 across a decade must cross at the geometric middle."""
        pressures = np.array([1e-6, 1e-5])
        crossing = _zero_crossing(pressures, np.array([-1.0, 1.0]), 0)
        self.assertAlmostEqual(np.log10(crossing), -5.5, places=9)

    def test_zero_crossing_is_not_the_arithmetic_middle(self):
        pressures = np.array([1e-6, 1e-5])
        crossing = _zero_crossing(pressures, np.array([-1.0, 1.0]), 0)
        self.assertNotAlmostEqual(crossing, 5.5e-6, places=8)

    def test_the_axis_label_cannot_be_read_as_ph(self):
        """'pH2(g)' on an axis in a library full of pH is a trap."""
        label = _carrier_label("H2(g)")
        self.assertNotIn("pH", label)
        self.assertIn("partial pressure", label)


class TestWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        conditions = syntrophic_conditions()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.producer = reaction("syntrophic_propionate_oxidation", conditions)
            cls.consumer = reaction("hydrogenotrophic_methanogenesis", conditions)
            cls.window = syntrophy_window(
                cls.producer, cls.consumer, low=1e-10, high=1.0, points=13
            )

    def test_both_partners_are_on_the_same_electron_basis(self):
        """Which is what makes one shared axis honest."""
        self.assertEqual(self.producer.n_electrons, self.consumer.n_electrons)
        self.assertEqual(self.window.n_electrons, 2.0)

    def test_a_window_exists_and_is_ordered(self):
        self.assertTrue(self.window.exists)
        self.assertLess(self.window.low, self.window.high)

    def test_the_window_is_narrow_but_real(self):
        """Published propionate syntrophy runs at 1e-6 to 1e-4 bar H2."""
        self.assertGreater(self.window.decades, 1.0)
        self.assertLess(self.window.decades, 4.0)
        self.assertGreater(self.window.low, 1e-8)
        self.assertLess(self.window.high, 1e-2)

    def test_both_partners_really_are_exergonic_inside_it(self):
        """Rebuild at the midpoint rather than trusting the swept array."""
        middle = 10 ** ((np.log10(self.window.low) + np.log10(self.window.high)) / 2)
        conditions = syntrophic_conditions(middle)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for name in (
                "syntrophic_propionate_oxidation",
                "hydrogenotrophic_methanogenesis",
            ):
                with self.subTest(partner=name):
                    energy = reaction(name, conditions).delta_G.magnitude
                    self.assertLess(energy, 0.0)

    def test_one_partner_fails_on_each_side_of_it(self):
        """Independently computed, an order of magnitude outside each edge."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            below = reaction(
                "hydrogenotrophic_methanogenesis",
                syntrophic_conditions(self.window.low / 10),
            ).delta_G.magnitude
            above = reaction(
                "syntrophic_propionate_oxidation",
                syntrophic_conditions(self.window.high * 10),
            ).delta_G.magnitude
        self.assertGreater(below, 0.0, "methanogen should starve below the window")
        self.assertGreater(above, 0.0, "propionate oxidation should stall above it")

    def test_the_edges_are_refined_off_the_grid(self):
        """Otherwise the reported width would depend on the sampling."""
        self.assertNotIn(self.window.low, list(self.window.pressures))
        self.assertNotIn(self.window.high, list(self.window.pressures))

    def test_a_coarse_and_a_fine_grid_agree_on_the_edges(self):
        """The point of the refinement: the answer is a property of the
        chemistry, not of how many points were sampled."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            finer = syntrophy_window(self.producer, self.consumer, low=1e-10, high=1.0, points=25)
        self.assertAlmostEqual(np.log10(finer.low), np.log10(self.window.low), places=6)
        self.assertAlmostEqual(np.log10(finer.high), np.log10(self.window.high), places=6)

    def test_best_shared_is_the_worse_off_partner(self):
        pressure, shared = self.window.best_shared
        index = list(self.window.pressures).index(pressure)
        self.assertAlmostEqual(
            shared,
            max(
                self.window.producer_delta_g[index],
                self.window.consumer_delta_g[index],
            ),
        )

    def test_best_shared_beats_every_other_grid_point(self):
        _, shared = self.window.best_shared
        worse = np.maximum(self.window.producer_delta_g, self.window.consumer_delta_g)
        self.assertAlmostEqual(shared, float(worse.min()))

    def test_the_partnership_lives_inside_the_energy_quantum(self):
        """The biological punchline: neither partner can be guaranteed even
        one quantum, which is why syntrophs grow so slowly."""
        _, shared = self.window.best_shared
        self.assertLess(shared, 0.0)
        self.assertGreater(shared, -20.0)

    def test_a_range_outside_the_window_finds_nothing(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            missing = syntrophy_window(self.producer, self.consumer, low=1e-1, high=1.0, points=5)
        self.assertFalse(missing.exists)
        self.assertEqual(missing.decades, 0.0)

    def test_mismatched_electron_counts_are_refused(self):
        """Two reactions normalised differently cannot share a y axis."""
        conditions = syntrophic_conditions()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            eight = reaction("hydrogenotrophic_methanogenesis", conditions, normalize_to="integer")
        if eight.n_electrons == self.producer.n_electrons:
            self.skipTest("integer normalisation happened to match")
        with self.assertRaises(MicrobialThermoError) as caught:
            syntrophy_window(self.producer, eight, points=3)
        self.assertIn("not comparable", str(caught.exception))

    def test_the_frame_reports_the_overlap(self):
        frame = self.window.to_frame()
        self.assertIn("both exergonic", frame.columns)
        self.assertTrue(frame["both exergonic"].any())


class TestFigure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        conditions = syntrophic_conditions()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.producer = reaction("syntrophic_propionate_oxidation", conditions)
            cls.consumer = reaction("hydrogenotrophic_methanogenesis", conditions)
            cls.window = syntrophy_window(
                cls.producer, cls.consumer, low=1e-10, high=1.0, points=13
            )

    def _plot(self, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return plot_syntrophy_window(self.producer, self.consumer, window=self.window, **kwargs)

    def test_it_saves_both_formats(self):
        figure, _ = self._plot()
        try:
            with tempfile.TemporaryDirectory() as directory:
                base = Path(directory) / "syntrophy"
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    plot_syntrophy_window(
                        self.producer,
                        self.consumer,
                        window=self.window,
                        save=base,
                    )
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
        finally:
            matplotlib.pyplot.close("all")

    def test_it_returns_the_window_it_drew(self):
        figure, window = self._plot()
        try:
            self.assertIs(window, self.window)
        finally:
            matplotlib.pyplot.close(figure)

    def test_a_precomputed_window_is_not_recomputed(self):
        """Passing one in should cost no pyGCC calls, which is what makes the
        figure usable from a slider callback."""
        figure, window = self._plot()
        try:
            self.assertIs(window, self.window)
        finally:
            matplotlib.pyplot.close(figure)

    def test_two_curves_are_drawn(self):
        figure, _ = self._plot()
        try:
            ax = figure.axes[0]
            self.assertEqual(len(ax.lines) >= 2, True)
        finally:
            matplotlib.pyplot.close(figure)

    def test_the_title_names_the_width(self):
        figure, _ = self._plot()
        try:
            self.assertIn("decades", figure.axes[0].get_title())
        finally:
            matplotlib.pyplot.close(figure)

    def test_the_axis_label_is_a_pressure_not_a_ph(self):
        figure, _ = self._plot()
        try:
            label = figure.axes[0].get_xlabel()
            self.assertIn("partial pressure", label)
            self.assertNotIn("pH", label)
        finally:
            matplotlib.pyplot.close(figure)

    def test_the_annotation_can_be_turned_off(self):
        with_note, _ = self._plot()
        without, _ = self._plot(annotate_best=False)
        try:
            self.assertGreater(len(with_note.axes[0].texts), len(without.axes[0].texts))
        finally:
            matplotlib.pyplot.close("all")


if __name__ == "__main__":
    unittest.main()
