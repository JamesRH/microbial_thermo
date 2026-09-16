"""The interactive redox tower.

The tower is not a fixed table, and the point of making it interactive is to
show that couples move at different rates and can change places. What is
tested here is that they move by the *right* amount.

Two checks carry most of the weight. The pH slope of every couple is turned
back into a proton count, which must come out an exact integer matching the
half reaction -- that catches a wrong sign or a dropped RT far more sharply
than comparing potentials to stored numbers. And the activity-ratio shift is
checked against the textbook 59.16/n mV per decade.

Building a grid costs a round of backend calls per temperature, so the grids
here are small and built once per class.
"""

import json
import re
import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

from microbial_thermo.figures.interactive_tower import (
    _draw_static,
    _grid_payload,
    _slider_script,
    _subtitle,
    plot_interactive_tower,
)
from microbial_thermo.tower import (
    DEFAULT_GRID_PH,
    DEFAULT_GRID_TEMPERATURES,
    TowerGrid,
    tower_grid,
)

#: 2.303 RT/F at 25 C, in millivolts. The number every textbook prints.
NERNST_MV_25C = 59.1596

#: Protons consumed by each reference couple's reduction half reaction, as
#: written at the electron count the library normalises it to. Written out
#: rather than derived, so the test is an independent statement.
EXPECTED_PROTONS = {
    "2H+/H2": 2,
    "O2/H2O": 2,
    "S0/HS-": 1,
    "SO4-2/SO3-2": 2,
    "SO4-2/HS-": 9,
    "CO2/CH4": 8,
    "CO2/acetate": 7,
    "NO3-/NO2-": 2,
    "NO2-/NH4+": 8,
    "NO3-/NH4+": 10,
    "NO3-/N2": 12,
    "Fe3+/Fe2+": 0,
}


def small_grid():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tower_grid(ph_values=(7.0, 8.0, 9.0), temperature_values=(25.0, 55.0))


class TestGridShape(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_it_is_indexed_ph_temperature_couple(self):
        self.assertEqual(
            self.grid.potentials.shape,
            (3, 2, len(self.grid.labels)),
        )

    def test_every_reference_couple_is_present(self):
        self.assertEqual(set(self.grid.labels), set(EXPECTED_PROTONS))

    def test_nothing_is_missing_at_these_conditions(self):
        self.assertFalse(np.isnan(self.grid.potentials).any())

    def test_labels_groups_and_electrons_line_up(self):
        self.assertEqual(len(self.grid.labels), len(self.grid.groups))
        self.assertEqual(len(self.grid.labels), len(self.grid.n_electrons))

    def test_nearest_index_snaps_to_the_grid(self):
        self.assertEqual(self.grid.nearest_index("pH", 7.9), 1)
        self.assertEqual(self.grid.nearest_index("T", 60.0), 1)

    def test_the_defaults_are_a_usable_grid(self):
        """Cheap on pH, sparse on temperature -- that asymmetry is deliberate
        and reversing it would make the grid far slower to build."""
        self.assertGreater(len(DEFAULT_GRID_PH), len(DEFAULT_GRID_TEMPERATURES))
        self.assertEqual(min(DEFAULT_GRID_PH), 3.0)
        self.assertEqual(max(DEFAULT_GRID_PH), 11.0)
        self.assertLessEqual(max(DEFAULT_GRID_TEMPERATURES), 100.0)
        self.assertGreater(min(DEFAULT_GRID_TEMPERATURES), 0.0)


class TestProtonSlopes(unittest.TestCase):
    """The pH slope of a couple is -59.16 mV per proton per electron.

    Reading the proton count back out of the computed potentials is a strong
    check: an integer is an unlikely thing to get by accident.
    """

    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_each_couple_implies_its_own_proton_count(self):
        at_seven = self.grid.at(7.0, 25.0)
        at_eight = self.grid.at(8.0, 25.0)
        for label, n, low, high in zip(
            self.grid.labels, self.grid.n_electrons, at_seven, at_eight, strict=True
        ):
            with self.subTest(couple=label):
                slope_mv = 1000.0 * (high - low)
                protons = slope_mv / -NERNST_MV_25C * n
                self.assertAlmostEqual(
                    protons,
                    EXPECTED_PROTONS[label],
                    places=2,
                    msg=f"{label} implies {protons:.3f} protons",
                )

    def test_a_couple_without_protons_does_not_move_with_ph(self):
        index = self.grid.labels.index("Fe3+/Fe2+")
        self.assertAlmostEqual(
            self.grid.at(7.0, 25.0)[index],
            self.grid.at(9.0, 25.0)[index],
            places=9,
        )

    def test_most_couples_do_move_with_ph(self):
        """Guard against a grid that silently ignored pH altogether."""
        moved = [
            label
            for label, low, high in zip(
                self.grid.labels,
                self.grid.at(7.0, 25.0),
                self.grid.at(9.0, 25.0),
                strict=True,
            )
            if abs(high - low) > 1e-4
        ]
        self.assertEqual(len(moved), len(self.grid.labels) - 1)  # all but Fe


class TestActivityRatio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_a_ratio_of_one_changes_nothing(self):
        np.testing.assert_array_equal(self.grid.at(7.0, 25.0), self.grid.at(7.0, 25.0, ratio=1.0))

    def test_the_shift_is_the_textbook_value_per_decade(self):
        base = self.grid.at(7.0, 25.0)
        shifted = self.grid.at(7.0, 25.0, ratio=10.0)
        for label, n, low, high in zip(
            self.grid.labels, self.grid.n_electrons, base, shifted, strict=True
        ):
            with self.subTest(couple=label):
                self.assertAlmostEqual(1000.0 * (high - low), NERNST_MV_25C / n, places=2)

    def test_a_one_electron_couple_swings_furthest(self):
        """Which is the reason the slider teaches anything: couples do not
        move together, so the order can change."""
        base = self.grid.at(7.0, 25.0)
        shifted = self.grid.at(7.0, 25.0, ratio=1e3)
        moves = {
            label: abs(b - s) for label, b, s in zip(self.grid.labels, base, shifted, strict=True)
        }
        self.assertEqual(max(moves, key=moves.get), "Fe3+/Fe2+")

    def test_the_shift_grows_with_temperature(self):
        """RT/nF, so a hotter tower spreads further for the same ratio."""
        index = self.grid.labels.index("Fe3+/Fe2+")
        cool = self.grid.at(7.0, 25.0, 1e3)[index] - self.grid.at(7.0, 25.0)[index]
        warm = self.grid.at(7.0, 55.0, 1e3)[index] - self.grid.at(7.0, 55.0)[index]
        self.assertGreater(warm, cool)

    def test_a_ratio_below_one_lowers_the_potential(self):
        base = self.grid.at(7.0, 25.0)
        lowered = self.grid.at(7.0, 25.0, ratio=1e-3)
        self.assertTrue((lowered < base).all())

    def test_the_order_really_can_change(self):
        """The whole claim the figure makes. At pH 7 and 25 C oxygen sits
        above iron; enrich the ferric:ferrous ratio and iron overtakes it,
        because a one-electron couple moves 59 mV per decade against the
        two-electron oxygen couple's 30."""
        at_parity = [r[0] for r in self.grid.entries_at(7.0, 25.0, 1.0)]
        enriched = [r[0] for r in self.grid.entries_at(7.0, 25.0, 1e3)]
        self.assertNotEqual(at_parity, enriched)
        self.assertEqual(at_parity[-1], "O2/H2O")
        self.assertEqual(enriched[-1], "Fe3+/Fe2+")

    def test_the_crossover_is_between_ten_and_a_hundred(self):
        """Pinning it down, so a change in the iron data or the shift is
        caught rather than merely moving the crossing somewhere else."""
        iron = self.grid.labels.index("Fe3+/Fe2+")
        oxygen = self.grid.labels.index("O2/H2O")
        below = self.grid.at(7.0, 25.0, 1e1)
        above = self.grid.at(7.0, 25.0, 1e2)
        self.assertLess(below[iron], below[oxygen])
        self.assertGreater(above[iron], above[oxygen])


class TestEntriesAndFrame(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_entries_come_back_ordered(self):
        volts = [r[3] for r in self.grid.entries_at(7.0, 25.0)]
        self.assertEqual(volts, sorted(volts))

    def test_unavailable_couples_are_dropped_not_zeroed(self):
        potentials = self.grid.potentials.copy()
        potentials[:, :, 0] = np.nan
        blanked = TowerGrid(
            ph_values=self.grid.ph_values,
            temperature_values=self.grid.temperature_values,
            labels=self.grid.labels,
            groups=self.grid.groups,
            n_electrons=self.grid.n_electrons,
            potentials=potentials,
        )
        rows = blanked.entries_at(7.0, 25.0)
        self.assertEqual(len(rows), len(self.grid.labels) - 1)
        self.assertNotIn(self.grid.labels[0], [r[0] for r in rows])

    def test_the_frame_reports_the_primed_potential(self):
        frame = self.grid.to_frame(7.0, 25.0)
        self.assertIn("E0' (V)", frame.columns)
        self.assertEqual(len(frame), len(self.grid.labels))

    def test_the_subtitle_states_all_three_settings(self):
        text = _subtitle(7.5, 40.0, 2.0)
        self.assertIn("pH 7.5", text)
        self.assertIn("40 °C", text)
        self.assertIn("10^2", text)

    def test_a_parity_ratio_is_written_as_one_to_one(self):
        self.assertIn("1:1", _subtitle(7.0, 25.0, 0.0))


class TestMatplotlibFrame(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_it_draws_one_rule_per_couple(self):
        import matplotlib.pyplot as plt

        figure = _draw_static(self.grid, 7.0, 25.0, 0.0)
        try:
            ax = figure.axes[0]
            # One LineCollection per hlines call, plus the zero line.
            self.assertEqual(len(ax.collections), len(self.grid.labels))
        finally:
            plt.close(figure)

    def test_the_title_states_the_conditions(self):
        import matplotlib.pyplot as plt

        figure = _draw_static(self.grid, 9.0, 55.0, 3.0)
        try:
            title = figure.axes[0].get_title()
            self.assertIn("pH 9", title)
            self.assertIn("55 °C", title)
        finally:
            plt.close(figure)


class TestHtmlExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = small_grid()

    def test_the_payload_is_valid_json(self):
        """NaN is not valid JSON. json.dumps will happily emit a bare NaN
        literal that JSON.parse then refuses, which would break the exported
        page silently -- so unavailable couples must serialise as null."""
        potentials = self.grid.potentials.copy()
        potentials[0, 0, 0] = np.nan
        blanked = TowerGrid(
            ph_values=self.grid.ph_values,
            temperature_values=self.grid.temperature_values,
            labels=self.grid.labels,
            groups=self.grid.groups,
            n_electrons=self.grid.n_electrons,
            potentials=potentials,
        )
        text = json.dumps(_grid_payload(blanked))
        self.assertNotIn("NaN", text)
        self.assertIsNone(json.loads(text)["potentials"][0][0][0])

    def test_the_payload_carries_the_computed_potentials(self):
        payload = _grid_payload(self.grid)
        np.testing.assert_allclose(
            np.array(payload["potentials"], dtype=float),
            np.asarray(self.grid.potentials),
        )

    def test_the_payload_carries_the_electron_counts(self):
        """The page redoes the Nernst shift itself, so it needs them."""
        payload = _grid_payload(self.grid)
        self.assertEqual(tuple(payload["nElectrons"]), self.grid.n_electrons)

    def test_the_script_embeds_the_grid_and_wires_three_sliders(self):
        script = _slider_script(self.grid, (-6.0, 6.0), "Redox tower")
        self.assertIn("var GRID =", script)
        self.assertEqual(script.count("control("), 4)  # one definition, three uses
        self.assertIn("Plotly.restyle", script)

    def test_the_script_uses_the_same_constants_as_the_library(self):
        """If these drift, the exported page and the notebook disagree."""
        from microbial_thermo.units import FARADAY, R

        script = _slider_script(self.grid, (-6.0, 6.0), "t")
        gas_constant = float(re.search(r"var rtOverF[^=]*= ([\d.]+)", script).group(1))
        faraday = float(re.search(r"/ (\d+\.\d+);", script).group(1))
        # R is defined as ``1 * ureg.molar_gas_constant``, so its magnitude is
        # 1 until it is converted.
        self.assertAlmostEqual(gas_constant, R.to("J/(mol*K)").magnitude, places=6)
        self.assertAlmostEqual(faraday, FARADAY.to("C/mol").magnitude, places=3)

    def test_it_writes_a_standalone_page(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "tower"
            plot_interactive_tower(self.grid, save_html=base)
            written = base.with_suffix(".html")
            self.assertTrue(written.exists())
            text = written.read_text()
            # Self-contained: plotly.js inlined, no network fetch.
            self.assertIn("tower-plot", text)
            self.assertNotIn('<script src="https://', text)
            self.assertIn("input", text)

    def test_the_figure_has_one_trace_per_group(self):
        figure = plot_interactive_tower(self.grid)
        self.assertEqual(len(figure.data), len(set(self.grid.groups)))


if __name__ == "__main__":
    unittest.main()
