"""The three element diagrams under sliders, and the exported page.

Two things are worth testing here and the rest is plumbing.

The first is that **only temperature costs anything**. pH and activity are
applied in closed form on lookup, so a series taken from the grid at a new pH
or activity must equal one built from scratch at those conditions -- not
approximately, exactly, because it is the same arithmetic in a different
order. If that ever stops holding, the sliders are lying.

The second is that the **exported page agrees with the library**. The HTML
file carries the decomposition and rebuilds the diagrams in JavaScript, which
means there are two implementations of the same equations and they can drift.
The check runs the real exported script in node and compares; it skips where
node is not installed rather than pretending to pass.
"""

import json
import shutil
import subprocess
import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from microbial_thermo.figures.basis import element_grid, element_series
from microbial_thermo.figures.frost import frost_diagram
from microbial_thermo.figures.interactive_element import (
    ELEMENT_DIV_ID,
    PANELS,
    _figure,
    _slider_script,
    plot_interactive_element,
)
from microbial_thermo.figures.latimer import latimer_diagram
from microbial_thermo.figures.pourbaix import pourbaix_field

NODE = shutil.which("node") or shutil.which("nodejs")

#: A browser stub thin enough to be obviously inert: the script only ever
#: creates elements, appends them and calls Plotly, none of which affect the
#: numbers being checked.
BROWSER_STUB = """
function El() { this.style = {}; this.children = []; }
El.prototype.appendChild = function (c) { this.children.push(c); };
El.prototype.addEventListener = function () {};
Object.defineProperty(El.prototype, 'textContent',
  {set: function (v) { this._t = v; }, get: function () { return this._t; }});
Object.defineProperty(El.prototype, 'innerHTML',
  {set: function (v) { this._h = v; }, get: function () { return this._h; }});
var root = new El();
root.parentNode = { insertBefore: function () {} };
root.nextSibling = null;
root.data = [];
global.document = {
  createElement: function () { return new El(); },
  getElementById: function () { return root; }
};
global.Plotly = { restyle: function () {}, relayout: function () {} };
"""

PROBE = """
var out = {};
var g = energies(TI, PH, LOGA);
var volts = g.map(function (v) { return v / GRID.faraday; });
var picked = rungs(g);
var stable = hull(picked, volts);
out.species = picked.map(function (i) { return GRID.species[i]; });
out.volts = picked.map(function (i) { return volts[i]; });
out.stable = stable.map(function (i) { return GRID.species[i]; });
out.steps = [];
for (var k = picked.length - 1; k > 0; k--) {
  var high = picked[k], low = picked[k - 1];
  var n = GRID.states[high] - GRID.states[low];
  out.steps.push([GRID.species[high], GRID.species[low],
                  (g[high] - g[low]) / (n * GRID.faraday)]);
}
var cell = field(TI, LOGA);
var m = Math.floor(cell.text.length / 2);
out.mid = cell.text[m][m];
out.midPh = cell.ph[m];
out.midEh = cell.eh[m];
out.corners = [cell.text[0][0], cell.text[cell.text.length - 1][cell.text.length - 1]];
console.log(JSON.stringify(out));
"""


def run_in_node(grid, temperature_index, pH, log_activity, points=60):
    script = _slider_script(grid, (-1.0, 1.4), points, (-9.0, 0.0), "t", 7.0, 25.0, -6.0)
    probe = (
        PROBE.replace("TI", str(temperature_index))
        .replace("PH", repr(pH))
        .replace("LOGA", repr(log_activity))
    )
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "check.js"
        path.write_text(BROWSER_STUB + script + probe)
        result = subprocess.run(
            [NODE, str(path)], capture_output=True, text=True, check=True, timeout=120
        )
    return json.loads(result.stdout)


class TestOnlyTemperatureCosts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("As", temperatures=(2.0, 25.0, 60.0))

    def test_a_rescaled_series_equals_one_built_from_scratch(self):
        taken = self.grid.at(25.0, pH=8.3, activity=1e-7)
        built = element_series("As", temperature_c=25.0, activity=1e-7, pH=8.3)
        self.assertEqual([e.backend for e in taken], [e.backend for e in built])
        for a, b in zip(taken, built, strict=True):
            with self.subTest(species=a.backend):
                self.assertAlmostEqual(a.base, b.base, places=9)
                self.assertAlmostEqual(a.protons, b.protons, places=12)
                self.assertAlmostEqual(a.electrons, b.electrons, places=12)

    def test_the_diagrams_agree_with_freshly_built_ones(self):
        taken = frost_diagram("As", series=self.grid.at(60.0, pH=5.0, activity=1e-3))
        built = frost_diagram("As", temperature_c=60.0, pH=5.0, activity=1e-3)
        self.assertEqual([p.backend for p in taken.points], [p.backend for p in built.points])
        for a, b in zip(taken.points, built.points, strict=True):
            self.assertAlmostEqual(a.volt_equivalent, b.volt_equivalent, places=9)

    def test_the_temperature_axis_snaps_to_the_grid(self):
        self.assertEqual(self.grid.at(30.0).temperature_c, 25.0)
        self.assertEqual(self.grid.at(70.0).temperature_c, 60.0)

    def test_a_prepared_series_sets_the_conditions_it_was_built_at(self):
        """A figure must not be able to state conditions it was not drawn at."""
        series = self.grid.at(60.0, pH=5.0, activity=1e-3)
        ladder = latimer_diagram("As", series=series, temperature_c=25.0, activity=1.0)
        self.assertEqual(ladder.temperature_c, 60.0)
        self.assertEqual(ladder.activity, 1e-3)
        self.assertEqual(ladder.pH, 5.0)


class TestSpeciesAvailableAtEveryTemperature(unittest.TestCase):
    def test_a_single_temperature_phase_is_dropped_from_the_whole_grid(self):
        """Manganite carries one log K, at 25 C.

        It belongs on the static 25 C diagram and cannot be on an interactive
        one, because a species that appears and vanishes as the slider moves
        is a different diagram each time rather than the same one at new
        conditions.
        """
        grid = element_grid("Mn", temperatures=(10.0, 25.0, 40.0))
        self.assertNotIn("Manganite", grid.species)
        self.assertIn("Manganite", dict(grid.dropped))
        self.assertIn("25", dict(grid.dropped)["Manganite"])

    def test_the_static_diagram_still_has_it(self):
        series = element_series("Mn", temperature_c=25.0)
        self.assertIn("Manganite", [entry.backend for entry in series])

    def test_every_temperature_gives_the_same_species_list(self):
        grid = element_grid("Mn", temperatures=(10.0, 25.0, 40.0))
        lists = [[entry.backend for entry in series] for series in grid.series]
        for other in lists[1:]:
            self.assertEqual(other, lists[0])

    def test_an_out_of_range_species_is_an_error_by_default(self):
        from microbial_thermo.exceptions import OutOfRangeError

        with self.assertRaises(OutOfRangeError):
            element_series("Mn", species=["Manganite", "Mn++"], temperature_c=2.0)


@unittest.skipIf(NODE is None, "node is not installed")
class TestTheExportedPageAgrees(unittest.TestCase):
    """Two implementations of the same equations, checked against each other."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("Mn", temperatures=(10.0, 25.0, 40.0))
        cls.index, cls.temperature = 2, 40.0
        cls.ph, cls.log_activity = 8.5, -6.0
        cls.js = run_in_node(cls.grid, cls.index, cls.ph, cls.log_activity)
        cls.series = cls.grid.at(cls.temperature, pH=cls.ph, activity=10.0**cls.log_activity)

    def test_the_same_species_are_on_the_ladder(self):
        frost = frost_diagram("Mn", series=self.series)
        self.assertEqual([p.backend for p in frost.points], self.js["species"])

    def test_the_volt_equivalents_agree(self):
        frost = frost_diagram("Mn", series=self.series)
        for point, value in zip(frost.points, self.js["volts"], strict=True):
            with self.subTest(species=point.backend):
                self.assertAlmostEqual(point.volt_equivalent, value, places=9)

    def test_the_hull_agrees(self):
        frost = frost_diagram("Mn", series=self.series)
        self.assertEqual([p.backend for p in frost.stable], self.js["stable"])

    def test_the_latimer_potentials_agree(self):
        ladder = latimer_diagram("Mn", series=self.series)
        self.assertEqual(len(ladder.steps), len(self.js["steps"]))
        for step, (high, low, value) in zip(ladder.steps, self.js["steps"], strict=True):
            with self.subTest(step=f"{high}->{low}"):
                self.assertEqual((step.oxidized.backend, step.reduced.backend), (high, low))
                self.assertAlmostEqual(step.potential, value, places=9)

    def test_the_eh_ph_field_agrees(self):
        field = pourbaix_field("Mn", series=self.series, points=60, eh_range=(-1.0, 1.4))
        self.assertEqual(field.at(self.js["midPh"], self.js["midEh"]), self.js["mid"])


class TestPanelSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("As", temperatures=(25.0,))

    def draw(self, **kwargs):
        return _figure(self.grid, 25.0, 7.0, 1e-6, (11.0, 7.0), (-1.0, 1.4), 40, **kwargs)

    def test_all_three_by_default(self):
        import matplotlib.pyplot as plt

        figure = self.draw()
        try:
            self.assertEqual(len(figure.axes), 3)
        finally:
            plt.close(figure)

    def test_the_cheap_pair_can_be_asked_for_alone(self):
        """Latimer and Frost need no expensive axis, so they are the light widget."""
        import matplotlib.pyplot as plt

        figure = self.draw(panels=("latimer", "frost"))
        try:
            self.assertEqual(len(figure.axes), 2)
            titles = " ".join(ax.get_title() for ax in figure.axes)
            self.assertIn("Latimer", titles)
            self.assertIn("Frost", titles)
            self.assertNotIn("predominance", titles)
        finally:
            plt.close(figure)

    def test_one_panel_alone(self):
        import matplotlib.pyplot as plt

        for panel in PANELS:
            with self.subTest(panel=panel):
                figure = self.draw(panels=(panel,))
                try:
                    self.assertEqual(len(figure.axes), 1)
                finally:
                    plt.close(figure)

    def test_the_order_is_fixed_whatever_order_is_asked_for(self):
        import matplotlib.pyplot as plt

        figure = self.draw(panels=("pourbaix", "latimer"))
        try:
            self.assertIn("Latimer", figure.axes[0].get_title())
        finally:
            plt.close(figure)

    def test_an_unknown_panel_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.draw(panels=("frost", "tower"))
        self.assertIn("tower", str(caught.exception))

    def test_no_panels_is_refused(self):
        with self.assertRaises(ValueError):
            self.draw(panels=())


class TestTheFrostControls(unittest.TestCase):
    """`show` and `label` reach the Frost panel of the interactive figure."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("As", temperatures=(25.0,))

    def names(self, **kwargs):
        import matplotlib.pyplot as plt

        figure = _figure(
            self.grid, 25.0, 7.0, 1e-6, (11.0, 5.0), (-1.0, 1.4), 40, panels=("frost",), **kwargs
        )
        try:
            return {
                text.get_text()
                for text in figure.axes[0].texts
                if text.get_text() and "most stable" not in text.get_text()
            }
        finally:
            plt.close(figure)

    def test_the_panel_is_built_with_every_form(self):
        """Otherwise there would be nothing for the control to show."""
        self.assertIn("$H_3AsO_4$", self.names(label="all"))

    def test_label_predominant_is_the_default(self):
        self.assertEqual(self.names(), self.names(label="predominant"))
        self.assertNotIn("$H_3AsO_4$", self.names())

    def test_show_predominant_drops_the_minority_markers(self):
        import matplotlib.pyplot as plt

        everything = _figure(
            self.grid, 25.0, 7.0, 1e-6, (11.0, 5.0), (-1.0, 1.4), 40, panels=("frost",)
        )
        fewer = _figure(
            self.grid,
            25.0,
            7.0,
            1e-6,
            (11.0, 5.0),
            (-1.0, 1.4),
            40,
            panels=("frost",),
            show="predominant",
        )
        try:
            self.assertGreater(len(everything.axes[0].lines), len(fewer.axes[0].lines))
        finally:
            plt.close("all")


class TestTheExportedControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("As", temperatures=(25.0,))

    def test_the_page_carries_both_dropdowns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "as"
            plot_interactive_element("As", grid=self.grid, save_html=path)
            written = path.with_suffix(".html").read_text()
        self.assertIn("createElement('select')", written)
        self.assertIn("showControl", written)
        self.assertIn("labelControl", written)
        self.assertIn("predominant only", written)

    def test_the_minority_trace_exists_to_be_toggled(self):
        figure = plot_interactive_element("As", grid=self.grid)
        self.assertEqual(figure.data[4].name, "other forms of the same state")
        self.assertTrue(len(figure.data[4].x) > 0)

    @unittest.skipIf(NODE is None, "node is not installed")
    def test_the_page_agrees_on_which_forms_are_minority(self):
        """The same split, computed twice in two languages."""
        javascript = run_in_node(self.grid, 0, 7.0, -6.0)
        series = self.grid.at(25.0, pH=7.0, activity=1e-6)
        frost = frost_diagram("As", series=series, predominant_only=False)
        self.assertEqual([p.backend for p in frost.predominant], javascript["species"])


class TestFigures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.grid = element_grid("As", temperatures=(25.0,))

    def test_the_three_panels_render(self):
        import matplotlib.pyplot as plt

        figure = _figure(self.grid, 25.0, 7.0, 1e-6, (12.0, 9.0), (-1.0, 1.4), 40)
        try:
            self.assertEqual(len(figure.axes), 3)
            titles = " ".join(ax.get_title() for ax in figure.axes)
            self.assertIn("Latimer", titles)
            self.assertIn("Frost", titles)
            self.assertIn("predominance", titles)
        finally:
            plt.close(figure)

    def test_the_ph_slider_is_drawn_on_the_field(self):
        import matplotlib.pyplot as plt

        figure = _figure(self.grid, 25.0, 9.0, 1e-6, (12.0, 9.0), (-1.0, 1.4), 40)
        try:
            field_ax = figure.axes[2]
            verticals = [
                line
                for line in field_ax.lines
                if len(set(line.get_xdata())) == 1 and set(line.get_xdata()) == {9.0}
            ]
            self.assertTrue(verticals, "the pH line should be on the Eh-pH panel")
        finally:
            plt.close(figure)

    def test_the_html_export_is_standalone_and_carries_the_sliders(self):
        import plotly.graph_objects as go

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "as"
            figure = plot_interactive_element("As", grid=self.grid, save_html=path)
            self.assertIsInstance(figure, go.Figure)
            written = path.with_suffix(".html").read_text()
        self.assertIn(ELEMENT_DIV_ID, written)
        self.assertIn("plotly.js", written.lower())
        self.assertIn("input.type = 'range'", written)
        # The decomposition travels, not the drawn diagram.
        self.assertIn('"protons"', written)
        self.assertIn('"electrons"', written)

    def test_the_export_pins_its_div_id(self):
        """So two exports of the same figure are byte-identical."""
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "a"
            second = Path(directory) / "b"
            plot_interactive_element("As", grid=self.grid, save_html=first)
            plot_interactive_element("As", grid=self.grid, save_html=second)
            self.assertEqual(
                first.with_suffix(".html").read_bytes(),
                second.with_suffix(".html").read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
