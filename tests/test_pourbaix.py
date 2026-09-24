"""Eh-pH predominance diagrams.

The construction is short enough to get subtly wrong, so most of these tests
check the *chemistry* of the result rather than the plumbing: does iron come
out ferric in acid and oxidising conditions, does arsenic's field structure
reproduce the pKa ladder, does lowering the dissolved activity grow the solid
fields.

The sharpest check is the last kind. The vertical boundaries of the arsenate
fields are set by acid-base equilibria, and the speciation layer computes
those pKa values by a completely different route. They agree to within a grid
cell.
"""

import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from microbial_thermo.figures.pourbaix import (
    DEFAULT_ACTIVITY,
    FIXED_DEFAULTS,
    _basis_coefficients,
    plot_pourbaix,
    pourbaix_field,
    water_stability,
)
from microbial_thermo.speciation import pKa_ladder
from microbial_thermo.species import resolve


def boundaries(field, row=-1):
    """pH values where the predominant species changes along one Eh row."""
    winner = field.winner[row]
    return [field.ph[i + 1] for i in range(len(winner) - 1) if winner[i] != winner[i + 1]]


class TestTheBasisDecomposition(unittest.TestCase):
    """n_E E + b H2O + c H+ + d e- (+ auxiliaries) -> species."""

    def test_arsenate_comes_out_as_a_five_electron_oxidation(self):
        n, b, c, d, extras = _basis_coefficients(resolve("HAsO4-2"), "As")
        self.assertEqual((n, b, c, d), (1, 4, -7, -5))
        self.assertEqual(extras, [])

    def test_hematite_counts_two_irons(self):
        n, b, c, d, _ = _basis_coefficients(resolve("Hematite"), "Fe")
        self.assertEqual((n, b, c, d), (2, 3, -6, -6))

    def test_selenate_is_a_six_electron_oxidation(self):
        n, _, _, d, _ = _basis_coefficients(resolve("SeO4-2"), "Se")
        self.assertEqual((n, d), (1, -6))

    def test_siderite_needs_a_carbonate_auxiliary(self):
        """Fe + HCO3- -> FeCO3 + H+ + 2 e-."""
        aux = [("C", resolve("HCO3-"))]
        n, b, c, d, extras = _basis_coefficients(resolve("Siderite"), "Fe", aux)
        self.assertEqual((n, b, c, d), (1, 0, -1, -2))
        self.assertEqual(len(extras), 1)
        self.assertEqual(extras[0][2], 1.0)

    def test_pyrite_needs_two_sulfates(self):
        aux = [("S", resolve("SO4-2"))]
        _, _, _, _, extras = _basis_coefficients(resolve("Pyrite"), "Fe", aux)
        self.assertEqual(extras[0][2], 2.0)

    def test_an_unaccounted_element_is_refused_not_mis_weighted(self):
        """The bug this replaced: without an auxiliary, siderite's carbon was
        simply ignored, its energy came out spuriously low, and it won the
        entire diagram."""
        with self.assertRaises(ValueError) as caught:
            _basis_coefficients(resolve("Siderite"), "Fe")
        self.assertIn("C", str(caught.exception))
        self.assertIn("auxiliary", str(caught.exception))


class TestIronDiagram(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.field = pourbaix_field("Fe", points=160)

    def test_ferric_iron_wins_in_acid_and_oxidising_conditions(self):
        self.assertEqual(self.field.at(pH=1.0, eh=0.9), "Fe+++")

    def test_ferrous_iron_wins_lower_down(self):
        self.assertEqual(self.field.at(pH=3.0, eh=0.3), "Fe++")

    def test_a_ferric_oxide_takes_the_neutral_oxidising_field(self):
        self.assertIn(self.field.at(pH=7.0, eh=0.6), {"Hematite", "Goethite", "Magnetite"})

    def test_pyrite_takes_the_reducing_field_when_sulfate_is_present(self):
        self.assertEqual(self.field.at(pH=7.0, eh=-0.4), "Pyrite")

    def test_the_sulfur_auxiliary_is_recorded_on_the_field(self):
        held = {element for element, _, _ in self.field.fixed}
        self.assertIn("S", held)

    def test_without_sulfate_there_is_no_pyrite_field(self):
        """Pyrite exists on the diagram only because sulfide is about. Drop
        the sulfur basis and the species is refused rather than drawn."""
        field = pourbaix_field("Fe", points=100, fixed={"C": ("HCO3-", 2e-3)})
        self.assertNotIn("Pyrite", field.present)


class TestArsenicReproducesItsOwnPkaLadder(unittest.TestCase):
    """The strongest check available: the vertical field boundaries are
    acid-base equilibria, and the speciation layer computes those pKa values
    by a completely different route."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        # A band high enough that only the As(V) species compete.
        cls.field = pourbaix_field("As", points=400, eh_range=(0.6, 1.2))

    def test_three_vertical_boundaries_appear(self):
        self.assertEqual(len(boundaries(self.field)), 3)

    def test_they_land_on_the_arsenate_pka_values(self):
        found = boundaries(self.field)
        expected = pKa_ladder("arsenate", temperature_c=25.0)
        # One grid cell is 14/400 = 0.035 pH units.
        for got, want in zip(found, expected, strict=True):
            with self.subTest(pKa=want):
                self.assertAlmostEqual(got, want, delta=0.06)

    def test_the_fields_run_in_the_right_order(self):
        self.assertEqual(self.field.at(pH=1.0, eh=1.0), "H3AsO4(aq)")
        self.assertEqual(self.field.at(pH=4.0, eh=1.0), "H2AsO4-")
        self.assertEqual(self.field.at(pH=9.0, eh=1.0), "HAsO4--")
        self.assertEqual(self.field.at(pH=13.0, eh=1.0), "AsO4---")

    def test_arsenite_holds_the_middle_and_the_element_the_bottom(self):
        full = pourbaix_field("As", points=200)
        self.assertEqual(full.at(pH=6.0, eh=0.1), "As(OH)3(aq)")
        self.assertEqual(full.at(pH=4.0, eh=-0.6), "As")


class TestActivityMatters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_lowering_the_activity_shrinks_the_solid_field(self):
        """Dilution favours the *dissolved* form, not the solid: at a lower
        total dissolved activity the solution is further from saturation, so
        the solid precipitates over less of the diagram. I asserted this
        backwards at first -- it is worth stating plainly, because "the solid
        is at unit activity so dilution favours it" is a tempting and wrong
        line of reasoning."""

        def solid_fraction(activity):
            field = pourbaix_field("As", points=140, activity=activity)
            index = field.species.index("As")
            return float((field.winner == index).mean())

        self.assertLess(solid_fraction(1e-8), solid_fraction(1e-4))

    def test_the_default_is_the_conventional_one(self):
        self.assertEqual(DEFAULT_ACTIVITY, 1e-6)

    def test_the_fixed_defaults_cover_sulfur_and_carbon(self):
        self.assertIn("S", FIXED_DEFAULTS)
        self.assertIn("C", FIXED_DEFAULTS)


class TestWaterStability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.ph, cls.upper, cls.lower = water_stability()

    def test_the_lower_line_passes_through_zero_at_ph_zero(self):
        """By definition of the standard hydrogen electrode."""
        self.assertAlmostEqual(self.lower[0], 0.0, delta=0.01)

    def test_the_upper_line_reproduces_the_published_1_229_volts(self):
        """Exactly, because it uses O2(g) at 1 bar -- the state the published
        value refers to. Against O2(aq) at unit activity it comes out 43 mV
        higher, which is a different reference state rather than an error."""
        self.assertAlmostEqual(self.upper[0], 1.229, delta=0.005)

    def test_both_fall_at_59_millivolts_per_ph_unit(self):
        for line in (self.upper, self.lower):
            slope = (line[-1] - line[0]) / (self.ph[-1] - self.ph[0])
            with self.subTest(slope=slope):
                self.assertAlmostEqual(slope, -0.0592, delta=0.002)

    def test_they_stay_the_right_way_up(self):
        self.assertTrue((self.upper > self.lower).all())


class TestOtherElements(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_every_default_element_draws_something(self):
        for element in ("As", "Fe", "Mn", "S", "Se", "N", "Cr", "U", "Cu"):
            with self.subTest(element=element):
                field = pourbaix_field(element, points=90)
                self.assertGreaterEqual(len(field.present), 2)

    def test_selenium_reduces_to_the_element_at_the_bottom(self):
        """The contrast with arsenic: Se(0) is where selenium ends up, which
        is why selenium is remediable by reduction and arsenic is not."""
        field = pourbaix_field("Se", points=140)
        self.assertEqual(field.at(pH=7.0, eh=-0.4), "Se")

    def test_sulfur_splits_into_sulfate_and_sulfide(self):
        field = pourbaix_field("S", points=140)
        self.assertEqual(field.at(pH=7.0, eh=0.6), "SO4--")
        self.assertIn(field.at(pH=7.0, eh=-0.4), {"H2S(aq)", "HS-"})

    def test_an_element_with_no_default_list_says_so(self):
        with self.assertRaises(ValueError) as caught:
            pourbaix_field("Zr", points=20)
        self.assertIn("species=", str(caught.exception))


class TestFigure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.field = pourbaix_field("Fe", points=120)

    def test_it_renders_and_saves(self):
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "fe"
            figure, _ = plot_pourbaix("Fe", field=self.field, save=base)
            try:
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
            finally:
                plt.close(figure)

    def test_the_title_states_the_activity_and_what_was_held_fixed(self):
        import matplotlib.pyplot as plt

        figure, _ = plot_pourbaix("Fe", field=self.field)
        try:
            title = figure.axes[0].get_title()
            self.assertIn("activity", title)
            self.assertIn("held fixed", title)
        finally:
            plt.close(figure)

    def test_the_water_lines_can_be_turned_off(self):
        import matplotlib.pyplot as plt

        with_lines, _ = plot_pourbaix("Fe", field=self.field)
        without, _ = plot_pourbaix("Fe", field=self.field, show_water=False)
        try:
            self.assertGreater(len(with_lines.axes[0].lines), len(without.axes[0].lines))
        finally:
            plt.close("all")


if __name__ == "__main__":
    unittest.main()
