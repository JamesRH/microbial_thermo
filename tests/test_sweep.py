"""Sweeps over environmental variables."""

import unittest

import numpy as np

import microbial_thermo as mt
from microbial_thermo.sweep import (
    concentration_axis,
    default_axes,
    partial_pressure_axis,
    ph_axis,
    sweep,
    temperature_axis,
)


class TestAxes(unittest.TestCase):
    def test_temperature_axis_avoids_absolute_zero_celsius(self):
        """Exactly 0 C breaks the IAPWS-95 water equation of state."""
        axis = temperature_axis(low=0.0, high=100.0)
        self.assertGreaterEqual(axis.values.min(), 0.01)
        self.assertLessEqual(axis.values.max(), 100.0)

    def test_concentration_axes_are_logarithmic(self):
        self.assertTrue(concentration_axis("SO4-2").log_scale)
        self.assertTrue(partial_pressure_axis("H2(g)").log_scale)
        self.assertFalse(ph_axis().log_scale)

    def test_apply_sets_the_right_field(self):
        base = mt.Conditions(temperature_c=25.0, pH=7.0)
        self.assertEqual(ph_axis().apply(base, 5.0).pH, 5.0)
        self.assertEqual(temperature_axis().apply(base, 60.0).temperature_c, 60.0)
        applied = partial_pressure_axis("H2(g)").apply(base, 1e-5)
        self.assertEqual(applied.partial_pressures["H2(g)"], 1e-5)

    def test_apply_does_not_mutate_the_original(self):
        base = mt.Conditions(temperature_c=25.0, pH=7.0)
        partial_pressure_axis("H2(g)").apply(base, 1e-5)
        self.assertEqual(base.partial_pressures, {})


class TestSweep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(
            temperature_c=25.0,
            pH=7.0,
            partial_pressures={"H2(g)": 1e-4},
            activity_model="ideal",
        )
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("methane", "CO2(aq)"),
            conditions=cls.conditions,
        )

    def test_sweep_returns_one_value_per_grid_point(self):
        result = sweep(self.reaction, ph_axis())
        self.assertEqual(len(result.delta_g), len(result.axis.values))
        self.assertEqual(result.failures, [])

    def test_hydrogen_sweep_crosses_zero_at_the_threshold(self):
        """The classic result: methanogenesis stops paying below ~1e-5 bar H2."""
        result = sweep(self.reaction, partial_pressure_axis("H2(g)"))
        self.assertLess(result.delta_g.min(), 0.0)
        self.assertGreater(result.delta_g.max(), 0.0)

    def test_free_energy_rises_monotonically_as_hydrogen_falls(self):
        result = sweep(self.reaction, partial_pressure_axis("H2(g)"))
        self.assertTrue(np.all(np.diff(result.delta_g) < 0))

    def test_ph_has_no_effect_when_protons_cancel(self):
        """H2/CO2 methanogenesis has no net H+, so it must be pH-independent."""
        result = sweep(self.reaction, ph_axis())
        self.assertLess(float(np.ptp(result.delta_g)), 1e-6)

    def test_ph_matters_when_protons_do_not_cancel(self):
        """Sulfate reduction consumes a net proton, so pH must move dG."""
        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            normalize_to="integer",
        )
        result = sweep(reaction, ph_axis())
        self.assertGreater(float(np.ptp(result.delta_g)), 1.0)

    def test_temperature_sweep_spans_the_supported_range(self):
        result = sweep(self.reaction, temperature_axis())
        self.assertEqual(result.failures, [])
        self.assertTrue(np.all(np.isfinite(result.delta_g)))

    def test_per_electron_is_dg_over_n(self):
        result = sweep(self.reaction, ph_axis())
        np.testing.assert_allclose(
            result.delta_g_per_electron,
            result.delta_g / result.n_electrons,
            rtol=1e-12,
        )

    def test_default_axes_cover_ph_temperature_and_participants(self):
        labels = [a.label for a in default_axes(self.reaction)]
        self.assertIn("pH", labels)
        self.assertIn("temperature (°C)", labels)
        self.assertTrue(any("H2(g)" in label for label in labels))
        # Water and protons are handled by pH and water activity, not swept.
        self.assertFalse(any("H2O" in label for label in labels))

    def test_frame_round_trip(self):
        frame = sweep(self.reaction, ph_axis()).to_frame()
        self.assertEqual(len(frame), 33)
        self.assertIn("dG (kJ/mol)", frame.columns)


class TestExplorer(unittest.TestCase):
    def test_explorer_builds_a_dropdown_per_axis(self):
        from microbial_thermo.figures.explorer import plot_energy_explorer

        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
            normalize_to="integer",
        )
        figure = plot_energy_explorer(reaction)
        buttons = figure.layout.updatemenus[0].buttons
        self.assertGreaterEqual(len(buttons), 3)
        self.assertEqual(buttons[0].label, "pH")
        self.assertIsNotNone(figure.layout.yaxis2)

    def test_svg_config_requests_svg(self):
        from microbial_thermo.figures.explorer import SVG_CONFIG

        self.assertEqual(SVG_CONFIG["toImageButtonOptions"]["format"], "svg")


if __name__ == "__main__":
    unittest.main()
