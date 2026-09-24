"""Light as an energy input, and the light-budget figure.

A phototroph runs a reaction that does not pay and makes up the shortfall
from absorbed photons. The photon arithmetic is simple enough that the risk
is not getting it wrong but *overclaiming* with it, so several of these tests
pin the caveats rather than the numbers.
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
from microbial_thermo.figures.light import (
    light_budget,
    light_budget_table,
    plot_light_budget,
)
from microbial_thermo.photons import (
    REACTION_CENTRES,
    light_available,
    net_with_light,
    photon_energy,
    photons_required,
)
from microbial_thermo.reaction import Couple, Reaction
from microbial_thermo.units import Quantity

#: Textbook einstein energies, kJ/mol, from E = N_A h c / lambda.
PUBLISHED = {680: 175.9, 700: 170.9, 870: 137.5, 960: 124.6}


def fixation(donor, pH=7.0):
    return Reaction.from_couples(
        donor=donor,
        acceptor=Couple.make("Biomass", "CO2"),
        conditions=mt.Conditions(temperature_c=25.0, pH=pH),
    )


IRON = Couple.make("Fe++", "Fe(OH)3", key_element="Fe")
ARSENITE = Couple.make("As(III)", "As(V)")
NITRITE = Couple.make("Nitrite", "Nitrate")


class TestPhotonEnergy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_it_matches_the_textbook_einstein(self):
        for wavelength, expected in PUBLISHED.items():
            with self.subTest(nm=wavelength):
                self.assertAlmostEqual(photon_energy(wavelength).magnitude, expected, delta=0.1)

    def test_reaction_centre_names_work_too(self):
        self.assertAlmostEqual(
            photon_energy("P870").magnitude, photon_energy(870).magnitude, places=9
        )

    def test_names_are_case_insensitive(self):
        self.assertAlmostEqual(
            photon_energy("p870").magnitude, photon_energy("P870").magnitude, places=9
        )

    def test_shorter_wavelengths_carry_more(self):
        ordered = [photon_energy(nm).magnitude for nm in sorted(REACTION_CENTRES.values())]
        self.assertEqual(ordered, sorted(ordered, reverse=True))

    def test_an_unknown_centre_says_what_is_known(self):
        with self.assertRaises(MicrobialThermoError) as caught:
            photon_energy("P123")
        self.assertIn("P870", str(caught.exception))

    def test_a_nonsense_wavelength_is_refused(self):
        for bad in (0, -500):
            with self.subTest(nm=bad), self.assertRaises(MicrobialThermoError):
                photon_energy(bad)

    def test_the_result_carries_units(self):
        self.assertEqual(str(photon_energy(870).units), "kilojoule / mole")


class TestLightAvailable(unittest.TestCase):
    def test_photons_add_linearly(self):
        self.assertAlmostEqual(
            light_available(3).magnitude, 3 * light_available(1).magnitude, places=9
        )

    def test_efficiency_scales_it(self):
        self.assertAlmostEqual(
            light_available(1, efficiency=0.5).magnitude,
            light_available(1, efficiency=1.0).magnitude / 2,
            places=9,
        )

    def test_efficiency_must_be_a_fraction(self):
        for bad in (0.0, -0.1, 1.5):
            with self.subTest(efficiency=bad), self.assertRaises(MicrobialThermoError):
                light_available(1, efficiency=bad)

    def test_negative_photons_are_refused(self):
        with self.assertRaises(MicrobialThermoError):
            light_available(-1)


class TestNetWithLight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_light_is_credited_not_charged(self):
        """An uphill reaction must become less positive, not more."""
        uphill = fixation(ARSENITE).delta_G_standard_prime
        self.assertLess(net_with_light(uphill, 1).magnitude, uphill.magnitude)

    def test_one_photon_covers_arsenite_fixation(self):
        uphill = fixation(ARSENITE).delta_G_standard_prime
        self.assertGreater(uphill.magnitude, 0.0)
        self.assertLess(net_with_light(uphill, 1).magnitude, 0.0)

    def test_one_photon_does_not_cover_nitrite_fixation(self):
        """The sharpest contrast among the three donors: nitrite is weak
        enough that a single 870 nm photon leaves it uphill."""
        uphill = fixation(NITRITE).delta_G_standard_prime
        self.assertGreater(net_with_light(uphill, 1).magnitude, 0.0)
        self.assertLess(net_with_light(uphill, 2).magnitude, 0.0)


class TestPhotonsRequired(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_a_reaction_that_already_pays_needs_none(self):
        self.assertEqual(photons_required(Quantity(-150.0, "kJ/mol")), 0.0)

    def test_it_measures_down_to_the_quantum_not_to_zero(self):
        """A reaction resting at exactly zero does no useful work, so the
        target is the energy quantum. That makes the count slightly larger
        than a naive divide."""
        naive = 100.0 / photon_energy(870).magnitude
        self.assertGreater(photons_required(Quantity(100.0, "kJ/mol")), naive)

    def test_the_three_donors_rank_as_expected(self):
        counts = {
            name: photons_required(fixation(couple).delta_G_standard_prime)
            for name, couple in (("Fe", IRON), ("As", ARSENITE), ("NO2", NITRITE))
        }
        self.assertLess(counts["Fe"], counts["As"])
        self.assertLess(counts["As"], counts["NO2"])

    def test_only_nitrite_needs_more_than_one(self):
        for name, couple, within_one in (
            ("Fe", IRON, True),
            ("As", ARSENITE, True),
            ("NO2", NITRITE, False),
        ):
            with self.subTest(donor=name):
                needed = photons_required(fixation(couple).delta_G_standard_prime)
                self.assertEqual(needed <= 1.0, within_one)

    def test_a_lower_efficiency_needs_more_photons(self):
        uphill = fixation(ARSENITE).delta_G_standard_prime
        self.assertGreater(
            photons_required(uphill, efficiency=0.4),
            photons_required(uphill, efficiency=1.0),
        )


class TestLightBudgetFigure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.reaction = fixation(ARSENITE)
        cls.result = light_budget(cls.reaction, ph_values=np.linspace(4.0, 10.0, 7))

    def test_it_returns_the_dark_curve_and_one_per_photon_count(self):
        for key in ("dark", 1, 2):
            with self.subTest(key=key):
                self.assertIn(key, self.result)

    def test_light_shifts_the_curve_by_exactly_the_photon_energy(self):
        shift = self.result["dark"] - self.result[1]
        np.testing.assert_allclose(shift, photon_energy(870).magnitude, rtol=1e-9)

    def test_two_photons_shift_twice_as_far(self):
        np.testing.assert_allclose(
            self.result["dark"] - self.result[2],
            2 * (self.result["dark"] - self.result[1]),
            rtol=1e-9,
        )

    def test_the_dark_curve_really_is_uphill(self):
        self.assertTrue((self.result["dark"] > 0).all())

    def test_higher_ph_is_kinder(self):
        """Protons are released, so raising pH lowers the cost."""
        self.assertLess(self.result["dark"][-1], self.result["dark"][0])

    def test_the_table_ranks_the_donors(self):
        frame = light_budget_table({"As": fixation(ARSENITE), "NO2": fixation(NITRITE)})
        self.assertEqual(len(frame), 2)
        self.assertIn("photons needed", frame.columns)
        self.assertTrue((frame["dG0' (kJ/mol)"] > 0).all())

    def test_it_renders_and_saves(self):
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "light"
            figure, _ = plot_light_budget(self.reaction, result=self.result, save=base)
            try:
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
            finally:
                plt.close(figure)

    def test_the_figure_is_footnoted_for_the_placeholder(self):
        import matplotlib.pyplot as plt

        figure, _ = plot_light_budget(self.reaction, result=self.result)
        try:
            texts = " ".join(t.get_text() for t in figure.texts)
            self.assertIn("Biomass(aq)", texts)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
