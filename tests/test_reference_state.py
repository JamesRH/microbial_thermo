"""The biochemical reference state, and what ionic strength does and does not do.

The motivating question was whether defaulting to I = 0.25 M would bring this
library's numbers closer to Alberty-convention biochemical tables. Measured,
the answer is no -- not because the correction is small, but because it does
not apply: dG0' is defined at unit activity, where no activity coefficient
appears at all. These tests pin that, so the claim cannot quietly come back.
"""

import unittest
import warnings

import microbial_thermo as mt
from microbial_thermo.conditions import (
    BIOCHEMICAL_IONIC_STRENGTH,
    SEAWATER_IONIC_STRENGTH,
)
from microbial_thermo.library import reaction

CONCENTRATIONS = {"HS-": 1e-4, "SO4--": 1e-3}
PRESSURES = {"H2(g)": 1e-5}
METABOLISM = "hydrogenotrophic_sulfate_reduction"


class TestBiochemicalConstructor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_the_defaults_are_the_biochemical_reference_state(self):
        conditions = mt.Conditions.biochemical()
        self.assertEqual(conditions.pH, 7.0)
        self.assertEqual(conditions.temperature_c, 25.0)
        self.assertEqual(conditions.ionic_strength, BIOCHEMICAL_IONIC_STRENGTH)
        self.assertEqual(conditions.activity_model, "bdot")

    def test_the_constant_is_the_published_one(self):
        self.assertEqual(BIOCHEMICAL_IONIC_STRENGTH, 0.25)

    def test_seawater_is_offered_too(self):
        conditions = mt.Conditions.biochemical(ionic_strength=SEAWATER_IONIC_STRENGTH)
        self.assertEqual(conditions.ionic_strength, 0.7)

    def test_anything_can_be_overridden(self):
        conditions = mt.Conditions.biochemical(
            pH=8.2, temperature_c=4.0, concentrations={"HS-": 1e-6}
        )
        self.assertEqual(conditions.pH, 8.2)
        self.assertEqual(conditions.temperature_c, 4.0)
        self.assertEqual(conditions.concentrations, {"HS-": 1e-6})

    def test_the_activity_model_can_be_forced_back_to_ideal(self):
        """Because b-dot needs a positive ionic strength, and someone may want
        the reference pH and temperature without the correction."""
        conditions = mt.Conditions.biochemical(activity_model="ideal")
        self.assertEqual(conditions.activity_model, "ideal")


class TestIonicStrengthDoesNotMoveTheStandardState(unittest.TestCase):
    """The finding that corrected the plan.

    An ionic strength scales activity *coefficients*. dG0' is defined at unit
    activity, so there is nothing for a coefficient to multiply, and the shift
    is exactly zero -- not small, zero.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_the_standard_prime_energy_is_untouched(self):
        at_zero = reaction(
            METABOLISM, mt.Conditions(temperature_c=25.0, pH=7.0)
        ).delta_G_standard_prime.magnitude
        at_quarter = reaction(
            METABOLISM, mt.Conditions.biochemical()
        ).delta_G_standard_prime.magnitude
        self.assertEqual(at_zero, at_quarter)

    def test_that_holds_at_seawater_strength_too(self):
        at_zero = reaction(
            METABOLISM, mt.Conditions(temperature_c=25.0, pH=7.0)
        ).delta_G_standard_prime.magnitude
        at_sea = reaction(
            METABOLISM,
            mt.Conditions.biochemical(ionic_strength=SEAWATER_IONIC_STRENGTH),
        ).delta_G_standard_prime.magnitude
        self.assertEqual(at_zero, at_sea)

    def test_but_it_does_move_a_concentration_dependent_energy(self):
        """Otherwise the machinery would simply be broken."""
        ideal = reaction(
            METABOLISM,
            mt.Conditions(
                temperature_c=25.0,
                pH=7.0,
                concentrations=CONCENTRATIONS,
                partial_pressures=PRESSURES,
            ),
        ).delta_G.magnitude
        corrected = reaction(
            METABOLISM,
            mt.Conditions.biochemical(concentrations=CONCENTRATIONS, partial_pressures=PRESSURES),
        ).delta_G.magnitude
        self.assertNotAlmostEqual(ideal, corrected, places=2)


class TestActivityCorrection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def _built(self, ionic_strength):
        return reaction(
            METABOLISM,
            mt.Conditions.biochemical(
                ionic_strength=ionic_strength,
                concentrations=CONCENTRATIONS,
                partial_pressures=PRESSURES,
            ),
        )

    def test_an_ideal_model_has_no_correction(self):
        built = reaction(
            METABOLISM,
            mt.Conditions(temperature_c=25.0, pH=7.0, concentrations=CONCENTRATIONS),
        )
        self.assertEqual(built.activity_correction().magnitude, 0.0)

    def test_it_equals_the_gap_to_the_ideal_calculation(self):
        built = self._built(0.25)
        ideal = reaction(
            METABOLISM,
            mt.Conditions(
                temperature_c=25.0,
                pH=7.0,
                concentrations=CONCENTRATIONS,
                partial_pressures=PRESSURES,
            ),
        )
        self.assertAlmostEqual(
            built.activity_correction().magnitude,
            built.delta_G.magnitude - ideal.delta_G.magnitude,
            places=9,
        )

    def test_it_grows_with_ionic_strength(self):
        weak = abs(self._built(0.1).activity_correction().magnitude)
        strong = abs(self._built(0.7).activity_correction().magnitude)
        self.assertGreater(strong, weak)

    def test_it_stays_small_for_this_reaction(self):
        """Under a kJ/mol between fresh water and seawater, which is the point
        worth making: ionic strength is a real correction but not the reason
        two sources disagree by ten."""
        self.assertLess(abs(self._built(0.7).activity_correction().magnitude), 2.0)

    def test_it_carries_units(self):
        self.assertEqual(str(self._built(0.25).activity_correction().units), "kilojoule / mole")


if __name__ == "__main__":
    unittest.main()
