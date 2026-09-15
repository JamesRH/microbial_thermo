"""Thermodynamics against the pyGCC backend and published reference values.

Two kinds of assertion appear here, deliberately distinguished:

*Tight* assertions are on quantities fixed by convention or by internal
consistency -- dGf(H+) = 0, E(SHE) = 0, the two-path cross-check. These must
hold to numerical precision.

*Loose* assertions compare against textbook tables. Published values differ by
tens of millivolts between sources because of differing standard states,
database vintages, and whether a gas or aqueous form was used. Asserting
equality against them would encode one textbook's conventions as truth, so
these are range checks with a stated tolerance.
"""

import unittest

import microbial_thermo as mt
from microbial_thermo.exceptions import OutOfRangeError, SpeciesNotFoundError
from microbial_thermo.reaction import Couple, HalfReactionResult

#: Tolerance for comparisons against published tables, in volts.
TEXTBOOK_TOLERANCE_V = 0.035


def half_result(reduced, oxidized, conditions, key_element=None):
    couple = Couple.make(reduced, oxidized, key_element)
    return HalfReactionResult(
        couple=couple,
        half=couple.half_reaction(),
        conditions=conditions,
        backend=mt.get_backend(),
    )


class TestFormationEnergies(unittest.TestCase):
    """Spot values against standard tables, in kJ/mol at 25 C."""

    REFERENCE = [
        ("H+", 0.0, 0.01),  # zero by convention, exactly
        ("H2O", -237.1, 0.5),
        ("SO4--", -744.5, 1.0),
        ("HS-", 12.0, 1.0),
        ("Acetate", -369.3, 1.5),
        ("CO2(aq)", -386.0, 1.0),
        ("O2(aq)", 16.4, 1.0),
        ("H2(aq)", 17.6, 1.0),
        ("CO2(g)", -394.4, 0.5),
        ("CH4(g)", -50.7, 1.0),
        ("H2(g)", 0.0, 0.01),  # element in its standard state
        ("O2(g)", 0.0, 0.01),
    ]

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()

    def test_reference_values(self):
        for name, expected, tolerance in self.REFERENCE:
            with self.subTest(species=name):
                value = self.backend.delta_Gf(name, 25.0).magnitude
                self.assertAlmostEqual(value, expected, delta=tolerance)

    def test_water_is_liquid_across_the_range(self):
        """Water must stay in the liquid field from 0.01 to 100 C."""
        for t in (0.01, 25.0, 60.0, 99.0, 100.0):
            with self.subTest(temperature=t):
                self.assertLess(self.backend.delta_Gf("H2O", t).magnitude, -230.0)

    def test_hundred_celsius_uses_saturation_pressure(self):
        """At 100 C, 1 bar puts water in the vapour field and yields NaN."""
        value = self.backend.delta_Gf("SO4--", 100.0).magnitude
        self.assertTrue(-760.0 < value < -730.0)

    def test_below_range_raises(self):
        with self.assertRaises(OutOfRangeError):
            self.backend.delta_Gf("SO4--", 0.0)

    def test_above_range_raises(self):
        with self.assertRaises(OutOfRangeError):
            self.backend.delta_Gf("SO4--", 150.0)

    def test_unknown_species_suggests_alternatives(self):
        with self.assertRaises(SpeciesNotFoundError) as caught:
            self.backend.delta_Gf("SO4", 25.0)
        self.assertTrue(caught.exception.suggestions)

    def test_caching_returns_identical_values(self):
        first = self.backend.delta_Gf("SO4--", 25.0).magnitude
        second = self.backend.delta_Gf("SO4--", 25.0).magnitude
        self.assertEqual(first, second)


class TestSolventProperties(unittest.TestCase):
    def test_debye_huckel_parameters_at_25C(self):
        props = mt.get_backend().solvent_properties(25.0)
        self.assertAlmostEqual(props["A"], 0.5114, delta=0.002)
        self.assertAlmostEqual(props["B"], 0.3288, delta=0.002)
        self.assertAlmostEqual(props["bdot"], 0.0374, delta=0.002)
        self.assertAlmostEqual(props["dielectric_constant"], 78.24, delta=0.2)

    def test_activity_coefficients_at_ionic_strength_0p1(self):
        """Values every geochemistry text tabulates for I = 0.1 molal."""
        backend = mt.get_backend()
        gamma_na = backend.activity_coefficient(1, 0.1, 25.0, species_name="Na+")
        gamma_so4 = backend.activity_coefficient(-2, 0.1, 25.0, species_name="SO4--")
        self.assertAlmostEqual(gamma_na, 0.77, delta=0.03)
        self.assertAlmostEqual(gamma_so4, 0.35, delta=0.05)

    def test_neutral_species_are_ideal(self):
        gamma = mt.get_backend().activity_coefficient(0, 0.5, 25.0)
        self.assertEqual(gamma, 1.0)

    def test_infinite_dilution_gives_unity(self):
        gamma = mt.get_backend().activity_coefficient(-2, 0.0, 25.0)
        self.assertEqual(gamma, 1.0)


class TestHalfReactionPotentials(unittest.TestCase):
    """E-standard-prime at pH 7 against published redox-tower values."""

    CASES = [
        # (reduced, oxidized, textbook E0' in volts)
        ("H2(g)", "H+", -0.414),
        ("H2O", "O2(g)", +0.816),
        ("NO2-", "NO3-", +0.430),
        ("N2(aq)", "NO3-", +0.740),
        ("NH4+", "NO3-", +0.360),
        ("HS-", "SO4-2", -0.217),
        ("methane", "CO2(aq)", -0.244),
        ("Fe+2", "Fe+3", +0.770),
        ("HS-", "Sulfur(s)", -0.270),
    ]

    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_against_published_tables(self):
        for reduced, oxidized, expected in self.CASES:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                result = half_result(reduced, oxidized, self.conditions)
                value = result.E_standard_prime.to("V").magnitude
                self.assertAlmostEqual(value, expected, delta=TEXTBOOK_TOLERANCE_V)

    def test_standard_hydrogen_electrode_is_zero(self):
        """Fixed by definition, so this is a tight assertion."""
        result = half_result("H2(g)", "H+", self.conditions)
        self.assertAlmostEqual(result.E_standard.to("V").magnitude, 0.0, places=6)

    def test_primed_hydrogen_potential_is_exactly_nernstian(self):
        """E0' = -2.303 RT/F * 7 for the proton couple, to the millivolt."""
        result = half_result("H2(g)", "H+", self.conditions)
        self.assertAlmostEqual(result.E_standard_prime.to("V").magnitude, -0.4139, places=3)

    def test_ph_independent_couple_is_unshifted(self):
        """Fe3+/Fe2+ involves no protons, so E0 and E0' must coincide."""
        result = half_result("Fe+2", "Fe+3", self.conditions)
        self.assertAlmostEqual(
            result.E_standard.to("V").magnitude,
            result.E_standard_prime.to("V").magnitude,
            places=9,
        )

    def test_ph_shifts_proton_coupled_potential_downward(self):
        """Raising pH must lower the potential of a proton-consuming couple."""
        low = half_result("HS-", "SO4-2", mt.Conditions(pH=5.0)).E_standard_prime
        high = half_result("HS-", "SO4-2", mt.Conditions(pH=9.0)).E_standard_prime
        self.assertLess(high.to("V").magnitude, low.to("V").magnitude)


class TestReactionConsistency(unittest.TestCase):
    """The two-path cross-check, over a range of metabolisms."""

    METABOLISMS = [
        ("hydrogenotrophic sulfate reduction", ("H2(g)", "H+"), ("HS-", "SO4-2")),
        ("hydrogenotrophic methanogenesis", ("H2(g)", "H+"), ("methane", "CO2(aq)")),
        ("aerobic acetate oxidation", ("acetate", "CO2(aq)"), ("H2O", "O2(aq)")),
        ("denitrification with H2", ("H2(g)", "H+"), ("N2(aq)", "NO3-")),
        ("aerobic ammonium oxidation", ("NH4+", "NO2-"), ("H2O", "O2(aq)")),
        ("iron oxidation with oxygen", ("Fe+2", "Fe+3"), ("H2O", "O2(aq)")),
        ("sulfide oxidation with nitrate", ("HS-", "SO4-2"), ("N2(aq)", "NO3-")),
    ]

    def test_cross_check_passes_for_every_metabolism(self):
        """Construction runs verify_consistency(), so this must not raise."""
        conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
        for label, donor, acceptor in self.METABOLISMS:
            with self.subTest(metabolism=label):
                reaction = mt.Reaction.from_couples(
                    donor=donor, acceptor=acceptor, conditions=conditions
                )
                reaction.verify_consistency()

    def test_cross_check_holds_away_from_standard_state(self):
        conditions = mt.Conditions(
            temperature_c=60.0,
            pH=8.2,
            concentrations={"SO4-2": 2.8e-2, "HS-": 1e-6},
            partial_pressures={"H2(g)": 1e-5},
            activity_model="ideal",
        )
        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"), acceptor=("HS-", "SO4-2"), conditions=conditions
        )
        reaction.verify_consistency()

    def test_reaction_is_balanced(self):
        reaction = mt.Reaction.from_couples(
            donor=("methane", "CO2(aq)"), acceptor=("H2O", "O2(aq)"), n_electrons=8
        )
        from microbial_thermo.balance import verify_conservation

        verify_conservation(reaction.coefficients)


class TestReactionEnergetics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_aerobic_acetate_oxidation_energy(self):
        """Published dG0' is near -850 to -895 kJ/mol per acetate."""
        reaction = mt.Reaction.from_couples(
            donor=("acetate", "CO2(aq)"),
            acceptor=("H2O", "O2(aq)"),
            conditions=self.conditions,
            normalize_to="donor",
        )
        value = reaction.delta_G_standard_prime.magnitude
        self.assertTrue(-920.0 < value < -820.0, f"dG0' was {value}")

    def test_hydrogenotrophic_sulfate_reduction_energy(self):
        """Published dG0' is near -150 kJ/mol for the 8-electron reaction."""
        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            n_electrons=8,
        )
        value = reaction.delta_G_standard_prime.magnitude
        self.assertTrue(-175.0 < value < -130.0, f"dG0' was {value}")

    def test_per_electron_normalization(self):
        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            n_electrons=8,
        )
        self.assertAlmostEqual(
            reaction.delta_G_per_electron.magnitude,
            reaction.delta_G.magnitude / 8.0,
            places=9,
        )

    def test_electron_pair_is_the_default_normalization(self):
        reaction = mt.Reaction.from_couples(donor=("H2(g)", "H+"), acceptor=("HS-", "SO4-2"))
        self.assertEqual(reaction.n_electrons, 2)

    def test_scaling_leaves_intensive_quantities_unchanged(self):
        """dG scales with the reaction; dE and dG-per-electron must not."""
        pair = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            n_electrons=2,
        )
        octet = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("HS-", "SO4-2"),
            conditions=self.conditions,
            n_electrons=8,
        )
        self.assertAlmostEqual(
            pair.delta_E_standard_prime.to("V").magnitude,
            octet.delta_E_standard_prime.to("V").magnitude,
            places=9,
        )
        self.assertAlmostEqual(
            pair.delta_G_per_electron.magnitude,
            octet.delta_G_per_electron.magnitude,
            places=6,
        )
        self.assertAlmostEqual(
            octet.delta_G_standard_prime.magnitude,
            4.0 * pair.delta_G_standard_prime.magnitude,
            places=6,
        )

    def test_low_hydrogen_makes_methanogenesis_less_favourable(self):
        """The syntrophy argument: dropping pH2 raises dG toward zero."""

        def energy(p_h2):
            conditions = mt.Conditions(
                temperature_c=25.0,
                pH=7.0,
                partial_pressures={"H2(g)": p_h2},
                activity_model="ideal",
            )
            return mt.Reaction.from_couples(
                donor=("H2(g)", "H+"),
                acceptor=("methane", "CO2(aq)"),
                conditions=conditions,
            ).delta_G.magnitude

        self.assertLess(energy(1.0), energy(1e-4))
        self.assertLess(energy(1e-4), energy(1e-7))


class TestVersions(unittest.TestCase):
    def test_versions_reports_pygcc(self):
        versions = mt.versions()
        self.assertIn("pygcc", versions)
        self.assertNotEqual(versions["pygcc"], "not installed")
        self.assertEqual(versions["microbial_thermo"], mt.__version__)


class TestHalfReactionConvenience(unittest.TestCase):
    """``half_reaction()`` exists because reaching one potential otherwise
    takes six lines of construction."""

    def test_matches_the_long_form_exactly(self):
        from microbial_thermo.reaction import Couple, HalfReactionResult

        short = mt.half_reaction("HS-", "SO4-2")
        couple = Couple.make("HS-", "SO4-2")
        long = HalfReactionResult(
            couple=couple,
            half=couple.half_reaction(),
            conditions=mt.Conditions(),
            backend=mt.get_backend(),
        )
        self.assertEqual(short.E_standard, long.E_standard)
        self.assertEqual(short.E_standard_prime, long.E_standard_prime)
        self.assertEqual(short.half.coefficients, long.half.coefficients)

    def test_standard_potential_ignores_the_ph_field(self):
        """E_standard holds every activity at 1, the proton included."""
        for ph in (0.0, 7.0, 11.0):
            with self.subTest(pH=ph):
                value = (
                    mt.half_reaction("HS-", "SO4-2", mt.Conditions(temperature_c=25.0, pH=ph))
                    .E_standard.to("V")
                    .magnitude
                )
                self.assertAlmostEqual(value, 0.2491, places=3)

    def test_primed_potential_follows_the_ph_field(self):
        low = mt.half_reaction("HS-", "SO4-2", mt.Conditions(pH=5.0))
        high = mt.half_reaction("HS-", "SO4-2", mt.Conditions(pH=9.0))
        self.assertGreater(
            low.E_standard_prime.to("V").magnitude,
            high.E_standard_prime.to("V").magnitude,
        )

    def test_temperature_is_honoured(self):
        cold = mt.half_reaction("HS-", "SO4-2", mt.Conditions(temperature_c=5.0))
        hot = mt.half_reaction("HS-", "SO4-2", mt.Conditions(temperature_c=95.0))
        self.assertNotAlmostEqual(
            cold.E_standard.to("V").magnitude,
            hot.E_standard.to("V").magnitude,
            places=3,
        )

    def test_rescaling_changes_the_equation_but_not_the_potential(self):
        """E is intensive, so n only affects how the half reaction reads."""
        default = mt.half_reaction("HS-", "SO4-2")
        pair = mt.half_reaction("HS-", "SO4-2", n_electrons=2)
        self.assertEqual(pair.half.n_electrons, 2)
        self.assertEqual(default.half.n_electrons, 8)
        self.assertAlmostEqual(
            default.E_standard.to("V").magnitude,
            pair.E_standard.to("V").magnitude,
            places=9,
        )

    def test_accepts_a_multi_species_side(self):
        result = mt.half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])
        self.assertEqual(result.half.n_electrons, 6)
        self.assertFalse(result.half.is_simple)

    def test_key_element_hint_is_passed_through(self):
        result = mt.half_reaction("acetate", "CO2(aq)", key_element="C")
        self.assertEqual(result.half.key_element, "C")

    def test_reproduces_published_potentials(self):
        for reduced, oxidized, expected in [
            ("H2(g)", "H+", -0.414),
            ("H2O", "O2(g)", 0.816),
            ("Fe+2", "Fe+3", 0.770),
        ]:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                value = (
                    mt.half_reaction(reduced, oxidized, mt.Conditions(temperature_c=25.0, pH=7.0))
                    .E_standard_prime.to("V")
                    .magnitude
                )
                self.assertAlmostEqual(value, expected, delta=0.035)


if __name__ == "__main__":
    unittest.main()
