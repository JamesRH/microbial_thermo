"""Iron and manganese: aqueous ions and the solid phases.

Aqueous Fe and Mn come from the HKF database. The solid phases are derived
from log K values in the GWB database (see backends/gwb.py), because no
direct-access database pyGCC ships contains the manganese oxides.

Goethite, hematite and magnetite appear in *both* the GWB database and
supcrtbl.dat, so they give two independent routes to the same number. Those
are asserted tightly against published values here; the manganese oxides,
which have only one route, are asserted more loosely.
"""

import unittest

import microbial_thermo as mt
from microbial_thermo.exceptions import OutOfRangeError
from microbial_thermo.species import resolve


class TestMineralFormationEnergies(unittest.TestCase):
    """Published dGf at 25 C, in kJ/mol."""

    IRON = [
        ("goethite", -490.6, 4.0),
        ("hematite", -742.7, 4.0),
        ("magnetite", -1012.6, 5.0),
        ("siderite", -679.0, 5.0),
        ("pyrite", -160.1, 3.0),
        ("ferrihydrite", -699.0, 6.0),
    ]
    MANGANESE = [
        ("pyrolusite", -465.1, 6.0),
        ("manganite", -557.7, 6.0),
        ("hausmannite", -1283.2, 12.0),
        ("rhodochrosite", -816.7, 5.0),
    ]

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()

    def _check(self, cases):
        for name, expected, tolerance in cases:
            with self.subTest(mineral=name):
                species = resolve(name)
                value = self.backend.delta_Gf(species.backend, 25.0).magnitude
                self.assertAlmostEqual(value, expected, delta=tolerance)

    def test_iron_minerals(self):
        self._check(self.IRON)

    def test_manganese_minerals(self):
        self._check(self.MANGANESE)

    def test_aqueous_metal_ions(self):
        for name, expected, tolerance in [
            ("Fe+2", -91.5, 3.0),
            ("Fe+3", -17.2, 3.0),
            ("Mn+2", -230.5, 3.0),
            ("Mn+3", -84.9, 5.0),
        ]:
            with self.subTest(species=name):
                species = resolve(name)
                value = self.backend.delta_Gf(species.backend, 25.0).magnitude
                self.assertAlmostEqual(value, expected, delta=tolerance)


class TestGwbRoute(unittest.TestCase):
    def test_database_parses_many_phases(self):
        self.assertGreater(len(mt.get_backend().minerals), 1000)

    def test_minerals_with_a_full_grid_vary_with_temperature(self):
        backend = mt.get_backend()
        cold = backend.delta_Gf("Goethite", 5.0).magnitude
        warm = backend.delta_Gf("Goethite", 95.0).magnitude
        self.assertNotAlmostEqual(cold, warm, places=1)

    def test_isothermal_minerals_are_flagged(self):
        backend = mt.get_backend()
        self.assertTrue(backend.minerals["Manganite"].is_isothermal)
        self.assertFalse(backend.minerals["Goethite"].is_isothermal)

    def test_isothermal_minerals_refuse_other_temperatures(self):
        """Better to fail than to invent a temperature dependence."""
        with self.assertRaises(OutOfRangeError):
            mt.get_backend().delta_Gf("Manganite", 60.0)

    def test_isothermal_minerals_work_at_their_own_temperature(self):
        value = mt.get_backend().delta_Gf("Manganite", 25.0).magnitude
        self.assertLess(value, 0.0)

    def test_minerals_appear_in_available_species(self):
        available = mt.get_backend().available_species()
        self.assertIn("Pyrolusite", available)
        self.assertIn("Goethite", available)


class TestMetalMetabolisms(unittest.TestCase):
    """Whole reactions, which also exercise the two-path cross-check."""

    METABOLISMS = [
        ("aqueous iron oxidation", ("Fe+2", "Fe+3"), ("H2O", "O2(aq)")),
        ("iron oxidation to goethite", ("Fe+2", "goethite"), ("H2O", "O2(aq)")),
        ("manganese oxidation to pyrolusite", ("Mn+2", "pyrolusite"), ("H2O", "O2(aq)")),
        ("goethite reduction with acetate", ("acetate", "CO2(aq)"), ("Fe+2", "goethite")),
        ("ferrihydrite reduction with H2", ("H2(g)", "H+"), ("Fe+2", "ferrihydrite")),
        ("pyrolusite reduction with acetate", ("acetate", "CO2(aq)"), ("Mn+2", "pyrolusite")),
    ]

    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_every_metal_metabolism_is_self_consistent(self):
        for label, donor, acceptor in self.METABOLISMS:
            with self.subTest(metabolism=label):
                reaction = mt.Reaction.from_couples(
                    donor=donor, acceptor=acceptor, conditions=self.conditions
                )
                reaction.verify_consistency()

    def test_solid_phases_are_held_at_unit_activity(self):
        """A pure solid has activity 1 whatever concentrations are set."""
        conditions = self.conditions.replace(
            concentrations={"goethite": 1e-6}, activity_model="ideal"
        )
        reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("Fe+2", "goethite"),
            conditions=conditions,
        )
        baseline = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("Fe+2", "goethite"),
            conditions=self.conditions.replace(activity_model="ideal"),
        )
        self.assertAlmostEqual(reaction.delta_G.magnitude, baseline.delta_G.magnitude, places=9)

    def test_ferrihydrite_is_easier_to_reduce_than_goethite(self):
        """Poorly crystalline Fe(III) is the more available electron acceptor,
        which is why it dominates microbial iron reduction in sediments."""

        def energy(acceptor):
            return mt.Reaction.from_couples(
                donor=("H2(g)", "H+"),
                acceptor=("Fe+2", acceptor),
                conditions=self.conditions,
            ).delta_G_standard_prime.magnitude

        self.assertLess(energy("ferrihydrite"), energy("goethite"))

    def test_oxidation_states_across_metal_couples(self):
        from fractions import Fraction

        from microbial_thermo.oxidation import mean_oxidation_state

        self.assertEqual(mean_oxidation_state("Fe", "FeOOH"), 3)
        self.assertEqual(mean_oxidation_state("Fe", "Fe3O4"), Fraction(8, 3))
        self.assertEqual(mean_oxidation_state("Mn", "MnO2"), 4)
        self.assertEqual(mean_oxidation_state("Mn", "MnOOH"), 3)
        self.assertEqual(mean_oxidation_state("Mn", "Mn3O4"), Fraction(8, 3))


class TestVanadium(unittest.TestCase):
    """Vanadium, all of it from pyGCC rather than the supplemental table."""

    @classmethod
    def setUpClass(cls):
        cls.backend = mt.get_backend()
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_formation_energies_against_published(self):
        for name, expected, tolerance in [
            ("vanadyl", -446.4, 2.0),  # VO2+ aqueous
            ("V(V)", -587.0, 3.0),  # VO2+ dioxovanadium
            ("vanadium dioxide", -1318.6, 4.0),  # V2O4
        ]:
            with self.subTest(species=name):
                species = resolve(name)
                value = self.backend.delta_Gf(species.backend, 25.0).magnitude
                self.assertAlmostEqual(value, expected, delta=tolerance)

    def test_vo2_means_the_solid_not_the_cation(self):
        """'VO2' is vanadium(IV) oxide; 'VO2+' is dioxovanadium(V). Two
        different oxidation states one character apart."""
        solid = resolve("VO2")
        cation = resolve("VO2+")
        self.assertEqual(solid.backend, "V2O4")
        self.assertEqual(solid.phase, "s")
        self.assertEqual(cation.backend, "VO2+")
        self.assertEqual(cation.phase, "aq")
        from microbial_thermo.oxidation import mean_oxidation_state

        self.assertEqual(mean_oxidation_state("V", solid.formula), 4)
        self.assertEqual(mean_oxidation_state("V", cation.formula), 5)

    def test_vanadyl_sulfate_resolves_to_the_ion_pair(self):
        species = resolve("vanadyl sulfate")
        self.assertEqual(species.backend, "VOSO4(aq)")
        self.assertEqual(species.parsed.elements, {"V": 1, "S": 1, "O": 5})

    def test_vanadyl_sulfate_is_tabulated_at_one_temperature(self):
        self.assertTrue(self.backend.minerals["VOSO4(aq)"].is_isothermal)
        with self.assertRaises(OutOfRangeError):
            self.backend.delta_Gf("VOSO4(aq)", 60.0)

    def test_vanadium_dioxide_has_a_full_temperature_grid(self):
        self.assertFalse(self.backend.minerals["V2O4"].is_isothermal)
        cold = self.backend.delta_Gf("V2O4", 5.0).magnitude
        warm = self.backend.delta_Gf("V2O4", 95.0).magnitude
        self.assertNotAlmostEqual(cold, warm, places=1)

    def test_the_two_reduction_couples(self):
        for reduced, oxidized, expected in [
            ("VO+2", "VO2+", 0.173),  # V(V)/V(IV)
            ("V+3", "VO+2", -0.486),  # V(IV)/V(III)
        ]:
            with self.subTest(couple=f"{oxidized}/{reduced}"):
                result = mt.half_reaction(reduced, oxidized, self.conditions)
                self.assertEqual(result.half.n_electrons, 1)
                self.assertAlmostEqual(
                    result.E_standard_prime.to("V").magnitude, expected, delta=0.01
                )

    def test_vanadate_reduction_is_exergonic_and_consistent(self):
        for equation in ("VO2+ + H2 -> VO++", "acetate + VO2+ -> CO2 + VO++"):
            with self.subTest(equation=equation):
                reaction = mt.Reaction.from_equation(
                    equation, conditions=self.conditions, normalize_to="integer"
                )
                reaction.verify_consistency()
                self.assertLess(reaction.delta_G_standard_prime.magnitude, 0.0)

    def test_the_solid_and_the_vanadyl_ion_are_not_a_redox_couple(self):
        """V2O4 and VO2+ are both V(IV); pairing them is a dissolution, and
        the library should say so rather than inventing a potential."""
        with self.assertRaises(ValueError) as caught:
            _ = mt.half_reaction("VO+2", "V2O4", self.conditions).E_standard
        self.assertIn("no electrons", str(caught.exception))

    def test_vanadyl_sulfate_cannot_stand_in_for_the_vanadyl_ion(self):
        """It carries sulfur, so a couple against VO2+ fails conservation --
        the right answer, since VOSO4 is an ion pair, not a redox form."""
        from microbial_thermo.exceptions import BalancingError

        with self.assertRaises(BalancingError):
            mt.half_reaction("VOSO4(aq)", "VO2+", self.conditions)


class TestVanadylSulfateAsACouple(unittest.TestCase):
    """VOSO4 can be a redox form, but only as a multi-species side.

    Sulfate has to travel with it or sulfur does not conserve, and its own
    vanadium oxidation state cannot be read from the formula, since V and S are
    both non-spectator elements.
    """

    @classmethod
    def setUpClass(cls):
        from microbial_thermo.reaction import Couple

        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
        cls.acceptor = Couple.make(["VOSO4(aq)"], ["VO2+", "SO4-2"], key_element="V")

    def test_it_balances_with_sulfate_carried_along(self):
        half = self.acceptor.half_reaction()
        self.assertEqual(half.n_electrons, 1)
        names = {s.backend for s in half.coefficients}
        self.assertIn("SO4--", names)
        self.assertIn("VOSO4(aq)", names)

    def test_the_whole_reaction_balances_and_cross_checks(self):
        from microbial_thermo.reaction import Couple, Reaction

        reaction = Reaction.from_couples(
            donor=Couple.make("H2(g)", "H+"),
            acceptor=self.acceptor,
            conditions=self.conditions,
            normalize_to="integer",
        )
        reaction.verify_consistency()
        self.assertLess(reaction.delta_G_standard_prime.magnitude, 0.0)

    def test_an_unassignable_oxidation_state_does_not_lose_the_figure(self):
        """VOSO4's vanadium cannot be assigned from the formula. That species
        goes unannotated; it must not take the whole diagram down."""
        import matplotlib

        matplotlib.use("Agg")
        from microbial_thermo.figures import plot_half_reactions
        from microbial_thermo.reaction import Couple, Reaction
        from microbial_thermo.typeset import build_equation_tokens

        reaction = Reaction.from_couples(
            donor=Couple.make("H2(g)", "H+"),
            acceptor=self.acceptor,
            conditions=self.conditions,
            normalize_to="integer",
        )
        plot_half_reactions(reaction)  # must not raise

        layout = build_equation_tokens(
            reaction.acceptor_half.half, "reduction", annotate_element="V"
        )
        annotations = {t.species.backend: t.oxidation_state for t in layout.species_tokens()}
        self.assertEqual(annotations["VO2+"], "+5")  # still labelled
        self.assertIsNone(annotations["VOSO4(aq)"])  # quietly skipped


if __name__ == "__main__":
    unittest.main()
