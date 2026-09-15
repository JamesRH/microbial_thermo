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


if __name__ == "__main__":
    unittest.main()
