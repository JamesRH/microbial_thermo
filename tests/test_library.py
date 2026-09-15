"""The curated metabolism library.

This doubles as the broadest regression test in the suite: every catalogued
metabolism is balanced, put through the two-path cross-check, and checked
against what is known about where it sits in the redox ladder.
"""

import unittest
import warnings

import microbial_thermo as mt
from microbial_thermo.exceptions import MicrobialThermoError
from microbial_thermo.library import (
    default_library,
    energy_table,
    metabolism,
    reaction,
)


class TestCatalogue(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = default_library()

    def test_catalogue_is_populated(self):
        self.assertGreaterEqual(len(self.library), 25)

    def test_names_are_unique(self):
        names = self.library.names()
        self.assertEqual(len(names), len(set(names)))

    def test_lookup_by_name(self):
        entry = metabolism("anammox")
        self.assertEqual(entry.group, "nitrogen")
        self.assertIn("anammox", entry.label.lower())

    def test_unknown_name_suggests_alternatives(self):
        with self.assertRaises(MicrobialThermoError) as caught:
            metabolism("anamox")
        self.assertIn("anammox", str(caught.exception))

    def test_groups_are_reported(self):
        groups = self.library.groups()
        for expected in ("methanogenesis", "sulfate reduction", "nitrification"):
            self.assertIn(expected, groups)

    def test_in_group_filters(self):
        entries = self.library.in_group("methanogenesis")
        self.assertGreaterEqual(len(entries), 3)
        self.assertTrue(all(e.group == "methanogenesis" for e in entries))


class TestEveryMetabolismBuilds(unittest.TestCase):
    """The corpus check: all of them, balanced and self-consistent."""

    @classmethod
    def setUpClass(cls):
        cls.library = default_library()
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_all_build_and_cross_check(self):
        for entry in self.library:
            with self.subTest(metabolism=entry.name):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    built = entry.build(self.conditions)
                built.verify_consistency()

    def test_all_conserve_atoms_and_charge(self):
        from microbial_thermo.balance import verify_conservation

        for entry in self.library:
            with self.subTest(metabolism=entry.name):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    built = entry.build(self.conditions)
                verify_conservation(built.coefficients, label=entry.name)

    def test_key_element_hint_applies_only_where_it_belongs(self):
        """A carbon hint must not be forced onto an inorganic acceptor.

        Applying it to both couples is what broke the organic-donor entries
        when this library was first written.
        """
        entry = metabolism("acetotrophic_sulfate_reduction")
        self.assertEqual(entry.key_element, "C")
        built = entry.build(self.conditions)
        self.assertEqual(built.donor_half.half.key_element, "C")
        self.assertEqual(built.acceptor_half.half.key_element, "S")

    def test_normalisation_is_honoured(self):
        built = reaction("knallgas", self.conditions, n_electrons=4)
        self.assertEqual(built.n_electrons, 4)


class TestRedoxLadder(unittest.TestCase):
    """The ordering the catalogue must reproduce.

    Redox zonation is the organising fact of sedimentary biogeochemistry: as
    electron acceptors are consumed in order of the energy they yield, the
    sequence oxygen, nitrate, manganese, iron, sulfate, CO2 emerges. If the
    library did not reproduce it, something would be wrong with the
    thermodynamics rather than with the ordering.
    """

    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.table = energy_table(mt.Conditions(temperature_c=25.0, pH=7.0))
        cls.by_name = {row["name"]: row["dG per e- (kJ/mol)"] for _, row in cls.table.iterrows()}

    def test_no_entry_failed_to_evaluate(self):
        failures = self.table[self.table["problem"] != ""]
        self.assertEqual(len(failures), 0, f"failed: {list(failures['name'])}")

    def test_zonation_order(self):
        order = [
            "aerobic_acetate_oxidation",
            "denitrification_acetate",
            "manganese_reduction_acetate",
            "ferrihydrite_reduction_acetate",
            "acetotrophic_sulfate_reduction",
            "acetoclastic_methanogenesis",
        ]
        energies = [self.by_name[name] for name in order]
        pairs = zip(order[:-1], order[1:], energies[:-1], energies[1:], strict=True)
        for earlier, later, a, b in pairs:
            with self.subTest(pair=f"{earlier} before {later}"):
                self.assertLess(a, b, f"{earlier} ({a:.1f}) should outrank {later} ({b:.1f})")

    def test_oxygen_is_the_best_acceptor(self):
        aerobic = self.by_name["aerobic_acetate_oxidation"]
        for name in ("denitrification_acetate", "acetotrophic_sulfate_reduction"):
            self.assertLess(aerobic, self.by_name[name])

    def test_marginal_metabolisms_are_marginal(self):
        """AOM and acetoclastic methanogenesis run near equilibrium, which is
        why they need such low product concentrations in situ."""
        for name in ("anaerobic_methane_oxidation", "acetoclastic_methanogenesis"):
            with self.subTest(metabolism=name):
                self.assertGreater(self.by_name[name], -25.0)
                self.assertLess(self.by_name[name], 0.0)

    def test_crystalline_iron_is_unfavourable_at_standard_state(self):
        """Goethite reduction only pays once Fe(II) is dilute, which is the
        point of the standard-state caveat."""
        self.assertGreater(self.by_name["goethite_reduction_hydrogen"], 0.0)

    def test_table_is_sorted_by_energy(self):
        values = [v for v in self.table["dG per e- (kJ/mol)"] if v == v]
        self.assertEqual(values, sorted(values))

    def test_group_filter(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            frame = energy_table(mt.Conditions(temperature_c=25.0, pH=7.0), group="methanogenesis")
        self.assertTrue(all(g == "methanogenesis" for g in frame["group"]))


class TestFailuresAreVisible(unittest.TestCase):
    def test_unevaluable_entries_are_reported_not_dropped(self):
        """At 60 C the isothermal manganese phases cannot be evaluated; those
        rows must still appear, carrying the reason."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            warm = energy_table(mt.Conditions(temperature_c=60.0, pH=7.0))
            cool = energy_table(mt.Conditions(temperature_c=25.0, pH=7.0))
        self.assertEqual(len(warm), len(cool))


if __name__ == "__main__":
    unittest.main()
