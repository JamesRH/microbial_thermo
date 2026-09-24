"""The layered thermodynamic data sources.

Three HKF-lineage sources are consulted in order -- the primary database, the
GWB log K minerals, then the supplementary files -- and the ordering is the
whole safety property. A later source can only ever *add* a species; it can
never change a number an earlier one already produced.

That matters because SUPCRTBL is not a newer edition of SUPCRT92. It revised
the mineral end-members against Holland & Powell (2011), and the two disagree
by real amounts. Preferring it silently would move published numbers.
"""

import unittest
import warnings

import numpy as np

import microbial_thermo as mt
from microbial_thermo.backends.pygcc_backend import (
    DEFAULT_HKF_DATABASE,
    SUPPLEMENTARY_HKF_DATABASES,
    PygccBackend,
    bundled_database,
)


class TestPrimaryDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_the_default_is_pinned_by_name(self):
        """Not left to pyGCC's own default, so the provenance record can say
        which file was read."""
        self.assertEqual(DEFAULT_HKF_DATABASE, "speq23.dat")
        self.assertEqual(self.backend.database_name, "speq23.dat")

    def test_speq23_is_a_strict_superset_of_speq21(self):
        """The reason the upgrade was free. Measured, not assumed."""
        older = set(PygccBackend(database=bundled_database("speq21.dat")).species_dict)
        newer = set(self.backend.species_dict)
        self.assertEqual(older - newer, set(), "speq23 dropped species speq21 had")
        self.assertEqual(len(newer - older), 3)

    def test_the_upgrade_moved_no_value_we_expose(self):
        """Every species in the registry must give the same formation energy
        under speq21 and speq23, or the upgrade was not free after all."""
        from microbial_thermo.species import default_registry

        older = PygccBackend(database=bundled_database("speq21.dat"))
        for name in sorted({s.backend for s in default_registry().all_species()}):
            with self.subTest(species=name):
                self.assertAlmostEqual(
                    older.delta_Gf(name, 25.0).magnitude,
                    self.backend.delta_Gf(name, 25.0).magnitude,
                    places=9,
                )


class TestSupplementarySources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_supcrtbl_is_declared_with_its_equation_of_state(self):
        """pyGCC defaults to SUPCRT and would misread a Holland & Powell
        record entirely, so the method has to travel with the file."""
        self.assertEqual(SUPPLEMENTARY_HKF_DATABASES["supcrtbl.dat"], "HP11")

    #: Species only supcrtbl carries. Realgar and orpiment are deliberately
    #: not here: they look like obvious candidates and are in fact already in
    #: the GWB file, so they come from there.
    UNIQUE = ("Scorodite", "Arsenopyrite-R", "Barium-As", "As2O5(s)")

    def test_it_supplies_species_nothing_else_has(self):
        for name in self.UNIQUE:
            with self.subTest(species=name):
                self.assertTrue(name not in self.backend.species_dict)
                self.assertTrue(name not in self.backend.minerals)
                self.assertTrue(name in self.backend.supplementary_species)

    def test_it_adds_a_useful_number_of_species(self):
        unique = (
            set(self.backend.supplementary_species)
            - set(self.backend.species_dict)
            - set(self.backend.minerals)
        )
        self.assertGreater(len(unique), 200)

    def test_scorodite_lands_near_the_published_value(self):
        """Ferric arsenate, the phase that controls arsenic solubility in
        oxidised sediments. Published dGf is near -1280 kJ/mol."""
        value = self.backend.delta_Gf("Scorodite", 25.0).to("kJ/mol").magnitude
        self.assertAlmostEqual(value, -1287.1, delta=15.0)

    def test_the_unit_is_not_double_converted(self):
        """A SUPCRTBL record stores kJ where SUPCRT92 stores calories, but
        heatcap converts internally and returns calories either way.
        Converting again lands 4.184x out on a number that still looks
        plausible -- which is exactly what happened while writing this.
        """
        value = self.backend.delta_Gf("Scorodite", 25.0).to("kJ/mol").magnitude
        self.assertGreater(value, -2000.0)
        self.assertLess(value, -500.0)

    def test_it_never_overrides_an_earlier_source(self):
        """The ordering guarantee. Every species supcrtbl shares with the
        primary database or the GWB minerals must still come from those."""
        shared = set(self.backend.supplementary_species) & (
            set(self.backend.species_dict) | set(self.backend.minerals)
        )
        self.assertGreater(len(shared), 100, "expected a substantial overlap")
        for name in sorted(shared):
            with self.subTest(species=name):
                self.assertNotIn("supcrtbl", self.backend.resolve(name).source)

    def test_the_overlap_disagreements_are_real_and_worth_knowing(self):
        """SUPCRTBL revised these against Holland & Powell (2011). Asserting
        the disagreement exists stops anyone assuming the sources are
        interchangeable."""
        from pygcc import heatcap

        for name, minimum in (("Hematite", 1.0), ("Magnetite", 2.0)):
            with self.subTest(species=name):
                entry, method, _ = self.backend.supplementary_species[name]
                raw = heatcap(T=25.0, P=1.0, Species_ppt=entry, Species=name, method=method).dG
                revised = float(np.asarray(raw).ravel()[0])
                current = self.backend._compute_gibbs_cal(name, 25.0, 1.0)
                gap_kj = abs(revised - current) * 4.184 / 1000.0
                self.assertGreater(gap_kj, minimum)

    def test_everything_available_still_resolves(self):
        unresolvable = []
        for name in self.backend.available_species():
            try:
                self.backend.resolve(name)
            except Exception:
                unresolvable.append(name)
        self.assertEqual(unresolvable, [])


if __name__ == "__main__":
    unittest.main()
