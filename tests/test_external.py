"""Imported external data, and the ordering that makes it safe.

An external source is admitted only if its standard state already matches --
aqueous, 1 molal, 25 C, dGf(H+) = 0, revised HKF. OBIGT qualifies, and that is
checkable rather than asserted: acetate's parameters are identical in OBIGT and
in speq23 to the last digit.

But matching conventions does not mean matching numbers. These sources
disagree substantially on some species, which is exactly why an import is
consulted *last* and can only ever add.
"""

import re
import unittest
import warnings

import numpy as np

import microbial_thermo as mt
from microbial_thermo.external import (
    external_sources,
    external_species,
    source_named,
)


def normalised(name):
    name = name.lower().replace("(aq)", "").strip()
    name = re.sub(r"-\d$", "", name)
    return name.replace("-", "").replace(" ", "").replace("_", "")


class TestTheManifest(unittest.TestCase):
    def test_obigt_is_declared_with_its_citation(self):
        source = source_named("OBIGT")
        self.assertIn("Dick", source.citation)
        self.assertEqual(source.energy_units, "cal")
        self.assertEqual(source.kind, "hkf_aqueous")

    def test_every_source_records_where_it_came_from(self):
        for source in external_sources():
            with self.subTest(source=source.name):
                self.assertTrue(source.url)
                self.assertTrue(source.citation)
                self.assertTrue(source.retrieved)
                self.assertTrue(source.license)

    def test_the_note_records_what_was_skipped(self):
        """A negative result is part of the provenance: 80 rows are in joules
        and were deliberately not converted."""
        note = source_named("OBIGT").note
        self.assertIn("joule", note)
        self.assertIn("536", note)


class TestImportedSpecies(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.species = external_species()
        cls.backend = mt.get_backend()

    def test_the_expected_number_arrived(self):
        self.assertEqual(len(self.species), 536)

    def test_each_carries_its_own_primary_citation(self):
        """Not just the database's citation -- the paper the number came
        from. That is what distinguishes an import from a transcription."""
        cited = [s for s in self.species.values() if s.citation]
        self.assertGreater(len(cited) / len(self.species), 0.9)

    def test_entries_are_in_pygccs_aqueous_layout(self):
        for name in ("citrate-3", "glutamate", "cysteine"):
            with self.subTest(species=name):
                self.assertEqual(len(self.species[name].entry), 13)

    def test_the_new_metabolites_resolve_and_are_sensible(self):
        for name, low, high in (
            ("citrate-3", -1200.0, -1100.0),
            ("isocitrate-3", -1200.0, -1100.0),
            ("oxaloacetate-2", -850.0, -750.0),
            ("cysteine", -400.0, -300.0),
        ):
            with self.subTest(species=name):
                value = self.backend.delta_Gf(name, 25.0).to("kJ/mol").magnitude
                self.assertGreater(value, low)
                self.assertLess(value, high)

    def test_they_carry_real_temperature_dependence(self):
        """The reason full HKF parameters were imported rather than a single
        25 C value: these have to work at the temperature asked for."""
        values = [
            self.backend.delta_Gf("citrate-3", t).to("kJ/mol").magnitude
            for t in (5.0, 25.0, 55.0, 85.0)
        ]
        self.assertEqual(len(set(round(v, 3) for v in values)), 4)

    def test_provenance_names_the_source_and_the_paper(self):
        record = self.backend.resolve("citrate-3")
        self.assertIn("OBIGT", record.source)
        self.assertIn("LaRowe", record.source)


class TestTheOrderingIsWhatMakesItSafe(unittest.TestCase):
    """The load-bearing property of the whole mechanism."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()
        cls.imported = external_species()

    def test_an_import_never_overrides_the_primary_database(self):
        shared = [n for n in cls_names(self.imported) if n in self.backend.species_dict]
        for name in shared:
            with self.subTest(species=name):
                self.assertNotIn("OBIGT", self.backend.resolve(name).source)

    def test_acetate_is_identical_in_both(self):
        """The agreement that licensed importing OBIGT at all."""
        from pygcc.species_eos import supcrtaq

        entry = list(self.imported["acetate"].entry)
        imported = float(np.asarray(supcrtaq(25.0, 1.0, entry), dtype=float).ravel()[0])
        own = self.backend._compute_gibbs_cal("Acetate", 25.0, 1.0)
        self.assertAlmostEqual(imported, own, places=6)

    def test_but_they_disagree_badly_on_some_amino_acids(self):
        """Recorded because it is the reason for the ordering. Methionine
        differs by about 180 kJ/mol between the two sources. Preferring OBIGT
        would have moved it by that much, silently."""
        from pygcc.species_eos import supcrtaq

        entry = list(self.imported["methionine"].entry)
        imported = float(np.asarray(supcrtaq(25.0, 1.0, entry), dtype=float).ravel()[0])
        own = self.backend._compute_gibbs_cal("Methionine(aq)", 25.0, 1.0)
        gap_kj = abs(imported - own) * 4.184 / 1000.0
        self.assertGreater(gap_kj, 100.0)

    def test_methionine_still_comes_from_the_primary_database(self):
        self.assertNotIn("OBIGT", self.backend.resolve("Methionine(aq)").source)

    def test_most_shared_species_do_agree(self):
        """The disagreements are real but a minority; if this ratio collapses,
        the import is no longer convention-matched and should be re-examined."""
        from pygcc.species_eos import supcrtaq

        primary = {}
        for name in self.backend.species_dict:
            primary.setdefault(normalised(name), name)

        close = total = 0
        for name, imported in self.imported.items():
            own = primary.get(normalised(name))
            if own is None:
                continue
            if imported.formula.replace(" ", "") != str(self.backend.species_dict[own][0]).replace(
                " ", ""
            ):
                continue
            try:
                value = float(
                    np.asarray(supcrtaq(25.0, 1.0, list(imported.entry)), dtype=float).ravel()[0]
                )
                mine = self.backend._compute_gibbs_cal(own, 25.0, 1.0)
            except Exception:
                continue
            total += 1
            if abs(value - mine) * 4.184 / 1000.0 < 0.01:
                close += 1
        self.assertGreater(total, 90)
        self.assertGreater(close / total, 0.7)

    def test_everything_available_still_resolves(self):
        unresolvable = [
            n for n in self.backend.available_species() if not _resolves(self.backend, n)
        ]
        self.assertEqual(unresolvable, [])


def cls_names(mapping):
    return list(mapping)


def _resolves(backend, name):
    try:
        backend.resolve(name)
        return True
    except Exception:
        return False


if __name__ == "__main__":
    unittest.main()
