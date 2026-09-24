"""Where a number came from, and whether an export is reproducible.

Two concerns that turn out to be the same concern. A result is only worth
quoting if you can say which database it came from and which version of the
software read it -- and an exported figure is only worth committing if
re-exporting it from unchanged inputs gives you the same bytes.
"""

import hashlib
import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import microbial_thermo as mt
from microbial_thermo.exceptions import SpeciesNotFoundError


class TestResolveCoversEveryRoute(unittest.TestCase):
    """``resolve`` must describe anything ``delta_Gf`` can evaluate.

    It used to check only the HKF dictionary, so it raised for Goethite and
    As(OH)3(aq) -- species that work perfectly well in reactions -- and then
    suggested the name just passed, because the suggester searched all three
    routes while the lookup searched one.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_the_hkf_route(self):
        record = self.backend.resolve("Fe++")
        self.assertIn("speq23", record.source)
        self.assertTrue(record.verified)

    def test_the_water_route(self):
        self.assertIn("IAPWS", self.backend.resolve("H2O").source)

    def test_the_mineral_log_k_route(self):
        """Goethite is the case that used to raise."""
        record = self.backend.resolve("Goethite")
        self.assertIn("log K", record.source)
        self.assertIn("thermo.com.dat", record.source)

    def test_an_aqueous_species_that_arrives_by_the_log_k_route(self):
        """As(OH)3(aq) is aqueous but reaches us through the GWB file, which
        is exactly the combination that made the old behaviour confusing."""
        self.assertIn("log K", self.backend.resolve("As(OH)3(aq)").source)

    def test_the_supplemental_route_reports_its_provenance(self):
        record = self.backend.resolve("Glucose(aq)")
        self.assertIn("hand-entered", record.source)
        self.assertIn("Amend", record.source)  # traced to its primary source
        self.assertTrue(record.verified)

    def test_a_still_unverified_entry_says_so(self):
        self.assertFalse(self.backend.resolve("NH2OH(aq)").verified)

    def test_the_placeholder_is_marked_unverified(self):
        self.assertFalse(self.backend.resolve("Biomass(aq)").verified)

    def test_everything_available_resolves(self):
        """The invariant that was broken: available_species() and resolve()
        must agree about what exists."""
        unresolvable = []
        for name in self.backend.available_species():
            try:
                self.backend.resolve(name)
            except SpeciesNotFoundError:
                unresolvable.append(name)
        self.assertEqual(unresolvable, [], f"{len(unresolvable)} names disagree")

    def test_the_supplementary_hkf_route(self):
        """Scorodite reaches us only through supcrtbl.dat, on the Holland &
        Powell equation of state rather than Maier-Kelley."""
        record = self.backend.resolve("Scorodite")
        self.assertIn("supcrtbl", record.source)
        self.assertIn("HP11", record.source)

    def test_a_genuinely_unknown_name_still_raises(self):
        with self.assertRaises(SpeciesNotFoundError):
            self.backend.resolve("Nonsense(aq)")

    def test_the_suggestion_is_not_the_name_itself(self):
        """The old symptom: 'did you mean: Goethite?' when you asked for
        Goethite. A suggestion identical to the query means the lookup and the
        suggester disagree about scope."""
        with self.assertRaises(SpeciesNotFoundError) as caught:
            self.backend.resolve("Nonsense(aq)")
        self.assertNotIn("Nonsense(aq)", str(caught.exception).split("did you mean:")[1])

    def test_provenance_summary_is_one_line(self):
        from microbial_thermo.supplemental import supplemental_species

        for name, entry in supplemental_species().items():
            with self.subTest(species=name):
                summary = entry.provenance_summary
                self.assertNotIn("\n", summary)
                self.assertLess(len(summary), len(entry.provenance) + 1)


class TestTracedSupplementalValues(unittest.TestCase):
    """The two entries traced to a primary source.

    Both came from the OBIGT database of CHNOSZ, which stores calories. The
    conversion is pinned here so the values cannot drift away from the source
    they claim -- which is the whole point of having verified them.
    """

    #: OBIGT's own numbers, cal/mol, and the citation each carries.
    OBIGT = {
        "Glucose(aq)": (-218100, -301600, "Amend"),
        "Pyruvate(aq)": (-113600, -137400, "Canovas"),
    }
    CAL_TO_J = 4.184

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        from microbial_thermo.supplemental import supplemental_species

        cls.table = supplemental_species()

    def test_the_gibbs_energies_match_obigt(self):
        for name, (gibbs_cal, _, _) in self.OBIGT.items():
            with self.subTest(species=name):
                self.assertAlmostEqual(
                    self.table[name].delta_Gf_kJ_mol,
                    gibbs_cal * self.CAL_TO_J / 1000.0,
                    places=3,
                )

    def test_the_enthalpies_match_obigt(self):
        for name, (_, enthalpy_cal, _) in self.OBIGT.items():
            with self.subTest(species=name):
                self.assertAlmostEqual(
                    self.table[name].delta_Hf_kJ_mol,
                    enthalpy_cal * self.CAL_TO_J / 1000.0,
                    places=3,
                )

    def test_each_names_its_primary_source(self):
        for name, (_, _, author) in self.OBIGT.items():
            with self.subTest(species=name):
                self.assertIn(author, self.table[name].provenance)
                self.assertIn("OBIGT", self.table[name].provenance)

    def test_they_are_marked_verified(self):
        for name in self.OBIGT:
            with self.subTest(species=name):
                self.assertTrue(self.table[name].verified)

    def test_having_an_enthalpy_now_allows_extrapolation(self):
        """Tracing them gained an enthalpy as well as a citation, which is
        what makes a van 't Hoff correction possible at all."""
        for name in self.OBIGT:
            with self.subTest(species=name):
                self.assertIsNotNone(self.table[name].delta_Hf_kJ_mol)

    def test_hydroxylamine_records_where_it_was_looked_for(self):
        """A negative result is worth writing down: the next person should
        not repeat the search."""
        provenance = self.table["NH2OH(aq)"].provenance
        self.assertIn("OBIGT", provenance)
        self.assertIn("Wagman", provenance)
        self.assertFalse(self.table["NH2OH(aq)"].verified)


class TestExportsAreReproducible(unittest.TestCase):
    """Re-exporting an unchanged figure must give unchanged bytes.

    Plotly stamps a random div id unless told otherwise, which left the two
    tracked HTML exports showing a one-line diff after every test run -- churn
    that hides a real content change when one happens.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.conditions = mt.Conditions(
            temperature_c=25.0,
            pH=7.0,
            activity_model="ideal",
            partial_pressures={"H2(g)": 1e-4},
        )
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("methane", "CO2(aq)"),
            conditions=cls.conditions,
        )

    @staticmethod
    def _digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def test_the_explorer_export_is_byte_stable(self):
        from microbial_thermo.figures import plot_energy_explorer
        from microbial_thermo.sweep import ph_axis

        axes = [ph_axis(4.0, 9.0, 4)]
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "a"
            second = Path(directory) / "b"
            plot_energy_explorer(self.reaction, axes=axes, save_html=first)
            plot_energy_explorer(self.reaction, axes=axes, save_html=second)
            self.assertEqual(
                self._digest(first.with_suffix(".html")),
                self._digest(second.with_suffix(".html")),
            )

    def test_the_tower_export_is_byte_stable(self):
        from microbial_thermo.figures import plot_interactive_tower
        from microbial_thermo.tower import tower_grid

        grid = tower_grid(ph_values=(7.0,), temperature_values=(25.0,))
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "a"
            second = Path(directory) / "b"
            plot_interactive_tower(grid, save_html=first)
            plot_interactive_tower(grid, save_html=second)
            self.assertEqual(
                self._digest(first.with_suffix(".html")),
                self._digest(second.with_suffix(".html")),
            )

    def test_the_div_ids_are_pinned_not_random(self):
        from microbial_thermo.figures import plot_energy_explorer
        from microbial_thermo.figures.explorer import EXPLORER_DIV_ID
        from microbial_thermo.sweep import ph_axis

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "explorer"
            plot_energy_explorer(self.reaction, axes=[ph_axis(4.0, 9.0, 4)], save_html=path)
            self.assertIn(EXPLORER_DIV_ID, path.with_suffix(".html").read_text())

    def test_a_real_change_still_changes_the_bytes(self):
        """Pinning the id must not make every export look identical."""
        from microbial_thermo.figures import plot_energy_explorer
        from microbial_thermo.sweep import ph_axis

        with tempfile.TemporaryDirectory() as directory:
            narrow = Path(directory) / "narrow"
            wide = Path(directory) / "wide"
            plot_energy_explorer(self.reaction, axes=[ph_axis(4.0, 9.0, 4)], save_html=narrow)
            plot_energy_explorer(self.reaction, axes=[ph_axis(4.0, 9.0, 8)], save_html=wide)
            self.assertNotEqual(
                self._digest(narrow.with_suffix(".html")),
                self._digest(wide.with_suffix(".html")),
            )


if __name__ == "__main__":
    unittest.main()
