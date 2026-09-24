"""The wider mineral and trace-metal set.

Phases a sediment actually contains, so a reaction can be written against the
solid that is really there. What is tested is mostly that they *work* --
resolve, balance, and give sensible energies -- plus the two places where
adding them exposed something.
"""

import unittest
import warnings
from fractions import Fraction

import microbial_thermo as mt
from microbial_thermo.balance import verify_conservation
from microbial_thermo.library import reaction
from microbial_thermo.reaction import Couple, Reaction
from microbial_thermo.species import default_registry, resolve

CONDITIONS = mt.Conditions(temperature_c=25.0, pH=7.0)

#: Every mineral added, with a bracket its formation energy must fall inside.
#: Loose brackets on purpose: this catches a unit slip or a wrong phase, not a
#: few kJ/mol of database revision.
MINERALS = {
    "Pyrrhotite": (-150.0, -50.0),
    "Troilite": (-150.0, -50.0),
    "Calcite": (-1200.0, -1050.0),
    "Aragonite": (-1200.0, -1050.0),
    "Dolomite": (-2250.0, -2100.0),
    "Magnesite": (-1100.0, -950.0),
    "Gypsum": (-1900.0, -1700.0),
    "Anhydrite": (-1400.0, -1250.0),
    "Barite": (-1450.0, -1300.0),
    "Jarosite": (-3400.0, -3200.0),
    "Scorodite": (-1350.0, -1200.0),
    "Orpiment": (-250.0, -100.0),
    "Claudetite": (-650.0, -500.0),
    "Sphalerite": (-250.0, -150.0),
    "Galena": (-150.0, -50.0),
    "Chalcopyrite": (-250.0, -120.0),
    "Covellite": (-100.0, -10.0),
    "Millerite": (-120.0, -30.0),
    "Cinnabar": (-90.0, -10.0),
    "Hydroxyapatite": (-6500.0, -6200.0),
    "Vivianite": (-4500.0, -4250.0),
    "Uraninite": (-1100.0, -950.0),
}


class TestTheMineralsResolveAndEvaluate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_each_is_in_the_registry(self):
        known = {s.backend for s in default_registry().all_species()}
        for name in MINERALS:
            with self.subTest(mineral=name):
                self.assertIn(name, known)

    def test_each_evaluates_into_a_sensible_bracket(self):
        for name, (low, high) in MINERALS.items():
            with self.subTest(mineral=name):
                value = self.backend.delta_Gf(name, 25.0).to("kJ/mol").magnitude
                self.assertGreater(value, low)
                self.assertLess(value, high)

    def test_the_hydrates_were_written_out_not_left_with_a_colon(self):
        """The database writes gypsum CaSO4:2H2O, which is not a formula."""
        gypsum = resolve("gypsum")
        counts = gypsum.parsed.elements
        self.assertEqual(counts["Ca"], 1)
        self.assertEqual(counts["S"], 1)
        self.assertEqual(counts["O"], 6)
        self.assertEqual(counts["H"], 4)

    def test_vivianite_carries_its_eight_waters(self):
        counts = resolve("vivianite").parsed.elements
        self.assertEqual(counts["Fe"], 3)
        self.assertEqual(counts["P"], 2)
        self.assertEqual(counts["O"], 16)
        self.assertEqual(counts["H"], 16)

    def test_each_is_uncharged(self):
        for name in MINERALS:
            with self.subTest(mineral=name):
                self.assertEqual(resolve(name).charge, 0)


class TestPolymorphsHadToBeDeclared(unittest.TestCase):
    """Adding these created two genuinely ambiguous formulas, which is what
    the ``canonical`` machinery exists for."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.registry = default_registry()

    def test_nothing_is_left_undeclared(self):
        self.assertEqual(self.registry.undeclared_ambiguities, {})

    def test_a_bare_caco3_is_calcite_not_aragonite(self):
        """Calcite is the polymorph stable at Earth-surface conditions. They
        differ by under a kJ/mol, which is why leaving this to file order
        would have been hard to notice."""
        self.assertEqual(resolve("CaCO3").backend, "Calcite")
        self.assertEqual(resolve("aragonite").backend, "Aragonite")

    def test_a_bare_fes_is_troilite_not_pyrrhotite(self):
        """Troilite *is* stoichiometric FeS; pyrrhotite is Fe(1-x)S and is
        only written FeS by the database's convention."""
        self.assertEqual(resolve("FeS").backend, "Troilite")
        self.assertEqual(resolve("pyrrhotite").backend, "Pyrrhotite")


class TestMineralsNeedTheirCounterIons(unittest.TestCase):
    """A metal sulfide is unusable without somewhere for the metal to go, and
    this library will not invent it. That refusal is the feature."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_sulfate_to_sphalerite_alone_is_refused(self):
        from microbial_thermo.exceptions import BalancingError

        with self.assertRaises(BalancingError) as caught:
            Reaction.from_couples(
                donor=Couple.make("H2(g)", "H+"),
                acceptor=Couple.make("Sphalerite", "SO4-2", key_element="S"),
                conditions=CONDITIONS,
            )
        self.assertIn("conserve", str(caught.exception))

    def test_naming_the_counter_ion_makes_it_work(self):
        built = Reaction.from_couples(
            donor=Couple.make("Sphalerite", ["SO4-2", "Zn++"], key_element="S"),
            acceptor=Couple.make("H2O", "O2(aq)"),
            conditions=CONDITIONS,
        )
        verify_conservation(built.coefficients)
        self.assertLess(built.delta_G_standard_prime.magnitude, 0.0)

    def test_pyrite_forms_from_sulfate_and_ferrous_iron(self):
        built = Reaction.from_couples(
            donor=Couple.make("H2(g)", "H+"),
            acceptor=Couple.make("Pyrite", [["SO4-2", 2], ["Fe+2", 1]], key_element="S"),
            conditions=CONDITIONS,
        )
        verify_conservation(built.coefficients)


class TestElementalForms(unittest.TestCase):
    """Needed as the origin of a Frost or Latimer diagram."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_they_sit_at_zero_as_an_element_must(self):
        for name in ("As", "Se", "Fe", "Cu", "Zn", "Ni"):
            with self.subTest(element=name):
                value = self.backend.delta_Gf(name, 25.0).to("kJ/mol").magnitude
                self.assertAlmostEqual(value, 0.0, delta=0.01)

    def test_manganese_does_not_and_that_is_recorded(self):
        """Elemental Mn evaluates to about -2.5 kJ/mol through the GWB route,
        which is a database problem rather than rounding. A manganese Frost
        diagram drawn against this species inherits the offset, so it is
        asserted here rather than discovered later."""
        value = self.backend.delta_Gf("Mn", 25.0).to("kJ/mol").magnitude
        self.assertLess(value, -1.0)
        self.assertGreater(value, -5.0)

    def test_they_resolve_by_their_obvious_names(self):
        self.assertEqual(resolve("zero-valent iron").backend, "Fe")
        self.assertEqual(resolve("elemental selenium").backend, "Se")


class TestNewMetalMetabolisms(unittest.TestCase):
    NAMES = (
        "selenate_respiration_hydrogen",
        "selenate_respiration_acetate",
        "selenite_reduction_hydrogen",
        "chromate_reduction_acetate",
        "uranium_reduction_acetate",
    )

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_all_are_catalogued_balanced_and_exergonic(self):
        for name in self.NAMES:
            with self.subTest(metabolism=name):
                built = reaction(name, CONDITIONS)
                verify_conservation(built.coefficients)
                built.verify_consistency()
                self.assertLess(built.delta_G_standard_prime.magnitude, -20.0)

    def test_selenium_reduces_in_two_steps_and_the_first_is_easier(self):
        """Se(VI) to Se(IV) is the more favourable step, which is why selenate
        reduction is the one that supports growth."""
        first = reaction("selenate_respiration_hydrogen", CONDITIONS)
        second = reaction("selenite_reduction_hydrogen", CONDITIONS)
        self.assertLess(
            first.delta_G_standard_prime.magnitude,
            second.delta_G_standard_prime.magnitude,
        )

    def test_selenium_ends_as_an_insoluble_element(self):
        """The contrast with arsenic that makes selenium remediable: the
        reduced product drops out of solution."""
        built = reaction("selenite_reduction_hydrogen", CONDITIONS)
        product = [s for s in built.coefficients if s.backend == "Se"]
        self.assertTrue(product)
        self.assertEqual(product[0].phase, "s")

    def test_uranium_ends_as_a_solid_too(self):
        built = reaction("uranium_reduction_acetate", CONDITIONS)
        product = [s for s in built.coefficients if s.backend == "Uraninite"]
        self.assertTrue(product)
        self.assertEqual(product[0].phase, "s")

    def test_the_new_couples_sit_between_sulfur_and_oxygen_on_the_tower(self):
        from microbial_thermo.tower import reference_couples

        potentials = {e.label: e.potential_v for e in reference_couples(CONDITIONS)}
        for label in ("SeO4-2/SeO3-2", "SeO3-2/Se0", "CrO4-2/Cr3+", "UO2+2/UO2(s)"):
            with self.subTest(couple=label):
                self.assertGreater(potentials[label], potentials["SO4-2/HS-"])
                self.assertLess(potentials[label], potentials["O2/H2O"])


class TestDeclaredInterDatabaseDisagreements(unittest.TestCase):
    """Three of the newly exposed species disagree between the databases.

    Recorded here as well as in KNOWN_DISAGREEMENTS because a user reaching
    for hydroxyapatite deserves to meet the caveat in the place they are
    working, not only in the audit.
    """

    def test_hydroxyapatite_is_the_worst_in_the_library(self):
        """51 kJ/mol. It is a solid-solution mineral whose composition varies
        and the databases have made different choices. Usable, but not for a
        quantitative saturation state without checking both routes."""
        from tests.test_database_consistency import KNOWN_DISAGREEMENTS

        self.assertIn("Hydroxyapatite", KNOWN_DISAGREEMENTS)

    def test_the_chromium_ions_are_declared_too(self):
        from tests.test_database_consistency import KNOWN_DISAGREEMENTS

        self.assertIn("Cr++", KNOWN_DISAGREEMENTS)
        self.assertIn("Cr+++", KNOWN_DISAGREEMENTS)


class TestArsenopyriteIsDeliberatelyAbsent(unittest.TestCase):
    """Three routes, three answers, none of them agreed on.

    GWB gives -50.0 kJ/mol from a single finite log K point; supcrtbl's two
    polymorphs give -109.6 and -125.5. Exposing any of them would be picking
    one at random, so the registry leaves it out. This test records that as a
    decision rather than an oversight.
    """

    def test_it_is_not_in_the_registry(self):
        known = {s.backend for s in default_registry().all_species()}
        for name in ("Arsenopyrite", "Arsenopyrite-R", "Arsenopyrite-B"):
            with self.subTest(name=name):
                self.assertNotIn(name, known)

    def test_the_routes_really_do_disagree(self):
        warnings.simplefilter("ignore")
        backend = mt.get_backend()
        values = [
            backend.delta_Gf(n, 25.0).to("kJ/mol").magnitude
            for n in ("Arsenopyrite", "Arsenopyrite-R", "Arsenopyrite-B")
        ]
        self.assertGreater(max(values) - min(values), 50.0)


class TestOxidationStatesOfTheNewMineralsAreRefusedNotGuessed(unittest.TestCase):
    """Where our oxidation-state method stops, and why that is right.

    ``mean_oxidation_state`` works by assigning conventional states to the
    elements it knows -- hydrogen +1, oxygen -2 -- and solving for the rest.
    A metal sulfide defeats it: sulfur has no single conventional state (-2 in
    sulfide, +6 in sulfate), and neither does a second metal. So for ZnS,
    orpiment and scorodite it **refuses**, and says why.

    That is the correct behaviour for this library, and it is the same gap
    item #27 found from the other direction: our method is weakest on
    inorganic solids because it has no structure to work from. pymatgen,
    which does this from ICSD statistics, handles them -- which is exactly why
    it is the cross-check in tests/test_oxidation_crosscheck.py.
    """

    def _refuses(self, element, formula):
        from microbial_thermo.exceptions import MicrobialThermoError
        from microbial_thermo.oxidation import mean_oxidation_state

        with self.assertRaises(MicrobialThermoError) as caught:
            mean_oxidation_state(element, formula)
        return str(caught.exception)

    def test_metal_sulfides_are_refused(self):
        for formula, element in (("ZnS", "Zn"), ("PbS", "Pb"), ("NiS", "Ni")):
            with self.subTest(formula=formula):
                message = self._refuses(element, formula)
                self.assertIn("no conventional state for S", message)

    def test_orpiment_is_refused_for_the_same_reason(self):
        self.assertIn("no conventional state for S", self._refuses("As", "As2S3"))

    def test_scorodite_is_refused_because_of_its_second_metal(self):
        message = self._refuses("As", "FeAsH4O6")
        self.assertIn("conventional state", message)

    def test_the_refusal_says_what_to_do_instead(self):
        self.assertIn("SMILES", self._refuses("Zn", "ZnS"))

    def test_it_works_whenever_every_other_element_is_conventional(self):
        """Which is the actual rule, and it is narrower than "oxides are
        fine": asking for carbon in calcite works, asking for calcium in the
        same mineral does not, because carbon has no conventional state.
        """
        from microbial_thermo.oxidation import mean_oxidation_state

        for element, formula, expected in (
            ("As", "As2O3", 3),
            ("U", "UO2", 4),
            ("Cr", "CrO4-2", 6),
            ("Se", "SeO3-2", 4),
            ("C", "CaCO3", 4),
            ("S", "CaSO4", 6),
        ):
            with self.subTest(formula=formula, element=element):
                self.assertEqual(mean_oxidation_state(element, formula), Fraction(expected))

    def test_but_calcium_in_the_same_calcite_is_refused(self):
        self.assertIn("no conventional state for C", self._refuses("Ca", "CaCO3"))


if __name__ == "__main__":
    unittest.main()
