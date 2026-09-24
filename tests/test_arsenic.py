"""Arsenic: As(III) and As(V).

Arsenic earns a place in a microbial bioenergetics library because the
As(V)/As(III) step is run in *both* directions by different organisms --
arsenate as a respiratory acceptor, arsenite as a lithotrophic donor -- and
the couple sits close to zero volts at pH 7, which is why both can pay.

The data get an unusually strong check here. SUPCRT-lineage databases carry
two parallel representations of As(III) that differ by one water, and the two
reach this library by *different routes*: As(OH)3(aq) comes through the GWB
log K path while HAsO2(aq) is direct HKF. They should agree exactly. That they
do, to a fifth of a kJ/mol across 0-100 C, is evidence about the data rather
than about the arithmetic.
"""

import math
import unittest
import warnings
from fractions import Fraction
from pathlib import Path

import microbial_thermo as mt
from microbial_thermo.balance import balance_half_reaction, verify_conservation
from microbial_thermo.library import reaction
from microbial_thermo.oxidation import mean_oxidation_state
from microbial_thermo.speciation import fractions, pKa_ladder
from microbial_thermo.species import resolve
from microbial_thermo.tower import reference_couples

#: 2.303 RT/F at 25 C, in volts.
NERNST_25C = 0.0591596


def gibbs(backend, name, temperature_c=25.0):
    return backend.delta_Gf(name, temperature_c).to("kJ/mol").magnitude


def by_backend(family, **kwargs):
    """``fractions`` is keyed by Species objects; key it by backend name."""
    return {s.backend: v for s, v in fractions(family, **kwargs).items()}


class TestNamesResolve(unittest.TestCase):
    """What a user is likely to type."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_the_oxidation_state_names_work(self):
        self.assertEqual(resolve("As(III)").backend, "As(OH)3(aq)")
        self.assertEqual(resolve("As(V)").backend, "HAsO4--")

    def test_the_common_chemical_names_work(self):
        for name, expected in [
            ("arsenite", "As(OH)3(aq)"),
            ("arsenate", "HAsO4--"),
            ("H3AsO3", "As(OH)3(aq)"),
            ("H2AsO4-", "H2AsO4-"),
            ("HAsO2", "HAsO2(aq)"),
            ("AsO4-3", "AsO4---"),
        ]:
            with self.subTest(name=name):
                self.assertEqual(resolve(name).backend, expected)

    def test_arsenic_adds_no_new_ambiguity(self):
        """Every contested name must still have a declared winner."""
        from microbial_thermo.species import default_registry

        self.assertEqual(default_registry().undeclared_ambiguities, {})

    def test_the_parser_does_not_split_the_element_symbol(self):
        """'As' must not be read as a hydrogen-free A plus sulfur, which is
        the failure mode a two-letter symbol invites."""
        self.assertEqual(mean_oxidation_state("As", "AsH3O3"), Fraction(3))
        self.assertEqual(mean_oxidation_state("As", "HAsO4-2"), Fraction(5))


class TestTheTwoRepresentationsAgree(unittest.TestCase):
    """As(OH)3 / H2AsO3- against HAsO2 / AsO2-, which differ by one water.

    These are the same chemistry written two ways, and they arrive by two
    different data paths, so agreement is a genuine cross-check.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.backend = mt.get_backend()

    def test_they_really_do_come_by_different_routes(self):
        """If this stops being true the cross-check below is weaker than it
        looks, so it is worth asserting rather than assuming."""
        self.assertIn("HAsO2(aq)", self.backend.species_dict)
        self.assertNotIn("As(OH)3(aq)", self.backend.species_dict)
        self.assertIn("As(OH)3(aq)", self.backend.minerals)

    def test_the_acids_agree_across_the_whole_range(self):
        for temperature in (0.01, 25.0, 55.0, 85.0, 100.0):
            with self.subTest(temperature_c=temperature):
                direct = gibbs(self.backend, "As(OH)3(aq)", temperature)
                via_water = gibbs(self.backend, "HAsO2(aq)", temperature) + gibbs(
                    self.backend, "H2O", temperature
                )
                self.assertAlmostEqual(direct, via_water, delta=0.25)

    def test_the_conjugate_bases_agree(self):
        for temperature in (5.0, 25.0, 55.0, 85.0):
            with self.subTest(temperature_c=temperature):
                direct = gibbs(self.backend, "H2AsO3-", temperature)
                via_water = gibbs(self.backend, "AsO2-", temperature) + gibbs(
                    self.backend, "H2O", temperature
                )
                self.assertAlmostEqual(direct, via_water, delta=0.1)

    def test_both_give_the_same_pka(self):
        from microbial_thermo.speciation import pKa

        for temperature in (5.0, 25.0, 55.0):
            with self.subTest(temperature_c=temperature):
                arsenous = pKa("As(OH)3(aq)", "H2AsO3-", temperature_c=temperature)
                metarsenous = pKa("HAsO2(aq)", "AsO2-", temperature_c=temperature)
                self.assertAlmostEqual(arsenous, metarsenous, delta=0.05)


class TestSpeciation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_the_arsenate_ladder_is_triprotic(self):
        self.assertEqual(len(pKa_ladder("arsenate", temperature_c=25.0)), 3)

    def test_the_arsenite_ladder_is_monoprotic_in_range(self):
        self.assertEqual(len(pKa_ladder("arsenite", temperature_c=25.0)), 1)

    def test_arsenate_is_a_mixture_at_physiological_ph(self):
        """pKa2 falls at 6.76, so neither form can be called *the* arsenate
        species at pH 7. This is the reason the family exists."""
        at_seven = by_backend("arsenate", pH=7.0, temperature_c=25.0)
        for name in ("H2AsO4-", "HAsO4--"):
            with self.subTest(species=name):
                self.assertGreater(at_seven[name], 0.2)

    def test_the_bare_name_picks_the_majority_species(self):
        """'As(V)' resolves to HAsO4--, which had better be the form actually
        in the majority at pH 7, or the default is misleading."""
        at_seven = by_backend("arsenate", pH=7.0, temperature_c=25.0)
        majority = max(at_seven, key=at_seven.get)
        self.assertEqual(majority, "HAsO4--")
        self.assertEqual(resolve("As(V)").backend, majority)

    def test_arsenite_is_uncharged_at_physiological_ph(self):
        """The reason As(III) is the mobile and more toxic state: at pH 7 it
        is a neutral molecule, so it crosses membranes."""
        at_seven = by_backend("arsenite", pH=7.0, temperature_c=25.0)
        self.assertGreater(at_seven["As(OH)3(aq)"], 0.99)
        self.assertEqual(resolve("As(III)").charge, 0)

    def test_arsenite_does_ionise_in_an_alkaline_lake(self):
        """Mono Lake sits near pH 10, which is above arsenite's pKa."""
        at_ten = by_backend("arsenite", pH=10.0, temperature_c=25.0)
        self.assertGreater(at_ten["H2AsO3-"], 0.5)


class TestPotentials(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.at_ph_7 = mt.Conditions(temperature_c=25.0, pH=7.0)

    def test_the_standard_couple_matches_the_published_value(self):
        """H3AsO4 + 2 H+ + 2 e- -> H3AsO3 + H2O, E0 = +0.560 V.

        The standard-state form is the one anchor here whose convention is
        unambiguous, so it is what the database is held to. Agreement is
        14 mV, about 2.7 kJ/mol over the electron pair.
        """
        half = mt.half_reaction("As(III)", "H3AsO4", self.at_ph_7)
        self.assertEqual(half.half.n_electrons, 2)
        self.assertAlmostEqual(half.E_standard.to("V").magnitude, 0.560, delta=0.02)

    def test_the_couple_is_a_two_electron_step(self):
        half = mt.half_reaction("As(III)", "As(V)", self.at_ph_7)
        self.assertEqual(half.half.n_electrons, 2)

    def test_the_oxidation_states_are_exactly_five_and_three(self):
        half = mt.half_reaction("As(III)", "As(V)", self.at_ph_7)
        oxidised, reduced = half.half.oxidation_states()
        self.assertEqual(oxidised, Fraction(5))
        self.assertEqual(reduced, Fraction(3))

    def test_the_primed_potential_sits_near_zero(self):
        """Which is the whole biological point: with the couple this close to
        zero, a modest donor can reduce As(V) and a modest acceptor can
        oxidise As(III), so both directions support life."""
        half = mt.half_reaction("As(III)", "As(V)", self.at_ph_7)
        volts = half.E_standard_prime.to("V").magnitude
        self.assertGreater(volts, -0.10)
        self.assertLess(volts, 0.10)

    def test_the_primed_potential_depends_on_which_arsenate_is_written(self):
        """A caution worth encoding. The three arsenate forms differ in proton
        count, so E0' at pH 7 differs by more than a hundred millivolts between
        them -- published values disagree largely for this reason."""
        volts = {
            name: mt.half_reaction("As(III)", name, self.at_ph_7).E_standard_prime.to("V").magnitude
            for name in ("H3AsO4", "H2AsO4-", "As(V)")
        }
        self.assertGreater(volts["H3AsO4"] - volts["As(V)"], 0.10)

    def test_the_ph_slope_implies_four_protons(self):
        """HAsO4-- + 4 H+ + 2 e- -> As(OH)3 + H2O. Reading the proton count
        back out of the potentials is a sharp check on the whole path."""
        entries = {
            e.label: e.potential_v
            for e in reference_couples(mt.Conditions(temperature_c=25.0, pH=7.0))
        }
        higher = {
            e.label: e.potential_v
            for e in reference_couples(mt.Conditions(temperature_c=25.0, pH=8.0))
        }
        slope = higher["As(V)/As(III)"] - entries["As(V)/As(III)"]
        protons = slope / -NERNST_25C * 2
        self.assertAlmostEqual(protons, 4.0, places=2)

    def test_it_sits_between_sulfur_and_nitrogen_on_the_tower(self):
        """Its position is what makes it interesting: below the nitrogen and
        oxygen acceptors, above the sulfur and carbon donors."""
        entries = reference_couples(mt.Conditions(temperature_c=25.0, pH=7.0))
        potentials = {e.label: e.potential_v for e in entries}
        self.assertGreater(potentials["As(V)/As(III)"], potentials["SO4-2/HS-"])
        self.assertLess(potentials["As(V)/As(III)"], potentials["NO3-/NO2-"])

    def test_the_tower_entry_is_grouped_as_a_metalloid(self):
        entries = {e.label: e.group for e in reference_couples()}
        self.assertEqual(entries["As(V)/As(III)"], "metalloid")

    def test_the_metalloid_group_has_its_own_colour(self):
        """An unrecognised group falls back to grey, which would silently
        lose the distinction rather than fail."""
        from microbial_thermo.figures.tower import _GROUP_COLOURS

        self.assertIn("metalloid", _GROUP_COLOURS)


class TestMetabolisms(unittest.TestCase):
    """Every catalogued arsenic metabolism, balanced and exergonic."""

    REDUCTION = [
        "arsenate_respiration_hydrogen",
        "arsenate_respiration_acetate",
        "arsenate_respiration_lactate",
    ]
    OXIDATION = ["arsenite_oxidation_oxygen", "arsenite_oxidation_nitrate"]

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

    def _built(self, name):
        return reaction(name, self.conditions)

    def test_all_of_them_are_catalogued(self):
        from microbial_thermo.library import default_library

        names = default_library().names()
        for name in self.REDUCTION + self.OXIDATION:
            with self.subTest(metabolism=name):
                self.assertIn(name, names)

    def test_all_of_them_balance(self):
        for name in self.REDUCTION + self.OXIDATION:
            with self.subTest(metabolism=name):
                verify_conservation(self._built(name).coefficients)

    def test_all_of_them_pay(self):
        """Each is a metabolism someone has isolated an organism for, so each
        must be exergonic under standard biochemical conditions."""
        for name in self.REDUCTION + self.OXIDATION:
            with self.subTest(metabolism=name):
                energy = self._built(name).delta_G_standard_prime.magnitude
                self.assertLess(energy, 0.0)

    def test_all_of_them_clear_the_energy_quantum(self):
        for name in self.REDUCTION + self.OXIDATION:
            with self.subTest(metabolism=name):
                energy = self._built(name).delta_G_standard_prime.magnitude
                self.assertLess(energy, -20.0)

    def test_oxygen_beats_nitrate_as_the_arsenite_acceptor(self):
        with_oxygen = self._built("arsenite_oxidation_oxygen")
        with_nitrate = self._built("arsenite_oxidation_nitrate")
        self.assertLess(
            with_oxygen.delta_G_standard_prime.magnitude,
            with_nitrate.delta_G_standard_prime.magnitude,
        )

    def test_the_lactate_donor_is_an_incomplete_oxidation(self):
        """Lactate stops at acetate, releasing four electrons rather than the
        twelve a complete oxidation would give. Conservation cannot infer
        where the third carbon goes, so the couple names both products."""
        half = balance_half_reaction("Lactate(aq)", ["Acetate", "CO2(aq)"])
        self.assertEqual(half.n_electrons, 4)
        self.assertFalse(half.is_simple)
        verify_conservation(half.coefficients)

    def test_both_directions_pay_at_the_same_ph(self):
        """The point of the whole entry. Arsenate reduction and arsenite
        oxidation are both exergonic at pH 7 -- with different partners --
        which is why the two halves of the arsenic cycle coexist."""
        reducing = self._built("arsenate_respiration_hydrogen")
        oxidising = self._built("arsenite_oxidation_oxygen")
        self.assertLess(reducing.delta_G_standard_prime.magnitude, 0.0)
        self.assertLess(oxidising.delta_G_standard_prime.magnitude, 0.0)

    def test_arsenate_respiration_tracks_the_donor_strength(self):
        """Hydrogen is the stronger donor, so it should yield more per
        electron pair than acetate does."""
        with_hydrogen = self._built("arsenate_respiration_hydrogen")
        with_acetate = self._built("arsenate_respiration_acetate")
        self.assertLess(
            with_hydrogen.delta_G_standard_prime.magnitude,
            with_acetate.delta_G_standard_prime.magnitude,
        )

    def test_the_two_path_cross_check_passes(self):
        """Construction runs verify_consistency(), so this must not raise --
        it is the guard against a redox sign error in the new couple."""
        for name in self.REDUCTION + self.OXIDATION:
            with self.subTest(metabolism=name):
                self._built(name).verify_consistency()


class TestArseniteOxidationUnderMonoLakeConditions(unittest.TestCase):
    """Alkalilimnicola ehrlichii oxidises arsenite with nitrate in Mono Lake,
    which is alkaline and hypersaline. Worth checking the chemistry survives
    conditions that far from the textbook."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_it_still_pays_at_ph_ten(self):
        alkaline = mt.Conditions(temperature_c=15.0, pH=9.8)
        energy = reaction("arsenite_oxidation_nitrate", alkaline).delta_G.magnitude
        self.assertLess(energy, 0.0)

    def test_the_arsenite_is_ionised_there(self):
        at_lake_ph = by_backend("arsenite", pH=9.8, temperature_c=15.0)
        self.assertGreater(at_lake_ph["H2AsO3-"], 0.3)


class TestBiomassPlaceholder(unittest.TestCase):
    """The <CH2O> stand-in for cell carbon.

    It is a modelling convention, not a measured compound, so what is tested
    is mostly that it announces itself as one.
    """

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")

    def test_it_resolves_under_the_obvious_names(self):
        for name in ("biomass", "CH2O", "<CH2O>", "Biomass(aq)"):
            with self.subTest(name=name):
                self.assertEqual(resolve(name).backend, "Biomass(aq)")

    def test_its_carbon_sits_at_oxidation_state_zero(self):
        """Which is the whole reason CH2O is the conventional proxy: real
        biomass averages near zero too."""
        self.assertEqual(mean_oxidation_state("C", "CH2O"), Fraction(0))

    def test_it_is_flagged_unverified(self):
        from microbial_thermo.supplemental import unverified_names

        self.assertIn("Biomass(aq)", unverified_names())

    def test_using_it_warns(self):
        """A placeholder that passed silently would be the dangerous kind.

        The formation energy is memoised per species and temperature, so the
        warning fires on the first use rather than literally every call --
        clear the cache to make that deterministic rather than dependent on
        which test ran first.
        """
        mt.get_backend()._gibbs_cache.clear()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            reaction(
                "arsenite_carbon_fixation",
                mt.Conditions(temperature_c=25.0, pH=7.0),
            )
        messages = " ".join(str(w.message) for w in caught)
        self.assertIn("Biomass(aq)", messages)
        self.assertIn("not been traced to a primary source", messages)

    def test_any_figure_using_it_is_footnoted(self):
        from microbial_thermo.figures.style import unverified_footnote

        built = reaction("arsenite_carbon_fixation", mt.Conditions(temperature_c=25.0, pH=7.0))
        self.assertIn("Biomass(aq)", unverified_footnote(built))

    def test_it_is_consistent_with_the_library_glucose(self):
        """The value is glucose / 6, so the two cannot drift apart silently."""
        from microbial_thermo.supplemental import supplemental_species

        table = supplemental_species()
        glucose = table["Glucose(aq)"].delta_Gf_kJ_mol
        biomass = table["Biomass(aq)"].delta_Gf_kJ_mol
        self.assertAlmostEqual(biomass, glucose / 6.0, delta=0.01)


class TestArseniteCarbonFixation(unittest.TestCase):
    """The anabolic half: As(III) -> As(V) driving CO2 into biomass."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
        cls.built = reaction("arsenite_carbon_fixation", cls.conditions)

    def test_it_balances(self):
        verify_conservation(self.built.coefficients)

    def test_it_passes_the_two_path_cross_check(self):
        self.built.verify_consistency()

    def test_it_does_not_pay(self):
        """Endergonic by design -- this is the cost side of the ledger. A
        negative number here would mean carbon fixation came free."""
        self.assertGreater(self.built.delta_G_standard_prime.magnitude, 0.0)

    def test_the_catabolism_more_than_covers_it_at_ph_7(self):
        """Which is what makes the organism possible."""
        income = reaction("arsenite_oxidation_oxygen", self.conditions).delta_G.magnitude
        cost = self.built.delta_G.magnitude
        self.assertLess(income + cost, 0.0)

    def test_the_margin_is_thin_in_acid(self):
        """At pH 4 a one-to-one budget does not even clear the biological
        energy quantum, which is a real constraint rather than a rounding
        detail."""
        acid = mt.Conditions(temperature_c=25.0, pH=4.0)
        net = (
            reaction("arsenite_carbon_fixation", acid).delta_G.magnitude
            + reaction("arsenite_oxidation_oxygen", acid).delta_G.magnitude
        )
        self.assertLess(net, 0.0)
        self.assertGreater(net, -20.0)

    def test_the_ph_slope_is_two_protons_worth(self):
        """Two protons leave per electron pair, so the cost must fall by
        2 x RT ln(10) = 11.4 kJ/mol per pH unit. Computed independently of
        the reaction's own proton bookkeeping."""
        rt_ln10 = 8.31446261815324e-3 * 298.15 * math.log(10)
        low = reaction(
            "arsenite_carbon_fixation", mt.Conditions(temperature_c=25.0, pH=6.0)
        ).delta_G.magnitude
        high = reaction(
            "arsenite_carbon_fixation", mt.Conditions(temperature_c=25.0, pH=7.0)
        ).delta_G.magnitude
        self.assertAlmostEqual(high - low, -2 * rt_ln10, places=2)

    def test_alkaline_water_is_kinder(self):
        alkaline = mt.Conditions(temperature_c=25.0, pH=9.0)
        self.assertLess(
            reaction("arsenite_carbon_fixation", alkaline).delta_G.magnitude,
            self.built.delta_G.magnitude,
        )


class TestTheExampleScript(unittest.TestCase):
    """examples/arsenite_carbon_fixation.py, which the README points at."""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        import sys

        root = Path(__file__).resolve().parent.parent
        sys.path.insert(0, str(root / "examples"))
        import arsenite_carbon_fixation as module

        cls.module = module

    def test_the_sweep_returns_all_three_curves(self):
        result = self.module.fixation_vs_ph(ph_values=[5.0, 7.0, 9.0])
        for key in ("pH", "anabolic", "catabolic", "net"):
            with self.subTest(key=key):
                self.assertEqual(len(result[key]), 3)

    def test_the_reported_slope_matches_the_prediction(self):
        """The script prints this as a self-check; it had better hold."""
        result = self.module.fixation_vs_ph(ph_values=[4.0, 7.0, 10.0])
        slope = self.module.slope_per_ph_unit(result)
        self.assertAlmostEqual(slope, -2 * self.module.RT_LN10_25C, places=2)

    def test_the_biomass_override_moves_the_answer(self):
        """Half a CH2O per electron pair, so a change in its formation energy
        must move dG by half that."""
        base = self.module.fixation_vs_ph(ph_values=[7.0], include_catabolic=False)
        shifted = self.module.fixation_vs_ph(
            ph_values=[7.0], include_catabolic=False, biomass_dgf=-130.0
        )
        moved = shifted["anabolic"][0] - base["anabolic"][0]
        from microbial_thermo.supplemental import supplemental_species

        default = supplemental_species()["Biomass(aq)"].delta_Gf_kJ_mol
        self.assertAlmostEqual(moved, 0.5 * (-130.0 - default), places=2)

    def test_the_override_is_restored_afterwards(self):
        """It mutates a cached table, so a leak would quietly corrupt every
        later calculation in the process."""
        before = self.module.fixation_vs_ph(ph_values=[7.0], include_catabolic=False)["anabolic"][0]
        self.module.fixation_vs_ph(ph_values=[7.0], include_catabolic=False, biomass_dgf=-130.0)
        after = self.module.fixation_vs_ph(ph_values=[7.0], include_catabolic=False)["anabolic"][0]
        self.assertAlmostEqual(before, after, places=9)

    def test_the_conclusion_survives_the_placeholder(self):
        """The point of the sensitivity option: across any plausible biomass
        energy, fixation stays uphill and the catabolism still covers it."""
        for value in (-120.0, -152.9, -185.0):
            with self.subTest(biomass_dgf=value):
                result = self.module.fixation_vs_ph(ph_values=[7.0], biomass_dgf=value)
                self.assertGreater(result["anabolic"][0], 0.0)
                self.assertLess(result["net"][0], 0.0)

    def test_it_renders_and_saves(self):
        import tempfile

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        result = self.module.fixation_vs_ph(ph_values=[5.0, 7.0, 9.0])
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "fixation"
            figure, _ = self.module.plot_fixation_vs_ph(result=result, save=base)
            try:
                self.assertTrue(base.with_suffix(".svg").exists())
                self.assertTrue(base.with_suffix(".png").exists())
            finally:
                plt.close(figure)


if __name__ == "__main__":
    unittest.main()
