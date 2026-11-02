"""Methanogenic substrates: the Frost diagram, and the energies against the literature.

Two independent checks. First, the Frost diagram and the reaction balancer
share the database and nothing else, so the height of a substrate above the
CH4-HCO3- line must equal the energy of its disproportionation computed as a
balanced reaction. Second, in the gas-phase convention the methanogenic
reactions must reproduce Thauer, Jungermann & Decker (1977, Bacteriol. Rev.
41:100), the table every textbook copies:

    4 H2 + CO2 -> CH4 + 2 H2O                 -131    kJ/mol CH4
    CH3COO- + H+ -> CH4 + CO2                  -36
    CH3OH + H2 -> CH4 + H2O                   -112.5
    2 (CH3)2S + 2 H2O -> 3 CH4 + CO2 + 2 H2S   -49 per CH4

DMS is checked with HS- rather than H2S, which is the form at pH 7; it lands
within a kJ regardless.
"""

import unittest
import warnings

import microbial_thermo as mt
from microbial_thermo import Couple
from microbial_thermo.figures.basis import FARADAY_KJ, _auxiliary_state_matches, element_series
from microbial_thermo.figures.frost import frost_diagram
from microbial_thermo.library import reaction

HELD = {"S": ("HS-", 1.0), "N": ("NH4+", 1.0)}
METHYL = ["Methanol(aq)", "Methanamine(aq)", "methanethiol", "dimethyl sulfide"]
STANDARD = mt.Conditions(temperature_c=25.0, pH=7.0)


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


class TestMethylCompoundsOnTheCarbonLadder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        species = ["HCO3-", "Acetate", "Methane(aq)", *METHYL]
        cls.held = _quiet(
            frost_diagram, "C", species=species, pH=7.0, predominant_only=False, fixed=HELD
        )
        cls.default = _quiet(frost_diagram, "C", species=species, pH=7.0, predominant_only=False)

    def test_against_bisulfide_and_ammonium_they_are_carbon_minus_two(self):
        states = {point.backend: point.oxidation_state for point in self.held.points}
        for name in METHYL:
            self.assertAlmostEqual(states[name], -2.0, places=9, msg=name)

    def test_against_sulfate_the_sulfur_compounds_are_refused(self):
        """S(-II) in DMS is not S(+VI) in sulfate, so the count is not carbon's."""
        excluded = dict(self.default.excluded)
        self.assertIn("dimethyl sulfide", excluded)
        self.assertIn("methanethiol", excluded)

    def test_no_smiles_means_no(self):
        """Minerals have no structure to check, so pyrite stays off the iron ladder."""
        pyrite = mt.resolve("Pyrite")
        self.assertFalse(_auxiliary_state_matches(pyrite, "S", "HS-"))

    def test_iron_still_refuses_pyrite_and_siderite(self):
        series = _quiet(element_series, "Fe", activity=1.0)
        kept = {entry.backend for entry in series if entry.state_is_conventional}
        self.assertNotIn("Pyrite", kept)
        self.assertNotIn("Siderite", kept)

    def test_methane_and_bicarbonate_are_the_hull(self):
        self.assertEqual({point.backend for point in self.held.stable}, {"Methane(aq)", "HCO3-"})


class TestTwoRoutesAgree(unittest.TestCase):
    def test_acetate_height_above_the_chord_is_its_disproportionation(self):
        diagram = _quiet(frost_diagram, "C", species=["HCO3-", "Acetate", "Methane(aq)"], pH=7.0)
        points = {point.backend: point for point in diagram.points}
        low, high, acetate = points["Methane(aq)"], points["HCO3-"], points["Acetate"]
        chord = low.volt_equivalent + (acetate.oxidation_state - low.oxidation_state) / (
            high.oxidation_state - low.oxidation_state
        ) * (high.volt_equivalent - low.volt_equivalent)
        from_frost = -2 * (acetate.volt_equivalent - chord) * FARADAY_KJ

        balanced = _quiet(
            mt.Reaction.from_couples,
            donor=Couple.make("acetate", "HCO3-", key_element="C"),
            acceptor=Couple.make("methane", "HCO3-", key_element="C"),
            conditions=STANDARD,
            normalize_to="methane",
        )
        self.assertAlmostEqual(from_frost, balanced.delta_G_standard_prime.magnitude, delta=0.05)


class TestAgainstThauer1977(unittest.TestCase):
    def _per_ch4(self, donor, acceptor):
        built = _quiet(
            mt.Reaction.from_couples,
            donor=donor,
            acceptor=acceptor,
            conditions=STANDARD,
            normalize_to="CH4(g)",
        )
        return built.delta_G_standard_prime.magnitude

    def test_hydrogenotrophic(self):
        value = self._per_ch4(("H2(g)", "H+"), ("CH4(g)", "CO2(g)"))
        self.assertAlmostEqual(value, -131.0, delta=1.0)

    def test_aceticlastic(self):
        carbon = Couple.make("CH4(g)", "CO2(g)", key_element="C")
        value = self._per_ch4(Couple.make("acetate", "CO2(g)", key_element="C"), carbon)
        self.assertAlmostEqual(value, -36.0, delta=1.0)

    def test_methyl_reduction(self):
        value = self._per_ch4(("H2(g)", "H+"), Couple.make("CH4(g)", "methanol", key_element="C"))
        self.assertAlmostEqual(value, -112.5, delta=1.0)

    def test_dms_disproportionation(self):
        donor = Couple.make("DMS", [["CO2(g)", 2], ["HS-", 1]], key_element="C")
        value = self._per_ch4(donor, Couple.make("CH4(g)", "CO2(g)", key_element="C"))
        self.assertAlmostEqual(value, -49.0, delta=1.0)


class TestTheLibraryEntries(unittest.TestCase):
    def test_methyl_reduction_beats_co2_reduction_per_h2(self):
        """The whole basis of the lower H2 floor: same donor, better acceptor."""
        methyl = _quiet(reaction, "methyl_reducing_methanogenesis", STANDARD)
        carbon = _quiet(reaction, "hydrogenotrophic_methanogenesis", STANDARD)
        gap = carbon.delta_G_standard_prime.magnitude - methyl.delta_G_standard_prime.magnitude
        self.assertAlmostEqual(gap, 65.0, delta=1.0)

    def test_disproportionations_agree_per_methane_however_written(self):
        """Through CO2 or through the CH3OH/CH4 couple: same energy per CH4."""
        direct = _quiet(reaction, "methylotrophic_methanogenesis", STANDARD, normalize_to="methane")
        through_co2 = _quiet(
            mt.Reaction.from_couples,
            donor=Couple.make("methanol", "CO2(aq)", key_element="C"),
            acceptor=Couple.make("methane", "CO2(aq)", key_element="C"),
            conditions=STANDARD,
            normalize_to="methane",
        )
        self.assertAlmostEqual(
            through_co2.delta_G_standard_prime.magnitude,
            direct.delta_G_standard_prime.magnitude,
            places=6,
        )
        self.assertEqual(direct.n_electrons * 4, through_co2.n_electrons)

    def test_methyl_disproportionations_count_the_electrons_the_methyls_move(self):
        """4 CH3X -> 3 CH4 + CO2 moves 6 e-, which is 2 per CH4 (README #30).

        Counted here from the balanced reaction and oxidation states written
        out by hand, not from the library's own electron bookkeeping. Every
        methyl carbon is C(-II). The one that ends as CO2 (C(+IV)) gives up
        6 e-, and each one that ends as CH4 (C(-IV)) takes 2. Written through
        CO2, as these entries were until November 2026, the library counted
        8 per CH4.
        """
        methyl, carbon_dioxide, methane = -2, 4, -4
        for name in (
            "methylotrophic_methanogenesis",
            "methylamine_methanogenesis",
            "dms_methanogenesis",
        ):
            with self.subTest(metabolism=name):
                built = _quiet(reaction, name, STANDARD, normalize_to="methane")
                made = {species.formula: nu for species, nu in built.coefficients.items()}
                given_up = made["CO2"] * (carbon_dioxide - methyl)
                taken = made["CH4"] * (methyl - methane)
                self.assertEqual(given_up, taken)
                self.assertEqual(built.n_electrons, given_up)
                self.assertEqual(built.n_electrons, 2)

    def test_the_rewrite_kept_the_energy_per_methane(self):
        """#30 changed only the electron count. These are the per-CH4 values
        the through-CO2 entries gave before it. Methylamine is not here: #31
        then moved it to methylammonium, which is tested separately."""
        before = {
            "methylotrophic_methanogenesis": -86.6,
            "dms_methanogenesis": -29.3,
        }
        for name, per_ch4 in before.items():
            with self.subTest(metabolism=name):
                built = _quiet(reaction, name, STANDARD, normalize_to="methane")
                self.assertAlmostEqual(built.delta_G_standard_prime.magnitude, per_ch4, delta=0.06)


class TestMethylammonium(unittest.TestCase):
    """README #31: methylamine methanogenesis is written for CH3NH3+."""

    def test_the_substrate_is_methylammonium(self):
        built = _quiet(reaction, "methylamine_methanogenesis", STANDARD, normalize_to="methane")
        formulas = {species.formula for species in built.coefficients}
        self.assertIn("CH6N+", formulas)
        self.assertNotIn("CH5N", formulas)

    def test_it_costs_the_deprotonation_energy_per_methylamine(self):
        """Two paths: the entry as built, and the neutral-amine version plus
        2.303 RT (pKa - 7) for each of the 4/3 methylamines per CH4."""
        ammonium = _quiet(reaction, "methylamine_methanogenesis", STANDARD, normalize_to="methane")
        neutral = _quiet(
            mt.Reaction.from_couples,
            donor=Couple.make("methylamine", ["CO2(aq)", "NH4+"], key_element="C"),
            acceptor=Couple.make(["methane", "NH4+"], "methylamine", key_element="C"),
            conditions=STANDARD,
            normalize_to="methane",
        )
        pka = float(mt.pKa_ladder("methylamine", temperature_c=25.0)[0])
        deprotonation = 2.302585 * 8.314462618e-3 * 298.15 * (pka - 7.0)
        self.assertAlmostEqual(
            ammonium.delta_G_standard_prime.magnitude,
            neutral.delta_G_standard_prime.magnitude + (4 / 3) * deprotonation,
            delta=0.01,
        )
        self.assertAlmostEqual(ammonium.delta_G_standard_prime.magnitude, -58.1, delta=0.1)


if __name__ == "__main__":
    unittest.main()
