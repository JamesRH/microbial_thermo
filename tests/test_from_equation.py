"""Building a reaction from a written equation, and the supplemental table.

``Reaction.from_equation`` reads an unbalanced equation, works out which
elements change oxidation state, forms the donor and acceptor couples from
that, and balances the rest. The awkward part is that the partner of an
elemental reactant is usually left implicit -- nobody writes the proton that
H2 oxidises to -- so it has to be supplied.
"""

import unittest
import warnings

import microbial_thermo as mt
from microbial_thermo.balance import infer_couples
from microbial_thermo.exceptions import AmbiguousReactionError, OutOfRangeError
from microbial_thermo.supplemental import supplemental_species, unverified_names


def coefficient_map(reaction):
    return {s.backend: v for s, v in reaction.coefficients.items()}


class TestCoupleInference(unittest.TestCase):
    def test_implicit_partner_is_supplied(self):
        """In 'NO3- + H2 -> NH2OH' the proton is never written, but H2 has to
        oxidise to something."""
        donor, acceptor = infer_couples("NO3- + H2 -> NH2OH")
        self.assertEqual(donor, ("H2(aq)", "H+"))
        self.assertEqual(acceptor, ("NH2OH(aq)", "NO3-"))

    def test_hydrogens_of_a_product_are_not_mistaken_for_oxidised_h2(self):
        """Hydroxylamine contains hydrogen, but those are not an oxidation
        product of H2 -- they are just hydrogens. Preferring non-auxiliary
        elements is what keeps this straight."""
        _, acceptor = infer_couples("NO3- + H2 -> NH2OH")
        self.assertEqual(acceptor[1], "NO3-")

    def test_oxygen_partner_is_supplied(self):
        donor, acceptor = infer_couples("CH4 + O2 -> CO2 + H2O")
        self.assertEqual(donor, ("Methane(aq)", "CO2(aq)"))
        self.assertEqual(acceptor, ("H2O", "O2(aq)"))

    def test_both_couples_written_out(self):
        donor, acceptor = infer_couples("H2 + SO4-2 -> HS-")
        self.assertEqual(donor, ("H2(aq)", "H+"))
        self.assertEqual(acceptor, ("HS-", "SO4--"))

    def test_one_element_can_supply_both_couples(self):
        """Anammox: ammonium is oxidised and nitrite reduced, both to N2."""
        donor, acceptor = infer_couples("NH4+ + NO2- -> N2")
        self.assertEqual(donor, ("NH4+", "N2(aq)"))
        self.assertEqual(acceptor, ("N2(aq)", "NO2-"))

    def test_disproportionation_is_inferred(self):
        """Thiosulfate is both donor and acceptor; the couples disambiguate a
        reaction that conservation alone leaves underdetermined."""
        donor, acceptor = infer_couples("S2O3-2 + H2O -> SO4-2 + HS-")
        # Both couples are returned as (reduced, oxidized): thiosulfate is the
        # reduced member of the donor couple and the oxidized member of the
        # acceptor couple, which is exactly what disproportionation means.
        self.assertEqual(donor, ("S2O3--", "SO4--"))
        self.assertEqual(acceptor, ("HS-", "S2O3--"))


class TestFromEquation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.7)

    def _build(self, equation, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mt.Reaction.from_equation(equation, conditions=self.conditions, **kwargs)

    def test_the_scratchbook_example(self):
        reaction = self._build("NO3- + H2 -> NH2OH", normalize_to="NO3-")
        self.assertEqual(
            coefficient_map(reaction),
            {"H+": -1, "H2(aq)": -3, "NO3-": -1, "NH2OH(aq)": 1, "H2O": 2},
        )
        self.assertEqual(reaction.n_electrons, 6)
        reaction.verify_consistency()

    def test_normalising_to_the_first_reactant(self):
        reaction = self._build("NO3- + H2 -> NH2OH", normalize_to="NO3-")
        self.assertEqual(coefficient_map(reaction)["NO3-"], -1)

    def test_renormalising_preserves_intensive_quantities(self):
        per_nitrate = self._build("NO3- + H2 -> NH2OH", normalize_to="NO3-")
        per_pair = per_nitrate.renormalized(2)
        self.assertEqual(per_pair.n_electrons, 2)
        self.assertAlmostEqual(
            per_nitrate.delta_G_per_electron.magnitude,
            per_pair.delta_G_per_electron.magnitude,
            places=6,
        )
        self.assertAlmostEqual(
            per_nitrate.delta_E_standard_prime.to("V").magnitude,
            per_pair.delta_E_standard_prime.to("V").magnitude,
            places=9,
        )
        self.assertAlmostEqual(
            per_nitrate.delta_G_standard_prime.magnitude,
            3.0 * per_pair.delta_G_standard_prime.magnitude,
            places=6,
        )

    def test_equations_balance_and_cross_check(self):
        for equation in (
            "NO3- + H2 -> NH2OH",
            "CH4 + O2 -> CO2 + H2O",
            "H2 + SO4-2 -> HS-",
            "NH4+ + NO2- -> N2",
            "Fe+2 + O2 -> Fe+3",
            "NH4+ + O2 -> NO2-",
            "Mn+2 + O2 -> pyrolusite",
        ):
            with self.subTest(equation=equation):
                self._build(equation, normalize_to="integer").verify_consistency()

    def test_charge_markers_are_not_read_as_separators(self):
        """'Fe+3' is one species; only a spaced '+' separates two."""
        reaction = self._build("Fe+2 + O2 -> Fe+3", normalize_to="integer")
        self.assertIn("Fe+++", coefficient_map(reaction))

    def test_unspaced_separator_still_reported(self):
        from microbial_thermo.exceptions import BalancingError

        with self.assertRaises(BalancingError) as caught:
            self._build("methane+O2 -> CO2 + H2O")
        self.assertIn("spaces", str(caught.exception))

    def test_no_redox_pair_is_reported(self):
        with self.assertRaises(AmbiguousReactionError):
            self._build("CO2(aq) + H2O -> HCO3- + H+")


class TestSupplementalTable(unittest.TestCase):
    def test_hydroxylamine_resolves_and_has_a_value(self):
        from microbial_thermo.species import resolve

        species = resolve("hydroxylamine")
        self.assertEqual(species.backend, "NH2OH(aq)")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value = mt.get_backend().delta_Gf("NH2OH(aq)", 25.0).magnitude
        self.assertAlmostEqual(value, -23.5, places=3)

    def test_nitrogen_in_hydroxylamine_is_minus_one(self):
        from fractions import Fraction

        self.assertEqual(mt.mean_oxidation_state("N", "NH3O"), Fraction(-1))

    def test_unverified_values_warn_every_time(self):
        backend = mt.get_backend()
        backend._gibbs_cache.clear()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            backend.delta_Gf("NH2OH(aq)", 25.0)
        messages = [str(w.message) for w in caught]
        self.assertTrue(
            any("not been traced to a primary source" in m for m in messages),
            messages,
        )

    def test_everything_shipped_is_currently_unverified(self):
        """If a value is ever verified, this test should be updated -- it
        exists so the flag cannot quietly drift."""
        self.assertEqual(unverified_names(), ["Glucose(aq)", "NH2OH(aq)", "Pyruvate(aq)"])

    def test_refuses_another_temperature_by_default(self):
        with self.assertRaises(OutOfRangeError):
            mt.get_backend().delta_Gf("NH2OH(aq)", 60.0)

    def test_extrapolates_only_when_asked_and_only_with_an_enthalpy(self):
        from microbial_thermo.backends import PygccBackend

        lenient = PygccBackend(allow_extrapolation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value = lenient.delta_Gf("NH2OH(aq)", 60.0).magnitude
        self.assertNotAlmostEqual(value, -23.5, places=2)

        # Glucose carries no enthalpy, so it cannot be extrapolated at all.
        with self.assertRaises(OutOfRangeError):
            lenient.delta_Gf("Glucose(aq)", 60.0)

    def test_provenance_is_recorded(self):
        for entry in supplemental_species().values():
            with self.subTest(species=entry.backend):
                self.assertTrue(entry.provenance.strip())

    def test_supplemental_species_appear_in_the_catalogue(self):
        available = mt.get_backend().available_species()
        for name in ("NH2OH(aq)", "Glucose(aq)", "Pyruvate(aq)"):
            self.assertIn(name, available)


class TestUnverifiedValuesAreVisibleOnFigures(unittest.TestCase):
    """A rendered figure has to say when it rests on a hand-entered value.

    The console warning is long gone by the time someone looks at a saved
    figure, so SPEC section 2.1 requires the figure itself to carry the notice.
    """

    @classmethod
    def setUpClass(cls):
        import matplotlib

        matplotlib.use("Agg")
        cls.conditions = mt.Conditions(temperature_c=25.0, pH=7.7)

    def _reaction(self, equation):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mt.Reaction.from_equation(equation, conditions=self.conditions)

    def test_footnote_names_the_unverified_species(self):
        from microbial_thermo.figures.style import unverified_footnote

        text = unverified_footnote(self._reaction("NO3- + H2 -> NH2OH"))
        self.assertIsNotNone(text)
        self.assertIn("NH2OH(aq)", text)

    def test_no_footnote_when_everything_is_from_the_database(self):
        from microbial_thermo.figures.style import unverified_footnote

        self.assertIsNone(unverified_footnote(self._reaction("NH4+ + O2 -> N2")))

    def test_the_figure_actually_carries_it(self):
        from microbial_thermo.figures import plot_half_reactions

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            figure = plot_half_reactions(self._reaction("NO3- + H2 -> NH2OH"))
        texts = [t.get_text() for t in figure.texts]
        self.assertTrue(
            any("unverified" in t for t in texts),
            f"no notice found among figure-level texts: {texts}",
        )

    def test_a_clean_figure_is_not_cluttered(self):
        from microbial_thermo.figures import plot_half_reactions

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            figure = plot_half_reactions(self._reaction("NH4+ + O2 -> N2"))
        texts = [t.get_text() for t in figure.texts]
        self.assertFalse(any("unverified" in t for t in texts))


if __name__ == "__main__":
    unittest.main()
