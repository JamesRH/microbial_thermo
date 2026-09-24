"""Our oxidation states against an independent method.

`oxidation.py` assigns oxidation states by partitioning bonds on
electronegativity. That is the right approach for organics, where the
structure is in the SMILES and RDKit can walk it. For **inorganic solids** it
is weakest: there is no structure, only a formula, so the assignment rests on
conventional states for the other elements.

pymatgen's ``Composition.oxi_state_guesses`` comes at the same problem from
statistics over the ICSD -- a genuinely independent route. Where the two agree
on a mineral, that is evidence; where they disagree, one of them is wrong and
it is worth knowing which.

pymatgen is a **test-only** dependency. Nothing at runtime imports it, and
these tests skip if it is absent.

Two things about the comparison, both learned by getting them wrong first:

* ``oxi_state_guesses`` returns *ranked candidates*, not an answer. For pyrite
  its first guess is Fe(III) with S at -1.5; the correct Fe(II)/S(-I) is
  second. So the assertion is that our value appears among the candidates,
  never that it equals the first.
* It rounds to four decimals, so an exact ``Fraction(8, 3)`` will not compare
  equal to ``2.6667`` at machine precision. Four of five apparent
  disagreements were this, not chemistry.
"""

import unittest
import warnings
from fractions import Fraction

try:
    from pymatgen.core import Composition

    HAVE_PYMATGEN = True
except ImportError:  # pragma: no cover - environment without the dev extra
    HAVE_PYMATGEN = False

from microbial_thermo.oxidation import mean_oxidation_state
from microbial_thermo.species import default_registry

#: pymatgen rounds its guesses to four decimals.
TOLERANCE = 1e-3

#: The solids we expose. Restricted deliberately: this cross-check is only
#: meaningful where our method has no structural information to work from.
MINERALS = (
    "Birnessite",
    "Bixbyite",
    "Fe(OH)3",
    "Goethite",
    "Hausmannite",
    "Hematite",
    "Magnetite",
    "Manganite",
    "Pyrolusite",
    "Sulfur(s)",
    "V2O4",
    "V3O5",
)


@unittest.skipUnless(HAVE_PYMATGEN, "pymatgen not installed")
class TestMineralsAgreeWithPymatgen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.registry = default_registry()

    def _formula(self, name):
        for entry in self.registry.all_species():
            if entry.backend == name:
                return entry.formula
        self.fail(f"{name} is no longer in the registry")

    def test_every_mineral_state_is_among_pymatgens_candidates(self):
        checked = 0
        for name in MINERALS:
            formula = self._formula(name)
            candidates = Composition(formula).oxi_state_guesses()
            self.assertTrue(candidates, f"pymatgen had no guess for {formula}")
            elements = {e for e in candidates[0]} - {"O", "H"}
            for element in sorted(elements):
                with self.subTest(mineral=name, element=element):
                    ours = float(mean_oxidation_state(element, formula))
                    options = [c[element] for c in candidates if element in c]
                    nearest = min(options, key=lambda o: abs(ours - o))
                    self.assertAlmostEqual(ours, nearest, delta=TOLERANCE)
                    checked += 1
        self.assertGreaterEqual(checked, 12)

    def test_the_mixed_valence_oxides_come_out_fractional(self):
        """Magnetite and hausmannite are the interesting cases: a single
        integer state is wrong for both, and both methods say so."""
        for name, expected in (("Magnetite", Fraction(8, 3)), ("Hausmannite", Fraction(8, 3))):
            with self.subTest(mineral=name):
                element = "Fe" if name == "Magnetite" else "Mn"
                self.assertEqual(mean_oxidation_state(element, self._formula(name)), expected)

    def test_birnessite_is_not_a_whole_number_either(self):
        value = mean_oxidation_state("Mn", self._formula("Birnessite"))
        self.assertEqual(value, Fraction(7, 2))
        self.assertNotEqual(value.denominator, 1)


@unittest.skipUnless(HAVE_PYMATGEN, "pymatgen not installed")
class TestPymatgenIsWrongForOrganics(unittest.TestCase):
    """The boundary of the cross-check, asserted so nobody widens it.

    ``oxi_state_guesses`` works from ICSD statistics over ionic solids. Given a
    covalent organic it returns assignments that are not chemically sensible --
    for formaldehyde it puts hydrogen at 0 and at -1 -- so extending this
    comparison to the organic species would produce failures that say nothing
    about our code.
    """

    def test_formaldehyde_carbon_is_zero(self):
        """H at +1 and O at -2 forces C to 0. This is not in doubt."""
        self.assertEqual(mean_oxidation_state("C", "CH2O"), Fraction(0))

    def test_and_pymatgen_disagrees_by_assigning_hydrogen_wrongly(self):
        candidates = Composition("CH2O").oxi_state_guesses()
        self.assertTrue(candidates)
        carbon = [c["C"] for c in candidates]
        self.assertNotIn(0.0, carbon)
        hydrogen = [c.get("H") for c in candidates]
        self.assertTrue(
            any(h is not None and h <= 0.0 for h in hydrogen),
            "expected pymatgen to put hydrogen at or below zero, which is the "
            "tell that it is not treating this as a covalent molecule",
        )

    def test_our_organic_states_come_from_structure_instead(self):
        """Acetate's two carbons differ, which a formula alone cannot say and
        pymatgen therefore cannot either."""
        from microbial_thermo.oxidation import per_atom_states

        self.assertEqual(per_atom_states("CC(=O)[O-]", "C"), [-3, 3])


if __name__ == "__main__":
    unittest.main()
