"""The derivation emitter: every number must be traceable."""

import re
import unittest

import microbial_thermo as mt


class TestDerivation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        conditions = mt.Conditions(
            temperature_c=25.0,
            pH=7.0,
            partial_pressures={"H2(g)": 1e-4},
            concentrations={"CO2(aq)": 1e-3, "methane": 1e-5},
            activity_model="ideal",
        )
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"),
            acceptor=("methane", "CO2(aq)"),
            conditions=conditions,
        )
        cls.work = cls.reaction.show_work()

    def test_has_all_steps(self):
        self.assertEqual(len(self.work.steps), 8)

    def test_markdown_contains_latex_and_headings(self):
        markdown = self.work.to_markdown()
        self.assertIn("## Derivation", markdown)
        self.assertIn("$$", markdown)
        self.assertIn("Step 1.", markdown)

    def test_plain_text_has_no_latex_commands(self):
        """Terminal output must not leak backslash commands."""
        text = self.work.to_text()
        leaked = sorted(set(re.findall(r"\\[a-zA-Z]+", text)))
        self.assertEqual(leaked, [], f"LaTeX leaked into plain text: {leaked}")

    def test_plain_text_has_no_stray_delimiters(self):
        text = self.work.to_text()
        self.assertNotIn("$$", text)

    def test_reports_the_computed_free_energy(self):
        """The dG quoted in the derivation must match the reaction's own value."""
        expected = f"{self.reaction.delta_G.magnitude:+.2f}"
        self.assertIn(expected, self.work.to_text())

    def test_lists_every_species_with_its_formation_energy(self):
        text = self.work.to_text()
        for species in self.reaction.coefficients:
            with self.subTest(species=species.backend):
                self.assertIn(species.backend, text)

    def test_cross_check_difference_is_reported(self):
        self.assertIn("difference", self.work.to_text())

    def test_notebook_repr_is_markdown(self):
        self.assertEqual(self.work._repr_markdown_(), self.work.to_markdown())


class TestCrossCheckStep(unittest.TestCase):
    """Step 8 works the Nernst correction out term by term.

    The tables are only useful if the numbers printed in them actually add up
    to the potentials printed underneath, so that is what is asserted here --
    parsed back out of the rendered text rather than recomputed, which would
    only prove the formula twice.
    """

    @classmethod
    def setUpClass(cls):
        conditions = mt.Conditions(
            temperature_c=25.0,
            pH=7.0,
            partial_pressures={"H2(g)": 1e-4},
            concentrations={"SO4-2": 2.8e-2, "HS-": 1e-6},
            activity_model="ideal",
        )
        cls.reaction = mt.Reaction.from_couples(
            donor=("H2(g)", "H+"), acceptor=("HS-", "SO4-2"), conditions=conditions
        )
        text = cls.reaction.show_work().to_text()
        cls.step8 = text[text.index("Step 8.") :]

    def _blocks(self):
        """Split step 8 into its two half-reaction sections."""
        parts = re.split(r"\n\s*(?:Donor|Acceptor): ", self.step8)
        return parts[1:]

    def test_both_half_reactions_are_tabulated(self):
        self.assertEqual(len(self._blocks()), 2)
        self.assertIn("ΔE term (V)", self.step8)

    def test_terms_sum_to_the_quoted_potential(self):
        """E0 plus every tabulated term must equal the quoted E."""
        for block in self._blocks():
            with self.subTest(block=block.splitlines()[0]):
                terms = [
                    float(row.split("|")[-2])
                    for row in block.splitlines()
                    if row.count("|") >= 5 and "species" not in row and "---" not in row
                ]
                self.assertTrue(terms)
                standard = float(re.search(r"E\^o = ([+-][\d.]+) V", block).group(1))
                actual = float(
                    re.search(r"E = [+-][\d.]+ ([+-][\d.]+) = ([+-][\d.]+) V", block).group(2)
                )
                self.assertAlmostEqual(standard + sum(terms), actual, places=3)

    def test_the_proton_term_alone_gives_the_primed_potential(self):
        for block in self._blocks():
            with self.subTest(block=block.splitlines()[0]):
                proton_row = next(r for r in block.splitlines() if r.strip().startswith("| H+ |"))
                proton_term = float(proton_row.split("|")[-2])
                standard = float(re.search(r"E\^o = ([+-][\d.]+) V", block).group(1))
                primed = float(
                    re.search(r"E\^o' = [+-][\d.]+ [+-][\d.]+ = ([+-][\d.]+) V", block).group(1)
                )
                self.assertAlmostEqual(standard + proton_term, primed, places=3)

    def test_quoted_primed_potentials_match_the_api(self):
        """The narrative must not drift from what the objects report."""
        for result in (self.reaction.donor_half, self.reaction.acceptor_half):
            value = result.E_standard_prime.to("V").magnitude
            self.assertIn(f"{value:+.4f}", self.step8)

    def test_the_difference_is_taken_in_the_right_order(self):
        donor = self.reaction.donor_half.E.to("V").magnitude
        acceptor = self.reaction.acceptor_half.E.to("V").magnitude
        self.assertIn(f"{acceptor:+.4f} - ({donor:+.4f})", self.step8)

    def test_both_routes_are_reported_and_agree(self):
        self.assertIn("sum of formation energies", self.step8)
        self.assertIn("−nFΔE", self.step8)
        # Scope to the comparison table; the Nernst tables above also have
        # signed numbers between pipes.
        table = self.step8[self.step8.index("| route |") :]
        quoted = re.findall(r"\| ([+-][\d.]+) \|", table)
        self.assertGreaterEqual(len(quoted), 2)
        self.assertAlmostEqual(float(quoted[0]), float(quoted[1]), places=4)

    def test_no_negative_zero_anywhere(self):
        self.assertNotIn("-0.0000", self.step8)


if __name__ == "__main__":
    unittest.main()
