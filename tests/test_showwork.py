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


if __name__ == "__main__":
    unittest.main()
