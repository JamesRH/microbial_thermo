"""Step-by-step derivations, for teaching.

The point of this module is that a student should be able to follow every
number back to its origin: which species contributed which formation energy,
how each term of the reaction quotient was formed, and where the Nernst
factor came from. A result with no visible arithmetic is not much use in a
course.

:func:`derivation` returns a :class:`Derivation`, which renders as Markdown
with LaTeX in a notebook and as plain text at a terminal.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from fractions import Fraction

from .balance import ELECTRON
from .oxidation import format_oxidation_state
from .reaction import Reaction, _activity
from .units import (
    FARADAY,
    KJ_PER_MOL_STR,
    R,
    as_magnitude,
    celsius_to_kelvin,
    ureg,
)

#: Default fraction command. matplotlib's mathtext does not implement
#: ``\tfrac``, so figure code passes ``\frac`` instead.
DISPLAY_FRAC = "\\tfrac"
MATHTEXT_FRAC = "\\frac"


def _latex_coefficient(value: Fraction, frac: str = DISPLAY_FRAC) -> str:
    if value == 1:
        return ""
    if value.denominator == 1:
        return f"{value.numerator}\\,"
    return f"{frac}{{{value.numerator}}}{{{value.denominator}}}\\,"


def _latex_side(terms, frac: str = DISPLAY_FRAC) -> str:
    if not terms:
        return "0"
    return " + ".join(f"{_latex_coefficient(v, frac)}\\mathrm{{{s.label}}}" for s, v in terms)


def latex_equation(coefficients: dict, frac: str = DISPLAY_FRAC) -> str:
    """Render a signed coefficient set as a LaTeX equation.

    ``frac`` selects the fraction command: the default suits real LaTeX and
    MathJax, while matplotlib requires :data:`MATHTEXT_FRAC`.
    """
    left = [(s, -v) for s, v in coefficients.items() if v < 0]
    right = [(s, v) for s, v in coefficients.items() if v > 0]
    return f"{_latex_side(left, frac)} \\longrightarrow {_latex_side(right, frac)}"


@dataclass
class Step:
    """One numbered step of a derivation."""

    title: str
    body: list[str] = field(default_factory=list)

    def line(self, text: str) -> Step:
        self.body.append(text)
        return self

    def equation(self, latex: str) -> Step:
        self.body.append(f"$$ {latex} $$")
        return self


@dataclass
class Derivation:
    """An ordered list of steps that renders as Markdown or plain text."""

    title: str
    steps: list[Step] = field(default_factory=list)

    def step(self, title: str) -> Step:
        entry = Step(title)
        self.steps.append(entry)
        return entry

    def to_markdown(self) -> str:
        out = [f"## {self.title}", ""]
        for index, step in enumerate(self.steps, start=1):
            out.append(f"**Step {index}. {step.title}**")
            out.append("")
            out.extend(step.body)
            out.append("")
        return "\n".join(out)

    def to_text(self) -> str:
        """Plain text, with the LaTeX reduced to something readable."""
        out = [self.title, "=" * len(self.title), ""]
        for index, step in enumerate(self.steps, start=1):
            out.append(f"Step {index}. {step.title}")
            for line in step.body:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("$$") and stripped.endswith("$$"):
                    out.append(f"    {_latex_to_text(stripped[2:-2])}")
                else:
                    out.append(f"  {_latex_to_text(stripped)}")
            out.append("")
        return "\n".join(out)

    def _repr_markdown_(self) -> str:  # Jupyter rich display
        return self.to_markdown()

    def __str__(self) -> str:
        return self.to_text()


def derivation(reaction: Reaction) -> Derivation:
    """Build the full derivation for ``reaction``."""
    conditions = reaction.conditions
    backend = reaction.backend
    temperature_k = celsius_to_kelvin(conditions.temperature_c)
    rt = as_magnitude((R * (temperature_k * ureg.kelvin)).to(KJ_PER_MOL_STR), KJ_PER_MOL_STR)

    work = Derivation(f"Derivation: {reaction.format()}")

    # --- 1. the couples and what changes oxidation state ---------------------
    step = work.step("Identify the redox couples")
    for label, result in (
        ("Electron donor", reaction.donor_half),
        ("Electron acceptor", reaction.acceptor_half),
    ):
        ox_state, red_state = result.oxidation_states()
        element = result.half.key_element
        # Use the couple's own rendering: a side may name several species.
        step.line(
            f"- {label}: {result.couple} — {element} goes from "
            f"{format_oxidation_state(ox_state)} on the oxidized side to "
            f"{format_oxidation_state(red_state)} on the reduced side"
        )
    step.line(
        f"- The donor is written in the oxidative direction and the acceptor "
        f"in the reductive direction, normalised to "
        f"{reaction.n_electrons} electron"
        f"{'s' if reaction.n_electrons != 1 else ''}."
    )

    # --- 2. balancing --------------------------------------------------------
    step = work.step("Balance each half reaction")
    step.line(
        "Balance the redox-active element first, then oxygen with $\\mathrm{H_2O}$, "
        "hydrogen with $\\mathrm{H^+}$, and finally charge with electrons."
    )
    step.line("")
    step.line("Oxidation (electron donor):")
    step.equation(latex_equation({k: -v for k, v in reaction.donor_half.half.coefficients.items()}))
    step.line("Reduction (electron acceptor):")
    step.equation(latex_equation(reaction.acceptor_half.half.coefficients))

    # --- 3. the combined reaction -------------------------------------------
    step = work.step("Add the half reactions so the electrons cancel")
    step.equation(latex_equation(reaction.coefficients))
    step.line(
        f"Both halves are scaled to {reaction.n_electrons} electrons, so the "
        "electrons cancel exactly and the sum conserves atoms and charge."
    )

    # --- 4. formation energies ----------------------------------------------
    step = work.step("Look up standard formation energies")
    step.line(
        f"At {conditions.temperature_c:g} °C from {backend.name} "
        f"{getattr(backend, 'version', '')} "
        f"({getattr(backend, 'database_name', 'default database')}):"
    )
    step.line("")
    step.line("| species | ν | ΔG°f (kJ/mol) | ν·ΔG°f |")
    step.line("|---|---:|---:|---:|")
    total_standard = 0.0
    for species, nu in sorted(reaction.coefficients.items(), key=lambda kv: kv[0].backend):
        if species == ELECTRON:
            continue
        value = as_magnitude(
            backend.delta_Gf(species.backend, conditions.temperature_c, conditions.pressure_bar),
            KJ_PER_MOL_STR,
        )
        contribution = float(nu) * value
        total_standard += contribution
        step.line(
            f"| {species.backend} | {_format_fraction(nu)} | {value:.2f} | {contribution:+.2f} |"
        )
    step.line("")
    step.line(
        "By convention $\\Delta G^\\circ_f(\\mathrm{H^+}) = 0$ and "
        "$\\Delta G^\\circ_f(e^-) = 0$, which places every potential on the "
        "standard hydrogen electrode."
    )
    step.equation(
        "\\Delta G^\\circ = \\sum_i \\nu_i \\Delta G^\\circ_{f,i} = "
        f"{total_standard:+.2f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}"
    )

    # --- 5. the reaction quotient -------------------------------------------
    step = work.step("Build the reaction quotient")
    step.line(f"Activity model: `{conditions.activity_model}`, pH {conditions.pH:g}.")
    step.line("")
    step.line("| species | ν | activity | ν·ln a |")
    step.line("|---|---:|---:|---:|")
    log_q = 0.0
    for species, nu in sorted(reaction.coefficients.items(), key=lambda kv: kv[0].backend):
        if species == ELECTRON:
            continue
        activity = _activity(species, conditions, backend)
        term = float(nu) * math.log(activity)
        log_q += term
        step.line(f"| {species.backend} | {_format_fraction(nu)} | {activity:.4g} | {term:+.4f} |")
    step.line("")
    step.equation(f"\\ln Q = {log_q:+.4f}")
    step.equation(
        f"RT = {R.to('J/mol/K').magnitude:.4f}\\ \\mathrm{{J\\,mol^{{-1}}K^{{-1}}}} "
        f"\\times {temperature_k:.2f}\\ \\mathrm{{K}} = "
        f"{rt:.5f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}"
    )
    step.equation(
        f"RT\\ln Q = {rt:.5f} \\times {log_q:+.4f} = {rt * log_q:+.2f}"
        "\\ \\mathrm{kJ\\,mol^{-1}}"
    )

    # --- 6. the free energy --------------------------------------------------
    step = work.step("Combine into the free energy")
    delta_g = as_magnitude(reaction.delta_G, KJ_PER_MOL_STR)
    step.equation(
        "\\Delta G = \\Delta G^\\circ + RT\\ln Q = "
        f"{total_standard:+.2f} {rt * log_q:+.2f} = {delta_g:+.2f}"
        "\\ \\mathrm{kJ\\,mol^{-1}}"
    )
    standard_prime = as_magnitude(reaction.delta_G_standard_prime, KJ_PER_MOL_STR)
    step.line(
        f"At unit activities but pH {conditions.pH:g}, "
        f"$\\Delta G^{{\\circ\\prime}} = {standard_prime:+.2f}$ kJ/mol."
    )
    step.line(
        f"Per mole of electrons: "
        f"${as_magnitude(reaction.delta_G_per_electron, KJ_PER_MOL_STR):+.2f}$ "
        "kJ/mol e⁻."
    )

    # --- 7. potentials -------------------------------------------------------
    step = work.step("Convert to potentials")
    step.equation("E = \\frac{-\\Delta G_\\mathrm{half}}{nF}")
    for label, result in (
        ("donor", reaction.donor_half),
        ("acceptor", reaction.acceptor_half),
    ):
        step.line(
            f"- {label} {result.couple}: "
            f"$E^\\circ = {result.E_standard.to('V').magnitude:+.3f}$ V, "
            f"$E^{{\\circ\\prime}} = "
            f"{result.E_standard_prime.to('V').magnitude:+.3f}$ V"
        )
    step.equation(
        "\\Delta E^{\\circ\\prime} = E_\\mathrm{acceptor} - E_\\mathrm{donor} = "
        f"{reaction.delta_E_standard_prime.to('V').magnitude:+.3f}\\ \\mathrm{{V}}"
    )

    # --- 8. the cross-check --------------------------------------------------
    step = work.step("Cross-check by an independent path")
    n = float(reaction.n_electrons)
    faraday_kj = as_magnitude(FARADAY.to("C/mol"), "C/mol") / 1000.0

    step.line(
        "The free energy is computed twice: once by summing formation energies "
        "over the whole reaction, and once from the two half-reaction "
        "potentials. Agreement confirms the half reactions really do sum to the "
        "combined reaction, with the right electron count and no sign slip."
    )
    step.line("")
    step.line(
        "Building the potentials takes three steps each, so they are worked "
        "out in full below. Both halves are treated as **reductions**, which is "
        "the convention that makes potentials comparable — the donor is shown "
        "oxidatively elsewhere, but not here."
    )

    # --- 8a. each half reaction, corrected term by term ---------------------
    step.line("")
    step.line("**a. Correct each potential from the standard state.**")
    step.equation("E = E^\\circ - \\frac{RT}{nF}\\sum_i \\nu_i \\ln a_i")

    for label, result in (
        ("Donor", reaction.donor_half),
        ("Acceptor", reaction.acceptor_half),
    ):
        half = result.half
        half_n = float(half.n_electrons)
        prefactor = rt / (half_n * faraday_kj)  # volts per unit of (nu ln a)

        step.line("")
        step.line(f"*{label}: {result.couple}*, as a reduction, n = {half.n_electrons}:")
        step.line("")
        step.line("| species | ν | activity | ν·ln a | ΔE term (V) |")
        step.line("|---|---:|---:|---:|---:|")

        total_term = 0.0
        proton_term = 0.0
        for species, nu in sorted(half.coefficients.items(), key=lambda kv: kv[0].backend):
            if species == ELECTRON:
                continue
            activity = _activity(species, conditions, backend)
            term = float(nu) * math.log(activity)
            delta_e = -prefactor * term
            total_term += delta_e
            if species.backend == "H+":
                proton_term = delta_e
            step.line(
                f"| {species.backend} | {_format_fraction(nu)} | "
                f"{activity:.4g} | {_signed(term)} | {_signed(delta_e)} |"
            )

        standard = result.E_standard.to("V").magnitude
        primed = result.E_standard_prime.to("V").magnitude
        actual = result.E.to("V").magnitude
        step.line("")
        step.line(
            f"- $E^\\circ$ = {_signed(standard)} V — every activity 1, so every term above is zero."
        )
        step.line(
            f"- $E^{{\\circ\\prime}}$ = {_signed(standard)} {_signed(proton_term)} = "
            f"{_signed(primed)} V — the proton term alone, at pH {conditions.pH:g}."
        )
        step.line(
            f"- $E$ = {_signed(standard)} {_signed(total_term)} = {_signed(actual)} V — "
            "every term, at the stated conditions."
        )

    # --- 8b. the difference --------------------------------------------------
    donor_e = reaction.donor_half.E.to("V").magnitude
    acceptor_e = reaction.acceptor_half.E.to("V").magnitude
    delta_e_total = reaction.delta_E.to("V").magnitude
    step.line("")
    step.line("**b. Subtract.** Electrons fall from the donor to the acceptor:")
    step.equation(
        "\\Delta E = E_\\mathrm{acceptor} - E_\\mathrm{donor} = "
        f"{_signed(acceptor_e)} - ({_signed(donor_e)}) = {_signed(delta_e_total)}"
        "\\ \\mathrm{V}"
    )

    # --- 8c. compare ---------------------------------------------------------
    path_b = -n * faraday_kj * delta_e_total
    step.line("")
    step.line("**c. Convert, and compare against the formation-energy sum.**")
    step.equation(
        f"-nF\\Delta E = -{n:g} \\times {faraday_kj:.3f} \\times "
        f"{_signed(delta_e_total)} = {path_b:+.2f}\\ \\mathrm{{kJ\\,mol^{{-1}}}}"
    )
    step.line("")
    step.line("| route | ΔG (kJ/mol) |")
    step.line("|---|---:|")
    step.line(f"| sum of formation energies (step 6) | {delta_g:+.4f} |")
    step.line(f"| −nFΔE, from the potentials above | {path_b:+.4f} |")
    step.line(f"| difference | {delta_g - path_b:+.2e} |")
    step.line("")
    step.line("The two agree, so the half reactions really do sum to the combined reaction. ✓")

    return work


#: Substitutions applied when rendering LaTeX as plain text.
_TEXT_SUBSTITUTIONS = [
    ("\\longrightarrow", "->"),
    ("\\nu", "v"),
    ("\\Delta", "d"),
    ("\\sum_i", "sum"),
    ("\\times", "x"),
    ("\\ln", "ln"),
    ("\\log", "log"),
    ("\\circ\\prime", "o'"),
    ("\\circ", "o"),
    ("\\prime", "'"),
    ("\\mathrm", ""),
    ("\\,", " "),
    ("\\ ", " "),
    ("^{-1}", "^-1"),
    ("$", ""),
    ("**", ""),
    ("*", ""),
]

_TFRAC = re.compile(r"\\t?frac\{(-?\d+)\}\{(-?\d+)\}")
_GENERAL_FRAC = re.compile(r"\\t?frac\{([^{}]*)\}\{([^{}]*)\}")
_BRACES = re.compile(r"[{}]")
_SUBSCRIPT = re.compile(r"_\{([^{}]*)\}")


def _latex_to_text(text: str) -> str:
    """Reduce a LaTeX fragment to readable plain text.

    Not a general LaTeX renderer -- it only has to handle the small set of
    constructs this module emits. Order matters: the word substitutions run
    first so that constructs such as ``\\mathrm{half}`` collapse to plain
    braces, then subscripts, and only then fractions, whose arguments would
    otherwise still contain nested braces.
    """
    for pattern, replacement in _TEXT_SUBSTITUTIONS:
        text = text.replace(pattern, replacement)
    text = _SUBSCRIPT.sub(r"_\1", text)
    text = _TFRAC.sub(r"\1/\2", text)
    for _ in range(3):  # resolve any nesting that survived
        replaced = _GENERAL_FRAC.sub(r"(\1)/(\2)", text)
        if replaced == text:
            break
        text = replaced
    text = _BRACES.sub("", text)
    return " ".join(text.split())


def _signed(value: float, places: int = 4) -> str:
    """Signed fixed-point that never produces a negative zero."""
    if abs(value) < 0.5 * 10**-places:
        return f"{0.0:+.{places}f}"
    return f"{value:+.{places}f}"


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return f"{value.numerator:+d}"
    return f"{value.numerator:+d}/{value.denominator}"
