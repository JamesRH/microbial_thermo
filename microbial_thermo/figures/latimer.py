"""Latimer diagrams: an element's redox ladder as a chain of potentials.

A Latimer diagram is the most compact statement of an element's redox
chemistry there is -- the oxidation states in a row, joined by arrows
carrying the potential of each one-step couple:

    HAsO4(2-)  --+0.60-->  As(OH)3  --+0.25-->  As

Nothing in Python draws these, which is the reason this module exists. Every
version a student meets is a static image drawn at pH 0 or pH 14, from one
table, at 25 C. This one is computed from the same data as everything else in
the library, at whatever pH and temperature are asked for, and that is the
whole point: the arsenic ladder at pH 7 is not the arsenic ladder at pH 0,
and being able to see the difference is worth more than the picture.

**How the potentials come out.** Each species is decomposed onto the common
basis (see :mod:`microbial_thermo.figures.basis`) giving its free energy per
mole of the element, ``g``, and its oxidation state ``z``. For two forms A
(oxidized) and B (reduced), the couple ``A + n e- -> B`` has ``n = z_A - z_B``
and

    E = (g_A - g_B) / (n F)

which is a two-line derivation from the formation reactions and needs no
table of standard potentials. It agrees with the published values where they
exist: the As(V)/As(III) step comes out at +0.574 V against a tabulated
+0.560, and As(III)/As(0) at +0.248 against +0.240.

**Which species sits at each oxidation state** depends on pH, and the diagram
picks the one that actually predominates -- the lowest free energy among the
forms of that state at the working pH. Arsenic(V) is H3AsO4 at pH 0, H2AsO4-
at pH 4 and HAsO4(2-) at pH 7, and the potentials differ accordingly. The
rejected forms are kept in ``alternatives`` rather than thrown away.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from .basis import FARADAY_KJ, element_series, ladder_entries, pretty
from .style import PALETTE, SIZES

#: Oxidation states within this of each other count as the same rung. Magnetite
#: is +8/3 and nothing else is near it, so the tolerance only ever collapses
#: float noise.
STATE_TOLERANCE = 1e-6


def _coefficient(value: float) -> str:
    """Format a per-element stoichiometric coefficient.

    These are fractions more often than not -- magnetite contributes 1/3 of
    itself per mole of iron -- and ``0.3333`` in a half reaction looks like a
    rounding error rather than a third.
    """
    fraction = Fraction(value).limit_denominator(12)
    if abs(float(fraction) - value) > 1e-6:
        return f"{value:.3g}"
    if fraction == 1:
        return ""
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


@dataclass(frozen=True)
class LatimerStep:
    """One arrow: the couple between two adjacent oxidation states."""

    oxidized: object  # SpeciesEnergy
    reduced: object  # SpeciesEnergy
    n_electrons: float
    potential: float  # volts vs SHE, at the diagram's pH and temperature
    delta_g: float  # kJ per mole of the element, for the step as written

    @property
    def label(self) -> str:
        return f"{self.potential:+.3f}"

    def half_reaction(self) -> str:
        """The step as a balanced half reaction, per mole of the element.

        Obtained by subtracting the two formation reactions, so the water and
        proton counts are whatever the decomposition already worked out:

            A + (b_B - b_A) H2O + (c_B - c_A) H+ + n e-  ->  B

        Coefficients are per mole of the element and are therefore fractional
        for a mixed-valence phase -- a third of a magnetite per iron.
        """
        left = [self.oxidized.backend]
        right = [self.reduced.backend]

        protons = self.reduced.protons - self.oxidized.protons
        if protons:
            (left if protons > 0 else right).append(f"{_coefficient(abs(protons))} H+".strip())

        # Oxygen is what is left once the protons are placed; recover it from
        # the species rather than carrying b through the dataclass.
        oxygen = (
            self.oxidized.species.parsed.element_count("O") / self.oxidized.n_element
            - self.reduced.species.parsed.element_count("O") / self.reduced.n_element
        )
        if oxygen:
            (right if oxygen > 0 else left).append(f"{_coefficient(abs(oxygen))} H2O".strip())

        left.append(f"{_coefficient(self.n_electrons)} e-".strip())
        return " + ".join(left) + " -> " + " + ".join(right)


@dataclass(frozen=True)
class LatimerDiagram:
    """An element's oxidation states in order, joined by couple potentials."""

    element: str
    rungs: tuple  # SpeciesEnergy, most oxidized first
    steps: tuple  # LatimerStep between consecutive rungs
    pH: float
    temperature_c: float
    activity: float
    series: object
    alternatives: tuple = ()  # (state, backend) not chosen at this pH
    excluded: tuple = ()  # (backend, why) -- no statable oxidation state
    fixed: tuple = ()

    @property
    def states(self) -> list:
        return [rung.oxidation_state for rung in self.rungs]

    @property
    def species(self) -> list:
        return [rung.backend for rung in self.rungs]

    def potential(self, oxidized: str, reduced: str) -> float:
        """The potential of any pair, adjacent or not.

        A skip-step potential is the electron-weighted mean of the steps it
        spans, never the arithmetic mean -- the classic Latimer exercise, and
        the reason this is a method rather than something a reader does by eye.
        """
        by_name = {rung.backend: rung for rung in self.rungs}
        try:
            high, low = by_name[oxidized], by_name[reduced]
        except KeyError as exc:
            raise KeyError(
                f"{exc.args[0]!r} is not on this diagram; it has {', '.join(self.species)}"
            ) from None
        n = high.oxidation_state - low.oxidation_state
        if n == 0:
            raise ValueError(f"{oxidized} and {reduced} are the same oxidation state")
        gap = self.series.energy(high, self.pH) - self.series.energy(low, self.pH)
        return float(gap / (n * FARADAY_KJ))

    @property
    def disproportionating(self) -> list:
        """Rungs unstable with respect to their two neighbours.

        A species disproportionates when the potential bringing it *down* to
        the next state exceeds the potential bringing it *in* from the state
        above -- the electrons it would take are worth more than the electrons
        it would give. Stated on this diagram rather than only on the Frost
        one because the comparison is right there on two adjacent arrows.
        """
        unstable = []
        for i, step_down in enumerate(self.steps):
            if i == 0:
                continue
            step_in = self.steps[i - 1]
            if step_down.potential > step_in.potential:
                unstable.append(self.rungs[i])
        return unstable

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(
            {
                "oxidized": [s.oxidized.backend for s in self.steps],
                "reduced": [s.reduced.backend for s in self.steps],
                "n electrons": [s.n_electrons for s in self.steps],
                "E (V)": [s.potential for s in self.steps],
                "dG (kJ/mol element)": [s.delta_g for s in self.steps],
                "half reaction": [s.half_reaction() for s in self.steps],
            }
        )


def latimer_diagram(
    element: str,
    species=None,
    pH: float | None = None,
    temperature_c: float = 25.0,
    activity: float = 1.0,
    fixed=None,
    backend=None,
    reference=None,
    series=None,
) -> LatimerDiagram:
    """Build the redox ladder of ``element`` at the given conditions.

    ``activity`` defaults to 1, the standard-state convention a published
    Latimer diagram is quoted at. Lowering it to an environmental value is
    legitimate and shifts every step involving a change in dissolved count.
    """
    if series is None:
        pH = 7.0 if pH is None else pH
        series = element_series(
            element,
            species=species,
            temperature_c=temperature_c,
            activity=activity,
            pH=pH,
            fixed=fixed,
            backend=backend,
            reference=reference,
        )
    else:
        # A prepared series carries the conditions its numbers were built at.
        # Taking temperature and activity from the arguments instead would let
        # a figure state conditions it was not drawn at, which is worse than
        # an error. pH is the exception: it is a coefficient, applied here.
        pH = series.pH if pH is None else pH
        temperature_c, activity = series.temperature_c, series.activity

    usable, excluded = ladder_entries(series)
    if not usable:
        raise ValueError(
            f"no species of {element} has an oxidation state the basis can state; "
            + "; ".join(f"{name}: {why}" for name, why in excluded)
        )
    energies = {entry: float(series.energy(entry, pH)) for entry in usable}

    # One species per oxidation state: the one that predominates at this pH.
    by_state: dict[float, list] = {}
    for entry in usable:
        state = round(entry.oxidation_state, 6)
        for existing in by_state:
            if abs(existing - state) < STATE_TOLERANCE:
                state = existing
                break
        by_state.setdefault(state, []).append(entry)

    rungs, alternatives = [], []
    for state in sorted(by_state, reverse=True):
        candidates = sorted(by_state[state], key=lambda e: energies[e])
        rungs.append(candidates[0])
        alternatives.extend((state, other.backend) for other in candidates[1:])

    steps = []
    for high, low in zip(rungs, rungs[1:], strict=False):
        n = high.oxidation_state - low.oxidation_state
        gap = energies[high] - energies[low]
        steps.append(
            LatimerStep(
                oxidized=high,
                reduced=low,
                n_electrons=n,
                potential=gap / (n * FARADAY_KJ),
                delta_g=-gap,
            )
        )

    return LatimerDiagram(
        element=element,
        rungs=tuple(rungs),
        steps=tuple(steps),
        pH=pH,
        temperature_c=temperature_c,
        activity=activity,
        series=series,
        alternatives=tuple(alternatives),
        excluded=tuple(excluded),
        fixed=series.fixed,
    )


def _state_label(state: float) -> str:
    """Oxidation state as a roman-numeral label where it is an integer."""
    numerals = {
        0: "0",
        1: "I",
        2: "II",
        3: "III",
        4: "IV",
        5: "V",
        6: "VI",
        7: "VII",
        8: "VIII",
    }
    rounded = round(state)
    if abs(state - rounded) < STATE_TOLERANCE and abs(rounded) in numerals:
        body = numerals[abs(rounded)]
        return f"-{body}" if rounded < 0 else body
    return _coefficient(state)


def plot_latimer(
    element: str,
    diagram: LatimerDiagram | None = None,
    figsize=None,
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    title: str | None = None,
    show_half_reactions: bool = False,
    ax=None,
    **kwargs,
):
    """Draw the ladder as a labelled chain. Returns ``(figure, diagram)``."""
    import matplotlib.pyplot as plt

    if diagram is None:
        diagram = latimer_diagram(element, **kwargs)

    n = len(diagram.rungs)
    if ax is not None:
        figure = ax.figure
        figsize = None
    elif figsize is None:
        figsize = (
            max(7.0, 2.6 * n),
            3.4 + (0.32 * len(diagram.steps) if show_half_reactions else 0),
        )
    if ax is None:
        figure, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-0.65, n - 0.35)
    ax.set_ylim(-1.0, 1.0)
    ax.axis("off")

    unstable = {rung.backend for rung in diagram.disproportionating}

    for i, rung in enumerate(diagram.rungs):
        edge = PALETTE["endergonic"] if rung.backend in unstable else PALETTE["annotation"]
        ax.text(
            i,
            0.0,
            pretty(rung.backend),
            ha="center",
            va="center",
            fontsize=SIZES["equation"] - 2,
            color=PALETTE["annotation"],
            bbox={
                "boxstyle": "round,pad=0.38",
                "facecolor": "white",
                "edgecolor": edge,
                "linewidth": 1.6 if rung.backend in unstable else 0.9,
            },
            zorder=5,
        )
        ax.text(
            i,
            -0.34,
            _state_label(rung.oxidation_state),
            ha="center",
            va="top",
            fontsize=SIZES["oxidation_state"],
            color=PALETTE["muted"],
        )

    for i, step in enumerate(diagram.steps):
        ax.annotate(
            "",
            xy=(i + 1 - 0.30, 0.0),
            xytext=(i + 0.30, 0.0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": PALETTE["electron"],
                "linewidth": 1.4,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )
        ax.text(
            i + 0.5,
            0.09,
            f"{step.potential:+.3f} V",
            ha="center",
            va="bottom",
            fontsize=SIZES["potential"],
            color=PALETTE["reduction"],
        )
        ax.text(
            i + 0.5,
            -0.06,
            f"{_coefficient(step.n_electrons) or '1'} e$^-$",
            ha="center",
            va="top",
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["muted"],
        )

    # The span: the electron-weighted potential across the whole ladder, which
    # is the number a reader most often wants and most often gets wrong by
    # averaging the arrows.
    if len(diagram.rungs) > 2:
        overall = diagram.potential(diagram.rungs[0].backend, diagram.rungs[-1].backend)
        ax.annotate(
            "",
            xy=(n - 1, -0.62),
            xytext=(0, -0.62),
            arrowprops={
                "arrowstyle": "-|>",
                "color": PALETTE["muted"],
                "linewidth": 1.0,
                "connectionstyle": "arc3,rad=-0.06",
            },
        )
        ax.text(
            (n - 1) / 2,
            -0.78,
            f"overall {overall:+.3f} V",
            ha="center",
            va="top",
            fontsize=SIZES["annotation"],
            color=PALETTE["muted"],
        )

    if unstable:
        ax.text(
            0.0,
            0.62,
            "outlined red: unstable to disproportionation at these conditions",
            ha="left",
            va="bottom",
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["endergonic"],
        )

    if show_half_reactions:
        lines = "\n".join(f"{s.half_reaction()}    E = {s.potential:+.3f} V" for s in diagram.steps)
        figure.text(
            0.5,
            0.02,
            lines,
            ha="center",
            va="bottom",
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["muted"],
            family="monospace",
        )

    if title is None:
        subtitle = (
            f"pH {diagram.pH:g} · {diagram.temperature_c:g} °C · "
            f"dissolved activity {diagram.activity:g}"
        )
        if diagram.fixed:
            held = ", ".join(f"{pretty(name)} {value:g}" for _, name, value in diagram.fixed)
            subtitle += f"\nheld fixed: {held}"
        title = f"{element} — Latimer diagram\n{subtitle}"
    ax.set_title(title, fontsize=SIZES["title"], pad=8)
    if len(figure.axes) == 1:
        figure.tight_layout()

    if save is not None:
        base = Path(save)
        if base.parent != Path(""):
            base.parent.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            figure.savefig(
                base.with_suffix(f".{suffix}"), format=suffix, dpi=dpi, bbox_inches="tight"
            )
    return figure, diagram
