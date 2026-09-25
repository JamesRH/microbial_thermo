"""Frost-Ebsworth diagrams: free energy against oxidation state.

A Frost diagram plots the *volt equivalent* of each form of an element
against its oxidation state. It is the same information as a Latimer diagram
rearranged so that the eye does the work:

* the **slope** between any two points is the potential of that couple, so a
  steep downward slope from left to right is a strong oxidant;
* a point **above** the line joining its neighbours is unstable and will
  disproportionate, and the vertical distance is the free energy released;
* the **lowest** point is the form the element ends up in.

**The convention, derived rather than looked up.** For a species S of element
E with oxidation state ``z``, write its formation from the element:

    E + b H2O + c H+ + (-z) e-  ->  S

and let ``g`` be that reaction's free energy per mole of E, at the working
pH, temperature and activity, with electrons at the SHE zero. Then

    volt equivalent  V = g / F

and that is the whole definition. The usual textbook statement -- "V is
N times E-standard for the couple X(N)/X(0)" -- is the same number, because
``E = (g_S - g_E)/(N F)`` and ``g_E`` is zero by construction.

Stating it this way settles the question that blocked this module: **negative
oxidation states need no special treatment.** A hydride simply has a positive
electron count in its formation reaction. Ammonium comes out at -0.823 V at
pH 0 against the textbook -0.82, sulfide at -0.289 against -0.288, with no
sign rule applied anywhere.

**What does need care is the reference.** On an Eh-pH diagram the element
reference cancels; here it is the zero of the y axis. Nitrogen drawn against
dissolved N2 rather than N2 gas puts every point 0.094 V low -- which is
exactly how the first version of this was found to be wrong. The reference
used is recorded on the diagram, and ``reference_is_element`` is False if no
elemental form was available, in which case the y axis has no meaning.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .basis import FARADAY_KJ, element_series, ladder_entries, pretty
from .latimer import STATE_TOLERANCE, _state_label
from .style import PALETTE, SIZES


@dataclass(frozen=True)
class FrostPoint:
    """One species, at its oxidation state and volt equivalent."""

    entry: object  # SpeciesEnergy
    gibbs: float  # kJ per mole of the element, relative to the element

    @property
    def backend(self) -> str:
        return self.entry.backend

    @property
    def oxidation_state(self) -> float:
        return self.entry.oxidation_state

    @property
    def volt_equivalent(self) -> float:
        return self.gibbs / FARADAY_KJ


@dataclass(frozen=True)
class Disproportionation:
    """A species that is unstable with respect to two others."""

    species: str
    oxidation_state: float
    into: tuple  # (more oxidized backend, more reduced backend)
    fractions: tuple  # moles of each, per mole of the species
    delta_g: float  # kJ per mole of the species, negative when it proceeds

    def __str__(self) -> str:
        high, low = self.into
        f_high, f_low = self.fractions
        return (
            f"{self.species} -> {f_high:.3g} {high} + {f_low:.3g} {low}  "
            f"({self.delta_g:+.1f} kJ/mol {self.species})"
        )


@dataclass(frozen=True)
class FrostDiagram:
    """Volt equivalents of one element's redox forms at fixed conditions."""

    element: str
    points: tuple  # FrostPoint, ordered by oxidation state
    pH: float
    temperature_c: float
    activity: float
    series: object
    excluded: tuple = ()
    fixed: tuple = ()

    @property
    def reference(self) -> str:
        return self.series.reference

    @property
    def reference_is_element(self) -> bool:
        return self.series.reference_is_element

    @property
    def states(self) -> list:
        return [point.oxidation_state for point in self.points]

    @property
    def volt_equivalents(self) -> list:
        return [point.volt_equivalent for point in self.points]

    @property
    def most_stable(self) -> FrostPoint:
        """The lowest point: where the element ends up, given time.

        Not the same as "what is present" -- that is a Pourbaix question and
        depends on Eh. This is the form with the lowest free energy at the
        stated conditions with no potential applied.
        """
        return min(self.points, key=lambda p: p.volt_equivalent)

    @property
    def stable(self) -> list:
        """Points on the lower convex hull: the forms that survive.

        Everything else is above a line joining two of these and will, given
        a mechanism, disproportionate into them.
        """
        ordered = sorted(self.points, key=lambda p: (p.oxidation_state, p.volt_equivalent))
        hull: list = []
        for point in ordered:
            while len(hull) >= 2 and _turn(hull[-2], hull[-1], point) <= 0:
                hull.pop()
            hull.append(point)
        return hull

    def slope(self, oxidized: str, reduced: str) -> float:
        """The couple potential between two points, which is the slope.

        This is the identity the diagram is built on, so it is worth being
        able to ask for: it must equal the corresponding Latimer arrow.
        """
        by_name = {point.backend: point for point in self.points}
        high, low = by_name[oxidized], by_name[reduced]
        span = high.oxidation_state - low.oxidation_state
        if span == 0:
            raise ValueError(f"{oxidized} and {reduced} are the same oxidation state")
        return (high.volt_equivalent - low.volt_equivalent) / span

    def disproportionation(self) -> list:
        """Every point above the hull, with the reaction it would run.

        The products are the hull points bracketing it, not its immediate
        neighbours: those are where the element actually ends up. The
        adjacent-pair test drawn on the Latimer diagram is the same test
        applied locally, so anything it flags appears here too.
        """
        hull = self.stable
        on_hull = {id(point) for point in hull}
        found = []
        for point in self.points:
            if id(point) in on_hull:
                continue
            segment = _bracketing_segment(hull, point.oxidation_state)
            if segment is None:
                continue
            low, high = segment
            span = high.oxidation_state - low.oxidation_state
            f_high = (point.oxidation_state - low.oxidation_state) / span
            f_low = 1.0 - f_high
            chord = f_high * high.volt_equivalent + f_low * low.volt_equivalent
            found.append(
                Disproportionation(
                    species=point.backend,
                    oxidation_state=point.oxidation_state,
                    into=(high.backend, low.backend),
                    fractions=(f_high, f_low),
                    delta_g=(chord - point.volt_equivalent) * FARADAY_KJ,
                )
            )
        return found

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(
            {
                "species": [p.backend for p in self.points],
                "oxidation state": self.states,
                "volt equivalent (V)": self.volt_equivalents,
                "dG from the element (kJ/mol)": [p.gibbs for p in self.points],
                "on hull": [any(p is q for q in self.stable) for p in self.points],
            }
        )


def _turn(a, b, c) -> float:
    """Cross product of ab x bc; <= 0 means b is not below the line ac."""
    return (b.oxidation_state - a.oxidation_state) * (c.volt_equivalent - a.volt_equivalent) - (
        c.oxidation_state - a.oxidation_state
    ) * (b.volt_equivalent - a.volt_equivalent)


def _bracketing_segment(hull, state):
    for low, high in zip(hull, hull[1:], strict=False):
        if low.oxidation_state <= state <= high.oxidation_state:
            return low, high
    return None


def frost_diagram(
    element: str,
    species=None,
    pH: float = 7.0,
    temperature_c: float = 25.0,
    activity: float = 1.0,
    fixed=None,
    backend=None,
    reference=None,
    series=None,
    predominant_only: bool = True,
) -> FrostDiagram:
    """Volt equivalents for ``element`` at the given conditions.

    ``predominant_only`` keeps one species per oxidation state -- the one
    with the lowest energy at this pH, which is the one that actually exists.
    Setting it False plots every form offered, which is how to see an acid
    dissociation as a vertical spread of points at one state.
    """
    if series is None:
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

    usable, excluded = ladder_entries(series)
    if not usable:
        raise ValueError(
            f"no species of {element} has an oxidation state the basis can state; "
            + "; ".join(f"{name}: {why}" for name, why in excluded)
        )

    points = [FrostPoint(entry=entry, gibbs=float(series.energy(entry, pH))) for entry in usable]

    if predominant_only:
        best: dict[float, FrostPoint] = {}
        for point in points:
            state = round(point.oxidation_state, 6)
            match = next(
                (key for key in best if abs(key - state) < STATE_TOLERANCE),
                state,
            )
            if match not in best or point.gibbs < best[match].gibbs:
                best[match] = point
        points = list(best.values())

    points.sort(key=lambda p: (p.oxidation_state, p.volt_equivalent))
    return FrostDiagram(
        element=element,
        points=tuple(points),
        pH=pH,
        temperature_c=temperature_c,
        activity=activity,
        series=series,
        excluded=tuple(excluded),
        fixed=series.fixed,
    )


def plot_frost(
    element: str,
    diagram: FrostDiagram | None = None,
    figsize=(8.0, 6.0),
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    title: str | None = None,
    annotate_slopes: bool = False,
    label_offset: float = 0.12,
    ax=None,
    **kwargs,
):
    """Draw the diagram. Returns ``(figure, diagram)``."""
    import matplotlib.pyplot as plt

    if diagram is None:
        diagram = frost_diagram(element, **kwargs)

    if ax is None:
        figure, ax = plt.subplots(figsize=figsize)
    else:
        figure = ax.figure

    states = np.array(diagram.states, dtype=float)
    values = np.array(diagram.volt_equivalents, dtype=float)

    hull = diagram.stable
    hull_states = np.array([p.oxidation_state for p in hull], dtype=float)
    hull_values = np.array([p.volt_equivalent for p in hull], dtype=float)
    on_hull = {id(p) for p in hull}

    # The chain through every point, then the hull over the top of it, so an
    # unstable species reads as a detour above a shortcut.
    ax.plot(
        states,
        values,
        color=PALETTE["muted"],
        linewidth=1.0,
        linestyle="--",
        zorder=2,
    )
    ax.plot(
        hull_states,
        hull_values,
        color=PALETTE["reduction"],
        linewidth=2.0,
        zorder=3,
        label="stable forms",
    )

    for point in diagram.points:
        stable = id(point) in on_hull
        ax.plot(
            [point.oxidation_state],
            [point.volt_equivalent],
            marker="o" if stable else "s",
            markersize=7 if stable else 6,
            color=PALETTE["reduction"] if stable else PALETTE["endergonic"],
            zorder=5,
        )
        ax.annotate(
            pretty(point.backend),
            xy=(point.oxidation_state, point.volt_equivalent),
            xytext=(0, 9 if stable else -16),
            textcoords="offset points",
            ha="center",
            fontsize=SIZES["annotation"],
            color=PALETTE["annotation"],
        )

    if annotate_slopes:
        for low, high in zip(hull, hull[1:], strict=False):
            span = high.oxidation_state - low.oxidation_state
            slope = (high.volt_equivalent - low.volt_equivalent) / span
            ax.annotate(
                f"{slope:+.2f} V",
                xy=(
                    (low.oxidation_state + high.oxidation_state) / 2,
                    (low.volt_equivalent + high.volt_equivalent) / 2,
                ),
                xytext=(0, -14),
                textcoords="offset points",
                ha="center",
                fontsize=SIZES["annotation"] - 1,
                color=PALETTE["reduction"],
            )

    best = diagram.most_stable
    ax.annotate(
        f"most stable: {pretty(best.backend)}",
        xy=(best.oxidation_state, best.volt_equivalent),
        xytext=(0, -34),
        textcoords="offset points",
        ha="center",
        fontsize=SIZES["annotation"],
        color=PALETTE["exergonic"],
        arrowprops={"arrowstyle": "-", "color": PALETTE["exergonic"], "linewidth": 0.8},
    )

    ax.axhline(0.0, color=PALETTE["annotation"], linewidth=0.8, zorder=1)
    ax.axvline(0.0, color=PALETTE["grid"], linewidth=0.8, zorder=1)
    ax.grid(color=PALETTE["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    ax.set_xlabel(f"oxidation state of {element}", fontsize=SIZES["potential"])
    ax.set_ylabel("volt equivalent  $N\\,E°$  (V)", fontsize=SIZES["potential"])
    ax.set_xticks(sorted({round(s, 3) for s in diagram.states}))
    ax.set_xticklabels([_state_label(s) for s in sorted({round(s, 3) for s in diagram.states})])

    span = float(values.max() - values.min()) or 1.0
    ax.set_ylim(values.min() - label_offset * span - 0.25, values.max() + label_offset * span + 0.2)

    if title is None:
        subtitle = (
            f"pH {diagram.pH:g} · {diagram.temperature_c:g} °C · "
            f"dissolved activity {diagram.activity:g} · zero at {pretty(diagram.reference)}"
        )
        if diagram.fixed:
            held = ", ".join(f"{pretty(name)} {value:g}" for _, name, value in diagram.fixed)
            subtitle += f"\nheld fixed: {held}"
        title = f"{element} — Frost-Ebsworth diagram\n{subtitle}"
    ax.set_title(title, fontsize=SIZES["title"], pad=10)

    if not diagram.reference_is_element:
        figure.text(
            0.5,
            0.005,
            f"No elemental form of {element} was available, so the zero of the "
            f"y axis is arbitrary ({diagram.reference}). Slopes are still valid; "
            "the heights are not.",
            ha="center",
            va="bottom",
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["endergonic"],
            wrap=True,
        )

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
