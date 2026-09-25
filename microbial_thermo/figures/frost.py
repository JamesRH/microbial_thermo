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

        **Computed over the predominant forms, one per oxidation state, not
        over every point.** Only one form of a state exists at a given pH, so
        a hull across all of them mixes two questions and gives a wrong
        answer to both: with every arsenate included it returned H3AsO4 as
        "stable" at pH 7 and drew a vertical hull segment down its own
        oxidation state, which is an acid dissociation drawn as a redox step.
        On the usual diagram, where there is one form per state already, this
        changes nothing.
        """
        ordered = sorted(self.predominant, key=lambda p: (p.oxidation_state, p.volt_equivalent))
        hull: list = []
        for point in ordered:
            while len(hull) >= 2 and _turn(hull[-2], hull[-1], point) <= 0:
                hull.pop()
            hull.append(point)
        return hull

    @property
    def predominant(self) -> list:
        """One point per oxidation state: the form that actually exists.

        The lowest free energy among the forms of that state at the working
        pH, which is the same rule the Latimer ladder uses to choose a rung.
        On a diagram built with ``predominant_only=True`` this is every point;
        on one built with it False it is the subset worth drawing a chain
        through.

        Distinct from :attr:`stable`, and the two are easy to confuse.
        *Predominant* is a comparison within one oxidation state -- which
        arsenate, of the four -- and is decided by pH. *Stable* is a
        comparison across states -- whether that arsenate survives at all --
        and is decided by the convex hull.
        """
        best: dict[float, FrostPoint] = {}
        for point in self.points:
            state = round(point.oxidation_state, 6)
            match = next(
                (key for key in best if abs(key - state) < STATE_TOLERANCE),
                state,
            )
            if match not in best or point.gibbs < best[match].gibbs:
                best[match] = point
        return sorted(best.values(), key=lambda p: p.oxidation_state)

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

        Only the **predominant** form of each state is considered. A minority
        acid form is not disproportionating -- it is simply not the form that
        dominates at this pH, and reporting "H3AsO4 -> HAsO4(2-)" as a
        disproportionation would be calling an acid dissociation a redox
        reaction.
        """
        hull = self.stable
        on_hull = {id(point) for point in hull}
        found = []
        for point in self.predominant:
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


def _selected(diagram: FrostDiagram, which: str) -> list:
    """``diagram.points`` or ``diagram.predominant``, by name."""
    if which in ("all", None):
        return list(diagram.points)
    if which in ("predominant", "prominent"):
        return list(diagram.predominant)
    raise ValueError(f"{which!r} is not a selection; use 'all' or 'predominant'")


def _draw_labels(ax, shown, labelled, values) -> None:
    """Name the points, nudging labels apart where forms pile up.

    Species of one oxidation state sit at the same x and are separated only
    by their energies, which for an acid and its conjugate base can be almost
    nothing -- H2AsO4- and HAsO4(2-) differ by 7 mV at pH 7, and their labels
    land squarely on top of each other. Within such a column the labels are
    staggered vertically by :func:`~microbial_thermo.figures.tower._stagger`,
    moved to the right so they clear the hull line, and given a leader back
    to the marker.

    A lone point at its oxidation state keeps the centred label above it, so
    the ordinary one-form-per-state diagram looks exactly as it did.
    """
    from .tower import _stagger

    span = float(values.max() - values.min()) or 1.0
    columns: dict[float, list] = {}
    for point in shown:
        columns.setdefault(round(point.oxidation_state, 6), []).append(point)

    spread = max(columns) - min(columns) if len(columns) > 1 else 1.0
    pad = spread * 0.04

    for state, members in columns.items():
        members = sorted(members, key=lambda p: p.volt_equivalent)
        if not any(id(point) in labelled for point in members):
            continue

        if len(members) == 1:
            point = members[0]
            ax.annotate(
                pretty(point.backend),
                xy=(state, point.volt_equivalent),
                xytext=(0, 9),
                textcoords="offset points",
                ha="center",
                fontsize=SIZES["annotation"],
                color=PALETTE["annotation"],
            )
            continue

        heights = _stagger([p.volt_equivalent for p in members], span * 0.08)
        for point, height in zip(members, heights, strict=True):
            if id(point) not in labelled:
                continue
            ax.annotate(
                pretty(point.backend),
                xy=(state, point.volt_equivalent),
                xytext=(state + pad, height),
                textcoords="data",
                ha="left",
                va="center",
                fontsize=SIZES["annotation"],
                color=PALETTE["annotation"],
                arrowprops={
                    "arrowstyle": "-",
                    "color": PALETTE["muted"],
                    "linewidth": 0.6,
                    "shrinkA": 1.0,
                    "shrinkB": 2.0,
                },
                zorder=6,
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
    pH: float | None = None,
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
    show: str = "all",
    label: str = "predominant",
    ax=None,
    **kwargs,
):
    """Draw the diagram. Returns ``(figure, diagram)``.

    ``show`` and ``label`` each take ``"all"`` or ``"predominant"``, and
    decide what is *drawn* from whatever the diagram holds. They matter only
    for a diagram built with ``predominant_only=False``, which keeps every
    form of every oxidation state -- all four arsenates, graphite and acetate
    both at C(0) -- and can therefore stack several points at one x.

    * ``show="predominant"`` draws only the form that exists at this pH. The
      same picture as ``predominant_only=True``, without rebuilding.
    * ``label="all"`` names every point drawn. Honest and crowded.
    * ``label="predominant"`` names one per oxidation state and leaves the
      rest as markers, which is the readable default.

    ``"prominent"`` is accepted as a spelling of ``"predominant"``.
    """
    import matplotlib.pyplot as plt

    if diagram is None:
        diagram = frost_diagram(element, **kwargs)

    if ax is None:
        figure, ax = plt.subplots(figsize=figsize)
    else:
        figure = ax.figure

    shown = _selected(diagram, show)
    labelled = {id(point) for point in _selected(diagram, label)} & {id(point) for point in shown}

    values = np.array([point.volt_equivalent for point in shown], dtype=float)

    hull = [point for point in diagram.stable if id(point) in {id(s) for s in shown}]
    on_hull = {id(p) for p in hull}

    # The dashed chain runs through the *predominant* forms only. Running it
    # through every point instead makes it double back vertically wherever an
    # element has several forms of one oxidation state, which reads as a
    # redox step and is not one.
    chain = [point for point in diagram.predominant if id(point) in {id(s) for s in shown}]
    ax.plot(
        [point.oxidation_state for point in chain],
        [point.volt_equivalent for point in chain],
        color=PALETTE["muted"],
        linewidth=1.0,
        linestyle="--",
        zorder=2,
    )
    ax.plot(
        [point.oxidation_state for point in hull],
        [point.volt_equivalent for point in hull],
        color=PALETTE["reduction"],
        linewidth=2.0,
        zorder=3,
        label="stable forms",
    )

    # Three classes, and the middle one is the reason this needs saying: a
    # minority form of an oxidation state is not unstable, it is just not the
    # form that dominates at this pH. Drawing it as "off the hull" -- which an
    # earlier version did -- says something false about it.
    predominant = {id(point) for point in diagram.predominant}
    seen_minor = seen_unstable = False
    for point in shown:
        if id(point) not in predominant:
            seen_minor = True
            ax.plot(
                [point.oxidation_state],
                [point.volt_equivalent],
                marker="s",
                markersize=5,
                markerfacecolor="white",
                markeredgecolor=PALETTE["muted"],
                markeredgewidth=1.2,
                zorder=4,
            )
            continue
        stable = id(point) in on_hull
        seen_unstable = seen_unstable or not stable
        ax.plot(
            [point.oxidation_state],
            [point.volt_equivalent],
            marker="o" if stable else "s",
            markersize=7 if stable else 6,
            color=PALETTE["reduction"] if stable else PALETTE["endergonic"],
            zorder=5,
        )
    if seen_minor:
        ax.plot(
            [],
            [],
            marker="s",
            linestyle="none",
            markersize=5,
            markerfacecolor="white",
            markeredgecolor=PALETTE["muted"],
            markeredgewidth=1.2,
            label="other forms of the same state",
        )
        if seen_unstable:
            ax.plot(
                [],
                [],
                marker="s",
                linestyle="none",
                markersize=6,
                color=PALETTE["endergonic"],
                label="disproportionates",
            )
        # Only drawn when there is something to explain, and only for the
        # classes actually present. A diagram with one form per state needs no
        # key and did not have one before.
        ax.legend(fontsize=SIZES["annotation"] - 1, frameon=False, loc="best")

    _draw_labels(ax, shown, labelled, values)

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

    owns_figure = len(figure.axes) == 1
    if not diagram.reference_is_element and owns_figure:
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

    if owns_figure:
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
