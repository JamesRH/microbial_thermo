"""Eh-pH (Pourbaix) diagrams: which form of an element wins, and where.

A predominance diagram over the two variables that actually decide aqueous
speciation. Every species of an element is formed from a common basis --
the element itself, water, protons and electrons -- and at each point in
(Eh, pH) the one with the lowest free energy per mole of element is drawn.

The construction, since it is short enough to state:

    n_E * E + b * H2O + c * H+ + d * e-  ->  S

with ``b`` set by the oxygen in S, ``c`` by its hydrogen after the water is
accounted for, and ``d`` by charge. Then at a given pH and Eh

    dG(S) = dGf(S) + RT ln(a) - n_E*dGf(E) - b*dGf(H2O)
            + c * RT ln(10) * pH + d * F * Eh

all divided by ``n_E``, because the comparison is per mole of the element.
Protons and electrons carry no formation energy on the SHE convention, which
is what makes those last two terms the whole pH and Eh dependence.

**Why this is built rather than borrowed.** pymatgen has a mature
``PourbaixDiagram``, but its pipeline wants a Materials Project API key and
builds entries from DFT solid energies plus experimental ion energies. Mixing
those with our HKF aqueous species would put two provenances inside one
diagram -- the same objection that rules out PHREEQC for speciation. The
geometry here is a few dozen lines; the data consistency is the part worth
protecting.

The dissolved activity matters and is not a detail. A Pourbaix diagram is
drawn at a chosen total dissolved activity, conventionally 1e-6, and the
solid fields grow as you lower it. The value is always printed on the figure
for that reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .basis import (
    DEFAULT_ACTIVITY,
    DEFAULT_SPECIES,
    ELEMENT_REFERENCE,
    FIXED_DEFAULTS,
    basis_coefficients,
    element_series,
    pretty,
)
from .style import PALETTE, SIZES

# The basis decomposition, the species lists and the element references are
# shared with the Latimer and Frost diagrams and live in `basis`. They are
# re-exported here because this module is where they were first written and
# where a reader looking for them will go.
__all__ = [
    "DEFAULT_ACTIVITY",
    "DEFAULT_SPECIES",
    "ELEMENT_REFERENCE",
    "FIXED_DEFAULTS",
    "PourbaixField",
    "plot_pourbaix",
    "pourbaix_field",
    "water_stability",
]

#: Kept as private aliases: the tests and any existing caller import these
#: names from this module.
_basis_coefficients = basis_coefficients
_pretty = pretty


@dataclass
class PourbaixField:
    """Which species predominates across an (Eh, pH) grid."""

    element: str
    ph: np.ndarray
    eh: np.ndarray
    winner: np.ndarray  # index into `species`, shape (len(eh), len(ph))
    species: tuple
    activity: float
    temperature_c: float
    fixed: tuple = ()

    @property
    def present(self) -> list:
        """Species that actually win somewhere. Often far fewer than offered."""
        return [self.species[i] for i in sorted(set(self.winner.ravel().tolist()))]

    def at(self, pH: float, eh: float) -> str:
        i = int(np.argmin(np.abs(self.eh - eh)))
        j = int(np.argmin(np.abs(self.ph - pH)))
        return self.species[int(self.winner[i, j])]


def pourbaix_field(
    element: str,
    species=None,
    temperature_c: float = 25.0,
    activity: float = DEFAULT_ACTIVITY,
    ph_range=(0.0, 14.0),
    eh_range=(-1.0, 1.4),
    points: int = 240,
    fixed=None,
    backend=None,
) -> PourbaixField:
    """Work out the predominant species of ``element`` across (Eh, pH).

    ``fixed`` supplies a basis species and activity for every element other
    than the one being drawn, as ``{"S": ("SO4-2", 1e-3)}``. Iron's diagram
    only has a siderite or pyrite field because carbonate and sulfide are
    present at some concentration, so that concentration is an input, not a
    detail -- and a species whose extra elements are not covered is dropped
    rather than silently mis-weighted.
    """
    series = element_series(
        element,
        species=species,
        temperature_c=temperature_c,
        activity=activity,
        fixed=fixed,
        backend=backend,
    )

    ph = np.linspace(*ph_range, points)
    eh = np.linspace(*eh_range, points)
    grid_ph, grid_eh = np.meshgrid(ph, eh)
    stacked = np.stack([series.energy(entry, grid_ph, grid_eh) for entry in series])

    return PourbaixField(
        element=element,
        ph=ph,
        eh=eh,
        winner=np.argmin(stacked, axis=0),
        species=tuple(entry.backend for entry in series),
        activity=activity,
        temperature_c=temperature_c,
        fixed=series.fixed,
    )


def water_stability(temperature_c: float = 25.0, ph=None, backend=None):
    """The Eh bounds within which liquid water is stable.

    Computed from this library's own O2/H2O and H+/H2 couples rather than
    quoted as 1.229 and 0, so the lines move correctly with temperature.
    Returns ``(pH, upper, lower)``.
    """
    from .. import get_backend
    from ..conditions import Conditions
    from ..reaction import Couple, HalfReactionResult

    backend = backend or get_backend()
    ph = np.linspace(0.0, 14.0, 57) if ph is None else np.asarray(ph)

    def line(reduced, oxidized):
        values = []
        for value in ph:
            conditions = Conditions(temperature_c=temperature_c, pH=float(value))
            couple = Couple.make(reduced, oxidized)
            result = HalfReactionResult(
                couple=couple,
                half=couple.half_reaction(),
                conditions=conditions,
                backend=backend,
            )
            values.append(result.E_standard_prime.to("V").magnitude)
        return np.asarray(values)

    # O2(g) and H2(g) at 1 bar, which is what the published 1.229 V refers to.
    # Against O2(aq) at unit activity the upper line sits 43 mV higher -- a
    # different reference state, not an error, but the wrong one to draw here.
    return ph, line("H2O", "O2(g)"), line("H2(g)", "H+")


def plot_pourbaix(
    element: str,
    field: PourbaixField | None = None,
    show_water: bool = True,
    figsize=(8.5, 6.5),
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    **kwargs,
):
    """Draw the predominance diagram. Returns ``(figure, field)``."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    if field is None:
        field = pourbaix_field(element, **kwargs)

    present = sorted(set(field.winner.ravel().tolist()))
    remap = {old: new for new, old in enumerate(present)}
    drawn = np.vectorize(remap.get)(field.winner)
    labels = [field.species[i] for i in present]

    palette = plt.cm.tab20(np.linspace(0, 1, max(len(labels), 2)))
    figure, ax = plt.subplots(figsize=figsize)
    ax.pcolormesh(
        field.ph,
        field.eh,
        drawn,
        cmap=ListedColormap(palette[: len(labels)]),
        shading="auto",
        alpha=0.75,
    )
    ax.contour(
        field.ph,
        field.eh,
        drawn,
        levels=np.arange(len(labels)) + 0.5,
        colors=PALETTE["annotation"],
        linewidths=0.8,
    )

    # Label each field at its centroid, which is where it is widest.
    for new, name in enumerate(labels):
        mask = drawn == new
        if mask.sum() < drawn.size * 0.004:
            continue
        rows, cols = np.nonzero(mask)
        ax.text(
            field.ph[int(np.median(cols))],
            field.eh[int(np.median(rows))],
            _pretty(name),
            ha="center",
            va="center",
            fontsize=SIZES["annotation"] + 1,
            color=PALETTE["annotation"],
            bbox={"facecolor": "white", "alpha": 0.55, "edgecolor": "none", "pad": 1.5},
        )

    if show_water:
        ph, upper, lower = water_stability(field.temperature_c)
        for values, label in ((upper, "O$_2$/H$_2$O"), (lower, "H$^+$/H$_2$")):
            ax.plot(ph, values, color=PALETTE["reduction"], linestyle="--", linewidth=1.3)
            ax.annotate(
                label,
                xy=(ph[-1], values[-1]),
                xytext=(-4, 4),
                textcoords="offset points",
                ha="right",
                fontsize=SIZES["annotation"] - 1,
                color=PALETTE["reduction"],
            )

    ax.set_xlabel("pH", fontsize=SIZES["potential"])
    ax.set_ylabel("Eh (V vs SHE)", fontsize=SIZES["potential"])
    ax.set_xlim(field.ph[0], field.ph[-1])
    ax.set_ylim(field.eh[0], field.eh[-1])
    subtitle = f"dissolved activity {field.activity:g} · {field.temperature_c:g} °C"
    if field.fixed:
        held = ", ".join(f"{_pretty(name)} {value:g}" for _, name, value in field.fixed)
        subtitle += f"\nheld fixed: {held}"
    ax.set_title(f"{element} predominance\n{subtitle}", fontsize=SIZES["title"], pad=10)
    figure.tight_layout()

    if save is not None:
        base = Path(save)
        if base.parent != Path(""):
            base.parent.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            figure.savefig(
                base.with_suffix(f".{suffix}"), format=suffix, dpi=dpi, bbox_inches="tight"
            )
    return figure, field
