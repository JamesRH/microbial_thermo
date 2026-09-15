"""Shared figure styling and the ATP / energy-quantum scale helpers."""

from __future__ import annotations

from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    DEFAULT_DELTA_G_ATP,
    KJ_PER_MOL_STR,
    as_magnitude,
)

#: Colour-blind-safe palette, consistent across every figure in the library.
PALETTE = {
    "oxidation": "#B4451F",  # donor, warm
    "reduction": "#1F6FB4",  # acceptor, cool
    "electron": "#5B5B5B",
    "annotation": "#333333",
    "muted": "#8A8A8A",
    "grid": "#DCDCDC",
    "exergonic": "#2E7D5B",
    "endergonic": "#A63A50",
    "quantum_band": "#E8B84B",
}

#: Font sizes, in points.
SIZES = {
    "equation": 15,
    "oxidation_state": 10,
    "role_label": 12,
    "potential": 11,
    "title": 13,
    "annotation": 9,
}


def atp_equivalents(delta_g, delta_g_atp=None) -> float:
    """How many ATP the given free energy could drive.

    A modelling convenience, not a measured quantity: it simply divides by the
    free energy of ATP hydrolysis under cellular conditions, which is itself an
    estimate. Reported as a positive count for an exergonic reaction.
    """
    delta_g_atp = delta_g_atp or DEFAULT_DELTA_G_ATP
    energy = as_magnitude(delta_g, KJ_PER_MOL_STR)
    per_atp = as_magnitude(delta_g_atp, KJ_PER_MOL_STR)
    return energy / per_atp


def add_atp_axis(ax, delta_g_atp=None, label: str = "ATP equivalents"):
    """Add a right-hand twin axis expressing the left axis in ATP units.

    Only valid where the left axis is a free energy in kJ/mol. On a redox tower
    the y axis is a potential, where the mapping to free energy depends on the
    electron count, so use :func:`atp_ruler_text` there instead.
    """
    delta_g_atp = delta_g_atp or DEFAULT_DELTA_G_ATP
    per_atp = as_magnitude(delta_g_atp, KJ_PER_MOL_STR)

    twin = ax.twinx()
    low, high = ax.get_ylim()
    twin.set_ylim(low / per_atp, high / per_atp)
    twin.set_ylabel(label, fontsize=SIZES["annotation"], color=PALETTE["muted"])
    twin.tick_params(axis="y", labelsize=SIZES["annotation"], colors=PALETTE["muted"])
    return twin


def add_energy_quantum_band(ax, quantum=None, orientation: str = "horizontal"):
    """Shade the region between zero and the biological energy quantum.

    Reactions falling inside the band are exergonic but probably too weakly so
    to conserve energy for a cell.
    """
    quantum = quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM
    value = as_magnitude(quantum, KJ_PER_MOL_STR)
    span = ax.axhspan if orientation == "horizontal" else ax.axvspan
    band = span(
        min(0.0, value),
        max(0.0, value),
        color=PALETTE["quantum_band"],
        alpha=0.18,
        zorder=0,
        label=f"below the biological energy quantum ({value:g} kJ/mol)",
    )
    line = ax.axhline if orientation == "horizontal" else ax.axvline
    line(0.0, color=PALETTE["annotation"], linewidth=0.8, zorder=1)
    return band


def atp_ruler_text(n_electrons, delta_g, delta_g_atp=None) -> str:
    """A short caption expressing a reaction's energy in ATP terms.

    Used where a twin axis would be misleading, such as beside a redox tower
    whose y axis is a potential.
    """
    count = atp_equivalents(delta_g, delta_g_atp)
    per_atp = as_magnitude(delta_g_atp or DEFAULT_DELTA_G_ATP, KJ_PER_MOL_STR)
    return f"{count:.1f} ATP equivalents per {n_electrons} e$^-$\n(at {per_atp:g} kJ/mol per ATP)"


def unverified_footnote(reaction) -> str | None:
    """Warning line for a figure that rests on a hand-entered value.

    SPEC section 2.1 requires that a figure depending on an unverified
    supplemental number say so. A reader looking at a rendered figure has no
    other way to know, since the console warning is long gone by then.
    """
    from ..supplemental import uses_unverified

    names = uses_unverified(s.backend for s in reaction.coefficients)
    if not names:
        return None
    return (
        "Depends on hand-entered, unverified formation energies for "
        + ", ".join(names)
        + " — check them before relying on this figure."
    )


def add_unverified_footnote(figure, reaction) -> None:
    """Draw :func:`unverified_footnote` along the bottom of a figure."""
    text = unverified_footnote(reaction)
    if text is None:
        return
    figure.text(
        0.5,
        0.005,
        text,
        fontsize=SIZES["annotation"] - 1,
        color=PALETTE["endergonic"],
        ha="center",
        va="bottom",
        wrap=True,
    )
