"""The redox tower.

Both half reactions of a reaction are marked on an axis of E-standard-prime,
with an arrow showing electrons falling from the donor to the acceptor. The
title carries the overall free energy.

The y axis is a potential, so a twin axis in ATP units would be wrong: the
mapping from potential to free energy depends on the electron count. The ATP
and biological-energy-quantum information therefore appears in the side panel,
labelled with this reaction's own n. (On the free-energy explorer, where the y
axis really is kJ/mol, a twin ATP axis is used instead.)

By convention the axis is inverted so that negative potentials sit at the top
and electrons visibly fall down the tower.
"""

from __future__ import annotations

from pathlib import Path

from ..showwork import MATHTEXT_FRAC, latex_equation
from ..tower import reference_couples
from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    DEFAULT_DELTA_G_ATP,
    FARADAY,
    KJ_PER_MOL_STR,
    as_magnitude,
)
from .style import PALETTE, SIZES, atp_equivalents

_GROUP_COLOURS = {
    "hydrogen": "#6A7BA2",
    "carbon": "#5B8C5A",
    "sulfur": "#C9A227",
    "nitrogen": "#7B5EA7",
    "metal": "#A8613C",
    "oxygen": "#3C7BA8",
    "other": PALETTE["muted"],
}


def plot_redox_tower(
    reaction,
    n_electrons: int | None = 2,
    figsize=(9.5, 8.0),
    primed: bool = True,
    show_reference: bool = True,
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    delta_g_atp=None,
    energy_quantum=None,
):
    """Draw a redox tower with ``reaction``'s two couples marked.

    Normalised to an electron pair by default, matching the half-reaction
    figure. Fixing n also makes the energy scale bar meaningful: with n held
    at 2, free energy and potential differ only by the constant 2F, so a
    span in volts converts directly to kJ/mol.

    Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    if n_electrons is not None and reaction.n_electrons != n_electrons:
        reaction = reaction.renormalized(n_electrons)

    conditions = reaction.conditions
    backend = reaction.backend

    donor_e = (
        (reaction.donor_half.E_standard_prime if primed else reaction.donor_half.E_standard)
        .to("V")
        .magnitude
    )
    acceptor_e = (
        (reaction.acceptor_half.E_standard_prime if primed else reaction.acceptor_half.E_standard)
        .to("V")
        .magnitude
    )

    entries = reference_couples(conditions, backend, primed=primed) if show_reference else []

    fig, ax = plt.subplots(figsize=figsize)
    values = [donor_e, acceptor_e] + [e.potential_v for e in entries]
    low, high = min(values), max(values)
    pad = max(0.12, 0.08 * (high - low))
    ax.set_ylim(high + pad, low - pad)  # inverted: negative at the top
    ax.set_xlim(0, 1)

    label = r"$E^{\circ\prime}$ (V)" if primed else r"$E^{\circ}$ (V)"
    ax.set_ylabel(label, fontsize=SIZES["potential"])
    ax.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)

    # Reference couples as a ladder of rungs on the left. Several couples sit
    # within a few millivolts of each other, so labels are nudged apart and
    # joined to their rung by a leader line where they had to move.
    span = (high + pad) - (low - pad)
    label_y = _stagger([e.potential_v for e in entries], span * 0.033)
    for entry, y in zip(entries, label_y, strict=True):
        colour = _GROUP_COLOURS.get(entry.group, PALETTE["muted"])
        ax.hlines(entry.potential_v, 0.20, 0.34, color=colour, linewidth=2.0, alpha=0.85)
        ax.text(
            0.185,
            y,
            entry.label,
            fontsize=SIZES["annotation"],
            color=colour,
            ha="right",
            va="center",
        )
        if abs(y - entry.potential_v) > span * 0.002:
            ax.plot(
                [0.175, 0.20],
                [y, entry.potential_v],
                color=colour,
                linewidth=0.7,
                alpha=0.6,
            )

    # The reaction's own couples.
    donor_x = 0.58
    for potential, colour, text in (
        (donor_e, PALETTE["oxidation"], f"donor  ${reaction.donor.label}$"),
        (acceptor_e, PALETTE["reduction"], f"acceptor  ${reaction.acceptor.label}$"),
    ):
        ax.hlines(potential, 0.44, 0.74, color=colour, linewidth=3.2, zorder=5)
        ax.plot([donor_x], [potential], "o", color=colour, markersize=8, zorder=6)
        ax.text(
            0.76,
            potential,
            f"{text}\n{potential:+.3f} V",
            fontsize=SIZES["annotation"],
            color=colour,
            ha="left",
            va="center",
        )

    # Electrons fall from the donor to the acceptor.
    ax.annotate(
        "",
        xy=(donor_x, acceptor_e),
        xytext=(donor_x, donor_e),
        arrowprops={
            "arrowstyle": "-|>",
            "color": PALETTE["electron"],
            "linewidth": 2.0,
            "shrinkA": 4,
            "shrinkB": 4,
        },
        zorder=7,
    )
    ax.text(
        donor_x - 0.03,
        (donor_e + acceptor_e) / 2.0,
        f"{reaction.n_electrons} e$^-$",
        fontsize=SIZES["annotation"],
        color=PALETTE["electron"],
        ha="right",
        va="center",
        fontweight="bold",
    )

    _add_energy_panel(ax, reaction, delta_g_atp, energy_quantum)
    _add_energy_scalebar(ax, reaction, delta_g_atp, energy_quantum)

    delta_g = as_magnitude(
        reaction.delta_G_standard_prime if primed else reaction.delta_G_standard,
        KJ_PER_MOL_STR,
    )
    delta_e = reaction.delta_E_standard_prime.to("V").magnitude if primed else None
    symbol = r"\Delta G^{\circ\prime}" if primed else r"\Delta G^{\circ}"
    title = (
        f"${latex_equation(reaction.coefficients, frac=MATHTEXT_FRAC)}$\n"
        f"${symbol}$ = {delta_g:+.1f} kJ/mol"
        + (f", $\\Delta E^{{\\circ\\prime}}$ = {delta_e:+.3f} V" if delta_e else "")
        + f"  ({reaction.n_electrons} e$^-$, pH {conditions.pH:g}, "
        f"{conditions.temperature_c:g} °C)"
    )
    ax.set_title(title, fontsize=SIZES["title"], pad=14)

    fig.tight_layout(rect=(0, 0.16, 1, 1))
    if save is not None:
        _save(fig, save, formats, dpi)
    return fig


def volts_per_kilojoule(n_electrons) -> float:
    """Volts corresponding to 1 kJ/mol at a fixed electron count.

    From :math:`\\Delta G = -nF\\Delta E`. Only meaningful once n is pinned,
    which is why the tower normalises to an electron pair by default.
    """
    faraday_kj = as_magnitude(FARADAY.to("C/mol"), "C/mol") / 1000.0
    return 1.0 / (float(n_electrons) * faraday_kj)


def _add_energy_scalebar(ax, reaction, delta_g_atp, energy_quantum):
    """A map-style scale bar converting potential difference to free energy.

    With n fixed, a *difference* in potential maps to free energy by the
    constant nF, so a scale bar is exact. An absolute kJ/mol axis would not
    be: it would imply each couple has an absolute free energy, which it does
    not -- only differences between couples carry energy.
    """
    from matplotlib.transforms import blended_transform_factory

    n = reaction.n_electrons
    per_volt = 1.0 / volts_per_kilojoule(n)  # kJ/mol per volt
    quantum = abs(as_magnitude(energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, KJ_PER_MOL_STR))
    atp = abs(as_magnitude(delta_g_atp or DEFAULT_DELTA_G_ATP, KJ_PER_MOL_STR))

    quantum_v = quantum / per_volt
    atp_v = atp / per_volt

    transform = blended_transform_factory(ax.transAxes, ax.transData)
    low, high = ax.get_ylim()  # inverted: low is the larger value
    span = low - high
    top = high + 0.42 * span  # anchor in the empty middle of the tower
    x = 0.90

    # The bar itself: zero to one ATP, with the energy quantum marked inside.
    ax.plot(
        [x, x],
        [top, top + atp_v],
        transform=transform,
        color=PALETTE["annotation"],
        linewidth=1.6,
        solid_capstyle="butt",
        clip_on=False,
        zorder=8,
    )
    ax.plot(
        [x, x],
        [top, top + quantum_v],
        transform=transform,
        color=PALETTE["quantum_band"],
        linewidth=5.0,
        alpha=0.75,
        solid_capstyle="butt",
        clip_on=False,
        zorder=7,
    )

    for offset, label in (
        (0.0, "0"),
        (quantum_v, f"{quantum_v:.3f} V  =  {quantum:g} kJ/mol\nenergy quantum"),
        (atp_v, f"{atp_v:.3f} V  =  {atp:g} kJ/mol\n1 ATP"),
    ):
        ax.plot(
            [x - 0.012, x + 0.012],
            [top + offset, top + offset],
            transform=transform,
            color=PALETTE["annotation"],
            linewidth=1.4,
            clip_on=False,
            zorder=8,
        )
        ax.text(
            x + 0.022,
            top + offset,
            label,
            transform=transform,
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["annotation"],
            ha="left",
            va="center",
            linespacing=1.3,
            clip_on=False,
        )

    ax.text(
        x,
        top - 0.035 * span,
        f"scale at {n} e$^-$",
        transform=transform,
        fontsize=SIZES["annotation"] - 1,
        color=PALETTE["muted"],
        ha="center",
        va="bottom",
        clip_on=False,
    )


def _stagger(values, minimum_separation):
    """Nudge sorted values apart so labels do not overlap.

    A single forward pass is enough here: ``values`` arrives sorted, and the
    reference set is small. Each label is pushed just far enough past its
    predecessor to clear it, so displacement stays minimal and ordering is
    preserved.
    """
    out = []
    previous = None
    for value in values:
        if previous is not None and value - previous < minimum_separation:
            value = previous + minimum_separation
        out.append(value)
        previous = value
    return out


def _add_energy_panel(ax, reaction, delta_g_atp, energy_quantum):
    """The ATP / energy-quantum ruler, in the side panel.

    Not a twin axis: the tower's y axis is a potential, and converting that to
    free energy depends on n, so the figure states n explicitly instead.
    """
    delta_g_atp = delta_g_atp or DEFAULT_DELTA_G_ATP
    energy_quantum = energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM

    delta_g = as_magnitude(reaction.delta_G_standard_prime, KJ_PER_MOL_STR)
    per_electron = delta_g / float(reaction.n_electrons)
    atp = atp_equivalents(delta_g, delta_g_atp)
    quantum = as_magnitude(energy_quantum, KJ_PER_MOL_STR)
    per_atp = as_magnitude(delta_g_atp, KJ_PER_MOL_STR)

    viable = delta_g <= quantum
    verdict = (
        "above the biological energy quantum" if viable else "below the biological energy quantum"
    )
    colour = PALETTE["exergonic"] if viable else PALETTE["endergonic"]

    lines = [
        "Energy yield",
        f"  ΔG°′        {delta_g:+.1f} kJ/mol",
        f"  per e⁻      {per_electron:+.1f} kJ/mol e⁻",
        f"  ATP equiv.  {atp:.1f} per {reaction.n_electrons} e⁻",
        f"              (at {per_atp:g} kJ/mol per ATP)",
        "",
        f"  {verdict}",
        f"  ({quantum:g} kJ/mol)",
    ]
    ax.get_figure().text(
        0.08,
        0.005,
        "\n".join(lines),
        fontsize=SIZES["annotation"],
        va="bottom",
        ha="left",
        family="monospace",
        color=PALETTE["annotation"],
        bbox={
            "boxstyle": "round,pad=0.5",
            "facecolor": "white",
            "edgecolor": colour,
            "linewidth": 1.2,
            "alpha": 0.95,
        },
    )


def _save(fig, save, formats, dpi):
    base = Path(save)
    base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in formats:
        fig.savefig(base.with_suffix(f".{suffix}"), format=suffix, dpi=dpi, bbox_inches="tight")
