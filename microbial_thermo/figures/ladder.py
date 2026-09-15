"""The affinity ladder.

Free energy per electron for many metabolisms at one set of conditions, sorted.
Where the redox tower shows what is *possible* from standard potentials, this
shows what actually pays under a particular measured water chemistry -- which is
how redox zonation arises in a sediment rather than being asserted.

The x axis is kJ/mol, so ATP equivalents are a genuine second axis here, and the
biological energy quantum is a shaded band. Bars crossing into the band are
exergonic but probably too weakly so to conserve energy.
"""

from __future__ import annotations

from pathlib import Path

from ..conditions import Conditions
from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    DEFAULT_DELTA_G_ATP,
    KJ_PER_MOL_STR,
    as_magnitude,
)
from .style import PALETTE, SIZES

#: One colour per metabolic group, so related metabolisms read together.
_GROUP_COLOURS = {
    "aerobic respiration": "#3C7BA8",
    "nitrification": "#7B5EA7",
    "denitrification": "#6A5ACD",
    "nitrogen": "#9370DB",
    "sulfur oxidation": "#C9A227",
    "sulfate reduction": "#B8860B",
    "sulfur reduction": "#DAA520",
    "metal oxidation": "#A8613C",
    "metal reduction": "#8B4513",
    "methanogenesis": "#5B8C5A",
    "methane oxidation": "#2E8B57",
    "acetogenesis": "#6B8E23",
}


def plot_affinity_ladder(
    conditions: Conditions | None = None,
    backend=None,
    n_electrons: int = 2,
    group: str | None = None,
    table=None,
    figsize=(10.0, 9.0),
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    delta_g_atp=None,
    energy_quantum=None,
    title: str | None = None,
):
    """Sorted bar chart of free energy per electron across the curated library.

    Pass ``table`` to supply a precomputed
    :func:`~microbial_thermo.library.energy_table`; otherwise one is built.
    Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    from ..library import energy_table

    conditions = conditions or Conditions()
    if table is None:
        table = energy_table(conditions, backend, n_electrons=n_electrons, group=group)

    usable = table[table["problem"] == ""].copy()
    skipped = table[table["problem"] != ""]
    if usable.empty:
        raise ValueError("no metabolism could be evaluated under these conditions")

    # Most negative at the top.
    usable = usable.sort_values("dG per e- (kJ/mol)", ascending=False)
    values = usable["dG per e- (kJ/mol)"].tolist()
    labels = usable["label"].tolist()
    groups = usable["group"].tolist()
    colours = [_GROUP_COLOURS.get(g, PALETTE["muted"]) for g in groups]

    quantum = as_magnitude(energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, KJ_PER_MOL_STR)
    per_atp = as_magnitude(delta_g_atp or DEFAULT_DELTA_G_ATP, KJ_PER_MOL_STR)

    fig, ax = plt.subplots(figsize=figsize)
    positions = range(len(values))
    ax.barh(list(positions), values, color=colours, height=0.72, zorder=3)

    ax.set_yticks(list(positions))
    ax.set_yticklabels(labels, fontsize=SIZES["annotation"])
    ax.set_xlabel(
        f"ΔG per electron (kJ/mol e$^-$), at {n_electrons} e$^-$",
        fontsize=SIZES["potential"],
    )
    ax.axvline(0.0, color=PALETTE["annotation"], linewidth=1.0, zorder=4)

    # The band between zero and the energy quantum: exergonic but too weak.
    ax.axvspan(
        min(0.0, quantum),
        max(0.0, quantum),
        color=PALETTE["quantum_band"],
        alpha=0.22,
        zorder=1,
    )
    # Label the band beneath the bars, where nothing can collide with it.
    ax.annotate(
        f"below the biological energy quantum ({quantum:g} kJ/mol)",
        xy=(quantum / 2.0, -0.9),
        xycoords=("data", "data"),
        fontsize=SIZES["annotation"] - 1,
        color=PALETTE["muted"],
        ha="center",
        va="top",
        annotation_clip=False,
    )
    ax.set_ylim(-1.6, len(values) - 0.3)

    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    _add_atp_axis(ax, per_atp)
    _add_group_legend(ax, groups)

    if title is None:
        title = f"What pays, at pH {conditions.pH:g} and {conditions.temperature_c:g} °C"
        if group:
            title = f"{group}: {title}"
    ax.set_title(title, fontsize=SIZES["title"], pad=14)

    if len(skipped):
        fig.text(
            0.01,
            0.005,
            "not evaluable here: " + ", ".join(skipped["label"].tolist()),
            fontsize=SIZES["annotation"] - 2,
            color=PALETTE["muted"],
            va="bottom",
        )

    fig.tight_layout()
    if save is not None:
        _save(fig, save, formats, dpi)
    return fig


def _add_atp_axis(ax, per_atp):
    """A true second axis: the primary axis is already kJ/mol."""
    twin = ax.twiny()
    low, high = ax.get_xlim()
    twin.set_xlim(low / per_atp, high / per_atp)
    twin.set_xlabel(
        f"ATP equivalents per electron (at {per_atp:g} kJ/mol per ATP)",
        fontsize=SIZES["annotation"],
        color=PALETTE["muted"],
    )
    twin.tick_params(axis="x", labelsize=SIZES["annotation"], colors=PALETTE["muted"])
    for side in ("top", "right", "left"):
        twin.spines[side].set_visible(False)
    return twin


def _add_group_legend(ax, groups):
    from matplotlib.patches import Patch

    seen: list[str] = []
    for group in groups:
        if group not in seen:
            seen.append(group)
    handles = [Patch(facecolor=_GROUP_COLOURS.get(g, PALETTE["muted"]), label=g) for g in seen]
    # Outside the axes: with a dozen groups a legend inside would sit on top
    # of the smallest bars, which are the interesting ones.
    ax.legend(
        handles=handles,
        fontsize=SIZES["annotation"] - 1,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )


def _save(fig, save, formats, dpi):
    base = Path(save)
    base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in formats:
        fig.savefig(base.with_suffix(f".{suffix}"), format=suffix, dpi=dpi, bbox_inches="tight")
