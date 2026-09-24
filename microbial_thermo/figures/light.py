"""The light budget: an uphill reaction, and the photons that pay for it.

A phototroph's defining trick is running a reaction that does not pay.
Photoferrotrophy, photoarsenotrophy and nitrite-driven photoautotrophy all
fix CO2 with a donor far too weak to do it in the dark, and make up the
shortfall from absorbed light.

This figure puts both halves on one axis: the reaction's free energy against
pH, and the same curve credited with one, two, three photons. Where a credited
curve drops below the biological energy quantum, that many photons' worth of
energy would cover the reaction.

What the figure does **not** claim is a quantum requirement. It answers the
thermodynamic question -- how many photons' worth of energy is the reaction
short by -- which is a floor. Real anoxygenic phototrophs run cyclic electron
flow and reverse electron transport and spend considerably more than the
floor. See :mod:`microbial_thermo.photons` for the rest of that caveat.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..conditions import Conditions
from ..photons import (
    DEFAULT_WAVELENGTH_NM,
    light_available,
    photon_energy,
    photons_required,
)
from ..showwork import MATHTEXT_FRAC, latex_equation
from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    KJ_PER_MOL_STR,
    as_magnitude,
)
from .style import PALETTE, SIZES, add_energy_quantum_band, unverified_footnote

#: Photon counts drawn by default. Two is enough for every reaction this was
#: built for -- the hungriest, Fe(II)/Fe(III) at pH 7, is short by 1.79 -- and
#: more lines than that just stretch the y axis.
DEFAULT_PHOTON_COUNTS = (1, 2)


def light_budget(
    reaction,
    ph_values=None,
    photons=DEFAULT_PHOTON_COUNTS,
    wavelength_nm=DEFAULT_WAVELENGTH_NM,
    efficiency: float = 1.0,
):
    """Free energy against pH, in the dark and with light credited.

    Returns ``{"pH": ..., "dark": ..., 1: ..., 2: ...}`` where the integer
    keys are photon counts. Energies are kJ/mol at the reaction's own electron
    normalisation.
    """
    from ..reaction import Reaction

    ph_values = np.asarray(
        ph_values if ph_values is not None else np.linspace(4.0, 10.0, 25),
        dtype=float,
    )

    dark = []
    for value in ph_values:
        conditions = _with_ph(reaction.conditions, float(value))
        rebuilt = Reaction.from_couples(
            donor=reaction.donor,
            acceptor=reaction.acceptor,
            conditions=conditions,
            backend=reaction.backend,
            n_electrons=reaction.n_electrons,
        )
        dark.append(as_magnitude(rebuilt.delta_G, KJ_PER_MOL_STR))
    dark = np.asarray(dark, dtype=float)

    out = {"pH": ph_values, "dark": dark}
    for count in photons:
        credit = light_available(count, wavelength_nm, efficiency).magnitude
        out[count] = dark - credit
    return out


def _with_ph(conditions: Conditions, pH: float) -> Conditions:
    return conditions.replace(pH=pH)


def plot_light_budget(
    reaction,
    ph_values=None,
    photons=DEFAULT_PHOTON_COUNTS,
    wavelength_nm=DEFAULT_WAVELENGTH_NM,
    efficiency: float = 1.0,
    figsize=(9.0, 5.75),
    label: str | None = None,
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    result=None,
):
    """Draw the uphill reaction and the photon credits against pH.

    Returns ``(figure, result)``.
    """
    import matplotlib.pyplot as plt

    if result is None:
        result = light_budget(reaction, ph_values, photons, wavelength_nm, efficiency)
    counts = [k for k in result if isinstance(k, int)]
    per_photon = light_available(1, wavelength_nm, efficiency).magnitude
    quantum = as_magnitude(DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, KJ_PER_MOL_STR)

    figure, ax = plt.subplots(figsize=figsize)
    add_energy_quantum_band(ax)

    ax.plot(
        result["pH"],
        result["dark"],
        color=PALETTE["endergonic"],
        linewidth=2.6,
        label=label or "in the dark — uphill",
        zorder=5,
    )

    # Successively lighter shades for successively more photons.
    shades = plt.cm.YlOrBr(np.linspace(0.75, 0.35, max(len(counts), 1)))
    for colour, count in zip(shades, counts, strict=True):
        pays = result[count] < quantum
        ax.plot(
            result["pH"],
            result[count],
            color=colour,
            linewidth=2.0,
            linestyle="--",
            label=f"+{count} photon{'s' if count > 1 else ''}" + ("  ✓ pays" if pays.any() else ""),
            zorder=4,
        )

    needed = photons_required(
        float(np.interp(7.0, result["pH"], result["dark"])),
        wavelength_nm,
        efficiency,
    )
    ax.annotate(
        f"at pH 7: short by {needed:.2f} photons\n"
        f"({per_photon:.0f} kJ/mol each at {_wavelength_label(wavelength_nm)})",
        xy=(7.0, float(np.interp(7.0, result["pH"], result["dark"]))),
        xytext=(12, 20),
        textcoords="offset points",
        fontsize=SIZES["annotation"],
        color=PALETTE["annotation"],
        arrowprops={"arrowstyle": "-", "color": PALETTE["muted"], "linewidth": 0.8},
    )

    ax.set_xlabel("pH", fontsize=SIZES["potential"])
    ax.set_ylabel(f"ΔG (kJ/mol per {reaction.n_electrons} e⁻)", fontsize=SIZES["potential"])
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(
        fontsize=SIZES["annotation"],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.12),
        ncol=2,
    )
    ax.set_title(
        f"${latex_equation(reaction.coefficients, frac=MATHTEXT_FRAC)}$\n"
        f"light budget at {_wavelength_label(wavelength_nm)}"
        + ("" if efficiency == 1.0 else f", {efficiency:.0%} efficiency"),
        fontsize=SIZES["title"],
        pad=12,
    )

    notice = unverified_footnote(reaction)
    if notice is not None:
        figure.text(
            0.5,
            0.005,
            notice,
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["endergonic"],
            ha="center",
            va="bottom",
            wrap=True,
        )

    figure.tight_layout()
    if save is not None:
        base = Path(save)
        if base.parent != Path(""):
            base.parent.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            figure.savefig(
                base.with_suffix(f".{suffix}"),
                format=suffix,
                dpi=dpi,
                bbox_inches="tight",
            )
    return figure, result


def _wavelength_label(wavelength_nm) -> str:
    if isinstance(wavelength_nm, str):
        from ..photons import REACTION_CENTRES

        return f"{wavelength_nm.upper()} ({REACTION_CENTRES[wavelength_nm.upper()]:g} nm)"
    return f"{float(wavelength_nm):g} nm"


def light_budget_table(reactions, wavelength_nm=DEFAULT_WAVELENGTH_NM, efficiency=1.0):
    """Compare several uphill reactions and what light each needs.

    ``reactions`` is a mapping of label to built reaction.
    """
    import pandas as pd

    rows = []
    for label, built in reactions.items():
        energy = built.delta_G_standard_prime
        rows.append(
            {
                "reaction": label,
                "dG0' (kJ/mol)": as_magnitude(energy, KJ_PER_MOL_STR),
                "photons needed": photons_required(energy, wavelength_nm, efficiency),
                "kJ/mol per photon": photon_energy(wavelength_nm).magnitude * efficiency,
            }
        )
    return pd.DataFrame(rows)
