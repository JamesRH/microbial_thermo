"""Arsenite-driven carbon fixation, and how pH moves it.

As(III) -> As(V) coupled to CO2 fixation into a biomass placeholder: the
anabolic half of chemolithoautotrophy in an arsenite oxidiser such as
*Rhizobium* sp. NT-26 or *Alkalilimnicola ehrlichii*.

    As(OH)3 + 1/2 CO2 + 1/2 H2O  ->  HAsO4(2-) + 2 H+ + 1/2 <CH2O>

The reaction is **endergonic** -- around +80 kJ/mol per electron pair at
pH 7 -- and that is the point. Fixing carbon costs energy. An autotroph pays
for it out of its catabolism, here arsenite oxidation with oxygen, so the
figure draws both on one axis: the cost, the income, and the balance.

Why pH matters. Arsenite oxidation releases two protons per electron pair, so
every term shifts by RT ln(10) x 2 = 11.4 kJ/mol per pH unit at 25 C. Both
halves move at that rate and in the same direction, so the *net* improves
twice as fast -- about 23 kJ/mol per pH unit. Alkaline water is
thermodynamically kinder to this organism, which is one reason the best
studied arsenite oxidisers come out of soda lakes.

The biomass term is a placeholder. <CH2O> is cell carbon at oxidation state
zero, and its formation energy is this library's glucose divided by six. Real
biomass is nearer CH1.8O0.5N0.2 and contains nitrogen this term ignores
entirely, so treat the numbers as sign-and-scale, not as quotable yields. Use
``--biomass-dgf`` to see how much your conclusion depends on that choice; the
answer, reassuringly, is that the sign does not.

Run it:

    python examples/arsenite_carbon_fixation.py --save figures/as_fixation
    python examples/arsenite_carbon_fixation.py --biomass-dgf -130 --no-catabolic

or import ``fixation_vs_ph`` and call it from a notebook.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import click
import numpy as np

import microbial_thermo as mt
from microbial_thermo.figures.style import (
    PALETTE,
    SIZES,
    add_energy_quantum_band,
    unverified_footnote,
)
from microbial_thermo.library import reaction

#: The anabolic reaction, and the catabolism that has to pay for it.
FIXATION = "arsenite_carbon_fixation"
CATABOLIC = "arsenite_oxidation_oxygen"

#: RT ln(10) at 25 C, kJ/mol -- one pH unit's worth of one proton.
RT_LN10_25C = 5.708


def _print_versions(context, parameter, value):
    """--version, reusing the library's own reporting so the two cannot drift."""
    if not value or context.resilient_parsing:
        return
    for name, version in mt.versions().items():
        click.echo(f"{name}: {version}")
    context.exit()


_version_option = click.option(
    "--version",
    is_flag=True,
    callback=_print_versions,
    expose_value=False,
    is_eager=True,
    help="Show the library and dependency versions, then exit.",
)


def _conditions(pH: float, temperature_c: float) -> mt.Conditions:
    return mt.Conditions(temperature_c=temperature_c, pH=pH)


def fixation_vs_ph(
    ph_values=None,
    temperature_c: float = 25.0,
    include_catabolic: bool = True,
    biomass_dgf: float | None = None,
):
    """Free energy of arsenite-driven carbon fixation across pH.

    Returns a dict of arrays: ``pH``, ``anabolic``, and -- when
    ``include_catabolic`` -- ``catabolic`` and ``net``. Energies are kJ/mol
    per electron pair, which is the basis both reactions are normalised to,
    so they may be added.

    ``biomass_dgf`` overrides the placeholder's formation energy for a
    sensitivity run; it is restored afterwards.
    """
    ph_values = np.asarray(
        ph_values if ph_values is not None else np.linspace(4.0, 10.0, 25),
        dtype=float,
    )

    with _biomass_energy(biomass_dgf), warnings.catch_warnings():
        # The placeholder warns on every use; once per run is plenty, and the
        # figure carries the notice anyway.
        warnings.simplefilter("ignore")
        anabolic = np.array(
            [reaction(FIXATION, _conditions(p, temperature_c)).delta_G.magnitude for p in ph_values]
        )
        out = {"pH": ph_values, "anabolic": anabolic}
        if include_catabolic:
            catabolic = np.array(
                [
                    reaction(CATABOLIC, _conditions(p, temperature_c)).delta_G.magnitude
                    for p in ph_values
                ]
            )
            out["catabolic"] = catabolic
            out["net"] = anabolic + catabolic
    return out


class _biomass_energy:
    """Temporarily override the biomass placeholder's formation energy.

    A context manager rather than a parameter because the value lives in the
    supplemental table, which the reaction layer reads through a cache.
    """

    def __init__(self, value):
        self.value = value
        self._saved = None

    def __enter__(self):
        if self.value is None:
            return self
        from microbial_thermo import supplemental

        supplemental.supplemental_species.cache_clear()
        entry = supplemental.supplemental_species()["Biomass(aq)"]
        self._saved = entry.delta_Gf_kJ_mol
        object.__setattr__(entry, "delta_Gf_kJ_mol", float(self.value))
        _clear_backend_cache()
        return self

    def __exit__(self, *exc):
        if self.value is None:
            return False
        from microbial_thermo import supplemental

        entry = supplemental.supplemental_species()["Biomass(aq)"]
        object.__setattr__(entry, "delta_Gf_kJ_mol", self._saved)
        _clear_backend_cache()
        return False


def _clear_backend_cache():
    """Formation energies are memoised per species and temperature."""
    backend = mt.get_backend()
    for attr in ("_gibbs_cache", "delta_Gf"):
        cache = getattr(backend, attr, None)
        if hasattr(cache, "cache_clear"):
            cache.cache_clear()
        elif isinstance(cache, dict):
            cache.clear()


def slope_per_ph_unit(result) -> float:
    """Observed kJ/mol per pH unit for the anabolic curve.

    Worth reporting because it is predictable: two protons leave per electron
    pair, so it should come out at 2 x RT ln(10).
    """
    fit = np.polyfit(result["pH"], result["anabolic"], 1)
    return float(fit[0])


def plot_fixation_vs_ph(
    result=None,
    temperature_c: float = 25.0,
    include_catabolic: bool = True,
    biomass_dgf: float | None = None,
    figsize=(9.0, 5.75),
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
):
    """Draw the cost, the income and the balance against pH.

    Returns ``(figure, result)``.
    """
    import matplotlib.pyplot as plt

    if result is None:
        result = fixation_vs_ph(
            temperature_c=temperature_c,
            include_catabolic=include_catabolic,
            biomass_dgf=biomass_dgf,
        )
    has_catabolic = "catabolic" in result

    figure, ax = plt.subplots(figsize=figsize)
    add_energy_quantum_band(ax)

    ax.plot(
        result["pH"],
        result["anabolic"],
        color=PALETTE["endergonic"],
        linewidth=2.4,
        label="CO₂ fixation to ⟨CH₂O⟩ — the cost",
        zorder=4,
    )
    if has_catabolic:
        ax.plot(
            result["pH"],
            result["catabolic"],
            color=PALETTE["reduction"],
            linewidth=2.4,
            label="As(III) + O₂ — the income",
            zorder=4,
        )
        ax.plot(
            result["pH"],
            result["net"],
            color=PALETTE["exergonic"],
            linewidth=2.4,
            linestyle="--",
            label="net, one fixation per oxidation",
            zorder=5,
        )

    slope = slope_per_ph_unit(result)
    ax.annotate(
        f"{slope:+.1f} kJ/mol per pH unit\n= 2 H⁺ × RT ln10 ({2 * RT_LN10_25C:.1f})",
        xy=(result["pH"][len(result["pH"]) // 2], result["anabolic"][len(result["pH"]) // 2]),
        xytext=(14, 26),
        textcoords="offset points",
        fontsize=SIZES["annotation"],
        color=PALETTE["annotation"],
        arrowprops={"arrowstyle": "-", "color": PALETTE["muted"], "linewidth": 0.8},
    )

    ax.set_xlabel("pH", fontsize=SIZES["potential"])
    ax.set_ylabel("ΔG (kJ/mol per 2 e⁻)", fontsize=SIZES["potential"])
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    # Below the axes: the three curves between them leave no clear corner,
    # and the net line runs straight through the middle.
    ax.legend(
        fontsize=SIZES["annotation"],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.12),
        ncol=2,
    )
    ax.set_title(
        "Arsenite-driven carbon fixation against pH\n"
        f"As(OH)₃ + ½CO₂ + ½H₂O → HAsO₄²⁻ + 2H⁺ + ½⟨CH₂O⟩   ({temperature_c:g} °C)",
        fontsize=SIZES["title"],
        pad=12,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        notice = unverified_footnote(reaction(FIXATION, _conditions(7.0, temperature_c)))
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


@click.command()
@click.option("--ph-low", default=4.0, show_default=True, help="Lowest pH.")
@click.option("--ph-high", default=10.0, show_default=True, help="Highest pH.")
@click.option("--points", default=25, show_default=True, help="Grid points.")
@click.option("--temperature", default=25.0, show_default=True, help="Temperature in °C.")
@click.option(
    "--biomass-dgf",
    default=None,
    type=float,
    help="Override ΔGf of the ⟨CH2O⟩ placeholder, kJ/mol, for a sensitivity run.",
)
@click.option(
    "--catabolic/--no-catabolic",
    default=True,
    show_default=True,
    help="Also draw the arsenite/O2 catabolism that pays for the fixation.",
)
@click.option("--save", default=None, help="Write SVG and PNG to this path stem.")
@click.option("--show/--no-show", default=False, help="Open an interactive window.")
@_version_option
def main(ph_low, ph_high, points, temperature, biomass_dgf, catabolic, save, show):
    """Plot arsenite-driven CO2 fixation against pH."""
    import matplotlib

    if not show:
        matplotlib.use("Agg")

    result = fixation_vs_ph(
        ph_values=np.linspace(ph_low, ph_high, points),
        temperature_c=temperature,
        include_catabolic=catabolic,
        biomass_dgf=biomass_dgf,
    )
    figure, _ = plot_fixation_vs_ph(result=result, temperature_c=temperature, save=save)

    click.echo(f"ΔG of CO2 fixation, kJ/mol per 2 e- ({temperature:g} °C)")
    for pH, value in zip(result["pH"], result["anabolic"], strict=True):
        if abs(pH - round(pH)) < 1e-9:
            line = f"  pH {pH:5.1f}   {value:+8.1f}"
            if "net" in result:
                index = int(np.argmin(np.abs(result["pH"] - pH)))
                line += f"   net {result['net'][index]:+8.1f}"
            click.echo(line)
    click.echo(f"\nslope: {slope_per_ph_unit(result):+.2f} kJ/mol per pH unit")
    click.echo(f"expected from 2 H+: {-2 * RT_LN10_25C:+.2f}")
    if save:
        click.echo(f"wrote {save}.svg and {save}.png")
    if show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
