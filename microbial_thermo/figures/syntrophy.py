"""The syntrophic window.

Two partners' free energies against hydrogen partial pressure, on one axis,
with the overlap shaded where both are exergonic.

The figure exists because syntrophy is otherwise hard to believe. A fermenter
oxidising propionate is *endergonic* under standard conditions and only pays
once hydrogen has been drawn down; its methanogenic partner only pays once
there is enough hydrogen to work with. Neither can run alone, and the range
where both can is narrow -- often a couple of decades of hydrogen pressure.
Drawing both curves on one axis makes that a picture rather than an assertion.

Both curves are normalised to the same electron count, which is what makes
putting them on one axis legitimate; the plot refuses if they are not.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..exceptions import MicrobialThermoError
from ..sweep import partial_pressure_axis, sweep
from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    KJ_PER_MOL_STR,
    as_magnitude,
)
from .style import (
    PALETTE,
    SIZES,
    add_energy_quantum_band,
    unverified_footnote,
)

#: Hydrogen is the usual currency, but formate works the same way.
DEFAULT_CARRIER = "H2(g)"


@dataclass
class SyntrophyWindow:
    """Both partners' free energies across a carrier pressure, and the
    overlap where both are exergonic."""

    carrier: str
    pressures: np.ndarray
    producer_delta_g: np.ndarray
    consumer_delta_g: np.ndarray
    n_electrons: float
    low: float | None = None
    high: float | None = None

    @property
    def exists(self) -> bool:
        return self.low is not None and self.high is not None

    @property
    def best_shared(self) -> tuple[float, float]:
        """The carrier pressure at which the *worse-off* partner does best,
        and that partner's free energy there.

        This is the honest ceiling on the partnership: at any pressure one of
        the two is doing at least this badly. It is usually a few kJ/mol --
        inside the biological energy quantum -- which is the whole reason
        syntrophs live so close to the thermodynamic edge.
        """
        worse = np.maximum(self.producer_delta_g, self.consumer_delta_g)
        index = int(np.nanargmin(worse))
        return float(self.pressures[index]), float(worse[index])

    @property
    def decades(self) -> float:
        """Width of the window in orders of magnitude, or 0 if there is none."""
        if not self.exists:
            return 0.0
        return float(np.log10(self.high / self.low))

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(
            {
                f"p{self.carrier} (bar)": self.pressures,
                "producer dG (kJ/mol)": self.producer_delta_g,
                "consumer dG (kJ/mol)": self.consumer_delta_g,
                "both exergonic": (self.producer_delta_g < 0) & (self.consumer_delta_g < 0),
            }
        )


def _carrier_label(carrier: str) -> str:
    """A pressure label that cannot be misread as pH.

    ``partial_pressure_axis`` labels its axis ``pH2(g) (bar)``, which in a
    library full of pH is a genuinely bad thing to write on an axis.
    """
    return f"{carrier} partial pressure (bar)"


def _zero_crossing(pressures, values, index):
    """Interpolate, in log pressure, where ``values`` crosses zero between
    ``index`` and ``index + 1``. Log because that is how the curve behaves:
    every term in RT ln Q is linear in log p."""
    before, after = values[index], values[index + 1]
    if before == after or not (np.isfinite(before) and np.isfinite(after)):
        return float(pressures[index])
    low, high = np.log10(pressures[index]), np.log10(pressures[index + 1])
    fraction = (0.0 - before) / (after - before)
    return float(10 ** (low + (high - low) * fraction))


def _longest_run(mask):
    """Start and end indices of the longest contiguous True run, or None."""
    best = current = None
    for i, flag in enumerate(mask):
        if flag:
            current = (i, i) if current is None else (current[0], i)
            if best is None or (current[1] - current[0]) > (best[1] - best[0]):
                best = current
        else:
            current = None
    return best


def _bounds(pressures, producer, consumer, mask):
    """Window edges, refined off the grid.

    At each edge one curve is the binding constraint -- whichever has just
    crossed zero -- so the edge is that curve's crossing rather than the grid
    point next to it. Without this the reported width depends on how finely
    the axis was sampled.
    """
    run = _longest_run(mask)
    if run is None:
        return None, None
    start, end = run

    low = float(pressures[start])
    if start > 0:
        # Whichever partner was failing just below the window sets this edge.
        blocking = consumer if consumer[start - 1] >= 0 else producer
        low = _zero_crossing(pressures, blocking, start - 1)

    high = float(pressures[end])
    if end < len(pressures) - 1:
        blocking = producer if producer[end + 1] >= 0 else consumer
        high = _zero_crossing(pressures, blocking, end)

    return low, high


def syntrophy_window(
    producer,
    consumer,
    carrier: str = DEFAULT_CARRIER,
    low: float = 1e-10,
    high: float = 1.0,
    points: int = 61,
) -> SyntrophyWindow:
    """Evaluate both partners across ``carrier`` partial pressure.

    ``producer`` releases the carrier and is inhibited by it; ``consumer``
    takes it up and needs enough of it. Both reactions are used at whatever
    conditions they already carry, with only the carrier pressure varied.
    """
    if producer.n_electrons != consumer.n_electrons:
        raise MicrobialThermoError(
            f"the partners are normalised differently -- producer is per "
            f"{producer.n_electrons} e- and consumer per {consumer.n_electrons} e-, "
            "so their free energies are not comparable on one axis. Rebuild both "
            "with the same n_electrons."
        )

    axis = partial_pressure_axis(carrier, low=low, high=high, n=points)
    produced = sweep(producer, axis)
    consumed = sweep(consumer, axis)

    mask = (produced.delta_g < 0) & (consumed.delta_g < 0)
    window_low, window_high = _bounds(axis.values, produced.delta_g, consumed.delta_g, mask)
    return SyntrophyWindow(
        carrier=carrier,
        pressures=axis.values,
        producer_delta_g=produced.delta_g,
        consumer_delta_g=consumed.delta_g,
        n_electrons=float(producer.n_electrons),
        low=window_low,
        high=window_high,
    )


def _footnotes(*reactions):
    """One combined unverified-data notice for both partners."""
    seen, lines = set(), []
    for built in reactions:
        text = unverified_footnote(built)
        if text is not None and text not in seen:
            seen.add(text)
            lines.append(text)
    return "\n".join(lines) if lines else None


def plot_syntrophy_window(
    producer,
    consumer,
    carrier: str = DEFAULT_CARRIER,
    producer_label: str = "producer — needs the carrier low",
    consumer_label: str = "consumer — needs the carrier high",
    low: float = 1e-10,
    high: float = 1.0,
    points: int = 61,
    figsize=(9.0, 5.5),
    energy_quantum=None,
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    title: str | None = None,
    annotate_best: bool = True,
    window: SyntrophyWindow | None = None,
):
    """Draw both partners against carrier pressure, shading where both pay.

    Returns ``(figure, window)`` so the numbers behind the picture are
    available without recomputing the sweep.
    """
    import matplotlib.pyplot as plt

    if window is None:
        window = syntrophy_window(producer, consumer, carrier, low, high, points)
    quantum = as_magnitude(energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, KJ_PER_MOL_STR)

    figure, ax = plt.subplots(figsize=figsize)

    if window.exists:
        ax.axvspan(
            window.low,
            window.high,
            color=PALETTE["exergonic"],
            alpha=0.13,
            zorder=0,
            label=f"both exergonic ({window.decades:.1f} decades)",
        )

    ax.semilogx(
        window.pressures,
        window.producer_delta_g,
        color=PALETTE["oxidation"],
        linewidth=2.2,
        label=producer_label,
        zorder=4,
    )
    ax.semilogx(
        window.pressures,
        window.consumer_delta_g,
        color=PALETTE["reduction"],
        linewidth=2.2,
        label=consumer_label,
        zorder=4,
    )

    # Zero line plus the band where a reaction is exergonic but too weakly so
    # to be worth anything to a cell.
    add_energy_quantum_band(ax, energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM)

    if window.exists:
        for edge in (window.low, window.high):
            ax.axvline(
                edge,
                color=PALETTE["exergonic"],
                linewidth=1.0,
                linestyle=":",
                zorder=3,
            )

    if annotate_best and window.exists:
        pressure, shared = window.best_shared
        ax.plot(
            [pressure],
            [shared],
            marker="o",
            markersize=6,
            color=PALETTE["annotation"],
            zorder=6,
        )
        inside_quantum = " — inside the quantum" if shared > quantum else ""
        ax.annotate(
            f"best either partner can be guaranteed:\n"
            f"{shared:.1f} kJ/mol at {pressure:.1e} bar{inside_quantum}",
            xy=(pressure, shared),
            # Above the crossing: on an X-shaped plot that wedge is the one
            # reliably empty, whatever the partners are.
            xytext=(0, 58),
            textcoords="offset points",
            fontsize=SIZES["annotation"],
            color=PALETTE["annotation"],
            ha="center",
            va="bottom",
            arrowprops={
                "arrowstyle": "-",
                "color": PALETTE["muted"],
                "linewidth": 0.8,
            },
            zorder=6,
        )

    ax.set_xlabel(_carrier_label(window.carrier), fontsize=SIZES["potential"])
    ax.set_ylabel(
        f"ΔG (kJ/mol per {window.n_electrons:g} e$^-$)",
        fontsize=SIZES["potential"],
    )
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(
        fontsize=SIZES["annotation"],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.13),
        ncol=2,
    )

    if title is None:
        title = "The syntrophic window"
        if window.exists:
            title += (
                f"\n{window.low:.1e} to {window.high:.1e} bar — {window.decades:.1f} decades wide"
            )
        else:
            title += " — none at these conditions"
    ax.set_title(title, fontsize=SIZES["title"], pad=12)

    notice = _footnotes(producer, consumer)
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
    return figure, window
