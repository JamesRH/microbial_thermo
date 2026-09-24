"""The interactive free-energy explorer.

Free energy on the y axis; the x-axis variable is chosen from a dropdown, in
the manner of a Gapminder chart. Built with Plotly because the whole figure --
dropdown included -- serialises into a single self-contained HTML file that
works with no live kernel, which is what makes it something you can hand to
students.

Two deliberate choices:

* Every curve is precomputed on a grid before the figure is built. pyGCC is
  far too slow to call per interaction, and a precomputed grid is also the only
  thing that survives export to static HTML.
* SVG export uses Plotly's own modebar button, which is handled client-side by
  plotly.js. The alternative, ``fig.write_image``, needs kaleido, and kaleido
  1.x requires a system Chrome installation -- an external, non-Python
  dependency that would break reproducibility on student machines.

Here the y axis really is kJ/mol, so ATP equivalents are shown as a genuine
second y axis (unlike the redox tower, where the axis is a potential).
"""

from __future__ import annotations

from pathlib import Path

from ..sweep import default_axes, sweep
from ..units import (
    DEFAULT_BIOLOGICAL_ENERGY_QUANTUM,
    DEFAULT_DELTA_G_ATP,
    KJ_PER_MOL_STR,
    as_magnitude,
)
from .style import PALETTE

#: Pinned so a re-export is byte-identical when nothing has changed. Plotly
#: stamps a fresh random div id otherwise, which left every test run showing a
#: one-line diff in a tracked file that meant nothing.
EXPLORER_DIV_ID = "energy-explorer"

#: Passed to Plotly so the modebar's download button yields SVG.
SVG_CONFIG = {
    "toImageButtonOptions": {"format": "svg", "filename": "free_energy"},
    "displaylogo": False,
}


def plot_energy_explorer(
    reaction,
    axes=None,
    per_electron: bool = False,
    delta_g_atp=None,
    energy_quantum=None,
    save_html: str | Path | None = None,
    title: str | None = None,
):
    """Build the interactive explorer for ``reaction``.

    ``axes`` defaults to pH, temperature, and every participating gas or
    aqueous species. Returns a Plotly figure; call ``fig.show(config=...)``
    with :data:`SVG_CONFIG` to get the SVG download button, or pass
    ``save_html`` to write a standalone file that already has it.
    """
    import plotly.graph_objects as go

    axes = axes or default_axes(reaction)
    delta_g_atp = delta_g_atp or DEFAULT_DELTA_G_ATP
    energy_quantum = energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM
    per_atp = as_magnitude(delta_g_atp, KJ_PER_MOL_STR)
    quantum = as_magnitude(energy_quantum, KJ_PER_MOL_STR)

    results = [sweep(reaction, axis) for axis in axes]

    figure = go.Figure()
    for index, result in enumerate(results):
        values = result.delta_g_per_electron if per_electron else result.delta_g
        figure.add_trace(
            go.Scatter(
                x=result.axis.values,
                y=values,
                mode="lines+markers",
                name=result.axis.label,
                visible=(index == 0),
                line={"color": PALETTE["reduction"], "width": 2.5},
                marker={"size": 5},
                hovertemplate=(
                    f"{result.axis.label}: %{{x}}<br>ΔG: %{{y:.1f}} kJ/mol<extra></extra>"
                ),
            )
        )

    unit = "kJ/mol e⁻" if per_electron else "kJ/mol"
    y_title = f"ΔG ({unit})"

    figure.update_layout(
        template="simple_white",
        title=title or _default_title(reaction),
        xaxis={
            "title": {"text": results[0].axis.label},
            "type": "log" if results[0].axis.log_scale else "linear",
        },
        yaxis={"title": {"text": y_title}},
        updatemenus=[_dropdown(results)],
        margin={"l": 70, "r": 90, "t": 110, "b": 60},
        hovermode="x unified",
        showlegend=False,
    )

    _add_reference_lines(figure, quantum, per_electron, results)
    _add_atp_axis(figure, results, per_electron, per_atp)

    if save_html is not None:
        path = Path(save_html).with_suffix(".html")
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.write_html(path, include_plotlyjs=True, config=SVG_CONFIG, div_id=EXPLORER_DIV_ID)
    return figure


def _dropdown(results):
    """Gapminder-style variable selector.

    Each entry toggles trace visibility and retitles the x axis, switching to
    a log scale for the concentration and partial-pressure axes.
    """
    buttons = []
    for index, result in enumerate(results):
        visibility = [i == index for i in range(len(results))]
        buttons.append(
            {
                "label": result.axis.label,
                "method": "update",
                "args": [
                    {"visible": visibility},
                    {
                        "xaxis.title.text": result.axis.label,
                        "xaxis.type": "log" if result.axis.log_scale else "linear",
                    },
                ],
            }
        )
    return {
        "buttons": buttons,
        "direction": "down",
        "showactive": True,
        "x": 0.0,
        "xanchor": "left",
        "y": 1.14,
        "yanchor": "top",
        "pad": {"r": 10, "t": 10},
    }


def _add_reference_lines(figure, quantum, per_electron, results):
    """Zero line and the biological-energy-quantum band."""
    limit = quantum / results[0].n_electrons if per_electron else quantum
    figure.add_hline(y=0.0, line={"color": PALETTE["annotation"], "width": 1})
    figure.add_hrect(
        y0=min(0.0, limit),
        y1=max(0.0, limit),
        fillcolor=PALETTE["quantum_band"],
        opacity=0.18,
        line_width=0,
        annotation_text="below the biological energy quantum",
        annotation_position="top left",
        annotation={"font": {"size": 10, "color": PALETTE["muted"]}},
    )


def _add_atp_axis(figure, results, per_electron, per_atp):
    """A genuine second y axis in ATP units.

    Valid here because the primary axis is a free energy, so the conversion is
    a single division and does not depend on the electron count.
    """
    import numpy as np

    values = np.concatenate(
        [r.delta_g_per_electron if per_electron else r.delta_g for r in results]
    )
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return
    low, high = float(np.min(finite)), float(np.max(finite))
    pad = max(5.0, 0.1 * (high - low))
    low, high = low - pad, high + pad

    figure.update_layout(
        yaxis={"range": [low, high]},
        yaxis2={
            "title": {"text": "ATP equivalents", "font": {"size": 11}},
            "overlaying": "y",
            "side": "right",
            "range": [low / per_atp, high / per_atp],
            "showgrid": False,
            "tickfont": {"size": 10},
        },
    )
    # An invisible trace anchors the second axis; Plotly will not render an
    # overlaying axis that nothing is attached to.
    import plotly.graph_objects as go

    figure.add_trace(
        go.Scatter(
            x=[results[0].axis.values[0]],
            y=[low / per_atp],
            yaxis="y2",
            mode="markers",
            marker={"opacity": 0},
            hoverinfo="skip",
            showlegend=False,
        )
    )


def _default_title(reaction) -> str:
    conditions = reaction.conditions
    return (
        f"{reaction.format()}<br>"
        f"<sub>{reaction.n_electrons} e⁻ · pH {conditions.pH:g} · "
        f"{conditions.temperature_c:g} °C · "
        f"activity model: {conditions.activity_model}</sub>"
    )
