"""The redox tower with live sliders.

The static tower in :mod:`~microbial_thermo.figures.tower` draws one set of
conditions. This one lets pH, temperature and the oxidised/reduced activity
ratio be moved, which is the only way to see the thing students most often get
wrong: the tower is not a fixed table. Couples move at different rates, and
they can cross. A tower drawn at pH 7 and one drawn at pH 9 do not always put
the acceptors in the same order.

Two front ends over the same precomputed
:class:`~microbial_thermo.tower.TowerGrid`:

* :func:`interactive_tower` -- ipywidgets, for use in a live notebook.
* :func:`plot_interactive_tower` -- Plotly, exportable to a standalone HTML
  file that keeps working with no kernel behind it.

Neither calls the backend while a slider moves. pH and temperature come from
the grid; the activity ratio is a closed-form Nernst shift applied on lookup,
so it is exact rather than interpolated and costs nothing.

Plotly's native sliders cannot express three *independent* dimensions -- each
slider's steps would have to know the other two sliders' positions, and the
slider state is not available to them. The HTML export therefore builds its own
sliders and wires them up with a short script, which is also why the whole grid
is serialised into the page.
"""

from __future__ import annotations

import json
from pathlib import Path

from ..tower import TowerGrid, tower_grid
from .style import PALETTE, SIZES
from .tower import _GROUP_COLOURS, _stagger

#: Ratio slider bounds, as powers of ten. Ten orders either way covers
#: everything from a pristine aquifer to a laboratory reagent bottle.
DEFAULT_LOG_RATIO_RANGE = (-6.0, 6.0)

#: Minimum vertical gap between two labels, in volts. Below this they are
#: nudged apart and given leader lines back to the rule.
LABEL_SEPARATION_V = 0.035

SVG_CONFIG = {
    "toImageButtonOptions": {"format": "svg", "filename": "redox_tower"},
    "displaylogo": False,
}


def _ordered(grid: TowerGrid, pH: float, temperature_c: float, ratio: float):
    """Couples at these conditions, weakest reductant first."""
    return grid.entries_at(pH, temperature_c, ratio)


def _subtitle(pH: float, temperature_c: float, log_ratio: float) -> str:
    kind = "E°′" if True else "E°"
    ratio_text = "1:1" if abs(log_ratio) < 1e-9 else f"10^{log_ratio:g} : 1"
    return f"pH {pH:g} · {temperature_c:g} °C · oxidised:reduced {ratio_text} · {kind} in volts"


# --------------------------------------------------------------------------
# ipywidgets front end
# --------------------------------------------------------------------------


def interactive_tower(
    grid: TowerGrid | None = None,
    figsize=(9.0, 7.5),
    log_ratio_range=DEFAULT_LOG_RATIO_RANGE,
    initial_ph: float = 7.0,
    initial_temperature: float = 25.0,
):
    """A matplotlib tower under three ipywidgets sliders.

    Returns the widget box; display it, or let a notebook cell return it.
    Builds the grid if one is not supplied, which costs a few seconds per
    temperature -- pass a prepared grid to avoid paying that twice.
    """
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    from IPython.display import display

    grid = grid if grid is not None else tower_grid()

    ph_slider = widgets.SelectionSlider(
        options=[(f"{v:g}", v) for v in grid.ph_values],
        value=min(grid.ph_values, key=lambda v: abs(v - initial_ph)),
        description="pH",
        continuous_update=False,
        style={"description_width": "initial"},
    )
    temperature_slider = widgets.SelectionSlider(
        options=[(f"{v:g} °C", v) for v in grid.temperature_values],
        value=min(grid.temperature_values, key=lambda v: abs(v - initial_temperature)),
        description="temperature",
        continuous_update=False,
        style={"description_width": "initial"},
    )
    ratio_slider = widgets.FloatSlider(
        value=0.0,
        min=log_ratio_range[0],
        max=log_ratio_range[1],
        step=0.5,
        description="log₁₀(ox/red)",
        continuous_update=False,
        readout_format=".1f",
        style={"description_width": "initial"},
    )

    output = widgets.Output()

    def redraw(*_):
        with output:
            output.clear_output(wait=True)
            figure = _draw_static(
                grid,
                ph_slider.value,
                temperature_slider.value,
                ratio_slider.value,
                figsize=figsize,
            )
            plt.show()
            plt.close(figure)

    for slider in (ph_slider, temperature_slider, ratio_slider):
        slider.observe(redraw, names="value")

    redraw()
    box = widgets.VBox(
        [
            widgets.HBox([ph_slider, temperature_slider]),
            ratio_slider,
            output,
        ]
    )
    display(box)
    return box


def _draw_static(grid, pH, temperature_c, log_ratio, figsize=(9.0, 7.5)):
    """One matplotlib frame of the tower at fixed conditions."""
    import matplotlib.pyplot as plt

    ratio = 10.0**log_ratio
    rows = _ordered(grid, pH, temperature_c, ratio)

    figure, ax = plt.subplots(figsize=figsize)

    # Rules sit at the true potential; labels are nudged apart and joined back
    # to their rule by a leader, so a crowded region stays readable without
    # misreporting where a couple actually is.
    label_positions = _stagger([r[3] for r in rows], LABEL_SEPARATION_V)
    for (label, group, n, volts), text_y in zip(rows, label_positions, strict=True):
        colour = _GROUP_COLOURS.get(group, PALETTE["muted"])
        ax.hlines(volts, 0.06, 0.52, color=colour, linewidth=2.6)
        if abs(text_y - volts) > 1e-9:
            ax.plot(
                [0.52, 0.545],
                [volts, text_y],
                color=colour,
                linewidth=0.7,
                alpha=0.6,
            )
        ax.text(
            0.55,
            text_y,
            f"{label}  ({n:g} e⁻)",
            va="center",
            fontsize=SIZES["annotation"] + 1,
            color=colour,
        )
        ax.text(
            0.04,
            text_y,
            f"{volts:+.3f}",
            va="center",
            ha="right",
            fontsize=SIZES["annotation"],
            color=PALETTE["annotation"],
        )

    ax.axhline(0.0, color=PALETTE["muted"], linewidth=0.8, linestyle=":")
    ax.set_xlim(-0.06, 1.25)
    ax.set_ylim(-0.75, 1.05)
    ax.set_ylabel("E (V) vs SHE", fontsize=SIZES["potential"])
    ax.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title(
        f"Redox tower\n{_subtitle(pH, temperature_c, log_ratio)}",
        fontsize=SIZES["title"],
        pad=10,
    )
    # Down the page is the direction electrons fall.
    ax.annotate(
        "electrons fall this way",
        xy=(1.16, 0.85),
        xytext=(1.16, -0.55),
        fontsize=SIZES["annotation"],
        color=PALETTE["muted"],
        rotation=90,
        ha="center",
        va="bottom",
        arrowprops={"arrowstyle": "<-", "color": PALETTE["muted"], "linewidth": 0.9},
    )
    figure.tight_layout()
    return figure


# --------------------------------------------------------------------------
# Plotly front end, exportable with no kernel
# --------------------------------------------------------------------------


def plot_interactive_tower(
    grid: TowerGrid | None = None,
    save_html: str | Path | None = None,
    log_ratio_range=DEFAULT_LOG_RATIO_RANGE,
    initial_ph: float = 7.0,
    initial_temperature: float = 25.0,
    title: str = "Redox tower",
):
    """Build the Plotly tower, with sliders that survive export to HTML.

    Returns the Plotly figure. Pass ``save_html`` to write the standalone
    page; the sliders are injected there by a post-script, so they exist in
    the written file but not in the returned figure object.
    """
    import plotly.graph_objects as go

    grid = grid if grid is not None else tower_grid()
    rows = _ordered(grid, initial_ph, initial_temperature, 1.0)

    figure = go.Figure()
    for group in sorted({r[1] for r in rows}):
        members = [r for r in rows if r[1] == group]
        figure.add_trace(
            go.Scatter(
                x=[r[3] for r in members],
                y=[r[0] for r in members],
                mode="markers",
                name=group,
                marker={
                    "size": 13,
                    "symbol": "line-ns",
                    "line": {
                        "width": 3,
                        "color": _GROUP_COLOURS.get(group, PALETTE["muted"]),
                    },
                },
                hovertemplate="%{y}<br>E = %{x:.3f} V<extra></extra>",
            )
        )

    figure.update_layout(
        template="simple_white",
        title={"text": f"{title}<br><sub>{_subtitle(initial_ph, initial_temperature, 0.0)}</sub>"},
        xaxis={"title": {"text": "E (V) vs SHE"}, "range": [-0.75, 1.05]},
        yaxis={
            "title": {"text": ""},
            "categoryorder": "array",
            "categoryarray": [r[0] for r in rows],
        },
        margin={"l": 190, "r": 40, "t": 110, "b": 60},
        legend={"orientation": "h", "y": -0.16},
        height=620,
    )
    figure.add_vline(x=0.0, line={"color": PALETTE["muted"], "width": 1, "dash": "dot"})

    if save_html is not None:
        path = Path(save_html).with_suffix(".html")
        if path.parent != Path(""):
            path.parent.mkdir(parents=True, exist_ok=True)
        figure.write_html(
            path,
            include_plotlyjs=True,
            config=SVG_CONFIG,
            div_id="tower-plot",
            post_script=_slider_script(grid, log_ratio_range, title),
        )
    return figure


def _grid_payload(grid: TowerGrid) -> dict:
    """The whole grid as plain JSON, for the exported page to index into."""
    import numpy as np

    potentials = np.asarray(grid.potentials)
    return {
        "ph": list(grid.ph_values),
        "temperatures": list(grid.temperature_values),
        "labels": list(grid.labels),
        "groups": list(grid.groups),
        "nElectrons": list(grid.n_electrons),
        # None rather than NaN: NaN is not valid JSON and json.dumps would
        # emit a literal that JSON.parse refuses.
        "potentials": [
            [[None if np.isnan(v) else float(v) for v in cell] for cell in row]
            for row in potentials
        ],
        "colours": [_GROUP_COLOURS.get(g, PALETTE["muted"]) for g in grid.groups],
    }


def _slider_script(grid: TowerGrid, log_ratio_range, title: str) -> str:
    """JavaScript injected into the exported page.

    Builds three range inputs, then on every move recomputes the tower from
    the embedded grid and restyles the plot. The Nernst shift is redone in
    JavaScript rather than precomputed, which is what keeps the ratio
    continuous instead of quantised to a handful of steps.
    """
    payload = json.dumps(_grid_payload(grid))
    low, high = log_ratio_range
    return f"""
var GRID = {payload};
var TITLE = {json.dumps(title)};
var LOG_LO = {low}, LOG_HI = {high};
var plot = document.getElementById('tower-plot');

function control(labelText, min, max, step, value, format) {{
  var wrap = document.createElement('div');
  wrap.style.cssText = 'display:flex;align-items:center;gap:10px;margin:6px 0;'
    + 'font:13px system-ui,sans-serif;';
  var name = document.createElement('label');
  name.textContent = labelText;
  name.style.cssText = 'width:130px;text-align:right;color:#333;';
  var input = document.createElement('input');
  input.type = 'range';
  input.min = min; input.max = max; input.step = step; input.value = value;
  input.style.cssText = 'flex:1;max-width:360px;';
  var readout = document.createElement('span');
  readout.style.cssText = 'width:90px;color:#555;font-variant-numeric:tabular-nums;';
  readout.textContent = format(parseFloat(value));
  input.addEventListener('input', function () {{
    readout.textContent = format(parseFloat(input.value));
    redraw();
  }});
  wrap.appendChild(name); wrap.appendChild(input); wrap.appendChild(readout);
  return {{ wrap: wrap, input: input }};
}}

var panel = document.createElement('div');
panel.style.cssText = 'max-width:760px;margin:12px auto 4px auto;';
plot.parentNode.insertBefore(panel, plot);

// pH and temperature are indices into the grid; the ratio is continuous.
var phControl = control('pH', 0, GRID.ph.length - 1, 1,
  nearest(GRID.ph, 7.0), function (i) {{ return GRID.ph[i].toFixed(1); }});
var tempControl = control('temperature', 0, GRID.temperatures.length - 1, 1,
  nearest(GRID.temperatures, 25.0),
  function (i) {{ return GRID.temperatures[i].toFixed(0) + ' \\u00B0C'; }});
var ratioControl = control('log\\u2081\\u2080(ox/red)', LOG_LO, LOG_HI, 0.1, 0,
  function (v) {{ return v.toFixed(1); }});

panel.appendChild(phControl.wrap);
panel.appendChild(tempControl.wrap);
panel.appendChild(ratioControl.wrap);

function nearest(values, target) {{
  var best = 0;
  for (var i = 1; i < values.length; i++) {{
    if (Math.abs(values[i] - target) < Math.abs(values[best] - target)) best = i;
  }}
  return best;
}}

function potentials(phIndex, tempIndex, logRatio) {{
  var cell = GRID.potentials[phIndex][tempIndex];
  // RT/F in volts, with T in kelvin; ln(10) converts the decadic slider.
  var rtOverF = 8.314462618 * (GRID.temperatures[tempIndex] + 273.15) / 96485.332;
  var shift = rtOverF * Math.LN10 * logRatio;
  var out = [];
  for (var i = 0; i < cell.length; i++) {{
    out.push(cell[i] === null ? null : cell[i] + shift / GRID.nElectrons[i]);
  }}
  return out;
}}

function redraw() {{
  var phIndex = parseInt(phControl.input.value, 10);
  var tempIndex = parseInt(tempControl.input.value, 10);
  var logRatio = parseFloat(ratioControl.input.value);
  var volts = potentials(phIndex, tempIndex, logRatio);

  var rows = [];
  for (var i = 0; i < GRID.labels.length; i++) {{
    if (volts[i] === null) continue;
    rows.push({{ label: GRID.labels[i], group: GRID.groups[i], v: volts[i] }});
  }}
  rows.sort(function (a, b) {{ return a.v - b.v; }});

  var groups = plot.data.map(function (trace) {{ return trace.name; }});
  var xs = [], ys = [];
  groups.forEach(function (group) {{
    var members = rows.filter(function (r) {{ return r.group === group; }});
    xs.push(members.map(function (r) {{ return r.v; }}));
    ys.push(members.map(function (r) {{ return r.label; }}));
  }});

  var ratioText = logRatio === 0 ? '1:1'
    : '10^' + logRatio.toFixed(1) + ' : 1';
  Plotly.restyle(plot, {{ x: xs, y: ys }});
  Plotly.relayout(plot, {{
    'yaxis.categoryarray': rows.map(function (r) {{ return r.label; }}),
    'title.text': TITLE + '<br><sub>pH ' + GRID.ph[phIndex].toFixed(1)
      + ' \\u00B7 ' + GRID.temperatures[tempIndex].toFixed(0) + ' \\u00B0C \\u00B7 '
      + 'oxidised:reduced ' + ratioText + ' \\u00B7 E\\u00B0\\u2032 in volts</sub>'
  }});
}}

redraw();
"""
