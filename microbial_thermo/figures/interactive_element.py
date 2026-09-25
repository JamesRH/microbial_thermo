"""One element's three diagrams, under live sliders.

The static Latimer, Frost and Eh-pH figures each draw one set of conditions.
This module puts all three on one canvas and lets temperature, pH and
dissolved activity be moved, which is the whole argument for computing these
rather than reproducing them from a textbook: a published Frost diagram is a
photograph of pH 0 and 25 C, and the interesting question is what happens
when you leave those conditions.

The three panels answer different questions about the same numbers, and
seeing them move together is the point:

* the **Latimer** chain says what each one-electron-at-a-time step costs;
* the **Frost** curve says which forms survive and which fall apart;
* the **Eh-pH** field says which one is actually present, and the pH slider
  is drawn on it as a vertical line -- so the other two panels are a slice
  through that line.

**What costs anything.** Only temperature: every other axis is closed form.
A pH change is already a coefficient on the decomposition, and an activity
change is ``RT ln(a2/a1)`` per mole of the element on the aqueous species
only. So the grid is built once over temperature (six backend passes) and
every slider move after that is arithmetic on floats already in memory --
including the Eh-pH panel, which is an argmin over a stack of planes.

Two front ends over the same :class:`~microbial_thermo.figures.basis.ElementGrid`:

* :func:`interactive_element` -- ipywidgets and matplotlib, for a live kernel.
* :func:`plot_interactive_element` -- Plotly, exportable to a standalone HTML
  file that keeps working with nothing behind it. The whole grid is serialised
  into the page and the diagrams are rebuilt in JavaScript, for the same
  reason the interactive tower does it: Plotly's own sliders cannot express
  three independent dimensions.
"""

from __future__ import annotations

import json
from pathlib import Path

from .basis import DEFAULT_TEMPERATURES, FARADAY_KJ, element_grid, pretty
from .frost import frost_diagram, plot_frost
from .latimer import latimer_diagram, plot_latimer
from .pourbaix import plot_pourbaix, pourbaix_field
from .style import PALETTE, SIZES

#: Activity slider bounds, as powers of ten. Unit activity is the
#: standard-state convention; 1e-6 is the Pourbaix convention; 1e-9 is a
#: trace metal in a real water.
DEFAULT_LOG_ACTIVITY_RANGE = (-9.0, 0.0)

SVG_CONFIG = {
    "toImageButtonOptions": {"format": "svg", "filename": "element_diagrams"},
    "displaylogo": False,
}

#: Fixed id so an exported page can be rebuilt byte-for-byte.
ELEMENT_DIV_ID = "element-plot"


def _figure(grid, temperature_c, pH, activity, figsize, eh_range, points):
    """Draw all three panels at one set of conditions."""
    import matplotlib.pyplot as plt

    series = grid.at(temperature_c, pH=pH, activity=activity)
    figure = plt.figure(figsize=figsize)
    layout = figure.add_gridspec(2, 2, height_ratios=(1.0, 1.9), hspace=0.42, wspace=0.26)

    ladder_ax = figure.add_subplot(layout[0, :])
    frost_ax = figure.add_subplot(layout[1, 0])
    field_ax = figure.add_subplot(layout[1, 1])

    ladder = latimer_diagram(grid.element, series=series)
    plot_latimer(grid.element, diagram=ladder, ax=ladder_ax)

    frost = frost_diagram(grid.element, series=series)
    plot_frost(grid.element, diagram=frost, ax=frost_ax, title=f"{grid.element} — Frost")

    field = pourbaix_field(grid.element, series=series, points=points, eh_range=eh_range)
    plot_pourbaix(grid.element, field=field, ax=field_ax)
    field_ax.axvline(pH, color=PALETTE["annotation"], linewidth=1.4, linestyle="-", alpha=0.7)
    field_ax.set_title(
        f"{grid.element} predominance — the line is the pH above",
        fontsize=SIZES["title"],
        pad=10,
    )

    # The panel-level warnings that plot_frost would have drawn on a figure
    # it owned: on a shared canvas they belong at the bottom, once.
    notes = ["  ·  ".join(str(event) for event in frost.disproportionation())]
    if not frost.reference_is_element:
        notes.append(
            f"no elemental {grid.element} available, so the volt-equivalent zero "
            f"({frost.reference}) is arbitrary — slopes hold, heights do not"
        )
    if getattr(grid, "dropped", ()):
        notes.append(
            "not on the interactive diagram, being unavailable at some temperature: "
            + ", ".join(name for name, _ in grid.dropped)
        )
    unstable = [note for note in notes if note]
    if unstable:
        figure.text(
            0.5,
            0.005,
            "\n".join(unstable),
            ha="center",
            va="bottom",
            fontsize=SIZES["annotation"] - 1,
            color=PALETTE["endergonic"],
            wrap=True,
        )
    return figure


def interactive_element(
    element: str,
    grid=None,
    species=None,
    temperatures=DEFAULT_TEMPERATURES,
    figsize=(12.5, 9.0),
    initial_ph: float = 7.0,
    initial_temperature: float = 25.0,
    initial_log_activity: float = -6.0,
    log_activity_range=DEFAULT_LOG_ACTIVITY_RANGE,
    eh_range=(-1.0, 1.4),
    points: int = 200,
    fixed=None,
    backend=None,
):
    """Three matplotlib panels under three ipywidgets sliders.

    Returns the widget box; display it, or let a notebook cell return it.
    Builds the grid if one is not supplied, which costs a backend pass per
    temperature -- pass a prepared grid to avoid paying that twice.
    """
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    from IPython.display import display

    grid = (
        grid
        if grid is not None
        else element_grid(
            element, species=species, temperatures=temperatures, fixed=fixed, backend=backend
        )
    )

    temperature_slider = widgets.SelectionSlider(
        options=[(f"{v:g} °C", v) for v in grid.temperatures],
        value=min(grid.temperatures, key=lambda v: abs(v - initial_temperature)),
        description="temperature",
        continuous_update=False,
        style={"description_width": "initial"},
    )
    ph_slider = widgets.FloatSlider(
        value=initial_ph,
        min=0.0,
        max=14.0,
        step=0.1,
        description="pH",
        continuous_update=False,
        readout_format=".1f",
        style={"description_width": "initial"},
    )
    activity_slider = widgets.FloatSlider(
        value=initial_log_activity,
        min=log_activity_range[0],
        max=log_activity_range[1],
        step=0.5,
        description="log₁₀(activity)",
        continuous_update=False,
        readout_format=".1f",
        style={"description_width": "initial"},
    )

    output = widgets.Output()

    def redraw(*_):
        with output:
            output.clear_output(wait=True)
            figure = _figure(
                grid,
                temperature_slider.value,
                ph_slider.value,
                10.0**activity_slider.value,
                figsize,
                eh_range,
                points,
            )
            display(figure)
            plt.close(figure)

    for slider in (temperature_slider, ph_slider, activity_slider):
        slider.observe(redraw, names="value")
    redraw()

    return widgets.VBox([widgets.HBox([temperature_slider, ph_slider, activity_slider]), output])


# --------------------------------------------------------------------------
# Plotly front end, exportable with no kernel
# --------------------------------------------------------------------------


def _grid_payload(grid, eh_range, points: int) -> dict:
    """The whole grid as plain JSON, for the exported page to rebuild from.

    What travels is the decomposition, not the diagrams: base energy, proton
    and electron coefficients per mole of the element, one set per
    temperature. Everything drawn in the page is derived from those four
    numbers per species, which is why the page needs no kernel and no
    thermodynamic data of its own.
    """
    first = grid.series[0]
    return {
        "element": grid.element,
        "temperatures": list(grid.temperatures),
        "activity": grid.activity,
        "reference": first.reference,
        "referenceIsElement": bool(first.reference_is_element),
        "species": [entry.backend for entry in first],
        "labels": [pretty(entry.backend) for entry in first],
        "aqueous": [bool(entry.species.is_aqueous) for entry in first],
        "nElement": [int(entry.n_element) for entry in first],
        "states": [float(entry.oxidation_state) for entry in first],
        "statable": [bool(entry.state_is_conventional) for entry in first],
        "base": [[float(entry.base) for entry in series] for series in grid.series],
        "protons": [float(entry.protons) for entry in first],
        "electrons": [float(entry.electrons) for entry in first],
        "rt": [float(series.rt) for series in grid.series],
        "faraday": FARADAY_KJ,
        "ehRange": list(eh_range),
        "points": int(points),
        "fixed": [[el, name, value] for el, name, value in first.fixed],
    }


def plot_interactive_element(
    element: str,
    grid=None,
    species=None,
    temperatures=DEFAULT_TEMPERATURES,
    save_html: str | Path | None = None,
    initial_ph: float = 7.0,
    initial_temperature: float = 25.0,
    initial_log_activity: float = -6.0,
    log_activity_range=DEFAULT_LOG_ACTIVITY_RANGE,
    eh_range=(-1.0, 1.4),
    points: int = 120,
    fixed=None,
    backend=None,
    title: str | None = None,
):
    """A Plotly Frost curve and Eh-pH field whose sliders survive export.

    Returns the Plotly figure. Pass ``save_html`` to write the standalone
    page; the sliders and the Latimer chain are injected there by a
    post-script, so they exist in the written file but not in the returned
    figure object.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    grid = (
        grid
        if grid is not None
        else element_grid(
            element, species=species, temperatures=temperatures, fixed=fixed, backend=backend
        )
    )
    title = title or f"{element} — redox diagrams"

    series = grid.at(initial_temperature, pH=initial_ph, activity=10.0**initial_log_activity)
    frost = frost_diagram(element, series=series)
    field = pourbaix_field(element, series=series, points=points, eh_range=eh_range)

    figure = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Frost–Ebsworth", "Eh–pH"),
        horizontal_spacing=0.12,
    )

    # Trace 0: every point. Trace 1: the hull. Trace 2: the field. Trace 3:
    # the pH line. The JavaScript restyles them by index, so the order here
    # is load-bearing.
    figure.add_trace(
        go.Scatter(
            x=frost.states,
            y=frost.volt_equivalents,
            mode="markers+text",
            text=[pretty(p.backend) for p in frost.points],
            textposition="top center",
            marker={"size": 9, "color": PALETTE["endergonic"]},
            name="all forms",
            hovertemplate="%{text}<br>state %{x}<br>%{y:.3f} V<extra></extra>",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=[p.oxidation_state for p in frost.stable],
            y=[p.volt_equivalent for p in frost.stable],
            mode="lines+markers",
            line={"color": PALETTE["reduction"], "width": 2},
            marker={"size": 9, "color": PALETTE["reduction"]},
            name="stable forms",
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Heatmap(
            x=field.ph,
            y=field.eh,
            z=field.winner,
            colorscale="Viridis",
            showscale=False,
            hovertemplate="pH %{x:.1f}<br>Eh %{y:.2f} V<br>%{text}<extra></extra>",
            text=[[field.species[int(i)] for i in row] for row in field.winner],
            name="",
        ),
        row=1,
        col=2,
    )
    figure.add_trace(
        go.Scatter(
            x=[initial_ph, initial_ph],
            y=list(eh_range),
            mode="lines",
            line={"color": "white", "width": 2},
            name="pH slider",
            hoverinfo="skip",
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    figure.update_xaxes(title_text=f"oxidation state of {element}", row=1, col=1)
    figure.update_yaxes(title_text="volt equivalent (V)", row=1, col=1)
    figure.update_xaxes(title_text="pH", row=1, col=2, range=[0.0, 14.0])
    figure.update_yaxes(title_text="Eh (V vs SHE)", row=1, col=2, range=list(eh_range))
    figure.update_layout(
        template="simple_white",
        title={"text": title},
        height=560,
        margin={"l": 70, "r": 40, "t": 150, "b": 60},
        legend={"orientation": "h", "y": -0.18},
    )

    if save_html is not None:
        path = Path(save_html).with_suffix(".html")
        if path.parent != Path(""):
            path.parent.mkdir(parents=True, exist_ok=True)
        figure.write_html(
            path,
            include_plotlyjs=True,
            config=SVG_CONFIG,
            div_id=ELEMENT_DIV_ID,
            post_script=_slider_script(
                grid,
                eh_range,
                points,
                log_activity_range,
                title,
                initial_ph,
                initial_temperature,
                initial_log_activity,
            ),
        )
    return figure


def _slider_script(
    grid,
    eh_range,
    points,
    log_activity_range,
    title,
    initial_ph,
    initial_temperature,
    initial_log_activity,
) -> str:
    """JavaScript injected into the exported page.

    Rebuilds all three diagrams from the decomposition on every slider move:
    the energies in closed form, the convex hull by monotone chain, and the
    Eh-pH field by an argmin over the species at each cell. That is a few
    hundred thousand floating-point operations for a 120x120 field, which is
    nothing, and it is what lets the exported page have a temperature slider
    at all.
    """
    payload = json.dumps(_grid_payload(grid, eh_range, points))
    low, high = log_activity_range
    return f"""
var GRID = {payload};
var TITLE = {json.dumps(title)};
var LOG_LO = {low}, LOG_HI = {high};
var LN10 = Math.LN10;
var plot = document.getElementById({json.dumps(ELEMENT_DIV_ID)});

function control(labelText, min, max, step, value, format) {{
  var wrap = document.createElement('div');
  wrap.style.cssText = 'display:flex;align-items:center;gap:10px;margin:6px 0;'
    + 'font:13px system-ui,sans-serif;';
  var name = document.createElement('label');
  name.textContent = labelText;
  name.style.cssText = 'width:140px;text-align:right;color:#333;';
  var input = document.createElement('input');
  input.type = 'range';
  input.min = min; input.max = max; input.step = step; input.value = value;
  input.style.cssText = 'flex:1;max-width:340px;';
  var readout = document.createElement('span');
  readout.style.cssText = 'width:96px;color:#555;font-variant-numeric:tabular-nums;';
  readout.textContent = format(parseFloat(value));
  input.addEventListener('input', function () {{
    readout.textContent = format(parseFloat(input.value));
    redraw();
  }});
  wrap.appendChild(name); wrap.appendChild(input); wrap.appendChild(readout);
  return {{ wrap: wrap, input: input }};
}}

function nearest(values, target) {{
  var best = 0;
  for (var i = 1; i < values.length; i++) {{
    if (Math.abs(values[i] - target) < Math.abs(values[best] - target)) best = i;
  }}
  return best;
}}

var panel = document.createElement('div');
panel.style.cssText = 'max-width:820px;margin:10px auto 2px auto;';
plot.parentNode.insertBefore(panel, plot);

var ladderBox = document.createElement('div');
ladderBox.style.cssText = 'max-width:960px;margin:2px auto 10px auto;text-align:center;'
  + 'font:13px system-ui,sans-serif;color:#333;line-height:1.9;';
plot.parentNode.insertBefore(ladderBox, plot.nextSibling);

var noteBox = document.createElement('div');
noteBox.style.cssText = 'max-width:960px;margin:2px auto 14px auto;text-align:center;'
  + 'font:12px system-ui,sans-serif;color:#A63A50;';
plot.parentNode.insertBefore(noteBox, ladderBox.nextSibling);

var tempControl = control('temperature', 0, GRID.temperatures.length - 1, 1,
  nearest(GRID.temperatures, {initial_temperature}),
  function (i) {{ return GRID.temperatures[i].toFixed(0) + ' \\u00B0C'; }});
var phControl = control('pH', 0, 14, 0.1, {initial_ph},
  function (v) {{ return v.toFixed(1); }});
var activityControl = control('log\\u2081\\u2080(activity)', LOG_LO, LOG_HI, 0.5,
  {initial_log_activity}, function (v) {{ return v.toFixed(1); }});

panel.appendChild(tempControl.wrap);
panel.appendChild(phControl.wrap);
panel.appendChild(activityControl.wrap);

// Free energy per mole of the element, at a temperature index, pH and Eh.
function energies(ti, pH, logActivity, eh) {{
  var rt = GRID.rt[ti];
  var base = GRID.base[ti];
  var shift = rt * LN10 * (logActivity - Math.log10(GRID.activity));
  var out = [];
  for (var i = 0; i < base.length; i++) {{
    var g = base[i] + GRID.protons[i] * rt * LN10 * pH;
    if (GRID.aqueous[i]) g += shift / GRID.nElement[i];
    if (eh !== undefined) g += GRID.electrons[i] * GRID.faraday * eh;
    out.push(g);
  }}
  return out;
}}

// One species per oxidation state: the lowest at this pH, as the static
// diagrams do. Only species whose electron count really is a state.
function rungs(g) {{
  var best = {{}};
  for (var i = 0; i < g.length; i++) {{
    if (!GRID.statable[i]) continue;
    var key = GRID.states[i].toFixed(6);
    if (!(key in best) || g[i] < g[best[key]]) best[key] = i;
  }}
  var picked = Object.keys(best).map(function (k) {{ return best[k]; }});
  picked.sort(function (a, b) {{ return GRID.states[a] - GRID.states[b]; }});
  return picked;
}}

function hull(indices, volts) {{
  var stack = [];
  for (var k = 0; k < indices.length; k++) {{
    var i = indices[k];
    while (stack.length >= 2) {{
      var a = stack[stack.length - 2], b = stack[stack.length - 1];
      var cross = (GRID.states[b] - GRID.states[a]) * (volts[i] - volts[a])
        - (GRID.states[i] - GRID.states[a]) * (volts[b] - volts[a]);
      if (cross <= 0) stack.pop(); else break;
    }}
    stack.push(i);
  }}
  return stack;
}}

function field(ti, logActivity) {{
  var n = GRID.points;
  var lo = GRID.ehRange[0], hi = GRID.ehRange[1];
  var z = [], text = [], phAxis = [], ehAxis = [];
  for (var j = 0; j < n; j++) phAxis.push(14.0 * j / (n - 1));
  for (var i = 0; i < n; i++) ehAxis.push(lo + (hi - lo) * i / (n - 1));
  for (var i = 0; i < n; i++) {{
    var row = [], names = [];
    for (var j = 0; j < n; j++) {{
      var g = energies(ti, phAxis[j], logActivity, ehAxis[i]);
      var best = 0;
      for (var s = 1; s < g.length; s++) if (g[s] < g[best]) best = s;
      row.push(best); names.push(GRID.species[best]);
    }}
    z.push(row); text.push(names);
  }}
  return {{ z: z, text: text, ph: phAxis, eh: ehAxis }};
}}

function redraw() {{
  var ti = parseInt(tempControl.input.value, 10);
  var pH = parseFloat(phControl.input.value);
  var logActivity = parseFloat(activityControl.input.value);

  var g = energies(ti, pH, logActivity);
  var volts = g.map(function (v) {{ return v / GRID.faraday; }});
  var picked = rungs(g);
  var stable = hull(picked, volts);

  var cell = field(ti, logActivity);

  Plotly.restyle(plot, {{
    x: [picked.map(function (i) {{ return GRID.states[i]; }})],
    y: [picked.map(function (i) {{ return volts[i]; }})],
    text: [picked.map(function (i) {{ return GRID.labels[i]; }})]
  }}, [0]);
  Plotly.restyle(plot, {{
    x: [stable.map(function (i) {{ return GRID.states[i]; }})],
    y: [stable.map(function (i) {{ return volts[i]; }})]
  }}, [1]);
  Plotly.restyle(plot, {{ z: [cell.z], text: [cell.text], x: [cell.ph], y: [cell.eh] }}, [2]);
  Plotly.restyle(plot, {{ x: [[pH, pH]], y: [GRID.ehRange] }}, [3]);

  // The Latimer chain, most oxidized first, with each arrow's potential.
  var chain = [];
  for (var k = picked.length - 1; k >= 0; k--) {{
    chain.push('<b>' + GRID.species[picked[k]] + '</b>');
    if (k > 0) {{
      var high = picked[k], low = picked[k - 1];
      var n = GRID.states[high] - GRID.states[low];
      var e = (g[high] - g[low]) / (n * GRID.faraday);
      chain.push('&nbsp;&mdash;<span style="color:#1F6FB4">' + (e >= 0 ? '+' : '')
        + e.toFixed(3) + ' V</span>&rarr;&nbsp;');
    }}
  }}
  ladderBox.innerHTML = chain.join('');

  var offHull = picked.filter(function (i) {{ return stable.indexOf(i) < 0; }});
  var notes = [];
  if (offHull.length) {{
    notes.push('unstable to disproportionation: '
      + offHull.map(function (i) {{ return GRID.species[i]; }}).join(', '));
  }}
  if (!GRID.referenceIsElement) {{
    notes.push('no elemental ' + GRID.element + ' available, so the volt-equivalent zero ('
      + GRID.reference + ') is arbitrary \\u2014 slopes hold, heights do not');
  }}
  if (GRID.fixed.length) {{
    notes.push('held fixed: ' + GRID.fixed.map(function (f) {{
      return f[1] + ' ' + f[2];
    }}).join(', '));
  }}
  noteBox.innerHTML = notes.join(' &nbsp;\\u00B7&nbsp; ');

  Plotly.relayout(plot, {{
    'title.text': TITLE + '<br><sub>pH ' + pH.toFixed(1) + ' \\u00B7 '
      + GRID.temperatures[ti].toFixed(0) + ' \\u00B0C \\u00B7 dissolved activity 10^'
      + logActivity.toFixed(1) + ' \\u00B7 volt-equivalent zero at ' + GRID.reference + '</sub>'
  }});
}}

redraw();
"""
