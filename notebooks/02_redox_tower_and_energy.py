# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: microbial-thermo
#     language: python
#     name: microbial-thermo
# ---

# %% [markdown]
# # 2. The redox tower, and how conditions change the energy
#
# Standard-state values say what is *possible*. This notebook is about the gap
# between that and what actually pays in a real environment.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt

warnings.filterwarnings("ignore", category=UserWarning)

conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

# %% [markdown]
# ## The reference couples
#
# Potentials are computed at the pH and temperature asked for, so this table is
# not a pH 7 snapshot — change the conditions and the numbers move.

# %%
from microbial_thermo.tower import tower_table

tower_table(conditions)

# %% [markdown]
# Compare against pH 5 to see the proton-coupled couples shift while
# Fe(III)/Fe(II) stays put:

# %%
tower_table(mt.Conditions(temperature_c=25.0, pH=5.0))

# %% [markdown]
# ## Drawing the tower
#
# The tower is normalised to an electron pair. That is what makes the energy
# scale bar exact: with *n* fixed, $\Delta G = -nF\Delta E$, so a *difference*
# in potential converts to kJ/mol by the constant $2F$ = 192.97 kJ mol⁻¹ V⁻¹.
#
# Read a reaction's yield straight off the figure by comparing the length of
# the donor-to-acceptor arrow against the scale bar.

# %%
from microbial_thermo.figures import plot_redox_tower

reaction = mt.Reaction.from_couples(
    donor=("H2(g)", "H+"), acceptor=("HS-", "SO4-2"), conditions=conditions
)
plot_redox_tower(reaction)
plt.show()

# %% [markdown]
# Note the scale bar is a *bar*, not an axis. An absolute kJ/mol axis would
# imply each couple has an absolute free energy; only differences between
# couples carry energy.

# %% [markdown]
# ## Standard state is a fiction
#
# $\Delta G^{\circ\prime}$ puts hydrogen at 1 bar. Real hydrogen in a sediment
# is nanomolar. Sweeping it shows where methanogenesis stops paying.

# %%
from microbial_thermo.sweep import partial_pressure_axis, sweep

realistic = mt.Conditions(
    temperature_c=25.0, pH=7.0, activity_model="ideal",
    partial_pressures={"H2(g)": 1e-4},
)
methanogenesis = mt.Reaction.from_couples(
    donor=("H2(g)", "H+"), acceptor=("methane", "CO2(aq)"), conditions=realistic
)

result = sweep(methanogenesis, partial_pressure_axis("H2(g)"))
frame = result.to_frame()
frame.iloc[::4]

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.semilogx(result.axis.values, result.delta_g, marker="o", markersize=4)
ax.axhline(0, color="black", linewidth=1)
ax.axhspan(-20, 0, color="#E8B84B", alpha=0.2)
ax.set_xlabel("$p_{H_2}$ (bar)")
ax.set_ylabel("ΔG (kJ/mol per 2 e$^-$)")
ax.set_title("Hydrogenotrophic methanogenesis against hydrogen partial pressure")
ax.annotate(
    "below the biological energy quantum",
    xy=(1e-9, -10), fontsize=9, color="#8A8A8A",
)
plt.show()

# %% [markdown]
# The curve crosses zero near $10^{-5}$ bar. That threshold is why
# methanogens and their syntrophic partners have to keep hydrogen so low, and
# why the two cannot both be comfortable at once.

# %% [markdown]
# ## Temperature
#
# pyGCC gives real temperature dependence rather than a van 't Hoff
# approximation, over 0.01 to 100 °C.

# %%
from microbial_thermo.sweep import temperature_axis

warm = sweep(methanogenesis, temperature_axis())
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(warm.axis.values, warm.delta_g)
ax.set_xlabel("temperature (°C)")
ax.set_ylabel("ΔG (kJ/mol per 2 e$^-$)")
ax.set_title("Temperature dependence at fixed composition")
plt.show()

# %% [markdown]
# ## The interactive explorer
#
# Pick the x-axis variable from the dropdown. Writing it to HTML gives a
# self-contained file — dropdown, hover and the SVG download button all work
# with no Python behind it, so it can be handed to students directly.

# %%
from microbial_thermo.figures import plot_energy_explorer
from microbial_thermo.figures.explorer import SVG_CONFIG

figure = plot_energy_explorer(methanogenesis)
figure.show(config=SVG_CONFIG)

# %%
plot_energy_explorer(methanogenesis, save_html="explorer_methanogenesis")
print("wrote explorer_methanogenesis.html")
