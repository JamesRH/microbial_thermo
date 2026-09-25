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
# # 5. Carbon
#
# Carbon is the element every other chapter is about, because it is the
# electron donor for almost everything. Its redox range is the widest of any
# biological element — methane at −IV to carbon dioxide at +IV — and the whole
# of heterotrophic metabolism lives somewhere on that span.
#
# This notebook is the template for the eight that follow: the ladder, what is
# stable, where each form lives, and which organisms make a living on each
# step. Every number is computed here, not quoted.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.figures import (
    frost_diagram,
    latimer_diagram,
    plot_frost,
    plot_latimer,
    plot_pourbaix,
)

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder
#
# A Latimer diagram is the element's oxidation states in a row with the
# potential of each step on the arrow. Ours is computed at whatever pH you ask
# for, which matters more for carbon than for most elements.

# %%
ladder = latimer_diagram("C", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("C", diagram=ladder)
plt.show()

# %% [markdown]
# Two things in that table are worth stopping on.
#
# The C(0) rung is **graphite**, and the species it beat is **acetate** —
# both are carbon at oxidation state zero. The library picks whichever has the
# lower free energy at the working pH and keeps the other in `.alternatives`:

# %%
ladder.alternatives

# %% [markdown]
# For a biologist acetate is the interesting C(0) species, and saying so is a
# one-line change. The zero of the diagram stays graphite either way, because
# the element reference is chosen independently of which species you plot:

# %%
biological = latimer_diagram(
    "C",
    species=["HCO3-", "CO3--", "Formate(aq)", "Acetate", "Methane(aq)"],
    pH=7.0,
)
biological.to_frame()[["oxidized", "reduced", "n electrons", "E (V)"]]

# %% [markdown]
# ## pH moves it, and not gently

# %%
for ph in (0.0, 4.0, 7.0, 9.0):
    steps = latimer_diagram("C", pH=ph).steps
    chain = "  ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps)
    print(f"pH {ph:4.1f}   {chain}")

# %% [markdown]
# The top rung changes species — CO₂(aq) below the first p*K*a, bicarbonate
# above it — and every potential drops as protons become scarcer. Carbon
# dioxide reduction is easier in acid, which is one reason methanogens are not
# fussy about it and acetogens are.

# %% [markdown]
# ## What is stable: the Frost diagram
#
# Volt equivalent against oxidation state. The slope between two points is the
# potential of that couple; a point above the line joining two others is
# unstable and will disproportionate; the lowest point is where carbon ends up.

# %%
plot_frost("C", pH=7.0, annotate_slopes=True)
plt.show()

# %%
frost = frost_diagram("C", pH=7.0)
print("most stable form:", frost.most_stable.backend)
for event in frost.disproportionation():
    print(event)

# %% [markdown]
# **Formate is thermodynamically doomed and completely ordinary.** The diagram
# says it should fall apart into bicarbonate and elemental carbon, releasing
# about 41 kJ per mole. Nothing of the sort happens in a bottle of formate on
# the bench, because there is no mechanism: graphite does not nucleate out of
# aqueous solution at 25 °C.
#
# That is the most useful thing a Frost diagram teaches, and it is worth
# stating plainly: **it is a statement about direction, not about rate.** Every
# metastable intermediate in biogeochemistry — formate, nitrite, elemental
# sulfur — is a species sitting above the hull with no way down.

# %% [markdown]
# ### Showing every form, not just the one that predominates
#
# By default each oxidation state contributes one point — the form with the
# lowest energy at the working pH, which is the one that actually exists.
# `predominant_only=False` keeps them all, and then two species can sit at the
# same x: graphite and acetate at C(0), and all three carbonate forms at
# C(+IV).

# %%
everything = frost_diagram("C", pH=7.0, predominant_only=False)
everything.to_frame()[["species", "oxidation state", "volt equivalent (V)", "on hull"]]

# %%
plot_frost("C", pH=7.0, predominant_only=False, label="all")
plt.show()

# %% [markdown]
# Two options control what is drawn, and they are independent:
#
# * `show="all"` or `show="predominant"` — which points get a marker.
# * `label="all"` or `label="predominant"` — which of them get named.
#
# The default is `show="all", label="predominant"`, which keeps the figure
# readable while hiding nothing. Compare:

# %%
figure, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
plot_frost("C", pH=7.0, predominant_only=False, label="all", ax=axes[0], title="label='all'")
plot_frost("C", pH=7.0, predominant_only=False, ax=axes[1], title="label='predominant' (default)")
plt.show()

# %% [markdown]
# **Predominant and stable are different words for different things**, and
# the open grey markers are where the difference shows.
#
# *Predominant* is a comparison **within** one oxidation state — which of the
# three carbonate forms, decided entirely by pH. *Stable* is a comparison
# **across** states — whether that carbonate survives at all, decided by the
# convex hull. An open grey square is a species that is neither: a real form
# of its state that some other form of the same state outcompetes at this pH.
# It is not falling apart, so it is not drawn in red.
#
# The distinction has teeth. Ask for the acetate/graphite pair directly:

# %%
by_name = {point.backend: point for point in everything.points}
for name in ("Graphite", "Acetate"):
    point = by_name[name]
    print(f"{name:10s} volt equivalent {point.volt_equivalent:+.3f} V, "
          f"ΔG {point.gibbs:+7.1f} kJ per mol C")

difference = by_name["Acetate"].gibbs - by_name["Graphite"].gibbs
print(f"\nacetate minus graphite: {difference:+.1f} kJ per mol C")

# %% [markdown]
# That difference is a real free energy and it is **not** a potential. Both
# species are C(0), so converting one into the other transfers *no* electrons,
# and `E = ΔG/nF` would divide by zero. The library refuses rather than
# printing something:

# %%
try:
    everything.slope("Acetate", "Graphite")
except ValueError as exc:
    print("refused:", exc)

# %% [markdown]
# So: distances **across** the diagram are potentials, distances **up and
# down within one oxidation state** are free energies, and the two are not
# interchangeable. That is worth saying out loud, because a Frost diagram
# invites you to read every gap as a voltage.

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("C", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# The carbonate boundaries are vertical, because they are acid–base equilibria
# and do not involve electrons at all; the methane boundary slopes, because it
# does. The dashed lines are the stability limits of water, computed from this
# library's own couples.

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
carbon = [
    "hydrogenotrophic_methanogenesis",
    "acetoclastic_methanogenesis",
    "acetogenesis",
    "aerobic_acetate_oxidation",
    "aerobic_methanotrophy",
    "anaerobic_methane_oxidation",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[carbon][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# Read that against the ladder. CO₂ → CH₄ is the bottom two steps run
# backwards, and it pays about as badly as anything a cell does — which is why
# methanogens are slow, live at the end of the electron chain, and are the
# organisms most often found sitting within a few kJ of the thermodynamic
# limit. Aerobic acetate oxidation runs the same carbon uphill in the other
# direction against oxygen and pays roughly twenty times better.

# %% [markdown]
# ## One thing worth taking away
#
# The nominal oxidation state of carbon in an organic molecule is the single
# most useful number about it, and the library computes it directly:

# %%
for name in ("Methane(aq)", "Acetate", "Formate(aq)", "HCO3-"):
    formula = mt.resolve(name).formula
    state = mt.mean_oxidation_state("C", formula)
    # An exact Fraction, not a float: acetate is 0 and the C(0) rung means it.
    print(f"{name:14s} {formula:8s}  mean C oxidation state {float(state):+.2f}  ({state})")

# %% [markdown]
# That is the x axis of the Frost diagram. A substrate's position on it tells
# you how many electrons it has to give, and the Frost curve tells you what
# they are worth.
