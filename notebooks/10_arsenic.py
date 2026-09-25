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
# # 10. Arsenic
#
# Arsenic poisons more people than any other element in groundwater, and the
# reason is thermodynamic rather than geological: the arsenic is already in
# the sediment, bound to iron oxides, and **reducing conditions let it go**.
#
# This notebook works out why, and then asks the question that makes arsenic
# worth teaching next to selenium: why does reduction immobilise selenium and
# mobilise arsenic, when the two sit next to each other in the periodic table
# and have superficially similar redox ladders?

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

import microbial_thermo as mt
from microbial_thermo.figures import frost_diagram, latimer_diagram, plot_frost, plot_latimer, plot_pourbaix
from microbial_thermo.figures.pourbaix import pourbaix_field, water_stability

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder, and why quoting it needs a pH

# %%
for ph in (0.0, 4.0, 7.0, 9.0, 12.0):
    ladder = latimer_diagram("As", pH=ph)
    chain = "  ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in ladder.steps)
    print(f"pH {ph:4.1f}   {chain}")

# %% [markdown]
# The As(V) species changes four times across that range — H₃AsO₄, H₂AsO₄⁻,
# HAsO₄²⁻, AsO₄³⁻ — and the potential changes with it. An arsenate potential
# quoted without naming the species is not a number, it is three numbers.

# %%
ladder = latimer_diagram("As", pH=7.0)
plot_latimer("As", diagram=ladder, show_half_reactions=True)
plt.show()

# %% [markdown]
# ## What is stable

# %%
plot_frost("As", pH=7.0, annotate_slopes=True)
plt.show()

# %%
frost = frost_diagram("As", pH=7.0)
print("most stable at pH 7:", frost.most_stable.backend)
print("most stable at pH 0:", frost_diagram("As", pH=0.0).most_stable.backend)
print("nothing disproportionates:", frost.disproportionation() == [])

# %% [markdown]
# Arsenic's curve is unremarkable — a clean convex chain with no unstable
# intermediate. The interesting thing is where the minimum sits: at pH 0 the
# element, at pH 7 **arsenite**, a dissolved species.

# %% [markdown]
# ### All four arsenates at once
#
# Arsenic is the crowded case: four As(V) forms and two As(III) forms, which
# is six species sharing two oxidation states. `predominant_only=False` keeps
# every one of them.

# %%
everything = frost_diagram("As", pH=7.0, activity=1e-6, predominant_only=False)
everything.to_frame()[["species", "oxidation state", "volt equivalent (V)", "on hull"]]

# %%
plot_frost("As", pH=7.0, activity=1e-6, predominant_only=False, label="all")
plt.show()

# %% [markdown]
# The labels are nudged apart and given leader lines back to their markers,
# because H₂AsO₄⁻ and HAsO₄²⁻ differ by only a few millivolts at pH 7 and
# would otherwise be printed on top of one another.
#
# The two options are independent — `show=` decides which points get a
# marker, `label=` which get named — and both take `"all"` or
# `"predominant"`. The default, `show="all", label="predominant"`, hides
# nothing and names only the form that exists:

# %%
figure, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
plot_frost("As", pH=7.0, activity=1e-6, predominant_only=False, ax=axes[0],
           title="show='all' (default)")
plot_frost("As", pH=7.0, activity=1e-6, predominant_only=False, show="predominant",
           ax=axes[1], title="show='predominant'")
plt.show()

# %% [markdown]
# **The open grey squares are the point of the figure.** They are the other
# arsenates — real species, at the right energies, simply not the ones that
# dominate at pH 7. Move the pH and they change places:

# %%
for ph in (0.0, 4.0, 7.0, 12.0):
    diagram = frost_diagram("As", pH=ph, activity=1e-6, predominant_only=False)
    winners = {p.oxidation_state: p.backend for p in diagram.predominant}
    print(f"pH {ph:5.1f}   As(V): {winners[5.0]:12s}  As(III): {winners[3.0]}")

# %% [markdown]
# Note what the grey squares are **not**: they are not drawn in red, because
# they are not disproportionating. A minority acid form is not falling apart,
# it is being outcompeted at the same oxidation state by a sibling that holds
# a different number of protons — and that is an acid–base question, not a
# redox one.
#
# The library keeps the two apart. Nothing at the same oxidation state can
# have a potential between it and its sibling, because no electrons move:

# %%
try:
    everything.slope("H3AsO4(aq)", "HAsO4--")
except ValueError as exc:
    print("refused:", exc)

# %% [markdown]
# and no minority form is ever reported as disproportionating:

# %%
minority = {p.backend for p in everything.points} - {p.backend for p in everything.predominant}
reported = {event.species for event in everything.disproportionation()}
print("minority forms:", sorted(minority))
print("reported as disproportionating:", sorted(reported) or "none")

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("As", figsize=(8.5, 6.5))
plt.show()

# %% [markdown]
# ## The measurement that matters
#
# How much of the range water can actually reach does each form hold? Only the
# band between the two water-stability lines is physically available.

# %%
field = pourbaix_field("As", points=241)
ph_axis, upper, lower = water_stability(25.0)

column = int(np.argmin(abs(field.ph - 7.0)))
window = (lower[int(np.argmin(abs(ph_axis - 7.0)))], upper[int(np.argmin(abs(ph_axis - 7.0)))])
print(f"water is stable from {window[0]:+.3f} to {window[1]:+.3f} V at pH 7\n")

previous = None
for row, winner in enumerate(field.winner[:, column]):
    name = field.species[int(winner)]
    if name != previous:
        inside = "inside" if window[0] <= field.eh[row] <= window[1] else "below"
        print(f"  Eh ≥ {field.eh[row]:+.3f} V : {name:14s} ({inside} the water window)")
        previous = name

# %% [markdown]
# Elemental arsenic only appears below about −0.28 V, and water itself gives
# out at −0.414 V. So As(0) has roughly 130 mV of the water window — a sliver
# — and everything above it is **dissolved**: arsenite from −0.28 to +0.02,
# arsenate above that.
#
# Now the crucial detail, which the Eh–pH diagram alone does not show:

# %%
for name in ("H3AsO4(aq)", "H2AsO4-", "HAsO4--", "As(OH)3(aq)", "H2AsO3-"):
    species = mt.resolve(name)
    print(f"{name:12s} charge {species.charge:+d}")

# %% [markdown]
# **Arsenite at pH 7 is As(OH)₃, and it is neutral.** Arsenate carries a
# charge and sorbs strongly to iron-oxide surfaces; uncharged arsenious acid
# does not. So reducing an aquifer does two things at once — it dissolves the
# iron oxide that was holding the arsenic, and it converts the arsenic to the
# form that would not stick to the oxide anyway.
#
# That is the Bengal-delta problem in two sentences, and both of them are on
# the diagrams above.

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
arsenic = [
    "arsenate_respiration_hydrogen",
    "arsenate_respiration_acetate",
    "arsenate_respiration_lactate",
    "arsenite_oxidation_oxygen",
    "arsenite_oxidation_nitrate",
    "arsenite_carbon_fixation",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[arsenic][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# Arsenate respiration is a genuine anaerobic metabolism, and it is the
# microbial engine of arsenic release: the organisms are not poisoning anyone
# deliberately, they are respiring a perfectly good electron acceptor.
#
# The last row is the striking one. Arsenite-driven carbon fixation is
# **endergonic** — it has to be, since arsenite is a poor donor and CO₂ a poor
# acceptor — and the organisms that do it are anoxygenic phototrophs paying
# the difference with light. Notebook `photoAs` works that budget out.

# %% [markdown]
# ## One thing worth taking away
#
# Arsenic has no insoluble sink at environmentally reachable potentials.
# Every form it takes between the water lines at pH 7 is dissolved. Selenium,
# next door, does have one — and notebook 11 measures the difference.
