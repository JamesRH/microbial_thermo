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
# # 7. Sulfur
#
# Sulfur is the element where this library's diagrams earn their keep, because
# sulfur's most interesting metabolism is one that a textbook Frost diagram —
# drawn at unit activity, as they all are — says is impossible.
#
# The claim to test: **elemental sulfur disproportionation**,
#
#     4 S⁰ + 4 H₂O → 3 H₂S + SO₄²⁻ + 2 H⁺
#
# a real metabolism run by *Desulfocapsa* and others, in which sulfur is both
# the donor and the acceptor. It is famously marginal, and famously requires a
# sulfide sink to run at all. We can compute exactly how marginal.

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

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

# %%
ladder = latimer_diagram("S", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("S", diagram=ladder)
plt.show()

# %% [markdown]
# Every rung changes species between pH 0 and pH 7 — bisulfate to sulfate,
# H₂S to HS⁻ — because sulfur's forms are acids with p*K*a values inside the
# environmental range. That is the reason a sulfur diagram quoted "at standard
# conditions" is nearly useless for a sediment:

# %%
for ph in (0.0, 2.0, 7.0, 9.0):
    steps = latimer_diagram("S", pH=ph).steps
    chain = "  ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps)
    print(f"pH {ph:4.1f}   {chain}")

# %% [markdown]
# ## Frost at unit activity: sulfur is stable

# %%
plot_frost("S", pH=7.0, activity=1.0, annotate_slopes=True)
plt.show()

# %%
textbook = frost_diagram("S", pH=7.0, activity=1.0)
print("on the hull:", [p.backend for p in textbook.stable])
for event in textbook.disproportionation():
    print(event)

# %% [markdown]
# Elemental sulfur is **on the hull**. At unit activity the diagram says
# disproportionation cannot happen, and every published sulfur Frost diagram
# says the same thing, because every published one is drawn at unit activity.
#
# Sulfite is above the hull, which is also true and also worth knowing:
# sulfite disproportionation is itself a microbial metabolism, and a much
# easier one.

# %% [markdown]
# ## Frost at a realistic activity: sulfur is not stable

# %%
plot_frost("S", pH=7.0, activity=1e-5, annotate_slopes=True)
plt.show()

# %%
real = frost_diagram("S", pH=7.0, activity=1e-5)
print("on the hull:", [p.backend for p in real.stable])
for event in real.disproportionation():
    print(event)

# %% [markdown]
# Elemental sulfur has left the hull. Drop the dissolved sulfur species from
# unit activity to 10 µM — an ordinary porewater — and disproportionation
# becomes exergonic by about 18 kJ per mole of sulfur.
#
# **Why lowering the activity changes the answer.** Sulfur is a solid, so its
# energy does not move; the dissolved products, sulfate and sulfide, both get
# cheaper as they get more dilute. The reaction makes four moles of dissolved
# product from a solid, so diluting the solution pulls it forward. This is Le
# Chatelier's principle, and it is the whole of the sulfur-disproportionation
# story.

# %% [markdown]
# ## Where the line is

# %%
activities = np.logspace(-6, 0, 25)
energies = []
for activity in activities:
    diagram = frost_diagram("S", pH=7.0, activity=activity)
    found = [e for e in diagram.disproportionation() if e.species == "Sulfur(s)"]
    energies.append(found[0].delta_g if found else 0.0)

figure, ax = plt.subplots(figsize=(8.0, 4.8))
ax.semilogx(activities, energies, linewidth=2.0, color="#1F6FB4")
ax.axhspan(-20.0, 0.0, color="#E8B84B", alpha=0.2, label="inside the biological energy quantum")
ax.axhline(0.0, color="#333333", linewidth=0.9)
ax.set_xlabel("dissolved sulfur activity")
ax.set_ylabel("ΔG of disproportionation\n(kJ per mol S)")
ax.set_title("4 S⁰ + 4 H₂O → 3 HS⁻ + SO₄²⁻ + 5 H⁺ at pH 7", fontsize=12)
ax.legend(frameon=False, fontsize=9)
plt.show()

# %% [markdown]
# The reaction turns exergonic at about 1.7 × 10⁻², but exergonic is not the
# same as *useful*: the curve stays inside the biological energy quantum — the
# shaded band, roughly the least free energy thought to sustain a cell — until
# the activity falls below about 5 × 10⁻⁶.
#
# That is the quantitative version of the textbook statement that sulfur
# disproportionators **need a sulfide sink**. Precipitating FeS keeps sulfide
# low; keeping sulfide low is what keeps this reaction worth running. The
# organism does not merely tolerate an iron-rich sediment, it requires one.

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("S", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# Note what is *not* on that diagram: elemental sulfur has no field at 10⁻⁶
# activity. It is metastable everywhere, which is the same result as the Frost
# curve in a different projection.

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
sulfur = [
    "hydrogenotrophic_sulfate_reduction",
    "acetotrophic_sulfate_reduction",
    "sulfur_reduction_hydrogen",
    "sulfide_oxidation_oxygen",
    "sulfur_oxidation_oxygen",
    "sulfide_oxidation_nitrate",
    "anaerobic_methane_oxidation",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[sulfur][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# ## One thing worth taking away
#
# A Frost diagram is a picture of a *chosen* state, not of the element. The
# activity it was drawn at is part of the claim it makes, and for sulfur the
# choice flips the conclusion completely: stable at 1 M, disproportionating at
# 10 µM, and a living organism in the gap.
#
# Every diagram in this library prints its conditions in the title for exactly
# this reason.
