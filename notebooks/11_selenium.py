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
# # 11. Selenium
#
# Selenium is the counter-example to arsenic, and the pair is worth teaching
# together. Both are metalloids, both are toxic as oxyanions, both are
# respired by anaerobes, and both have a Se(VI)/Se(IV)/Se(0) — As(V)/As(III)/
# As(0) — ladder that looks the same on paper.
#
# The difference is that **selenium's reduced element is an environmental
# sink and arsenic's is not**, and that single fact is why selenium
# bioremediation works and arsenic remediation is hard.

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

import microbial_thermo as mt
from microbial_thermo.figures import frost_diagram, latimer_diagram, plot_frost, plot_latimer, plot_pourbaix
from microbial_thermo.figures.pourbaix import pourbaix_field, water_stability

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder

# %%
ladder = latimer_diagram("Se", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("Se", diagram=ladder)
plt.show()

# %% [markdown]
# Both steps are strongly favourable at pH 7 — selenate to selenite at
# +0.45 V and selenite to the element at +0.26 V. Compare arsenic, where the
# first step is worth +0.013 V and the second is *negative*:

# %%
for element in ("Se", "As"):
    steps = latimer_diagram(element, pH=7.0).steps
    chain = "   ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps)
    print(f"{element}:  {chain}")

# %% [markdown]
# Reducing selenite all the way to elemental selenium pays. Reducing arsenite
# to elemental arsenic costs. That is the whole story, and everything below is
# a consequence of it.

# %% [markdown]
# ## What is stable

# %%
plot_frost("Se", pH=7.0, annotate_slopes=True)
plt.show()

# %%
frost = frost_diagram("Se", pH=7.0)
print("most stable at pH 7:", frost.most_stable.backend)
print("most stable at pH 0:", frost_diagram("Se", pH=0.0).most_stable.backend)

# %% [markdown]
# The element is the minimum at both pH values — unlike arsenic, where the
# minimum moves to the dissolved arsenite by pH 7.

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("Se", figsize=(8.5, 6.5))
plt.show()

# %% [markdown]
# ## Measuring the difference
#
# The question is how much of the *reachable* range — between the water
# stability lines — each element spends as an insoluble solid.

# %%
ph_axis, upper, lower = water_stability(25.0)
index = int(np.argmin(abs(ph_axis - 7.0)))
window = (lower[index], upper[index])


def solid_span(element, solid, pH=7.0):
    """Volts of the water window at this pH held by the elemental solid."""
    field = pourbaix_field(element, points=241)
    column = int(np.argmin(abs(field.ph - pH)))
    wanted = field.species.index(solid)
    inside = [
        eh
        for eh, winner in zip(field.eh, field.winner[:, column])
        if window[0] <= eh <= window[1] and int(winner) == wanted
    ]
    return (max(inside) - min(inside)) if inside else 0.0


print(f"water is stable over {window[1] - window[0]:.3f} V at pH 7\n")
for element, solid in (("Se", "Se"), ("As", "As")):
    span = solid_span(element, solid)
    print(f"  elemental {element}: {span:.3f} V of that window ({span / (window[1] - window[0]):.0%})")

# %% [markdown]
# Elemental selenium holds 0.58 V — essentially the whole reducing half of
# the window. Elemental arsenic holds 0.12 V, a sliver at the very bottom. An anaerobic sediment
# sits comfortably inside selenium's solid field and comfortably inside
# **arsenite's dissolved** field.
#
# So the same treatment — make it anoxic, let the anaerobes respire — takes
# selenium out of solution and puts arsenic into it.

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
selenium = [
    "selenate_respiration_hydrogen",
    "selenate_respiration_acetate",
    "selenite_reduction_hydrogen",
    "selenite_reduction_acetate",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[selenium][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# Selenate respiration on hydrogen pays about −83 kJ/mol e⁻ here: roughly
# three quarters of what nitrate pays and four times what sulfate pays. That
# is a comfortable living, which is why selenate-respiring bacteria are easy
# to enrich and why the reaction is used industrially — feed an anaerobic
# bioreactor a carbon source, and the organisms precipitate red elemental
# selenium you can filter out.
#
# There is no arsenic equivalent, because there is no arsenic solid to make.

# %% [markdown]
# Both steps have an entry on both donors, which is the point: **selenium
# reduction is a complete two-step pathway that ends in a solid**, and either
# step can be run on hydrogen or on organic carbon. That is what makes a
# bioreactor design possible — feed it either, and the selenium comes out as
# a filterable red precipitate.
#
# It is worth noticing what the second step is *not* doing:

# %%
for name in ("selenite_reduction_acetate", "selenate_respiration_acetate"):
    print(f"{table.loc[name, 'label']:38s} {table.loc[name, 'dG per e- (kJ/mol)']:+7.1f} kJ/mol e-")

# %% [markdown]
# Reducing selenite to the element pays *less* per electron than reducing
# selenate to selenite. An organism optimising purely for energy would stop at
# selenite — and some do. The ones used industrially go the whole way, which
# means the useful step for remediation is the one the organism has least
# incentive to take. Enrichment conditions matter as much as thermodynamics.
#
# ### The environmental reason anyone cares
#
# Selenium is the textbook case of a nutrient with a narrow window: essential
# in selenocysteine and selenoproteins at trace level, toxic not far above it,
# with the gap between the two among the smallest of any element. Irrigation
# of seleniferous soils concentrates selenate in drainage water, and selenate
# is the mobile, bioavailable form — the Kesterson Reservoir deformities in
# the 1980s are what put selenium into environmental regulation.
#
# The Eh–pH diagram above is the remediation strategy in one picture: get the
# water into the lower half of it and the problem precipitates.

# %% [markdown]
# ## One thing worth taking away
#
# "Reduction immobilises metals" is one of those rules that is true often
# enough to be dangerous. It is true for selenium, uranium and chromium, and
# false for arsenic and iron — and which way it goes is written on the Eh–pH
# diagram, in whether the reduced form has a solid field inside the water
# window.
#
# It takes about ten lines to check, as above, and it is worth checking every
# time.
