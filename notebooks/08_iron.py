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
# # 8. Iron
#
# Iron is the most abundant redox-active metal on Earth and the one where the
# solid phases matter most. Its aqueous chemistry is a footnote to its mineral
# chemistry: what a cell actually respires is usually a surface, and which
# surface depends on pH, on how much carbonate and sulfide are about, and on
# how much iron is dissolved.
#
# That makes iron the element where the three diagrams disagree most
# interestingly, and where knowing what each one is *for* matters.

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
from microbial_thermo.figures.basis import element_series
from microbial_thermo.figures.pourbaix import pourbaix_field

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder, and what is not allowed on it

# %%
ladder = latimer_diagram("Fe", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("Fe", diagram=ladder)
plt.show()

# %% [markdown]
# Magnetite sits between Fe(III) and Fe(II) at **+8/3**, which is not a
# typo: Fe₃O₄ is one ferrous iron and two ferric, and the mean state is
# exactly 8/3. The fractional electron counts on the arrows follow from it.
#
# Two iron minerals were offered and refused a rung, and the reason is worth
# reading:

# %%
for name, why in ladder.excluded:
    print(f"{name}: {why}\n")

# %% [markdown]
# **Pyrite is the case that forces the rule.** FeS₂ is iron(II) with sulfur at
# −I, a disulfide. Our basis accounts for sulfur as sulfate, so forming pyrite
# from Fe(0) takes twelve electrons per iron — a real number about a real
# reaction, and emphatically not iron's oxidation state. A rung on a Latimer
# or Frost diagram must be a species whose electron count *is* a state, so
# pyrite is refused.
#
# Siderite is refused with it, although FeCO₃ really is Fe(II) and would have
# been right. Telling those two apart mechanically needs an oxidation-state
# assignment independent of the decomposition, and the one this library has
# refuses metal sulfides outright — which is precisely the case in question. A
# rule that is never wrong beats a rule that is sometimes wrong and silent
# about it. Nothing is lost: siderite is Fe(II), and Fe²⁺ already holds that
# rung.
#
# Both are back on the Eh–pH diagram below, which is where a siderite or
# pyrite field belongs.

# %% [markdown]
# ## What is stable, and why unit activity is the wrong question for iron

# %%
plot_frost("Fe", pH=7.0, activity=1.0)
plt.show()

# %%
unit = frost_diagram("Fe", pH=7.0, activity=1.0)
for event in unit.disproportionation():
    print(event)

# %% [markdown]
# At unit activity the diagram claims Fe²⁺ disproportionates into magnetite
# and iron metal. That is arithmetically correct and physically absurd: a 1 M
# ferrous solution at pH 7 is supersaturated with respect to essentially every
# iron mineral there is. It is not a condition; it is a convention.
#
# Draw it at an activity a groundwater might actually have:

# %%
plot_frost("Fe", pH=7.0, activity=1e-6)
plt.show()

# %%
real = frost_diagram("Fe", pH=7.0, activity=1e-6)
print("on the hull:", [p.backend for p in real.stable])
for event in real.disproportionation():
    print(event)

# %% [markdown]
# Now Fe²⁺ is stable and **magnetite** is the marginal phase, unstable to
# hematite plus Fe²⁺ by a few kJ per mole. That is about right: magnetite's
# stability field is genuinely narrow and genuinely sensitive to how much
# ferrous iron is in solution. Watch it move:

# %%
for activity in (1e-2, 1e-4, 1e-6, 1e-8):
    diagram = frost_diagram("Fe", pH=7.0, activity=activity)
    off = [f"{e.species} ({e.delta_g:+.1f} kJ/mol)" for e in diagram.disproportionation()]
    print(f"activity {activity:.0e}: hull {[p.backend for p in diagram.stable]}"
          f"{'  off-hull: ' + ', '.join(off) if off else ''}")

# %% [markdown]
# ## Where each form lives
#
# This is the diagram iron actually needs, because it is the one that can hold
# a carbonate and a sulfide field at the same time.

# %%
plot_pourbaix("Fe", figsize=(8.5, 6.5))
plt.show()

# %% [markdown]
# The pyrite field at the bottom exists **only because there is sulfide**, at
# a concentration stated in the figure title. Change that number and the field
# changes size; set it to zero and the field disappears. It is a modelling
# input, not a property of iron:

# %%
plot_pourbaix("Fe", fixed={"S": ("SO4-2", 0), "C": ("HCO3-", 2e-3)}, figsize=(8.5, 6.5))
plt.show()

# %% [markdown]
# Pyrite is gone, and the subtitle no longer claims sulfate is held fixed —
# nothing on the figure contains sulfur, so there is nothing to hold.
#
# What moved into the space it held is worth measuring rather than guessing:

# %%
import collections

with_sulfide = pourbaix_field("Fe", fixed={"S": ("SO4-2", 1e-6), "C": ("HCO3-", 2e-3)})
without = pourbaix_field("Fe", fixed={"S": ("SO4-2", 0), "C": ("HCO3-", 2e-3)})

was_pyrite = with_sulfide.winner == with_sulfide.species.index("Pyrite")
print(f"pyrite held {was_pyrite.mean():.1%} of the panel; without sulfide that area is now")
counts = collections.Counter(without.species[int(i)] for i in without.winner[was_pyrite])
for name, count in counts.most_common():
    print(f"   {name:12s} {count / was_pyrite.sum():5.1%}")

# %% [markdown]
# **Siderite appears, and takes almost none of it.** Nearly half goes back to
# dissolved Fe²⁺ and most of the rest to iron metal below the water line;
# carbonate picks up about two percent.
#
# That is the useful result, and it is the opposite of what "remove the
# sulfide and the carbonate takes over" would predict. Sulfide is a far
# stronger sink for ferrous iron than carbonate is — which is why pyrite, not
# siderite, is the mineral that ends up holding iron in anoxic marine
# sediment, and why siderite is a freshwater mineral, where there is
# carbonate about and very little sulfate to reduce.

# %% [markdown]
# An activity of exactly zero means *the element is absent*, and every species
# needing it is dropped by name rather than evaluated at `log(0)`:

# %%
series = element_series("Fe", fixed={"S": ("SO4-2", 0), "C": ("HCO3-", 2e-3)})
for name, why in series.skipped:
    print(f"{name}: {why}")

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
iron = ["iron_oxidation_oxygen", "ferrihydrite_reduction_acetate", "goethite_reduction_hydrogen"]
table = mt.energy_table(conditions).set_index("name")
table.loc[iron][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# ## One thing worth taking away
#
# The three diagrams are not three views of one truth; they answer three
# questions, and iron is where the difference shows.
#
# * The **Latimer** ladder asks what a one-step reduction is worth. It needs
#   species whose oxidation state is meaningful, so it excludes pyrite.
# * The **Frost** curve asks which forms survive being left alone. It depends
#   on activity so strongly that the unit-activity version — the one in every
#   textbook — describes a solution that cannot exist.
# * The **Eh–pH** field asks what is actually present. It is the only one that
#   can hold other elements at stated concentrations, and therefore the only
#   one that can say anything about siderite or pyrite.
#
# Asking the wrong one of the three is the most common way to get a confident
# wrong answer about iron.
