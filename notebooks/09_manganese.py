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
# # 9. Manganese
#
# Manganese is iron's awkward neighbour. It has more accessible oxidation
# states, its oxides are stronger oxidants, and — unlike iron — its
# intermediate states are genuinely unstable with respect to the ends. Mn(III)
# is the textbook disproportionating species, and this notebook watches it
# stop disproportionating as the pH rises.
#
# It is also the element where our data has a visible flaw, so it is where we
# look at that honestly.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.figures import frost_diagram, latimer_diagram, plot_frost, plot_latimer, plot_pourbaix
from microbial_thermo.figures.basis import element_grid

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder

# %%
ladder = latimer_diagram("Mn", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("Mn", diagram=ladder)
plt.show()

# %% [markdown]
# Five rungs, two of them mixed-valence oxides. Hausmannite is Mn₃O₄ at
# **+8/3** — the manganese analogue of magnetite — and bixbyite is Mn₂O₃ at
# +III, with manganite (MnOOH) offered at the same state and beaten to the
# rung:

# %%
ladder.alternatives

# %% [markdown]
# ## Mn(III) disproportionates, until it does not

# %%
for ph in (0.0, 2.0, 4.0, 6.0, 7.0, 8.0):
    diagram = frost_diagram("Mn", pH=ph, activity=1e-6)
    off = ", ".join(f"{e.species} {e.delta_g:+.1f}" for e in diagram.disproportionation())
    print(f"pH {ph:4.1f}   {off or 'nothing above the hull'}")

# %% [markdown]
# At pH 2 hausmannite is unstable by 45 kJ per mole and bixbyite by 30. By pH
# 7 those have fallen to 7 and 1.5 — inside the biological energy quantum,
# which is to say effectively at equilibrium.
#
# That is the whole reason Mn(III) is a laboratory curiosity in acid and a
# **major environmental species** at circumneutral pH: manganese(III) complexes
# are now known to be a large fraction of dissolved manganese in porewaters,
# and this is the thermodynamic room they occupy.

# %%
plot_frost("Mn", pH=2.0, activity=1e-6)
plt.show()

# %%
plot_frost("Mn", pH=7.0, activity=1e-6)
plt.show()

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("Mn", figsize=(8.5, 6.5))
plt.show()

# %% [markdown]
# Mn²⁺ holds nearly the whole field, with the oxides confined to the top right
# — high pH *and* strongly oxidising. Compare the iron diagram, where the
# ferric oxides take most of the space above the water line. Manganese oxides
# are harder to make and easier to reduce, which is why manganese reduction
# sits above iron reduction on the redox ladder and why manganese oxides are
# the strongest oxidant a sediment normally has after oxygen and nitrate.

# %% [markdown]
# ## Two things about the data, stated plainly
#
# ### Elemental manganese is not zero

# %%
backend = mt.get_backend()
for name in ("Mn", "Fe", "Cu", "Zn", "Ni", "As", "Se"):
    print(f"{name:3s} {backend.delta_Gf(name, 25.0).to('kJ/mol').magnitude:+.4f} kJ/mol")

# %% [markdown]
# Every elemental form should be exactly zero by definition, and every one of
# them comes out within 0.004 kJ/mol of it — except manganese, at −2.5. That
# is a database problem, not a rounding one, and **any manganese Frost or
# Latimer diagram inherits it**: 2.5 kJ/mol over the two-electron Mn²⁺/Mn step
# is about 13 mV.
#
# It affects only the metal itself, which no sediment contains, so every
# environmentally interesting number above is unaffected. It is recorded here
# rather than quietly corrected, because a silently patched database is worse
# than a documented flaw.

# %% [markdown]
# ### Manganite has one temperature

# %%
grid = element_grid("Mn", temperatures=(10.0, 25.0, 40.0))
print("on the interactive grid:", grid.species)
for name, why in grid.dropped:
    print(f"\ndropped: {name}\n  {why}")

# %% [markdown]
# Manganite carries a single log K, at 25 °C, so it can be on the static 25 °C
# diagram and cannot be on one with a temperature slider. A phase that appears
# and vanishes as a slider moves is a different diagram each time, not the
# same diagram at new conditions — so it is dropped from the whole grid and
# named on the figure.
#
# Ask for the default manganese ladder away from 25 °C and the library
# refuses rather than dropping a phase quietly:

# %%
from microbial_thermo.exceptions import OutOfRangeError

try:
    latimer_diagram("Mn", pH=7.0, temperature_c=2.0)
except OutOfRangeError as exc:
    print("refused:", exc)

# %% [markdown]
# The fix is to say which species you mean — the grid has already worked out
# which ones survive the whole temperature range:

# %%
for temperature in (2.0, 25.0, 60.0):
    steps = latimer_diagram(
        "Mn", species=list(grid.species), pH=7.0, temperature_c=temperature
    ).steps
    chain = "  ".join(f"{s.potential:+.3f}" for s in steps)
    print(f"{temperature:5.1f} °C   {chain}")

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
table = mt.energy_table(conditions).set_index("name")
table.loc[["manganese_oxidation_oxygen", "manganese_reduction_acetate"]][
    ["label", "dG per e- (kJ/mol)", "reaction"]
]

# %% [markdown]
# ## One thing worth taking away
#
# Manganese oxidation by oxygen is exergonic but notoriously slow abiotically,
# and manganese-oxidising bacteria are catalysing a reaction that
# thermodynamics has already approved. Manganese reduction is the reverse and
# pays well. Between them sits Mn(III), too unstable to accumulate in acid and
# stable enough to matter at pH 7 — which is why the manganese cycle has a
# genuine three-species character that the iron cycle does not.
