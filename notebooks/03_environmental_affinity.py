# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: microbial-thermo
#     language: python
#     name: microbial-thermo
# ---

# %% [markdown]
# # 3. Speciation, and what pays in a real porewater
#
# Redox zonation is usually presented as a list to memorise. It is not: it
# falls out of the thermodynamics once you put real concentrations in.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## pH speciation
#
# Measurements come as totals — total sulfide, DIC, total ammonia — while a
# reaction is written with one specific form. At pH 7 sulfide is almost an even
# split, so treating a measured total as though it were all HS⁻ is wrong by
# about a factor of two.

# %%
mt.speciation_table("sulfide", 7.0)

# %%
for ph in (5.0, 6.0, 7.0, 8.0, 9.0):
    member, fraction = mt.dominant("sulfide", ph)
    print(f"pH {ph}: {member.backend:8s} holds {fraction:.1%}")

# %% [markdown]
# p*K*a values are computed from the backend at the working temperature, not
# looked up, so they move with temperature:

# %%
for temperature in (5.0, 25.0, 60.0, 95.0):
    ladder = mt.pKa_ladder("sulfide", temperature_c=temperature)
    print(f"{temperature:5.1f} °C   pKa = {ladder[0]:.2f}")

# %% [markdown]
# Give a family total and the split is handled for you:

# %%
anoxic = mt.Conditions(
    temperature_c=12.0,
    pH=7.4,
    activity_model="ideal",
    total_concentrations={"sulfide": 1e-5, "DIC": 3e-3},
    concentrations={
        "SO4-2": 2.0e-2, "NO3-": 5e-6, "Fe+2": 5e-5, "Mn+2": 2e-5,
        "acetate": 1e-5, "NH4+": 1e-4, "NO2-": 1e-7, "lactate": 1e-6,
        "methanol": 1e-6, "O2(aq)": 1e-6,
    },
    partial_pressures={"H2(g)": 5e-6, "CH4(g)": 1e-3},
)

# %% [markdown]
# ## The curated library
#
# Named metabolisms, so you do not have to remember which couples to pair.

# %%
from microbial_thermo.library import default_library, energy_table

library = default_library()
print(f"{len(library)} metabolisms in {len(library.groups())} groups")
print(library.groups())

# %%
mt.metabolism("anammox")

# %% [markdown]
# ## What actually pays here

# %%
table = energy_table(anoxic)
table[["label", "group", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# ## The affinity ladder
#
# The same data as a figure. Bars reaching into the shaded band are exergonic
# but probably too weakly so to conserve energy for a cell.

# %%
from microbial_thermo.figures import plot_affinity_ladder

plot_affinity_ladder(anoxic, figsize=(11.5, 9))
plt.show()

# %% [markdown]
# Notice what has happened. At standard state the aerobic metabolisms win by a
# wide margin and sulfate reduction is comfortably exergonic. Under this
# porewater — oxygen nearly gone, hydrogen at 5 nM — the sulfate and methane
# metabolisms have collapsed into the energy-quantum band, which is exactly the
# regime where syntrophy and very low product concentrations become necessary.
#
# ### Exercise
#
# Raise the hydrogen partial pressure to $10^{-3}$ bar and redraw. Which
# metabolisms move out of the band, and why those ones?

# %% [markdown]
# ## Comparing the same metabolism across conditions

# %%
from microbial_thermo.library import reaction as build

standard = mt.Conditions(temperature_c=25.0, pH=7.0)
for name in ("acetotrophic_sulfate_reduction", "hydrogenotrophic_methanogenesis"):
    at_standard = build(name, standard).delta_G_per_electron.magnitude
    in_situ = build(name, anoxic).delta_G_per_electron.magnitude
    print(f"{name:38s} standard {at_standard:+7.1f}   in situ {in_situ:+7.1f} kJ/mol e-")
