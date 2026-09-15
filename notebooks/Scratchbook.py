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
# # Scratchbook: from a written reaction to the full analysis
#
# Start with nothing but an unbalanced reaction typed as a string, and let the
# library do the rest: work out the couples, balance both half reactions, and
# draw everything.
#
# The worked example is nitrate reduction to hydroxylamine with hydrogen, at
# **pH 7.7**.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.balance import infer_couples
from microbial_thermo.figures import plot_energy_explorer, plot_half_reactions, plot_redox_tower

EQUATION = "NO3- + H2 -> NH2OH"
conditions = mt.Conditions(temperature_c=25.0, pH=7.7)

# %% [markdown]
# > **A caution before anything else.** Hydroxylamine is in none of pyGCC's
# > databases, so its formation energy comes from the hand-entered supplemental
# > table and has **not** been traced to a primary source. Every number below
# > inherits that uncertainty, and the library warns each time the value is
# > used. Check it before relying on any of this.

# %%
from microbial_thermo.supplemental import supplemental_species

entry = supplemental_species()["NH2OH(aq)"]
print(f"NH2OH(aq)  dGf = {entry.delta_Gf_kJ_mol} kJ/mol at {entry.temperature_c} C")
print(f"verified: {entry.verified}")
print(entry.provenance.strip())

# %% [markdown]
# ## 1. What the string implies
#
# `infer_couples` reads the equation, finds which elements change oxidation
# state, and forms the couples from that. Note that H₂'s partner — the proton —
# never appears in the written equation; it is supplied.

# %%
donor, acceptor = infer_couples(EQUATION)
print("equation:", EQUATION)
print("donor    (reduced, oxidized):", donor)
print("acceptor (reduced, oxidized):", acceptor)

# %% [markdown]
# Nitrogen goes from +5 in nitrate to −1 in hydroxylamine, so nitrate is the
# acceptor. Hydrogen goes 0 → +1, so H₂ is the donor.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    reaction = mt.Reaction.from_equation(EQUATION, conditions=conditions)

print(reaction.summary())

# %% [markdown]
# ## 2. The two half reactions

# %%
print("oxidation (donor):  ", reaction.donor_half.format())
print("reduction (acceptor):", reaction.acceptor_half.format())
print()
for label, result in (("donor", reaction.donor_half), ("acceptor", reaction.acceptor_half)):
    oxidized, reduced = result.oxidation_states()
    print(
        f"{label:9s} {result.couple}  "
        f"{result.half.key_element}: {oxidized} -> {reduced}   "
        f"E0' = {result.E_standard_prime.to('V').magnitude:+.3f} V"
    )

# %% [markdown]
# ## 3. Normalised to one mole of the first reactant
#
# The first reactant in the written equation is nitrate, so this is the
# reaction per mole of nitrate reduced.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    per_nitrate = mt.Reaction.from_equation(
        EQUATION, conditions=conditions, normalize_to="NO3-"
    )

print(per_nitrate.format())
print(f"n electrons = {per_nitrate.n_electrons}")
plot_half_reactions(per_nitrate, n_electrons=None)
plt.show()

# %% [markdown]
# ## 4. Normalised to an electron pair
#
# The same chemistry, scaled so exactly two electrons move. This is the
# library's default for figures, and it is what makes couples comparable
# against each other.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    per_pair = per_nitrate.renormalized(2)

print(per_pair.format())
plot_half_reactions(per_pair)
plt.show()

# %% [markdown]
# A detail worth noticing before the tower: the donor here is `H2(aq)`, because
# a bare `H2` resolves to dissolved hydrogen — the physiologically relevant
# form. The tower's `2H+/H2` reference rung uses H₂ **gas**, so the donor marker
# will not sit exactly on that rung:

# %%
from microbial_thermo.reaction import Couple, HalfReactionResult

for name in ("H2(aq)", "H2(g)"):
    pair = Couple.make(name, "H+")
    result = HalfReactionResult(
        couple=pair, half=pair.half_reaction(),
        conditions=conditions, backend=mt.get_backend(),
    )
    print(f"{name:8s} E0' = {result.E_standard_prime.to('V').magnitude:+.3f} V")

print("\nWrite H2(g) explicitly if you want the gas-phase convention.")

# %% [markdown]
# ## 5. The electron tower, both ways
#
# The couple potentials are **intensive**: they do not change with
# normalisation, so the rungs sit in exactly the same places on both towers.
# What changes is the total free energy and therefore the energy scale bar,
# which is drawn for the reaction's own electron count.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_redox_tower(per_nitrate, n_electrons=None)
plt.show()

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_redox_tower(per_pair)
plt.show()

# %% [markdown]
# Compare the two energy panels: ΔG°′ scales with the reaction while ΔG per
# electron and ΔE°′ do not. That is the whole reason the library normalises
# figures to an electron pair by default.

# %%
for label, built in (("per nitrate", per_nitrate), ("per electron pair", per_pair)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        print(
            f"{label:18s} n={built.n_electrons!s:>4s}  "
            f"dG0' = {built.delta_G_standard_prime.magnitude:+8.1f} kJ/mol   "
            f"per e- = {built.delta_G_per_electron.magnitude:+7.2f}   "
            f"dE0' = {built.delta_E_standard_prime.to('V').magnitude:+.3f} V"
        )

# %% [markdown]
# ## 6. Free energy against pH
#
# pH is the first entry in the dropdown; the other axes are there too. The
# temperature axis is deliberately left out: hydroxylamine's supplemental value
# is tabulated at 25 °C only and refuses to extrapolate, so a temperature sweep
# would be all gaps.

# %%
from microbial_thermo.sweep import concentration_axis, partial_pressure_axis, ph_axis

axes = [
    ph_axis(low=4.0, high=10.0, n=31),
    concentration_axis("NO3-"),
    concentration_axis("NH2OH(aq)"),
    partial_pressure_axis("H2(aq)"),
]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    figure = plot_energy_explorer(per_pair, axes=axes)

figure.show()

# %% [markdown]
# ## 7. Show the work
#
# Every number above, traced back to where it came from.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    work = per_pair.show_work()
work

# %% [markdown]
# ### Still to do
#
# - NH₃ + MnO₂ → N₂ + Mn²⁺ at pH 10
# - NH₄⁺ → N₂ under aerobic conditions
#
# Both should work through `Reaction.from_equation` the same way. The first
# needs pyrolusite, which is tabulated across temperature; several other
# manganese oxides are 25 °C only.
