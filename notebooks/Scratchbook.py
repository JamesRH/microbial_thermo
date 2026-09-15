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
# # Scratchbook

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.balance import infer_couples
from microbial_thermo.figures import plot_energy_explorer, plot_half_reactions, plot_redox_tower

# %% [markdown]
# # Nitrate to Hydroxylamine by H2 oxidation

# %%
EQUATION = "NO3- + H2 -> NH2OH"
conditions = mt.Conditions(temperature_c=25.0, pH=4)

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

# %%
mt.half_reaction("NO3-", "NH2OH(aq)").E_standard
mt.half_reaction("H2(aq)","H+").E_standard

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
#

# %% [markdown]
# # Ammonia oxidation by manganese oxide, pH 10
#
# pH 10 is chosen deliberately: ammonia's p*K*a is near 9.2, so above it the
# neutral NH₃ is the dominant form rather than ammonium. Worth checking before
# writing the reaction with one form or the other.

# %%
mt.speciation_table("ammonia", 10.0)

# %%
alkaline = mt.Conditions(temperature_c=25.0, pH=10.0)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    manganese = mt.Reaction.from_equation(
        "NH3 + MnO2 -> N2 + Mn+2", conditions=alkaline, normalize_to="integer"
    )

print(manganese.summary())

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_half_reactions(manganese)
plt.show()

# %% [markdown]
# Manganese(IV) is a strong enough oxidant to take ammonia all the way to N₂,
# and strongly exergonic at that. Note this is written with pyrolusite, which
# has a full temperature grid; several other manganese oxides in the database
# are tabulated at 25 °C only and will refuse other temperatures.

# %% [markdown]
# The tower, at an electron pair so it can be set against the other two
# questions directly. Note the acceptor rung sits where pH 10 puts it, not
# where a pH 7 table would.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_redox_tower(manganese)
plt.show()

# %% [markdown]
# Free energy against pH. Written for an electron pair this reaction **consumes**
# two protons, one per electron, so raising the pH makes it *less* favourable:
# −156 kJ/mol at pH 5, −88 at pH 11. The slope is
# $\nu_{H^+} RT\ln 10$ = 11.4 kJ/mol per pH unit.
#
# Temperature is in the dropdown here, unlike the first question: pyrolusite has
# a full log K grid, where hydroxylamine did not.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    manganese_explorer = plot_energy_explorer(manganese.renormalized(2))

manganese_explorer.show()

# %% [markdown]
# # Ammonium oxidation to N₂ under oxygen, pH 7

# %%
neutral = mt.Conditions(temperature_c=25.0, pH=7.0)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    aerobic = mt.Reaction.from_equation(
        "NH4+ + O2 -> N2", conditions=neutral, normalize_to="integer"
    )

print(aerobic.summary())

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_half_reactions(aerobic)
plt.show()

# %% [markdown]
# The tower again at an electron pair. The gap here is the widest of the three,
# which is what the free energy per electron already said.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_redox_tower(aerobic)
plt.show()

# %% [markdown]
# Free energy against pH — and compare it with the manganese case, because the
# two go in *opposite directions*. This reaction **releases** two thirds of a
# proton per electron pair, so raising the pH makes it *more* favourable:
# −205 kJ/mol at pH 5, −228 at pH 11.
#
# The slope is 3.8 kJ/mol per pH unit against manganese's 11.4 — a factor of
# three, which is exactly the ratio of their proton coefficients (1 per electron
# against 1/3). Nothing about pH sensitivity is arbitrary: it is the proton
# stoichiometry, read straight off the balanced half reactions.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    aerobic_explorer = plot_energy_explorer(aerobic.renormalized(2))

aerobic_explorer.show()

# %% [markdown]
# ### The three side by side
#
# Per electron, so they are comparable despite very different electron counts.

# %%
rows = [
    ("NO3- + H2 -> NH2OH", per_pair, conditions),
    ("NH3 + MnO2 -> N2 + Mn+2", manganese, alkaline),
    ("NH4+ + O2 -> N2", aerobic, neutral),
]
print(f"{'equation':28s} {'pH':>5s} {'n':>3s} {'dG0prime':>11s} {'per e-':>9s} {'dE0prime':>9s}")
for label, built, where in rows:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        print(
            f"{label:28s} {where.pH:5.1f} {built.n_electrons!s:>3s} "
            f"{built.delta_G_standard_prime.magnitude:11.1f} "
            f"{built.delta_G_per_electron.magnitude:9.2f} "
            f"{built.delta_E_standard_prime.to('V').magnitude:+9.3f}"
        )

# %% [markdown]
# A caution about reading that table. The three rows are **not** a fair
# comparison of acceptor strength: each uses a different electron donor, and
# each sits at a different pH. Oxygen looks best partly because it is, and
# partly because its donor is ammonium rather than ammonia.
#
# To compare acceptors properly, hold the donor and the pH fixed and look at the
# acceptor half reactions alone — which is exactly what the redox tower does:

# %%
from microbial_thermo.reaction import Couple, HalfReactionResult

same_conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
for reduced, oxidized in [
    ("H2O", "O2(aq)"),
    ("Mn+2", "pyrolusite"),
    ("NH2OH(aq)", "NO3-"),
]:
    pair = Couple.make(reduced, oxidized)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = HalfReactionResult(
            couple=pair, half=pair.half_reaction(),
            conditions=same_conditions, backend=mt.get_backend(),
        )
        value = result.E_standard_prime.to("V").magnitude
    print(f"{str(pair):26s} E0' = {value:+.3f} V  (pH 7)")

# %% [markdown]
# On a common footing the order is oxygen, then manganese(IV), then nitrate to
# hydroxylamine. Pyrolusite looked weak in the table above only because pH 10
# penalises a couple that consumes four protons — the same reaction is far more
# favourable at pH 7.
#
# Only the nitrate row depends on the unverified hydroxylamine value; the other
# two rest entirely on pyGCC data.
