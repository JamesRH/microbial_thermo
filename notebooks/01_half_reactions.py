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
# # 1. Half reactions, balancing, and oxidation states
#
# This notebook covers the basics: naming a redox couple, seeing it balanced,
# reading off the potentials, and producing the stacked half-reaction figure.
#
# Every number here is computed from pyGCC at the stated conditions. Nothing is
# read from a stored table.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt

warnings.filterwarnings("ignore", category=UserWarning)
mt.versions()

# %% [markdown]
# ## A single couple
#
# Couples are always written `(reduced, oxidized)`. The library works out which
# element changes oxidation state and balances the half reaction: the redox
# element first, then oxygen with water, hydrogen with protons, and charge with
# electrons.

# %%
from microbial_thermo.reaction import Couple, HalfReactionResult

conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

couple = Couple.make("HS-", "SO4-2")
half = HalfReactionResult(
    couple=couple,
    half=couple.half_reaction(),
    conditions=conditions,
    backend=mt.get_backend(),
)

print("as a reduction:", half.half.format("reduction"))
print("as an oxidation:", half.half.format("oxidation"))
print("electrons transferred:", half.n_electrons)
print(f"E0  = {half.E_standard.to('V').magnitude:+.3f} V")
print(f"E0' = {half.E_standard_prime.to('V').magnitude:+.3f} V")

# %% [markdown]
# `E0` is at unit activity of everything, including $[\mathrm{H^+}] = 1$ M.
# `E0'` holds everything at unit activity but puts protons at the pH given in
# the conditions. The gap between them is the proton stoichiometry at work.
#
# A couple with no protons in it shows no gap at all:

# %%
for reduced, oxidized in [("HS-", "SO4-2"), ("Fe+2", "Fe+3"), ("H2(g)", "H+")]:
    pair = Couple.make(reduced, oxidized)
    result = HalfReactionResult(
        couple=pair, half=pair.half_reaction(),
        conditions=conditions, backend=mt.get_backend(),
    )
    print(
        f"{str(pair):22s} E0 = {result.E_standard.to('V').magnitude:+.3f} V   "
        f"E0' = {result.E_standard_prime.to('V').magnitude:+.3f} V"
    )

# %% [markdown]
# ## Oxidation states
#
# The figures report the *mean* oxidation state of the redox-active element,
# as exact fractions rather than rounded decimals.

# %%
from microbial_thermo.oxidation import format_oxidation_state

for element, formula in [
    ("S", "SO4-2"), ("S", "HS-"), ("N", "NO3-"), ("C", "CH4"),
    ("C", "C2H3O2-"), ("Fe", "Fe3O4"), ("Mn", "MnO2"),
]:
    state = mt.mean_oxidation_state(element, formula)
    print(f"{element} in {formula:10s} {format_oxidation_state(state):>6s}")

# %% [markdown]
# For organic molecules the per-atom states are available too, from the bond
# graph. Acetate's two carbons are chemically distinct — the methyl carbon is
# at −3 and the carboxyl at +3 — even though their mean is 0.

# %%
print("acetate per-atom:", mt.atom_oxidation_states("CC(=O)[O-]"))
print("acetate NOSC (the mean):", mt.nosc("C2H3O2-"))

# %% [markdown]
# ## A whole reaction
#
# Pair a donor with an acceptor. The donor runs in the oxidative direction, the
# acceptor in the reductive one, and both are scaled so the electrons cancel.

# %%
reaction = mt.Reaction.from_couples(
    donor=("H2(g)", "H+"),
    acceptor=("HS-", "SO4-2"),
    conditions=conditions,
    normalize_to="integer",
)
print(reaction.summary())

# %% [markdown]
# ## Show your work
#
# Every number can be traced back to its origin: which species contributed
# which formation energy, how each term of the reaction quotient was formed,
# and where the Nernst factor came from.

# %%
reaction.show_work()

# %% [markdown]
# ## The half-reaction figure
#
# Normalised to an electron pair by default, whatever the reaction was built
# with. Oxidation states sit above the top equation and below the bottom one;
# potentials are on the right.

# %%
from microbial_thermo.figures import plot_half_reactions

figure = plot_half_reactions(reaction)
plt.show()

# %% [markdown]
# ### Exercise
#
# Draw the half-reaction diagram for anaerobic methane oxidation coupled to
# sulfate reduction, then work out from the figure how many electrons move per
# molecule of methane.

# %%
aom = mt.Reaction.from_couples(
    donor=("methane", "CO2(aq)"),
    acceptor=("HS-", "SO4-2"),
    conditions=conditions,
)
plot_half_reactions(aom)
plt.show()
