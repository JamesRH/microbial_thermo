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
# # 4. When the library refuses to guess
#
# Three places where the obvious guess is wrong often enough to matter, and
# where this library stops and asks instead. Each one would be easy to paper
# over, and papering over it would produce a number that looks entirely
# reasonable and is not.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt

warnings.filterwarnings("ignore", category=UserWarning)

conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

# %% [markdown]
# ## 1. A name that two species both claim
#
# `H2` is not one thing. There is dissolved hydrogen, which is what a cell
# actually sees, and hydrogen gas, which is what a headspace measurement
# reports. Both are in the database, and both answer to `H2`.
#
# The registry knows which names are contested:

# %%
from microbial_thermo.species import default_registry, resolve

registry = default_registry()
for name, (winner, claimants, declared) in sorted(registry.ambiguities().items()):
    others = ", ".join(c.backend for c in claimants if c is not winner)
    mark = "declared" if declared else "BY FILE ORDER"
    print(f"{name:5s} -> {winner.backend:12s} ({mark})   also claimed by: {others}")

# %% [markdown]
# Rather than let file ordering decide, each contested name has a declared
# winner in `species.yaml`. A bare name means the dissolved form:

# %%
print("H2   ->", resolve("H2").backend)
print("CH4  ->", resolve("CH4").backend)
print("H2 as a gas ->", resolve("H2", phase="g").backend)

# %% [markdown]
# This is not a bookkeeping nicety. The two forms are far enough apart that
# picking the wrong one would change a conclusion:

# %%
aqueous = mt.half_reaction("H2(aq)", "H+", conditions)
gas = mt.half_reaction("H2(g)", "H+", conditions)

gap_mv = 1000 * (
    aqueous.E_standard_prime.to("V").magnitude - gas.E_standard_prime.to("V").magnitude
)
print(f"H2(aq)/H+  E°' = {aqueous.E_standard_prime.to('V').magnitude:+.4f} V")
print(f"H2(g)/H+   E°' = {gas.E_standard_prime.to('V').magnitude:+.4f} V")
print(f"difference       {gap_mv:+.1f} mV")

# %% [markdown]
# A test asserts that every contested name has a declared answer, so adding a
# new aqueous/gas pair to the database fails loudly rather than quietly
# depending on where it was pasted into the file.

# %% [markdown]
# ## 2. An equation conservation cannot balance
#
# Most equations have exactly one balanced form, and the library finds it:

# %%
from microbial_thermo.balance import balance_equation, verify_conservation

coefficients = balance_equation("methane + O2(aq) -> CO2(aq) + H2O")
{s.backend: v for s, v in coefficients.items()}

# %% [markdown]
# But partial sulfide oxidation does not. Sulfide can go to sulfate, or to
# elemental sulfur, or to any mixture of the two — and conservation of mass
# and charge has nothing to say about which. Written with both products, the
# equation has a whole family of solutions:

# %%
from microbial_thermo.exceptions import AmbiguousReactionError

AMBIGUOUS = "HS- + O2 -> SO4-2 + Sulfur(s) + H2O + H+"

try:
    balance_equation(AMBIGUOUS)
except AmbiguousReactionError as error:
    print(error)

# %% [markdown]
# The error carries the solution basis — the two independent vectors whose
# combinations all conserve mass and charge. Read them as bookkeeping rather
# than as reactions: the signs say which way each species would have to move,
# and only some combinations put everything on a physically sensible side. The
# second vector has no elemental sulfur at all, and is complete oxidation.
#
# Note also that the backend names are legacy SUPCRT spellings — `SO4--`, not
# the `SO4-2` that was typed in. The registry accepts either.

# %%
try:
    balance_equation(AMBIGUOUS)
except AmbiguousReactionError as error:
    for index, option in enumerate(error.basis, start=1):
        verify_conservation(option)
        print(f"option {index}: " + str({s.backend: v for s, v in option.items()}))

# %% [markdown]
# Supplying the missing information resolves it. `fix=` states moles as
# positive numbers; the sign comes from the side the species was written on.

# %%
to_sulfur = balance_equation(AMBIGUOUS, fix={"HS-": 4, "Sulfur(s)": 2})
{s.backend: str(v) for s, v in to_sulfur.items()}

# %%
to_sulfate = balance_equation(AMBIGUOUS, fix={"HS-": 4, "Sulfur(s)": 1})
{s.backend: str(v) for s, v in to_sulfate.items()}

# %% [markdown]
# Both conserve every element and the charge — they are different *reactions*,
# not a right and a wrong answer. Compare how much oxygen each consumes per
# sulfide, which is the quantity a sulfide oxidiser is actually managing.
# Coefficients stay exact fractions throughout, never floats:

# %%
import pandas as pd

rows = []
for label, solution in [("half to S⁰", to_sulfur), ("mostly to SO₄²⁻", to_sulfate)]:
    verify_conservation(solution)
    by_name = {s.backend: v for s, v in solution.items()}
    rows.append(
        {
            "outcome": label,
            "HS⁻": -by_name["HS-"],
            "O₂ consumed": -by_name["O2(aq)"],
            "S⁰": by_name["Sulfur(s)"],
            "SO₄²⁻": by_name["SO4--"],
            "O₂ per HS⁻": by_name["O2(aq)"] / by_name["HS-"],
        }
    )
pd.DataFrame(rows)

# %% [markdown]
# Asking for too little still fails, and says how much is missing rather than
# picking something:

# %%
try:
    balance_equation(AMBIGUOUS, fix={"O2": 2})
except AmbiguousReactionError as error:
    print(error)

# %% [markdown]
# ## 3. An average that hides the structure
#
# The mean carbon oxidation state is the number that goes into the
# thermodynamics, and it is genuinely what matters for the energetics. It is
# also, read as chemistry, misleading.
#
# Acetate averages to zero. It contains no carbon at zero:

# %%
from fractions import Fraction

from microbial_thermo.oxidation import (
    format_per_atom_states,
    mean_oxidation_state,
    per_atom_states,
)

for name in ["acetate", "propionate", "butyrate", "glucose", "methanol"]:
    species = resolve(name)
    mean = mean_oxidation_state("C", species.formula)
    detail = (
        format_per_atom_states(species.smiles, "C") if species.smiles else "no SMILES"
    )
    print(f"{name:10s} {species.formula:8s} mean C {str(mean):>6s}   per atom: {detail}")

# %% [markdown]
# Acetate's two carbons sit at −3 and +3: a fully reduced methyl group and a
# fully oxidised carboxyl. Butyrate's four average to −1 while spanning −3 to
# +3. The mean is the right number for the energy balance and the wrong
# picture of the molecule.
#
# This matters biologically because it is the *methyl* carbon that becomes
# methane in acetoclastic methanogenesis and the carboxyl that becomes CO₂ —
# the molecule splits along exactly the division the average erases.

# %%
acetate = resolve("acetate")
print("acetate per-carbon states:", per_atom_states(acetate.smiles, "C"))
print("mean                     :", mean_oxidation_state("C", acetate.formula))
print(
    "check — the per-atom states average to the mean:",
    Fraction(sum(per_atom_states(acetate.smiles, "C")), 2),
)

# %% [markdown]
# The half-reaction figure will show either. `per_atom=True` annotates each
# species with its distinct carbon states instead of the average:

# %%
from microbial_thermo.figures import plot_half_reactions

acetate_oxidation = mt.Reaction.from_couples(
    donor=("Acetate", "CO2(aq)"),
    acceptor=("H2O", "O2(aq)"),
    conditions=conditions,
)
plot_half_reactions(acetate_oxidation, per_atom=True)
plt.show()

# %% [markdown]
# Species with no SMILES in the registry are simply left un-annotated rather
# than blocking the figure.

# %% [markdown]
# ## What these have in common
#
# In all three cases the library could have picked something. A bare `H2`
# could resolve by whichever entry came first; an underdetermined equation
# could return one arbitrary member of its solution family; an average could
# be reported without saying it is an average.
#
# Each of those would produce output that looks correct. The cost of asking is
# one extra argument; the cost of guessing is a plausible wrong number that
# nothing downstream can detect.
