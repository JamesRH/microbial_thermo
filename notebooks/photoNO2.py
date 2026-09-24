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
# # Nitrite-driven photoautotrophy
#
# Nitrite as the electron donor for photosynthetic CO₂ fixation. The
# hardest of the three: nitrate/nitrite sits at +0.407 V, so the pull
# against CO₂ is the longest.
#
# The pattern is the same in all three of the `photo*` notebooks: build the
# CO₂-fixation reaction, confirm it does **not** pay in the dark, then add the
# energy of absorbed photons and see what it takes to close the gap.
#
# Two conventions to keep straight:
#
# * A couple is always written **(reduced, oxidized)**. So the CO₂/biomass
#   acceptor is `Couple.make("Biomass", "CO2")` — biomass is the reduced form.
#   Writing it the other way round still gives the right energy, because the
#   direction comes from which couple is the donor and which the acceptor, but
#   it swaps the `reduced`/`oxidized` labels and the oxidation-state
#   annotations.
# * The reaction direction is CO₂ fixation because `acceptor=` is the couple
#   being reduced. Nothing else needs doing to get that direction.

# %%
import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.figures import (
    light_budget_table,
    plot_energy_explorer,
    plot_light_budget,
    plot_redox_tower,
)
from microbial_thermo.photons import photon_energy, photons_required
from microbial_thermo.reaction import Couple, Reaction
from microbial_thermo.sweep import concentration_axis, ph_axis

conditions = mt.Conditions(temperature_c=25.0, pH=7)

# %% [markdown]
# ## The biomass placeholder
#
# `<CH2O>` is cell carbon at oxidation state zero, and it is a **placeholder**,
# not a measured compound — it lives in the hand-entered supplemental table and
# warns on every use. Its value is this library's glucose divided by six, so
# writing the acceptor as `glucose` instead gives the same energy per carbon.

# %%
from microbial_thermo.supplemental import supplemental_species

supplemental_species()["Biomass(aq)"]

# %% [markdown]
# ## The reaction
#
# NO₂⁻ → NO₃⁻ is a two-electron step at +0.407 V at pH 7 — a much
# weaker donor than either arsenite or Fe(II).

# %%
reaction = Reaction.from_couples(
    donor=Couple.make("Nitrite", "Nitrate"),
    acceptor=Couple.make("Biomass", "CO2"),
    conditions=conditions,
)
print(reaction.summary())

# %% [markdown]
# ## It does not pay in the dark
#
# A positive ΔG is the whole point — this is the reaction a phototroph needs
# light for.

# %%
energy = reaction.delta_G_standard_prime
print(f"dG0' = {energy.magnitude:+.1f} kJ/mol per {reaction.n_electrons} e-")
print(f"uphill: {energy.magnitude > 0}")

# %%
plot_redox_tower(reaction)
plt.show()

# %% [markdown]
# ## Which species the explorer should sweep
#
# `Couple` itself has no `backend` — the backend name lives on the `Species`
# at `couple.oxidized` and `couple.reduced`. For a couple whose side names
# several species, `oxidized_side` is a tuple of `(Species, coefficient)`
# pairs, so index into that instead.

# %%
Elist = [
    reaction.acceptor_half.couple.oxidized.backend,
    reaction.acceptor_half.couple.reduced.backend,
    reaction.donor_half.couple.oxidized.backend,
    reaction.donor_half.couple.reduced.backend,
]
Elist

# %% [markdown]
# The general form, which also survives a multi-species side and drops the
# things there is no point sweeping:

# %%
sweepable = [
    s.backend for s in reaction.coefficients if s.backend not in ("H+", "H2O", "e-")
]
sweepable

# %%
axes = [ph_axis(low=0.0, high=14.0, n=29)] + [
    concentration_axis(name, n=13) for name in sweepable if name != "H2O"
]
figure = plot_energy_explorer(reaction, axes=axes)
figure.show()

# %% [markdown]
# ## Adding the energy from phototrophy
#
# A mole of photons carries $N_A hc/\lambda$. Purple bacteria use
# bacteriochlorophyll *a*, absorbing near 870 nm, which is 137.5 kJ/mol.

# %%
print(f"P870 photon: {photon_energy('P870').magnitude:.1f} kJ/mol")
print(f"photons needed here: {photons_required(energy):.2f}")

# %% [markdown]
# `photons_required` asks how many photons' worth of energy the reaction is
# short by, measured down to the biological energy quantum rather than to
# zero. It is a **floor**, not a quantum requirement: real anoxygenic
# phototrophs run cyclic electron flow and reverse electron transport and
# spend several photons per electron.

# %%
plot_light_budget(reaction)
plt.show()

# %% [markdown]
# ## All three donors side by side

# %%
others = {
    "photoarsenotrophy": Couple.make("As(III)", "As(V)"),
    "photoferrotrophy": Couple.make("Fe++", "Fe(OH)3", key_element="Fe"),
    "nitrite photoautotrophy": Couple.make("Nitrite", "Nitrate"),
}
built = {
    name: Reaction.from_couples(
        donor=couple, acceptor=Couple.make("Biomass", "CO2"), conditions=conditions
    )
    for name, couple in others.items()
}
light_budget_table(built)

# %% [markdown]
# This one needs 1.28 photons' worth at pH 7, so a single 870 nm photon
# per electron pair does *not* cover it. That is the sharpest contrast
# among the three donors.
