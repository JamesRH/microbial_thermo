# ---
# jupyter:
#   jupytext:
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

# %%
import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.balance import infer_couples
from microbial_thermo.figures import plot_energy_explorer, plot_half_reactions, plot_redox_tower
from microbial_thermo.supplemental import supplemental_species
entry = supplemental_species()["NH2OH(aq)"]

from microbial_thermo.reaction import Couple, Reaction



# %%
from microbial_thermo.reaction import Couple, Reaction

conditions = mt.Conditions(temperature_c=25.0, pH=7)

reaction = Reaction.from_couples(
    donor=Couple.make("H2(aq)", "H+"),
    acceptor=Couple.make(["VOSO4(aq)"], ["VO2+", "SO4-2"], key_element="V"),
    conditions=conditions,
    normalize_to="integer",
)
print(reaction.summary())

# %%
p = plot_redox_tower(reaction)

# %%
# What the acceptor couple's oxidized side holds, which the next cell indexes.
reaction.acceptor_half.couple.oxidized_side

# %%
Elist = [reaction.acceptor_half.couple.oxidized_side[0][0].backend, 
reaction.acceptor_half.couple.oxidized_side[1][0].backend, 
reaction.acceptor_half.couple.reduced.backend, 
reaction.donor_half.couple.reduced.backend]

# %%
from microbial_thermo.sweep import concentration_axis, partial_pressure_axis, ph_axis

axes = [ph_axis(low=1.0, high=13.0, n=31)] + [concentration_axis(E) for E in Elist]
figure = plot_energy_explorer(reaction, axes=axes)

figure.show()

# %%
