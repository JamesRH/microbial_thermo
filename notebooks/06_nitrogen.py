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
# # 6. Nitrogen
#
# Nitrogen has the widest span of oxidation states any element uses
# biologically — eight electrons, from ammonium at −III to nitrate at +V — and
# unlike carbon it has a *gas* sitting in the middle of the range. That one
# fact organises the whole nitrogen cycle: every path down from nitrate and
# every path up from ammonium runs into N₂, and N₂ is where nitrogen stays.

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

import microbial_thermo as mt
from microbial_thermo.figures import (
    frost_diagram,
    interactive_element,
    latimer_diagram,
    plot_latimer,
    plot_pourbaix,
    plot_redox_tower,
    pourbaix_field,
)
from microbial_thermo.figures.basis import element_series
from microbial_thermo.speciation import pKa_ladder

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder
#
# Every form of nitrogen this notebook uses, including both members of each
# acid-base pair; at a given pH each ladder and diagram keeps whichever one
# predominates.

# %%
NITROGEN_SPECIES = [
    "N2(aq)",
    "NO3-",
    "NO2-",
    "HNO2(aq)",
    "NO(aq)",
    "N2O(aq)",
    "NH4+",
    "NH3(aq)",
    "NH3OH+",
    "NH2OH(aq)",
    "N2H5+",
    "N2H4(aq)",
]
ladder = latimer_diagram("N", species=NITROGEN_SPECIES, pH=7.0)
ladder.to_frame()

# %%
plot_latimer("N", diagram=ladder)
plt.show()

# %% [markdown]
# **Every intermediate is on it, and they come from four places.** Not all of
# them are equally trustworthy, so here is each one's source:

# %%
for name in ("NO(aq)", "N2O(aq)", "N2H5+", "N2H4(aq)", "NH2OH(aq)", "NH3OH+"):
    record = mt.get_backend().resolve(name)
    print(f"{name:10s} {record.source}")

# %% [markdown]
# **NO and N₂O**, the two denitrification intermediates, are in none of the
# databases pyGCC bundles. They come from CHNOSZ's OBIGT database, and through
# it from Schulte, Shock and Wood (2001). Until they were imported, this
# ladder ran straight from nitrite to N₂ in one three-electron arrow.
# **Hydrazinium**, N₂H₅⁺, is the protonated hydrazine that anammox bacteria
# make, and the form that dominates below pH ≈ 8. It is in pyGCC's own HKF
# data, and was all along: an earlier version of this notebook reported it
# missing, because the library had never registered the name.
# **Hydroxylamine**, the intermediate of every ammonia oxidiser, is in no
# database at all; its value is typed in by hand and unverified.
#
# Hydroxylamine and hydrazine each also have a second form that takes over at
# one end of the pH scale, and neither of those is in any database: neutral
# **N₂H₄** above pH 8, and **NH₃OH⁺** below pH 6. Both are derived here from
# their partner and a pKa. Without them a pH slider would draw the minority
# form at each end and present it as the real one. All three hand-made values
# are unverified, and every figure that uses one says so.
#
# **Denitrification's four steps are wildly unequal.** Nitrite to NO is the
# weakest of the four, +0.24 V. NO to N₂O and N₂O to N₂ are
# both above +1.2 V. That is higher than O₂/H₂O at +0.815 V, so per electron
# N₂O is a better acceptor than oxygen.
#
# **Published values for these couples are usually for the gases**, and the
# ladder uses the dissolved forms. Dissolved NO sits 15.5 kJ/mol above the
# gas, which moves the nitrite step by 0.16 V:

# %%
for no in ("NO(aq)", "NO(g)"):
    couple = mt.half_reaction(no, "NO2-", conditions=mt.Conditions(pH=7.0))
    print(f"NO2-/{no:6s} E°' = {couple.E_standard_prime.to('V').magnitude:+.3f} V")

# %% [markdown]
# The red outlines are a *local* test: each rung against its two neighbours.
# Nitrite and hydrazinium pass it and still fall off the Frost hull below,
# because the species they actually fall apart into are not their neighbours.

# %% [markdown]
# ## What is stable

# %% [markdown]
# The Frost diagram below is interactive. Move pH and dissolved activity; the
# filled markers are the form that predominates at each oxidation state, and
# the open grey squares are the other forms of that state. It is fixed at
# 25 °C because hydroxylamine and the two derived forms have no temperature
# dependence. Spanning temperatures would mean dropping them. The
# temperature-capable version, without them, is in the *Interactive* section
# further down.

# %%
nitrogen_forms = element_series("N", species=NITROGEN_SPECIES)
interactive_element(
    "N",
    species=NITROGEN_SPECIES,
    temperatures=(25.0,),
    panels=("frost",),
    initial_ph=7.0,
    initial_log_activity=0.0,
    figsize=(9.0, 7.0),
)

# %%
frost = frost_diagram("N", species=NITROGEN_SPECIES, pH=7.0)
print("most stable form:", frost.most_stable.backend)
for event in frost.disproportionation():
    print(event)

# %% [markdown]
# Three results there, and all three matter.
#
# **N₂ is the bottom of the curve at pH 7.** Dinitrogen is where nitrogen goes
# and stays; that is the thermodynamic statement behind "fixed nitrogen is
# scarce". Nitrogenase spends sixteen ATP per N₂ because it is climbing out of
# this minimum, not because the enzyme is badly designed.
#
# **Nitrite is above the hull by a wide margin** — about 62 kJ per mole of
# nitrite at pH 7 — and would fall apart into nitrate and N₂ on its own. It is
# not the least stable intermediate, though. With every intermediate on the
# diagram it is one of the mildest; only N₂O, at about 56, sits closer to the
# hull. The organisms that oxidise or reduce nitrite are still exploiting a
# species that had nowhere good to sit.
#
# **Hydroxylamine, NO and hydrazinium sit higher above the hull than nitrite
# does**: hydroxylamine by about 221 kJ per mole of nitrogen, NO by 110,
# hydrazinium by 85. Hydroxylamine and hydrazinium fall to N₂ and ammonium, NO
# and N₂O to N₂ and nitrate. Those energies,
# like everything on a Frost diagram, are per mole of N; hydrazinium carries
# two, so per mole of N₂H₅⁺ its fall is about 169 kJ. Hydroxylamine and
# hydrazine are intermediates that nitrogen-cycle organisms make on purpose
# and never let accumulate: hydroxylamine in ammonia oxidisers, hydrazine in
# anammox. Their
# height is the point. A species far above the hull has a long way to fall,
# and anammox takes hydrazine's fall as four low-potential electrons, from a
# couple well below H⁺/H₂:

# %%
for reduced, oxidized in (("N2H5+", "N2(g)"), ("H2(g)", "H+")):
    couple = mt.half_reaction(reduced, oxidized, conditions=mt.Conditions(pH=7.0))
    print(f"{oxidized}/{reduced:6s} E°' = {couple.E_standard_prime.to('V').magnitude:+.3f} V")

# %% [markdown]
# **Two of the three results hold at any pH; the first does not.** Slide from
# 0 to 14 and the *form* at each oxidation state changes, at each family's
# pKa, but the hull is ammonia, N₂ and nitrate throughout and every
# intermediate stays above it. What moves is the bottom of the curve:

# %%
for family in ("nitrite", "hydroxylamine", "hydrazine", "ammonia"):
    print(f"{family:13s} pKa {pKa_ladder(family)[0]:5.2f}")
print()
for ph in (0.0, 3.0, 7.0, 11.0, 14.0):
    diagram = frost_diagram("N", series=nitrogen_forms, pH=ph)
    hull = ", ".join(point.backend for point in diagram.stable)
    falls = ", ".join(f"{e.species} {e.delta_g:+.0f}" for e in diagram.disproportionation())
    print(
        f"pH {ph:4.1f}  lowest: {diagram.most_stable.backend:7s} hull: {hull:24s} "
        f"above it (kJ/mol N): {falls}"
    )

# Each energy is base + protons * RT ln10 * pH, so two of them cross where
# those straight lines meet.
ammonium, dinitrogen = (
    next(entry for entry in nitrogen_forms if entry.backend == name) for name in ("NH4+", "N2(aq)")
)
crossover = (dinitrogen.base - ammonium.base) / (
    (ammonium.protons - dinitrogen.protons) * nitrogen_forms.rt_ln10
)
print(f"\nNH4+ and N2(aq) are equally low at pH {crossover:.2f}")

# %% [markdown]
# Below pH 3.9 the lowest point is ammonium, not N₂. That is the nitrogen
# Frost diagram as textbooks print it, at pH 0, with NH₄⁺ at −0.82 V in the
# basin. Forming ammonium from N₂ consumes four protons per nitrogen, so acid
# pushes it forward: each pH unit lifts ammonium by about 0.24 V on this axis
# (4 × 0.059), while N₂ does not move. The "N₂ is where nitrogen stays" argument above is
# a statement about neutral water, which is where it is needed.

# %% [markdown]
# The zero of that diagram is **N₂ gas at 1 bar**, and the species plotted at
# state 0 is **dissolved** N₂, which sits a little above it:

# %%
print("reference:", frost.reference, "— genuinely elemental:", frost.reference_is_element)
dissolved = [p for p in frost.points if p.backend == "N2(aq)"][0]
print(f"N2(aq) volt equivalent: {dissolved.volt_equivalent:+.3f} V")

# %% [markdown]
# That 0.094 V is the free energy of dissolving nitrogen, in volt-equivalent
# units. It is small, it is real, and getting it backwards — using dissolved
# N₂ as the zero — puts every point on the nitrogen diagram 0.094 V low. The
# figure always states which zero it used.

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("N", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# Almost the whole diagram is N₂, with nitrate only at the top and ammonium
# only at the bottom. Compare that to the real ocean, which is full of nitrate
# in the presence of ammonium — two species that this diagram says cannot
# coexist. They coexist because the reaction between them is slow without a
# catalyst, and **anammox bacteria are that catalyst**: they make a living on
# exactly the disequilibrium this figure says should not exist.

# %% [markdown]
# ## Interactive

# %% [markdown]
# All three panels, with temperature as well. Hydroxylamine and the two
# derived forms are left out here because they exist only at 25 °C. NO, N₂O
# and nitrous acid are kept, since they carry a full temperature dependence.

# %%
interactive_element(
    "N",
    species=["N2(aq)", "NO3-", "NO2-", "HNO2(aq)", "NO(aq)", "N2O(aq)", "NH4+", "NH3(aq)"],
    initial_ph=7.0,
    initial_log_activity=-5.0,
)

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
nitrogen = [
    "ammonia_oxidation",
    "nitrite_oxidation",
    "comammox",
    "denitrification_step_nitrate_nitrite",
    "denitrification_hydrogen",
    "dnra",
    "anammox",
    "n_damo",
    "nitrate_dependent_iron_oxidation",
    "sulfide_oxidation_nitrate",
    "arsenite_oxidation_nitrate",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[nitrogen][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# The last four rows are the ones worth dwelling on, because they are all the
# same trick: **nitrate and nitrite are good enough oxidants to be worth using
# on almost anything.** Methane, ferrous iron, sulfide and arsenite are all
# oxidised against them by organisms that do nothing else, and each of those
# couplings was discovered long after the nitrogen cycle was supposedly
# understood.
#
# `n_damo` is the sharpest example. Methane oxidation against **sulfate**
# pays about −3.7 kJ/mol e⁻ — inside the biological energy quantum, a
# metabolism that barely exists. Against **nitrite** the same methane pays
# −116:

# %%
pair = table.loc[["anaerobic_methane_oxidation", "n_damo"]]
pair[["label", "dG per e- (kJ/mol)"]]

# %% [markdown]
# Thirty times the yield, from changing nothing but the acceptor. Both
# organisms exist; the AOM consortia are slow, live in consortia, and double
# on a timescale of months, while *Methylomirabilis* grows on its own. The
# ladder explains both facts at once.
#
# *Methylomirabilis* has a further trick: it is thought to dismutate NO into
# N₂ and **O₂**, then use the oxygen in an ordinary methane monooxygenase —
# making its own oxidant inside an anoxic cell. The ladder says why that is
# possible at all. NO sits above the O₂/H₂O line, so it can oxidise water's
# oxygen, and splitting it runs downhill:

# %%
dismutation = mt.Reaction.from_couples(
    donor=("H2O", "O2(aq)"), acceptor=("N2(aq)", "NO(aq)"), conditions=conditions, normalize_to="O2(aq)"
)
print(f"{dismutation.delta_G_standard_prime.to('kJ/mol').magnitude:+.1f} kJ/mol   {dismutation}")

# %% [markdown]
# What the thermodynamics cannot see is that any enzyme does it. Worth
# remembering whenever a diagram says a reaction is impossible, or possible:
# the diagram constrains the overall energy, never the route.

# %% [markdown]
# ## One thing worth taking away
#
# DNRA and denitrification compete for the same nitrate, and the two diagrams
# say different things about them.
#
# Per electron, **denitrification wins outright** — about −110 kJ/mol e⁻
# against −75 for DNRA in the table above — and the Frost diagram says why:
# the slope from nitrate down to N₂ is +0.729 V, while the slope from nitrate
# down to ammonium is only +0.363 V. Ammonium is *not* the deeper point on the
# curve at pH 7; N₂ is.

# %%
print(f"NO3- to N2   {frost.slope('NO3-', 'N2(aq)'):+.3f} V")
print(f"NO3- to NH4+ {frost.slope('NO3-', 'NH4+'):+.3f} V")

# %% [markdown]
# So why does DNRA ever happen? Because energy per electron is not what is
# scarce in a carbon-rich sediment — **nitrate** is. Denitrification disposes
# of five electrons per nitrate; DNRA disposes of eight. When electrons are
# abundant and nitrate is not, the pathway that consumes more electrons per
# nitrate wins the competition despite paying less for each one.
#
# That is worth keeping in mind whenever a ladder is used to predict what an
# ecosystem does: the ladder ranks *yield*, and organisms compete over
# whatever is *limiting*. They are the same question only when the limiting
# resource is the electron donor.

# %% [markdown]
# ## Dirammox, and a first step outside water's reach
#
# *Alcaligenes ammonioxydans* oxidises ammonia straight to N₂ by way of
# hydroxylamine, a pathway called **Dirammox** (direct ammonia oxidation).
# With oxygen as the acceptor, the whole pathway and its second step run
# steeply downhill. The first step, in the middle line, does not:

# %%
overall = mt.Reaction.from_couples(
    donor=("NH4+", "N2(g)"), acceptor=("H2O", "O2(g)"), conditions=conditions, normalize_to="NH4+"
)
first = mt.Reaction.from_couples(
    donor=("NH4+","NH2OH"), acceptor=("H2O", "O2(g)"), conditions=conditions, normalize_to="NH4+"
)
second = mt.Reaction.from_couples(
    donor=("NH2OH", "N2(g)"), acceptor=("H2O", "O2(g)"), conditions=conditions, normalize_to="NH2OH"
)
for reaction in (overall, first, second):
    energy = reaction.delta_G_standard_prime.to("kJ/mol").magnitude
    print(f"{energy:+7.1f} kJ/mol   {reaction}")

# %% [markdown]
# All of the difficulty is in the first step, taking two electrons from
# ammonia to make hydroxylamine:
#
# $$\text{NH}_4^+ + \text{H}_2\text{O} \rightarrow \text{NH}_2\text{OH} + 3\,\text{H}^+ + 2e^-$$
#
# Any acceptor for those two electrons has to sit at a higher potential than
# this couple does. Here is the couple, beside the strongest oxidant that is
# stable in water:

# %%
for reduced, oxidized in (("NH4+", "NH2OH"), ("NH3(aq)", "NH2OH"), ("H2O", "O2(g)")):
    couple = mt.half_reaction(reduced, oxidized, conditions=conditions)
    print(f"E°' = {couple.E_standard_prime.to('V').magnitude:+.3f} V   {couple}")

# %% [markdown]
# **+0.898 V against +0.815 V, at pH 7.** An oxidant above the O₂/H₂O
# line would oxidise water itself, so no electron acceptor that is stable in
# water can take ammonia's electrons and leave hydroxylamine behind. Not
# nitrate, not Fe(III), not oxygen. Stronger oxidants do exist, peroxide
# among them, but only because they react with water slowly; none is a
# standing pool a respiratory chain can draw on.
#
# Mind which species the potential is written for. The figure usually quoted,
# +0.899 V, is NH₄⁺/NH₂OH, with three protons. Written for NH₃, the neutral form
# that ammonia monooxygenase uses, the same couple is +0.832 V. Both are right for the
# species they name; the difference is the energy of taking a proton off
# ammonium, shared over two electrons. At pH 7, where more than 99% of the
# ammonia is NH₄⁺, the first is the one that describes the solution. The
# second still sits above O₂/H₂O, but by only 17 mV.
#
# On the redox tower the problem is plain: the donor rung sits *below* the
# bottom of the tower, and the electron arrow points up. The tower here is
# interactive. At each pH it uses whichever forms of ammonia and of
# hydroxylamine predominate there, and redraws the reference couples to
# match. Under it, the Eh–pH diagram shows the whole of water's stability
# window, between the dashed lines. The black line is the ammonia /
# hydroxylamine couple drawn across pH, and the thin vertical line marks the
# slider.

# %%
import ipywidgets as widgets
from IPython.display import display


def first_step_at(ph):
    """Ammonia to hydroxylamine with oxygen, in the forms that exist at ``ph``."""
    predominant = frost_diagram("N", series=nitrogen_forms, pH=ph).predominant
    by_state = {round(point.oxidation_state): point.backend for point in predominant}
    return mt.Reaction.from_couples(
        donor=(by_state[-3], by_state[-1]),
        acceptor=("H2O", "O2(g)"),
        conditions=mt.Conditions(pH=ph),
    )


ph_values = np.arange(0.0, 14.01, 0.25)
hydroxylamine_line = [
    first_step_at(ph).donor_half.E_standard_prime.to("V").magnitude for ph in ph_values
]
nitrogen_field = pourbaix_field("N")


def draw_barrier(ph):
    tower = plot_redox_tower(first_step_at(ph))
    display(tower)
    plt.close(tower)
    figure, _ = plot_pourbaix("N", field=nitrogen_field, figsize=(8.0, 6.0))
    ax = figure.axes[0]
    ax.plot(ph_values, hydroxylamine_line, color="black", linewidth=2.2, label="ammonia / hydroxylamine")
    ax.axvline(ph, color="black", linewidth=1.0, alpha=0.6)
    ax.legend(loc="upper right", fontsize=9)
    display(figure)
    plt.close(figure)


ph_slider = widgets.FloatSlider(
    value=7.0, min=0.0, max=14.0, step=0.5, description="pH", continuous_update=False
)
barrier_output = widgets.Output()


def redraw_barrier(*_):
    with barrier_output:
        barrier_output.clear_output(wait=True)
        draw_barrier(ph_slider.value)


ph_slider.observe(redraw_barrier, names="value")
redraw_barrier()
widgets.VBox([ph_slider, barrier_output])

# %% [markdown]
# At pH 7, even with oxygen, the best acceptor biology has, the step costs
# about +16 kJ/mol. Every other acceptor on the tower is further up and costs
# more.
#
# The Frost diagram makes the same statement as a slope. The line from NH₄⁺
# up to NH₂OH is steeper than any oxidant stable in water can climb:

# %%
print(f"slope NH4+ to NH2OH  {frost.slope('NH2OH(aq)', 'NH4+'):+.3f} V")

# %% [markdown]
# ### Across pH: three plateaus, set by protonation
#
# Slide the pH and read the gap off the tower. The same numbers, in one place:

# %%
for ph in (0.0, 3.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0):
    step = first_step_at(ph)
    couple = step.donor_half.E_standard_prime.to("V").magnitude
    water = step.acceptor_half.E_standard_prime.to("V").magnitude
    energy = step.delta_G_standard_prime.to("kJ/mol").magnitude
    forms = f"{step.donor.reduced.backend} -> {step.donor.oxidized.backend}"
    print(
        f"pH {ph:4.1f}  {forms:20s} couple {couple:+.3f} V  O2/H2O {water:+.3f} V  "
        f"gap {1000 * (couple - water):+4.0f} mV  dG {energy:+5.1f} kJ/mol"
    )

# %% [markdown]
# Below pH 6 both sides of the couple are protonated, NH₄⁺ and NH₃OH⁺. The
# couple then moves one proton per electron, exactly as water does, so the
# two lines run parallel 114 mV apart: +22 kJ/mol to make hydroxylamine with
# oxygen. Between pH 6 and 9 hydroxylamine has lost its proton and ammonium
# has not, so the couple falls faster than water's and the gap narrows to
# 24 mV at pH 9. Above ammonium's pKa both are neutral, the lines are
# parallel again, and the gap is 17 mV, or +3.2 kJ/mol.
#
# **So the story does change, and it changes at the alkaline end.** At unit
# activity the first step lies outside water's window at every pH from 0 to
# 14. In acid it is outside by more than a tenth of a volt; in base, by less
# than two hundredths. Strong acid changes nothing: the acid plateau already
# begins at pH 6.
#
# ### Can a cell pull hydroxylamine low enough?
#
# Partly. Keeping the product scarce lowers the couple's potential by about
# 30 mV for every tenfold, while oxygen at air saturation, 0.21 bar, sits
# 10 mV below its standard line. With 1 mM ammonium:

# %%
air = {"O2(g)": 0.21}
gaps = {}
for hydroxylamine in (1e-3, 1e-4, 1e-5, 1e-6):
    at = mt.Conditions(
        pH=7.0, concentrations={"NH4+": 1e-3, "NH2OH(aq)": hydroxylamine}, partial_pressures=air
    )
    e_nitrogen = mt.half_reaction("NH4+", "NH2OH", conditions=at).E.to("V").magnitude
    e_oxygen = mt.half_reaction("H2O", "O2(g)", conditions=at).E.to("V").magnitude
    gaps[hydroxylamine] = e_nitrogen - e_oxygen
    print(
        f"NH2OH {hydroxylamine:.0e} M   E(NH2OH/NH4+) {e_nitrogen:+.3f} V   "
        f"E(O2/H2O) {e_oxygen:+.3f} V   gap {1000 * gaps[hydroxylamine]:+3.0f} mV"
    )

# At 1e-3 M the two are equal, so that gap is the one at a ratio of one.
per_decade = gaps[1e-3] - gaps[1e-4]
print(f"oxygen breaks even at NH4+/NH2OH = {10 ** (gaps[1e-3] / per_decade):,.0f}")

# %% [markdown]
# Only the ratio enters, not the ammonium level. Hydroxylamine has to be held
# more than a thousandfold below ammonium before oxygen even breaks even, and
# breaking even means zero energy, with nothing left to conserve.
#
# **How much of this rests on the unverified hydroxylamine value?** Each
# kJ/mol of error in it moves the couple by 5 mV. So trust only the order of
# magnitude of that ratio. The pH 7 conclusion is sturdier: the 83 mV gap
# would close only if the hand-entered value were wrong by 16 kJ/mol, and the
# acid plateau's 114 mV only if it were wrong by 22. The alkaline plateau's
# 17 mV is a different matter: about 3 kJ/mol, well inside what an untraced
# number could be off by. Above pH 9 or so, this notebook cannot tell you
# whether the first step is outside water's window at all.
#
# ### So how does anything make hydroxylamine?
#
# Not by handing ammonia's electrons to an acceptor. Ammonia monooxygenase
# in nitrifiers splits O₂ instead. One oxygen atom goes onto the nitrogen; the
# other is reduced to water by two electrons from the quinone pool, which
# come in turn from oxidising the hydroxylamine further down the pathway:
#
# $$\text{NH}_3 + \text{O}_2 + 2\text{H}^+ + 2e^- \rightarrow \text{NH}_2\text{OH} + \text{H}_2\text{O}$$
#
# The enzyme takes two electrons *in* rather than giving any out. Formally
# the nitrogen is still oxidised, by the oxygen atom it receives, and that
# transfer is still uphill. It is paid for inside the same reaction by the
# second oxygen atom, which falls all the way to water with two electrons
# from a far more negative donor. That is how the monooxygenase gets round
# the tower: not by finding a stronger acceptor, which water will not keep,
# but by spending a large downhill step to buy a small uphill one.
#
# Whatever enzyme makes hydroxylamine in a Dirammox organism faces the same
# arithmetic. At physiological pH it cannot be a dehydrogenase handing
# ammonia's electrons to a respiratory chain; it has to couple the step to
# something downhill, as the monooxygenase does. (The enzyme characterised
# so far, DnfA, works on the *second* step, hydroxylamine to N₂, and is
# itself a diiron oxygenase.) As with *Methylomirabilis* above, the
# diagram constrains the energy and not the route. Here, though, the energy
# rules out an entire class of route.
