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
# # 7. Sulfur
#
# Sulfur is the element where this library's diagrams earn their keep, because
# sulfur's most interesting metabolism is one that a textbook Frost diagram —
# drawn at unit activity, as they all are — says is impossible.
#
# The claim to test: **elemental sulfur disproportionation**,
#
#     4 S⁰ + 4 H₂O → 3 H₂S + SO₄²⁻ + 2 H⁺
#
# a real metabolism run by *Desulfocapsa* and others, in which sulfur is both
# the donor and the acceptor. It is famously marginal, and famously requires a
# sulfide sink to run at all. We can compute exactly how marginal.

# %%
import warnings

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np

import microbial_thermo as mt
from microbial_thermo.figures import (
    frost_diagram,
    latimer_diagram,
    plot_frost,
    plot_latimer,
    plot_pourbaix,
)

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## The ladder

# %%
ladder = latimer_diagram("S", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("S", diagram=ladder)
plt.show()

# %% [markdown]
# Every dissolved rung changes species between pH 0 and pH 7 — bisulfate to
# sulfate, SO₂(aq) to bisulfite, thiosulfuric acid to thiosulfate, H₂S to HS⁻
# — because sulfur's forms are acids with p*K*a values inside the
# environmental range. That is the reason a sulfur diagram quoted "at standard
# conditions" is nearly useless for a sediment:

# %%
for ph in (0.0, 2.0, 7.0, 9.0):
    steps = latimer_diagram("S", pH=ph).steps
    chain = "  ".join(
        f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps
    )
    print(f"pH {ph:4.1f}   {chain}")

# %% [markdown]
# ## Frost at unit activity: sulfur is stable

# %%
plot_frost("S", pH=7.0, activity=1.0, annotate_slopes=True)
plt.show()

# %%
textbook = frost_diagram("S", pH=7.0, activity=1.0)
print("on the hull:", [p.backend for p in textbook.stable])
for event in textbook.disproportionation():
    print(event)

# %% [markdown]
# Elemental sulfur is **on the hull**. At unit activity the diagram says
# disproportionation cannot happen, and every published sulfur Frost diagram
# says the same thing, because every published one is drawn at unit activity.
#
# Sulfite and thiosulfate are above the hull, which is also true and also
# worth knowing: both disproportionations are microbial metabolisms, and much
# easier ones. Thiosulfate sits at +2 because that is the *mean* of its two
# sulfur atoms; the diagram works per mole of sulfur, so a species with two
# of them is plotted at half its free energy.

# %% [markdown]
# ## Frost at a realistic activity: sulfur is not stable

# %%
plot_frost("S", pH=7.0, activity=1e-5, annotate_slopes=True)
plt.show()

# %%
real = frost_diagram("S", pH=7.0, activity=1e-5)
print("on the hull:", [p.backend for p in real.stable])
for event in real.disproportionation():
    print(event)

# %% [markdown]
# Elemental sulfur has left the hull. Drop the dissolved sulfur species from
# unit activity to 10 µM — an ordinary porewater — and disproportionation
# becomes exergonic by about 18 kJ per mole of sulfur.
#
# **Why lowering the activity changes the answer.** Sulfur is a solid, so its
# energy does not move; the dissolved products, sulfate and sulfide, both get
# cheaper as they get more dilute. The reaction makes four moles of dissolved
# product from a solid, so diluting the solution pulls it forward. This is Le
# Chatelier's principle, and it is the whole of the sulfur-disproportionation
# story.

# %% [markdown]
# ## Where the line is

# %%
activities = np.logspace(-6, 0, 25)
energies = []
for activity in activities:
    diagram = frost_diagram("S", pH=7.0, activity=activity)
    found = [e for e in diagram.disproportionation() if e.species == "Sulfur(s)"]
    energies.append(found[0].delta_g if found else 0.0)

figure, ax = plt.subplots(figsize=(8.0, 4.8))
ax.semilogx(activities, energies, linewidth=2.0, color="#1F6FB4")
ax.axhspan(-20.0, 0.0, color="#E8B84B", alpha=0.2, label="inside the biological energy quantum")
ax.axhline(0.0, color="#333333", linewidth=0.9)
ax.set_xlabel("dissolved sulfur activity")
ax.set_ylabel("ΔG of disproportionation\n(kJ per mol S)")
ax.set_title("4 S⁰ + 4 H₂O → 3 HS⁻ + SO₄²⁻ + 5 H⁺ at pH 7", fontsize=12)
ax.legend(frameon=False, fontsize=9)
plt.show()

# %% [markdown]
# The reaction turns exergonic at about 1.7 × 10⁻², but exergonic is not the
# same as *useful*: the curve stays inside the biological energy quantum — the
# shaded band, roughly the least free energy thought to sustain a cell — until
# the activity falls below about 5 × 10⁻⁶.
#
# That is the quantitative version of the textbook statement that sulfur
# disproportionators **need a sulfide sink**. Precipitating FeS keeps sulfide
# low; keeping sulfide low is what keeps this reaction worth running. The
# organism does not merely tolerate an iron-rich sediment, it requires one.

# %% [markdown]
# ## Where each form lives

# %%
plot_pourbaix("S", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# Note what is *not* on that diagram: elemental sulfur has no field at 10⁻⁶
# activity. It is metastable everywhere, which is the same result as the Frost
# curve in a different projection. Thiosulfate and sulfite have no field
# either, at this activity or at unit activity: they are intermediates, never
# the resting state.

# %% [markdown]
# ## The same thing, as a slider
#
# Everything above is one figure with one control. Drag **activity** and watch
# elemental sulfur leave the convex hull — the moment it does is the moment
# disproportionation starts paying.
#
# The diagram carries all three of sulfur's intermediate states: S⁰ at 0,
# thiosulfate at +2 and sulfite at +4, each in whichever acid–base form the pH
# selects.

# %%
from microbial_thermo.figures import interactive_element

interactive_element("S", initial_ph=7.0, initial_log_activity=-5.0)

# %% [markdown]
# Start at `log₁₀(activity) = 0` — the textbook diagram, sulfur on the hull —
# and drag left. Somewhere below −1.8 the S⁰ marker turns red and the note at
# the bottom of the figure names the reaction it would run.
#
# Then drag **pH** up: the disproportionation gets stronger, because the
# products include protons and a rising pH pulls them away. That is the same
# −59 mV per proton per electron that moves every other diagram in this
# series, showing up as a metabolism becoming viable.
#
# **Reading the energy off the picture.** A point above the hull can fall to
# the hull on either side of it — some of its sulfur oxidised, the rest
# reduced. How far it sits above the chord joining its two hull neighbours,
# times the Faraday constant, is the free energy of that disproportionation
# per mole of sulfur. That is the number printed in the note under the
# figure. Thiosulfate and sulfite are above the hull at every pH, activity and
# temperature the slider reaches; only S⁰ needs the activity pulled down
# before it pays:

# %%
for log_activity in (0, -3, -5, -7):
    diagram = frost_diagram("S", pH=7.0, activity=10.0**log_activity)
    print(f"activity 1e{log_activity:+d}")
    for event in diagram.disproportionation():
        print("   ", event)

# %% [markdown]
# Sulfite pays most at every activity, and thiosulfate always pays: these are
# the two that disproportionators were first grown on. Dilution helps S⁰ and
# thiosulfate, by the same Le Chatelier argument as before, because each makes
# more dissolved sulfur than it uses up. Sulfite's −57.6 kJ/mol S does not move
# at all once S⁰ is off the hull, since four dissolved sulfites make four
# dissolved products,
# $$4\,HSO_3^- \rightarrow 3\,SO_4^{2-} + HS^- + 3\,H^+$$
# and the activity cancels.

# %% [markdown]
# ## The step into S⁰ and the step out
#
# The slope of a Frost diagram between two forms *is* the potential of the
# couple joining them. That makes the diagram a map of how much each step of an
# oxidation pays, and it explains a second thing about sulfur: why some
# sulfide oxidisers — *Beggiatoa*, and purple sulfur bacteria such as
# *Allochromatium* — oxidise sulfide only as far as S⁰ at first and store it as
# globules, coming back for the rest later.
#
# A lower potential is a better electron donor, whatever the acceptor. With
# O₂, nitrate or a photosynthetic reaction centre, two steps differ in what
# they pay per electron by $F\,\Delta E$: 10 mV is 0.96 kJ/mol e⁻. So the
# question "which step pays more?" needs only the potentials. The steps of the
# two routes into a globule, and out of it:
#
# * **sulfide route**: HS⁻ → S⁰ (sulfide:quinone oxidoreductase), then
#   S⁰ → sulfite (reverse dissimilatory sulfite reductase, rDsr), then sulfite
#   → sulfate;
# * **thiosulfate route**: the Sox system without SoxCD releases one sulfur as
#   sulfate and leaves the other as S⁰,
#   $S_2O_3^{2-} + H_2O \rightarrow S^0 + SO_4^{2-} + 2H^+ + 2e^-$,
#   and the S⁰ goes out the same way as above.
#
# Concentrations matter, so they are stated: 1 mM total sulfide and sulfate, as
# in a sulfide-rich mat, and 1 µM sulfite and thiosulfate, which never
# accumulate. They are totals, so each couple is evaluated for whichever form
# the pH makes.

# %%
STEPS = {
    "HS⁻ → S⁰ (step in)": ("HS-", "Sulfur(s)"),
    "S⁰ → sulfite (rDsr, step out)": ("Sulfur(s)", "HSO3-"),
    "sulfite → sulfate": ("HSO3-", "SO4-2"),
    "S⁰ → sulfate (the rest, per e⁻)": ("Sulfur(s)", "SO4-2"),
    "S₂O₃²⁻ → S⁰ + SO₄²⁻ (Sox, step in)": (["S2O3-2"], ["Sulfur(s)", "SO4-2"]),
}
MAT = {"sulfide": 1e-3, "sulfate": 1e-3, "sulfite": 1e-6, "thiosulfate": 1e-6}


def step_potentials(pH, totals=MAT, temperature_c=25.0):
    """Potential (V) of each step in STEPS, at the given total concentrations."""
    conditions = mt.Conditions(
        temperature_c=temperature_c, pH=pH, total_concentrations=dict(totals)
    )
    return {
        name: mt.half_reaction(reduced, oxidized, conditions).E.to("V").magnitude
        for name, (reduced, oxidized) in STEPS.items()
    }


at_7 = step_potentials(7.0)
for name, potential in at_7.items():
    print(f"{name:36s} {potential:+.3f} V")
gap = at_7["S⁰ → sulfite (rDsr, step out)"] - at_7["HS⁻ → S⁰ (step in)"]
print(f"\nstep in beats step out by {1000 * gap:.0f} mV = {96.485 * gap:.1f} kJ/mol e⁻")

# %% [markdown]
# At pH 7 the step into S⁰ pays more per electron than the next step out of
# it, by a few kJ/mol e⁻. The cell takes the better step while sulfide lasts,
# and the step it puts off, S⁰ → sulfite, is the worst-paying step on the
# route. On a Frost diagram drawn at these concentrations, this is S⁰ sitting
# *below* the chord from HS⁻ to sulfite: a corner the curve bends up out of,
# which is the same geometry as
# "does not disproportionate into sulfide and sulfite". The disproportionation
# story above and the globule story are one picture read from either side. A
# point the curve bends *down* out of falls apart; a point it bends *up* out
# of is where electrons pile up.
#
# Two things qualify "the step to S⁰ pays more", and both are on the figure
# below:
#
# * **It is only true of the next step.** At these concentrations sulfite →
#   sulfate is the best-paying step of the three at every pH, and taken as a
#   whole, S⁰ → sulfate pays *more* per electron than HS⁻ → S⁰ above pH 4.6. The stored
#   sulfur still holds three quarters of the electrons and the most profitable
#   step. That is what makes a globule an energy reserve rather than waste.
# * **It is only true below about pH 7.8.** Above that the order of the
#   first two steps reverses, and lowering the sulfite moves the crossover
#   down (to about pH 7.2 at 10 nM).
#
# The thiosulfate route behaves the same way: the Sox step into S⁰ beats the
# rest of the way to sulfate below about pH 7.9 at 1 µM thiosulfate, and up
# to about pH 12.4 at 1 mM.
#
# Drag the concentrations and watch the crossovers move. All of this is at
# 25 °C.

# %%
PH_AXIS = np.arange(0.0, 14.01, 0.5)
STYLES = {
    "HS⁻ → S⁰ (step in)": {"color": "#1F6FB4", "linewidth": 2.4},
    "S⁰ → sulfite (rDsr, step out)": {"color": "#C0392B", "linewidth": 2.4},
    "sulfite → sulfate": {"color": "#7D7D7D", "linewidth": 1.6},
    "S⁰ → sulfate (the rest, per e⁻)": {"color": "#C0392B", "linewidth": 1.4, "linestyle": "--"},
    "S₂O₃²⁻ → S⁰ + SO₄²⁻ (Sox, step in)": {"color": "#2E8B57", "linewidth": 1.8},
}


def crossover(ph, first, second):
    """Every pH at which two potential curves cross, by linear interpolation."""
    difference = np.asarray(first) - np.asarray(second)
    found = []
    for i in np.nonzero(np.sign(difference[:-1]) != np.sign(difference[1:]))[0]:
        found.append(
            ph[i] - difference[i] * (ph[i + 1] - ph[i]) / (difference[i + 1] - difference[i])
        )
    return found


def plot_steps(log_sulfide=-3.0, log_sulfite=-6.0, log_thiosulfate=-6.0, log_sulfate=-3.0):
    totals = {
        "sulfide": 10.0**log_sulfide,
        "sulfite": 10.0**log_sulfite,
        "thiosulfate": 10.0**log_thiosulfate,
        "sulfate": 10.0**log_sulfate,
    }
    rows = [step_potentials(ph, totals) for ph in PH_AXIS]
    curves = {name: [row[name] for row in rows] for name in STEPS}

    figure, ax = plt.subplots(figsize=(9.0, 5.4))
    for name, values in curves.items():
        ax.plot(PH_AXIS, values, label=name, **STYLES[name])
    for ph in crossover(
        PH_AXIS, curves["HS⁻ → S⁰ (step in)"], curves["S⁰ → sulfite (rDsr, step out)"]
    ):
        ax.axvline(ph, color="#C0392B", linewidth=0.9, linestyle=":")
        ax.annotate(
            f"pH {ph:.1f}: step in = step out",
            xy=(ph, ax.get_ylim()[1]),
            xytext=(4, -12),
            textcoords="offset points",
            fontsize=9,
            color="#C0392B",
        )
    ax.set_xlabel("pH")
    ax.set_ylabel("E of the step (V)\nlower = pays more per electron")
    ax.set_title(
        f"Sulfide oxidation step by step, 25 °C — sulfide 1e{log_sulfide:g}, "
        f"sulfite 1e{log_sulfite:g}, thiosulfate 1e{log_thiosulfate:g}, "
        f"sulfate 1e{log_sulfate:g} M",
        fontsize=10,
    )
    ax.grid(color="#E5E5E5", linewidth=0.6)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    plt.show()


def concentration_slider(name, value):
    return widgets.FloatSlider(
        value=value,
        min=-9.0,
        max=-1.0,
        step=0.5,
        description=f"log₁₀ {name}",
        continuous_update=False,
        readout_format=".1f",
        style={"description_width": "initial"},
    )


widgets.interact(
    plot_steps,
    log_sulfide=concentration_slider("sulfide", -3.0),
    log_sulfite=concentration_slider("sulfite", -6.0),
    log_thiosulfate=concentration_slider("thiosulfate", -6.0),
    log_sulfate=concentration_slider("sulfate", -3.0),
)
# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
sulfur = [
    "hydrogenotrophic_sulfate_reduction",
    "acetotrophic_sulfate_reduction",
    "sulfur_reduction_hydrogen",
    "sulfide_oxidation_oxygen",
    "sulfur_oxidation_oxygen",
    "sulfide_oxidation_nitrate",
    "anaerobic_methane_oxidation",
]
table = mt.energy_table(conditions).set_index("name")
table.loc[sulfur][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# ### The disproportionations, as named metabolisms
#
# The three reactions this notebook has been circling are in the curated
# library, and they are the only entries in it where the donor and the
# acceptor are the same element:

# %%
table.loc[
    [
        "sulfur_disproportionation",
        "sulfite_disproportionation",
        "thiosulfate_disproportionation",
    ]
][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# **The +6.7 on the first row is the same result as the Frost diagram**, by a
# completely different route: the library balanced two half reactions and
# looked up formation energies, while the diagram compared a point to a chord
# on a convex hull. They agree because they are the same thermodynamics, and
# the agreement is worth checking rather than assuming:

# %%
standard = frost_diagram("S", pH=7.0, activity=1.0)
chord_energy = [e for e in standard.disproportionation() if e.species == "Sulfur(s)"]
print("Frost diagram says:", chord_energy or "S(0) is on the hull — nothing to disproportionate")
print(
    f"library says:       {table.loc['sulfur_disproportionation', 'dG per e- (kJ/mol)']:+.1f} kJ/mol e-"
)

# %% [markdown]
# The Frost diagram says nothing, because at unit activity elemental sulfur is
# *on* the hull; the library says +6.7 kJ/mol e⁻, which is the same statement
# with a sign: the reaction is uphill and will not run. Take the sulfide away
# and both change together.
#
# **Sulfite is the easy one** at about −39 kJ/mol e⁻, and thiosulfate sits
# between them at −3. That ordering is the reason sulfur disproportionators
# were found on sulfite and thiosulfate first: those two pay without needing
# anything else to happen, and elemental sulfur does not.
#
# ### Why an iron-rich sediment is the habitat
#
# Put the pieces together. Disproportionation needs the sulfide drawn below
# about 10⁻⁵ to clear the biological energy quantum. The thing that draws
# sulfide down in a real sediment is **iron**, precipitating FeS and
# ultimately pyrite. So the organism does not merely tolerate iron — it
# requires an iron supply, and the rate of its metabolism is set by a mineral
# reaction it has no control over.
#
# This is the general shape of a lot of biogeochemistry: an organism living on
# a reaction that is only viable because a *different* process is removing one
# of its products. Syntrophy (notebook 02) is the same argument with hydrogen
# instead of sulfide.

# %% [markdown]
# ## One thing worth taking away
#
# A Frost diagram is a picture of a *chosen* state, not of the element. The
# activity it was drawn at is part of the claim it makes, and for sulfur the
# choice flips the conclusion completely: stable at 1 M, disproportionating at
# 10 µM, and a living organism in the gap.
#
# Every diagram in this library prints its conditions in the title for exactly
# this reason.
