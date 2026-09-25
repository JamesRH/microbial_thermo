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
# # 12. Chromium and uranium
#
# Two short cases, chosen because they are the metals bioremediation is
# actually asked to deal with, and because each one shows a limit of these
# diagrams that the earlier notebooks did not.
#
# * **Chromium** has no elemental form in any database we carry, so its Frost
#   diagram has no honest zero — the slopes still mean something and the
#   heights do not.
# * **Uranium** has one step, and that step does not move with pH at all,
#   which is unusual enough to be worth explaining.

# %%
import warnings

import matplotlib.pyplot as plt

import microbial_thermo as mt
from microbial_thermo.figures import frost_diagram, latimer_diagram, plot_frost, plot_latimer, plot_pourbaix

warnings.filterwarnings("ignore", category=UserWarning)

# %% [markdown]
# ## Chromium: the diagram that says what it does not know

# %%
ladder = latimer_diagram("Cr", pH=7.0)
ladder.to_frame()

# %%
plot_latimer("Cr", diagram=ladder)
plt.show()

# %%
frost = frost_diagram("Cr", pH=7.0)
print("zero point used:", frost.reference)
print("is it a genuine elemental form?", frost.reference_is_element)

# %% [markdown]
# `False`. Chromium metal is in none of our databases, so the reference falls
# back to the first species offered — Cr²⁺ — and the y axis of the Frost
# diagram becomes an arbitrary offset. The figure prints that in red rather
# than quietly drawing a plausible-looking diagram:

# %%
plot_frost("Cr", pH=7.0, annotate_slopes=True)
plt.show()

# %% [markdown]
# **What survives and what does not.** A volt equivalent is a height, so every
# height on that figure is meaningless. A couple potential is a *difference*
# of heights, so every slope is exactly right — the arbitrary offset cancels.
# The Latimer diagram, which only ever shows differences, is unaffected.
#
# This is worth internalising because it is the general rule: on an Eh–pH
# diagram the element reference cancels everywhere, on a Frost diagram it sets
# the zero, and on a Latimer diagram it never appears.

# %% [markdown]
# ### The chemistry that matters

# %%
for ph in (2.0, 5.0, 7.0, 9.0):
    steps = latimer_diagram("Cr", pH=ph).steps
    chain = "  ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps)
    print(f"pH {ph:4.1f}   {chain}")

# %%
plot_pourbaix("Cr", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# Chromate — Cr(VI) — is a strong oxidant, soluble at every pH, and the
# carcinogenic form. Cr(III) is the benign one and is nearly insoluble as the
# hydroxide, though our species list carries only the aqueous ion. The
# remediation strategy is simply to reduce it, and the potential says that is
# easy: +0.36 V at pH 7 against almost any donor a sediment has.
#
# The pH dependence is steep, though. In acid the same couple is worth over a
# volt, which is why chromate is far more aggressive in acidic waste.

# %% [markdown]
# ## Uranium: one step, and it ignores pH

# %%
for ph in (2.0, 5.0, 7.0, 9.0):
    steps = latimer_diagram("U", pH=ph).steps
    chain = "  ".join(f"{s.oxidized.backend} –({s.potential:+.3f})→ {s.reduced.backend}" for s in steps)
    print(f"pH {ph:4.1f}   {chain}")

# %% [markdown]
# Identical at every pH, which is rare enough to be worth checking rather than
# assuming. The reason is in the half reaction:

# %%
ladder = latimer_diagram("U", pH=7.0)
for step in ladder.steps:
    print(step.half_reaction())

# %% [markdown]
# **No protons appear.** Uranyl carries its two oxygens into uraninite
# unchanged, so there is no H⁺ term, and with no H⁺ term there is no pH
# dependence — the −59 mV per proton per electron that moves every other
# diagram in this series has nothing to act on.
#
# That has a practical consequence: uranium bioremediation works over the
# whole pH range of natural waters, which is not true of chromium.

# %%
plot_latimer("U", diagram=ladder)
plt.show()

# %%
plot_pourbaix("U", figsize=(8.0, 6.0))
plt.show()

# %% [markdown]
# Uraninite — UO₂, the ore mineral — takes the entire lower half of the
# diagram. Uranium is the clean case of "reduction immobilises": U(VI) is
# soluble and mobile, U(IV) is a solid, and the boundary sits right in the
# middle of the range a sediment can reach.

# %% [markdown]
# ## Who lives on each step

# %%
conditions = mt.Conditions(temperature_c=25.0, pH=7.0)
table = mt.energy_table(conditions).set_index("name")
table.loc[
    [
        "chromate_reduction_acetate",
        "chromate_reduction_hydrogen",
        "uranium_reduction_acetate",
    ]
][["label", "dG per e- (kJ/mol)", "reaction"]]

# %% [markdown]
# Chromate on either donor pays like a respiration, because that is what it
# is: *Shewanella* and *Desulfovibrio* reduce Cr(VI) as an electron acceptor,
# not as a detoxification. The engineering problem is delivery, not energy.
#
# There is a catch that the thermodynamics cannot show, and it is worth
# stating because it is where this kind of analysis stops being sufficient.
# Cr(VI) crosses cell membranes on the **sulfate transporter** — chromate and
# sulfate are near-identical in size and charge — so it gets inside, and once
# inside, partial reduction produces Cr(V) and Cr(IV) radicals that damage
# DNA. The reaction that makes chromium safe in a sediment is the same
# reaction that makes it carcinogenic in a cell.
#
# Nothing on a Frost or Eh–pH diagram will ever tell you that. They constrain
# what is possible; toxicity lives in the kinetics and the transport.

# %% [markdown]
# ## One thing worth taking away
#
# Both of these are real bioremediation strategies and both are thermodynamic
# gifts: the organisms are respiring a contaminant that happens to become
# insoluble when reduced, and they would do it whether or not anyone wanted
# them to. The engineering is entirely about delivering electrons.
#
# Compare arsenic in notebook 10, where exactly the same microbial process
# does exactly the wrong thing.
