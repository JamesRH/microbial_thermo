# microbial_thermo

Thermodynamic calculations and figures for microbial physiology and
biogeochemistry. Given a redox couple or a whole reaction, it balances the
chemistry, computes $E^\circ$, $E^{\circ\prime}$ and $\Delta G$ from a real
thermodynamic database, and draws the three figures a course needs: a stacked
half-reaction diagram, a redox tower, and an interactive free-energy explorer.

Built for teaching first. Every result can be asked to show its work.

**Backend:** [pyGCC](https://bitbucket.org/Tutolo-RTG/pygcc/) only.
**Range:** 0.01–100 °C, near-surface pressure.

---

## Installation

```bash
mamba env create -f environment.yaml
mamba activate microbial-thermo
pip install -e .
```

To build the environment by hand instead:

```bash
mamba create -n microbial-thermo python=3.11
mamba activate microbial-thermo
mamba install -c conda-forge numpy pandas scipy chempy sympy pint matplotlib seaborn plotly rdkit click jupytext pyyaml ipykernel ipywidgets ruff
pip install pygcc==1.5.3
```

`pygcc` is installed with `pip` because it is not published on conda-forge;
PyPI is its only distribution channel. Everything else comes from conda-forge.

### Jupyter

Register the environment as a kernel:

```bash
mamba activate microbial-thermo
python -m ipykernel install --user --name microbial-thermo --display-name "microbial-thermo"
```

Then pick **microbial-thermo** as the notebook kernel. Notebooks in
`notebooks/` are paired to `.py` scripts with jupytext in `percent` format;
edit the `.py` and run `jupytext --sync <notebook>.ipynb`.

| notebook | covers |
|---|---|
| `01_half_reactions` | couples, balancing, oxidation states, show-your-work, the half-reaction figure |
| `02_redox_tower_and_energy` | the tower and its scale bar, hydrogen and temperature sweeps, the interactive explorer |
| `03_environmental_affinity` | pH speciation, the curated library, and the affinity ladder for a real porewater |

All three are committed with their outputs and execute clean end to end.

---

## Library usage

```python
import microbial_thermo as mt

conditions = mt.Conditions(temperature_c=25.0, pH=7.0)

reaction = mt.Reaction.from_couples(
    donor=("H2(g)", "H+"),        # couples are always (reduced, oxidized)
    acceptor=("HS-", "SO4-2"),
    conditions=conditions,
    normalize_to="integer",
)
print(reaction.summary())
```

```
Reaction:  H+ + 4 H2(g) + SO4-- -> HS- + 4 H2O
  donor     H+/H2(g)  (oxidative)  E0' = -0.414 V
  acceptor  SO4--/HS-  (reductive)  E0' = -0.217 V
  n electrons        8
  dE0'               +0.197 V
  dG0'               -152.2 kJ/mol
  dG (conditions)    -152.2 kJ/mol
  dG per electron    -19.0 kJ/mol e-
```

### Realistic conditions

Anything not named is held at unit activity, so the standard state is the
default and you only specify what you actually measured.

```python
conditions = mt.Conditions(
    temperature_c=60.0,
    pH=8.2,
    concentrations={"SO4-2": 2.8e-2, "HS-": 1e-6},
    partial_pressures={"H2(g)": 1e-5},
    ionic_strength=0.7,
    activity_model="bdot",     # extended Debye-Huckel; "ideal" or "unit" also
)
```

### Show your work

```python
reaction.show_work()           # renders as Markdown + LaTeX in a notebook
print(reaction.show_work().to_text())
```

Eight numbered steps: identifying the couples and their oxidation-state
changes, balancing each half reaction, cancelling electrons, a table of every
formation energy with its source, the reaction quotient term by term, the
$RT\ln Q$ arithmetic, conversion to potentials, and the independent
cross-check.

### Normalisation

Default is an electron pair. Override with `normalize_to=`:

| value | meaning |
|---|---|
| `"electron_pair"` | 2 electrons (default) |
| `"electron"` | 1 electron |
| `"donor"` / `"acceptor"` | 1 mol of that half reaction |
| `"integer"` | smallest count giving whole-number coefficients |
| any species name | 1 mol of that species |

`"integer"` is usually what you want for a figure: it turns
`H2 + ¼CO2 → ¼CH4 + ½H2O` into `4 H2 + CO2 → CH4 + 2 H2O`.

### Oxidation states

```python
mt.mean_oxidation_state("S", "SO4-2")     # Fraction(6)   -- what figures show
mt.nosc("C6H12O6")                        # Fraction(0)   -- LaRowe & Van Cappellen
mt.atom_oxidation_states("CC(=O)[O-]")    # per-atom, acetate: C at -3 and +3
```

Exact `Fraction`s, so magnetite iron reports `+8/3` rather than `2.6666667`.

### pH speciation

Measurements are reported as totals — total sulfide, DIC, total ammonia — while
a reaction is written with one specific form. Give the total and the split is
handled for you:

```python
conditions = mt.Conditions(
    temperature_c=25.0, pH=7.0, activity_model="ideal",
    total_concentrations={"sulfide": 1e-6, "DIC": 2.2e-3},
)
```

At pH 6 only 9% of total sulfide is HS⁻; at pH 8 it is 91%. Naming one form and
treating a measured total as all of it is wrong by roughly a factor of two near
a p*K*a, so the library warns when you do that.

```python
mt.speciation_table("sulfide", 7.0)
#  species  fraction   percent  pKa to next
#  H2S(aq)  0.492942   49.294   6.99
#  HS-      0.507058   50.706
```

p*K*a values are **computed from the backend at the working temperature**, not
tabulated, so sulfide's p*K*a falls from 6.99 at 25 °C to 6.65 at 60 °C and the
distribution shifts with it. At 25 °C the computed values reproduce the
literature: carbonate 6.34 / 10.33, ammonium 9.24, phosphate 2.17 / 7.21 /
12.32, acetate 4.76.

The carbonate step is balanced with water — CO₂(aq) + H₂O ⇌ HCO₃⁻ + H⁺ — which
a proton-only treatment gets wrong by tens of p*K* units.

Twelve families ship: sulfide, carbonate, ammonia, phosphate, acetate, lactate,
formate, propanoate, butanoate, sulfite, nitrite, sulfate.

### Iron and manganese

Aqueous ions come from the HKF database; solid phases are derived from log K
values in pyGCC's GWB database, because no direct-access database it ships
contains the manganese oxides.

```python
mt.Reaction.from_couples(donor=("Fe+2", "goethite"), acceptor=("H2O", "O2(aq)"))
mt.Reaction.from_couples(donor=("acetate", "CO2(aq)"), acceptor=("Mn+2", "pyrolusite"))
```

Available: `Fe+2`, `Fe+3`, `Mn+2`, `Mn+3`, `MnO4-`, permanganate/manganate;
goethite, hematite, magnetite, siderite, pyrite, ferrihydrite; pyrolusite,
manganite, hausmannite, bixbyite, birnessite, rhodochrosite.

Goethite, hematite and magnetite exist in two independent databases, giving two
routes to the same number. They agree with each other and with published values
to within 1–3 kJ/mol, which is what licenses trusting the manganese oxides,
where only one route exists.

### The curated metabolism library

Twenty-eight named metabolisms, so you need not remember which couples to pair.
Nothing thermodynamic is stored — energies are computed at whatever conditions
you ask for.

```python
from microbial_thermo.library import default_library, energy_table, reaction

reaction("anammox", conditions)
energy_table(conditions)          # every metabolism, sorted by energy per electron
default_library().groups()        # methanogenesis, sulfate reduction, nitrification, ...
```

At pH 7 and 25 °C the table reproduces canonical redox zonation unprompted:
oxygen (−109 kJ/mol e⁻) > denitrification (−96) > Mn(IV) (−65) > Fe(III) (−13)
> sulfate (−5) > methanogenesis (−1.4). Acetoclastic methanogenesis and AOM land
at the famously marginal −1.4 and −3.7 kJ/mol e⁻. The zonation order is asserted
in the test suite: if it broke, the thermodynamics would be wrong.

*Limitation:* a couple names one reduced and one oxidized species, so
metabolisms whose oxidation yields two carbon products cannot be catalogued.
Syntrophic propionate oxidation (propionate → acetate + CO₂ + H₂) is the case
that matters. Complete oxidations to CO₂ are fine.

### Figures

```python
from microbial_thermo.figures import (
    plot_half_reactions, plot_redox_tower, plot_energy_explorer,
)

plot_half_reactions(reaction, save="figures/sulfate")     # SVG + PNG
plot_redox_tower(reaction, save="figures/tower")          # SVG + PNG
plot_affinity_ladder(conditions, save="figures/ladder")   # SVG + PNG
fig = plot_energy_explorer(reaction, save_html="figures/explorer")
```

The **affinity ladder** is the counterpart to the tower: where the tower shows
what is possible from standard potentials, the ladder shows what actually pays
under one measured water chemistry, with the energy-quantum band drawn across
it. Its x axis is kJ/mol, so ATP equivalents are a genuine second axis there.

Both static figures normalise to **an electron pair** by default, whatever
electron count the reaction was built with; pass `n_electrons=None` to draw it
as given, or another integer to rescale.

Pinning n is also what makes the tower's energy scale bar exact. Since
$\Delta G = -nF\Delta E$, a *difference* in potential converts to free energy
by the constant $2F = 192.97$ kJ·mol⁻¹·V⁻¹, so the bar marks 0.104 V as the
−20 kJ/mol energy quantum and 0.259 V as one ATP. You can read a reaction's
yield straight off the tower by comparing its donor–acceptor gap to the bar. An
absolute kJ/mol *axis* would be wrong — it would imply each couple has an
absolute free energy, when only differences carry energy.

The explorer writes a self-contained HTML file: dropdown, hover, and the SVG
download button all work with no Python process behind it, so you can hand the
file to students directly.

### Sweeps without plotting

```python
from microbial_thermo.sweep import sweep, partial_pressure_axis

result = sweep(reaction, partial_pressure_axis("H2(g)"))
result.to_frame()          # a pandas DataFrame
```

---

## Command line

```
$ mthermo --help
Usage: mthermo [OPTIONS] COMMAND [ARGS]...

  Thermodynamics for microbial physiology and biogeochemistry.

Options:
  --version        Show the library and dependency versions, then exit.
  --database TEXT  Path to an alternative pyGCC database file.
  --help           Show this message and exit.

Commands:
  couple      Balance one redox couple and report its potentials.
  figure      Render a figure to SVG and PNG, or to standalone HTML for...
  reaction    Compute a whole reaction from a donor and an acceptor couple.
  speciation  Show acid-base speciation at a given pH and temperature.
  species     List known species, optionally filtered by a substring.
  tower       Print the reference redox couples, ordered by potential.
```

Examples:

```bash
mthermo couple HS- SO4-2
mthermo reaction -d "H2(g)" "H+" -a HS- SO4-2 --normalize-to integer --show-work
mthermo reaction -d "H2(g)" "H+" -a methane CO2 -p "H2(g)=1e-5" --ph 6.8
mthermo figure -d acetate CO2 -a H2O O2 --kind tower -o figures/acetate
mthermo figure -d "H2(g)" "H+" -a HS- SO4-2 --kind explorer -o figures/explorer
mthermo tower --ph 7 -T 25
mthermo speciation                      # every family and its pKa ladder
mthermo speciation sulfide --ph 7 -T 60
mthermo species sulf
```

```
$ mthermo couple HS- SO4-2
reduction   SO4-- + 9 H+ + 8 e- -> HS- + 4 H2O
oxidation   HS- + 4 H2O -> SO4-- + 9 H+ + 8 e-
n electrons 8
S oxidation state 6 -> -2
E0          +0.249 V
E0'         -0.217 V
```

### Version tracking

```
$ mthermo --version
microbial_thermo 0.1.0
computational dependencies:
  pygcc        1.5.3
  numpy        2.4.6
  pandas       2.3.3
  scipy        1.17.1
  chempy       0.10.1
  sympy        1.14.0
  pint         0.25.3
  matplotlib   3.11.1
  plotly       7.0.0
```

---

## How it works, and why

### The two-path cross-check

Every `Reaction` computes its free energy twice on construction: once by
summing formation energies over the whole reaction, and once as $-nF\Delta E$
from the two half-reaction potentials. Disagreement raises
`ThermodynamicConsistencyError`. Redox sign errors are the characteristic
failure in this field — textbooks included — so the check runs on every
calculation, not only in the test suite.

### Working with pyGCC

pyGCC exposes no arbitrary-reaction API: `calcRxnlogK` only evaluates a
species' *database-defined* formation reaction. So this library assembles
$\Delta G_\mathrm{rxn} = \sum_i \nu_i \Delta G_{f,i}$ itself from per-species
values. Three consequences worth knowing:

- **Water is not a database species.** HKF parameterises solutes, not the
  solvent, so H₂O comes from `iapws95` instead.
- **Gases and minerals are not HKF solutes.** They carry Maier–Kelley
  coefficients and are routed through `heatcap`.
- **There is a cliff at 100 °C.** At exactly 100 °C, 1 bar is below the
  saturation pressure and water would be steam, so `supcrtaq` returns `NaN`.
  Pressure follows the saturation curve at and above 99 °C.

### Are the databases mutually comparable?

Numbers come from four places: `speq21.dat` (HKF aqueous solutes), `speq21.dat`
again via `heatcap` (gases and some minerals), `thermo.com.dat` (minerals the
others lack, through log K), and `iapws95` (water). The mineral route *mixes*
sources — log K from one file, basis-species energies from another — so it has
to be justified rather than assumed.

It was checked directly. For the 900-odd species present in **both** databases,
the tabulated log K was compared against the log K predicted from the HKF
parameters. The median discrepancy is 0.0000 log K and 809/939 agree within
0.1, so the two files are the same SUPCRT lineage and may be combined. The
basis species the mineral route actually leans on (Fe³⁺, HS⁻, CO₂(aq), OH⁻)
agree to better than 0.15 kJ/mol.

The audit found one real inconsistency, now corrected. IAPWS-95 puts liquid
water at −237.140 kJ/mol while the rest of the data assumes the SUPCRT
convention of −237.18 — a 9 cal/mol reference-state difference, not an error
(−237.14 is modern CODATA). Back-calculating the water energy implied by the
OH⁻, Fe³⁺ and CO₂(aq) reactions gives −56687.7 cal/mol every time. Water is
therefore shifted by a constant onto the SUPCRT datum, which preserves the
IAPWS-95 temperature dependence exactly and removes a systematic 0.041 kJ/mol
error *per mole of water* from every reaction. Pass
`PygccBackend(water_convention="iapws")` to opt out.

Three genuine inter-database disagreements survive and are documented rather
than hidden: rhodochrosite (2.5 kJ/mol — we use the `speq21` value, which is
nearer the published one), FeOH⁺ (1.1 kJ/mol), and the Mn(III)/Mn(VII) set
(0.4 kJ/mol). All are held to declared tolerances by
`tests/test_database_consistency.py`, so a database change that widened them
would fail the suite.

### Validation

Checked against independent published values:

| quantity | computed | published |
|---|---|---|
| $\Delta G_f$(H₂O) | −237.14 kJ/mol | −237.1 |
| $\Delta G_f$(SO₄²⁻) | −744.46 kJ/mol | −744.5 |
| $\Delta G_f$(acetate) | −369.32 kJ/mol | −369.3 |
| $E^{\circ\prime}$(H⁺/H₂) | −0.414 V | −0.414 |
| $E^{\circ\prime}$(O₂/H₂O) | +0.815 V | +0.816 |
| $E^{\circ\prime}$(SO₄²⁻/HS⁻) | −0.217 V | −0.217 |
| $E^{\circ\prime}$(Fe³⁺/Fe²⁺) | +0.770 V | +0.770 |
| $\gamma$(Na⁺) at I = 0.1 | 0.775 | ~0.77 |

Textbook potentials are asserted as ranges rather than equalities: published
tables differ by tens of millivolts between sources because of differing
standard states and database vintages, and asserting equality would encode one
textbook's conventions as truth.

---

## Development

```bash
python -m unittest discover -s tests     # 175 tests
ruff format microbial_thermo tests
ruff check microbial_thermo tests
```

---

## Known limitations

- **Glucose and pyruvate are absent** from every database pyGCC ships. Acetate,
  formate, lactate, propanoate, butanoate, methanol and ethanol are present.
  The supplemental-data mechanism for filling such gaps is specified in
  `SPEC.md` §2.1 but is **not yet implemented**.
- **Isothermal minerals.** A few phases — manganite and birnessite among them —
  are tabulated at 25 °C only in the source database. They evaluate at 25 °C and
  raise `OutOfRangeError` elsewhere rather than inventing a temperature
  dependence the data does not contain.
- **N₂O is absent**, so the full denitrification pathway cannot be stepped
  through species by species.
- **Neutral species are treated as ideal** ($\gamma = 1$) under the B-dot model,
  the usual Helgeson convention. Dissolved-gas activities are therefore not
  salted out.
- **Speciation acts on activities, not on a transformed standard state.**
  Giving a family total yields the correct activity of each form, which is what
  in-situ ΔG needs. It is not the same as eQuilibrator's Legendre-transformed
  ΔG′ over pseudoisomer groups, where protons are implicit and reactions are
  balanced without H⁺. Here H⁺ is always explicit and balanced. See Future work
  item 1.

---

## Future work

Ordered by logical dependency first, then by effort within each tier. Earlier
items unblock later ones.

### Tier 1 — small, and builds directly on what exists

1. **Ionic-strength convention for literature comparison** *(small, worthwhile)*.
   Alberty-convention biochemical tables are quoted at I = 0.25 M rather than
   I = 0. Defaulting `Conditions` to that when someone is comparing against such
   a table would close most of the remaining gap, and the B-dot machinery for it
   already exists.

   **Not** recommended: the full transformed (Legendre) standard state. It was
   evaluated and deliberately declined. The reasoning, so it does not get
   relitigated:

   - There are two conventions, not one. The *microbial bioenergetics*
     literature — Thauer, Jungermann & Decker; Amend & Shock; LaRowe & Amend —
     writes reactions with explicit species and explicit protons. **That is
     already what this library computes**, and our −152.2 kJ/mol for
     hydrogenotrophic sulfate reduction sits on the ≈−152 kJ/mol tabulated
     there. The *biochemical* convention (Alberty, IUBMB, eQuilibrator) is the
     one we do not match.
   - The two agree more than expected. For a reaction whose every reactant is a
     single species, hydrogen conservation forces the transformed result to
     equal our species-level ΔG°′ *exactly*; the transform is pure rebookkeeping
     there. They diverge only through pseudoisomer grouping, bounded by the
     mixing entropy $RT\ln(\text{populated forms})$ — at most ≈1.7 kJ/mol, for a
     diprotic reactant sitting exactly on its p$K_a$.
   - Ionic strength is the larger discrepancy, and item 1 above addresses it
     without touching anything structural.
   - The cost is not small: a second balancing path (dropping both the hydrogen
     and charge rows from the conservation matrix), a transformed analogue of
     the two-path cross-check, group identities throughout the reaction layer,
     and figures that lose their H⁺ terms — which for a biogeochemistry course
     is a regression, since proton stoichiometry is precisely what teaches why
     pH moves the energetics.

   Revisit only to interoperate with eQuilibrator values or to do
   metabolic-pathway thermodynamics, where everything upstream is in Alberty's
   convention and mixing conventions silently corrupts a pathway sum.
2. **Supplemental formation-energy table** for species pyGCC lacks, chiefly
   glucose and pyruvate. Every entry carries a provenance string and a
   `verified` flag, and unverified values warn on use. Designed in `SPEC.md` §2.1.
3. **Multi-product couples**, so incomplete oxidations can be catalogued —
   syntrophic propionate and butyrate oxidation, which yield acetate *and* CO₂.
   The `Couple` abstraction currently names one species each side, which is what
   blocks them. This is the main thing missing from the metabolism library.
4. **Provenance export**: per-result record of pyGCC version, database file
   hash, and per-species source, dumpable as BibTeX.
5. **Sulfur disproportionation** and other reactions whose balancing is
   genuinely underdetermined — currently detected and refused, rather than
   solved by asking which products are intended.

### Tier 2 — moderate, mostly new figures over existing machinery

7. **Syntrophy window plot**: both partners' $\Delta G$ against $p\mathrm{H_2}$
   on one axis, showing the narrow overlap where both are exergonic.
8. **Environmental gradient profiles**: read a CSV of depth, T, pH and
   concentrations and plot the affinity of many metabolisms against depth — the
   figure that makes redox zonation fall out of thermodynamics.
9. **Two-dimensional contours** over pairs of variables ($p\mathrm{H_2}$ × pH,
   T × pH) with the $\Delta G = 0$ and energy-quantum contours drawn.
10. **Interactive redox tower** with live pH, temperature and concentration
    sliders, alongside the current static version.
11. **Per-atom oxidation states on the half-reaction figure**, so a species like
    acetate shows its two chemically distinct carbons rather than their mean.
    The per-atom engine already exists; only the figure work remains.
12. **Problem generator and grader**: randomised conditions with worked
    solutions, built on show-your-work.

### Tier 3 — larger, or needing data the current backend lacks

13. **Eh–pH (Pourbaix) diagrams** with water stability lines and the couples
    overlaid.
14. **Wider mineral support**: sulfides beyond pyrite, carbonates, and clays,
    all of which the GWB route already reaches — they need only registry
    entries and validation.
15. **Uncertainty propagation** through formation-energy uncertainties, with a
    tornado plot showing which variable dominates. Note that SUPCRT-lineage
    databases mostly do not carry uncertainties, so this likely needs
    user-supplied values.
16. **PHREEQC export**, emitting a reaction set as PHREEQC input so students can
    move from hand calculation to full speciation modelling.
17. **Pressure beyond near-surface**, opening up the hydrothermal and deep
    subsurface range pyGCC is actually built for.

---

## Files

| path | purpose |
|---|---|
| `SPEC.md` | design specification, including verified backend findings |
| `PROMPT.md` | the original project brief |
| `AGENTS.md` | repository conventions for AI agents |
| `microbial_thermo/` | the library |
| `tests/` | unittest suite |
| `notebooks/` | jupytext-paired notebooks (not yet populated) |
