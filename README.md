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
source setup.sh
```

That one command builds the `microbial-thermo` environment if it is missing,
activates it, installs this package in editable mode, and registers the Jupyter
kernel. Pass `--force` to rebuild from scratch.

It must be **sourced**, not executed — activating an environment changes the
current shell, and a subshell would discard that. Running it directly prints a
reminder and exits rather than half-working.

Then launch the notebooks:

```bash
./run_notebooks.sh              # JupyterLab, notebooks/ directory, right kernel
./run_notebooks.sh --classic    # the classic Notebook interface
```

`run_notebooks.sh` is executed rather than sourced: it activates the environment
for its own process and leaves your shell alone.

<details>
<summary>Doing it by hand instead</summary>

```bash
mamba env create -f environment.yaml
mamba activate microbial-thermo
pip install -e .
```

Or piece by piece:

```bash
mamba create -n microbial-thermo python=3.11
mamba activate microbial-thermo
mamba install -c conda-forge numpy pandas scipy chempy sympy pint matplotlib seaborn plotly rdkit click jupytext pyyaml ipykernel ipywidgets ruff
pip install pygcc==1.5.3
```

```bash
mamba activate microbial-thermo
python -m ipykernel install --user --name microbial-thermo --display-name "microbial-thermo"
```

</details>

`pygcc` is installed with `pip` because it is not published on conda-forge;
PyPI is its only distribution channel. Everything else comes from conda-forge.

### Jupyter

Pick **microbial-thermo** as the notebook kernel. Notebooks in
`notebooks/` are paired to `.py` scripts with jupytext in `percent` format;
edit the `.py` and run `jupytext --sync <notebook>.ipynb`.

| notebook | covers |
|---|---|
| `01_half_reactions` | couples, balancing, oxidation states, show-your-work, the half-reaction figure |
| `02_redox_tower_and_energy` | the tower and its scale bar, hydrogen and temperature sweeps, the interactive explorer |
| `03_environmental_affinity` | pH speciation, the curated library, and the affinity ladder for a real porewater |
| `Scratchbook` | a worked problem end to end, starting from nothing but an equation string |

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

### From a written equation

Type the reaction and let the library work out the rest:

```python
mt.Reaction.from_equation("NO3- + H2 -> NH2OH", conditions=conditions)
# H+ + 3 H2(aq) + NO3- -> NH2OH(aq) + 2 H2O
```

It finds which elements change oxidation state, forms the donor and acceptor
couples from that, and balances water, protons and electrons as usual. Species
the writer left implicit are supplied — nobody writes the proton that H₂
oxidises to.

Hydrogen and oxygen are considered last, since they double as the auxiliary
balancing species. In the example above the hydrogens of hydroxylamine are not
an oxidation product of H₂; they are just hydrogens, and preferring nitrogen
keeps that straight.

Disproportionation works too — `"S2O3-2 + H2O -> SO4-2 + HS-"` resolves to
thiosulfate as both donor and acceptor, which is a reaction conservation alone
leaves underdetermined. Where the equation does *not* resolve to exactly one
donor and one acceptor it raises rather than guessing; name the couples with
`from_couples` in that case.

### One couple, one line

```python
mt.half_reaction("HS-", "SO4-2").E_standard          # +0.2491 V
mt.half_reaction("HS-", "SO4-2").E_standard_prime    # -0.2168 V at pH 7
```

| property | meaning |
|---|---|
| `.E_standard` | E° — every activity 1, **H⁺ included** (pH 0) |
| `.E_standard_prime` | E°′ — every activity 1, H⁺ at the conditions' pH |
| `.E` | the actual potential at your concentrations |

`E_standard` ignores the `pH` field, since it holds the proton at 1 M too; only
temperature reaches it. Pass conditions, a `key_element` hint, or `n_electrons`
to rescale the displayed equation — the last changes how the half reaction
reads, not the potential, since E is intensive.

Either side may name several species, as couples do generally:
`mt.half_reaction("Propanoate(aq)", ["Acetate", "HCO3-"])`.

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

Step 8 works the potentials out in full rather than asserting them. For each
half reaction it tabulates the Nernst correction species by species —

$$E = E^\circ - \frac{RT}{nF}\sum_i \nu_i \ln a_i$$

— then shows $E^\circ$ becoming $E^{\circ\prime}$ by the proton term alone,
and $E$ by every term, before taking $\Delta E$ and comparing $-nF\Delta E$
against the formation-energy sum. The tests parse the rendered numbers back out
and check they add up, so the narrative cannot drift from the arithmetic.

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

### Species pyGCC does not have

Hydroxylamine, glucose and pyruvate are in none of pyGCC's databases, so they
come from a small hand-entered table instead.

**These values are of a different kind from everything else here.** They are
literature figures typed in by hand, not computed from an equation of state and
not cross-checked against a second source the way the mineral data was. Each
carries a `verified` flag, an unverified one **warns every time it is used**, and
any figure resting on one is footnoted in red. All three currently ship unverified — trace them to a primary source and
set the flag before relying on them.

They are single-temperature values, honoured at 25 °C and refused elsewhere
unless an enthalpy is available, in which case a van 't Hoff correction is
applied and a warning issued — the same stance taken for the isothermal
manganese oxides.

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

Thirty-one named metabolisms, so you need not remember which couples to pair.
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

### Multi-product couples and syntrophy

Either side of a couple may name several species, which is what lets
*incomplete* oxidations be expressed:

```python
mt.Reaction.from_couples(
    donor=("Propanoate(aq)", ["Acetate", "HCO3-"]),   # two products
    acceptor=("H2(g)", "H+"),
    conditions=conditions,
    normalize_to="integer",
)
# Propanoate(aq) + 3 H2O -> Acetate + HCO3- + H+ + 3 H2(g)
```

Give explicit proportions where they are not one-to-one: `[("Acetate", 2)]` for
butyrate. Proportions *within* a side are yours to state, because conservation
cannot supply them — propionate could in principle go to 1.5 acetate, or to
3 bicarbonate, and which happens is biochemistry. Conservation then fixes the
scale *between* the sides, and any element that fails to balance is caught.

This is what makes syntrophy expressible. Both reactions come out **endergonic
at standard state** — +78 kJ/mol for propionate against a published +76, +50 for
butyrate against +48 — which is exactly why they need a partner to draw the
hydrogen down. The test suite asserts that a syntrophic window exists: a range
of $p_{H_2}$ where propionate oxidation and hydrogenotrophic methanogenesis are
*both* exergonic, with neither working outside it.

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
by the constant $nF$, so at an electron pair the bar marks 0.104 V as the
−20 kJ/mol energy quantum and 0.259 V as one ATP. You can read a reaction's
yield straight off the tower by comparing its donor–acceptor gap to the bar. An
absolute kJ/mol *axis* would be wrong — it would imply each couple has an
absolute free energy, when only differences carry energy.

The bar carries **both conventions**, because they answer different questions
and only one of them depends on n:

| | at n = 2 | at n = 6 |
|---|---|---|
| per reaction — 1 ATP | 0.259 V | 0.086 V |
| per electron — 1 ATP per e⁻ | 0.518 V | 0.518 V |

The potential axis is intensive: every rung stays exactly where it is under
renormalisation. The reference quantities (one ATP, the energy quantum) are
per-reaction, so *their* position does move with n. Giving both readings makes
that visible instead of leaving it to a caption — and the gap between the two
sets is n itself.

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
python -m unittest discover -s tests     # 243 tests
ruff format microbial_thermo tests
ruff check microbial_thermo tests
```

---

## Known limitations

- **Glucose, pyruvate and hydroxylamine are absent** from every database pyGCC
  ships. Acetate, formate, lactate, propanoate, butanoate, methanol and ethanol
  are present. The three missing ones come from the hand-entered supplemental
  table instead, and all three are flagged `verified: false` — they warn on every
  use and footnote any figure that depends on them. Tracing them to primary
  sources is the one real correctness debt outstanding.
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
2. **Verify the supplemental values.** Hydroxylamine, glucose and pyruvate are
   hand-entered and flagged unverified. Each needs tracing to a primary source,
   confirming against its standard state, and the flag setting.
3. **Syntrophy window figure**: both partners' ΔG against $p_{H_2}$ on one
   axis, shading the overlap where each is exergonic. The energetics are in
   place and asserted in the tests; only the plotting remains.
4. **Provenance export**: per-result record of pyGCC version, database file
   hash, and per-species source, dumpable as BibTeX.
5. **Underdetermined full-equation balancing.** `Reaction.from_equation` now
   handles disproportionation, because inferring the couples supplies the
   information conservation cannot — thiosulfate resolves to donor and acceptor
   both being S₂O₃²⁻. But `balance_equation` alone still refuses such cases. It
   could offer the nullspace basis and ask which combination is meant rather
   than only raising.
6. **Execute the notebooks in the test suite.** They are currently executed by
   hand, so a library change can silently rot them. `nbconvert --execute` over
   `notebooks/` would catch it; the cost is roughly a minute.
7. **Phase-ambiguous species names.** A bare `H2` resolves to `H2(aq)` purely
   because that entry comes first in `species.yaml`, and the gas is 91 mV away
   at pH 7. The resolution is deliberate but undeclared — either make it
   explicit in the registry or warn when an ambiguous name is used.

### Tier 2 — moderate, mostly new figures over existing machinery

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
| `NOTES.md` | operational knowledge for anyone picking this up: pyGCC's quirks, the conventions that must not drift, and the traps |
| `SPEC.md` | design specification, including verified backend findings |
| `PROMPT.md` | the original project brief |
| `AGENTS.md` | repository conventions for AI agents |
| `microbial_thermo/` | the library |
| `tests/` | unittest suite |
| `notebooks/` | jupytext-paired teaching notebooks, committed with outputs |
