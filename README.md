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

### Figures

```python
from microbial_thermo.figures import (
    plot_half_reactions, plot_redox_tower, plot_energy_explorer,
)

plot_half_reactions(reaction, save="figures/sulfate")     # SVG + PNG
plot_redox_tower(reaction, save="figures/tower")          # SVG + PNG
fig = plot_energy_explorer(reaction, save_html="figures/explorer")
```

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
  couple    Balance one redox couple and report its potentials.
  figure    Render a figure to SVG and PNG, or to standalone HTML for the...
  reaction  Compute a whole reaction from a donor and an acceptor couple.
  species   List known species, optionally filtered by a substring.
  tower     Print the reference redox couples, ordered by potential.
```

Examples:

```bash
mthermo couple HS- SO4-2
mthermo reaction -d "H2(g)" "H+" -a HS- SO4-2 --normalize-to integer --show-work
mthermo reaction -d "H2(g)" "H+" -a methane CO2 -p "H2(g)=1e-5" --ph 6.8
mthermo figure -d acetate CO2 -a H2O O2 --kind tower -o figures/acetate
mthermo figure -d "H2(g)" "H+" -a HS- SO4-2 --kind explorer -o figures/explorer
mthermo tower --ph 7 -T 25
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
python -m unittest discover -s tests     # 116 tests
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
- **pH speciation is not yet implemented.** You currently name the species you
  want (`H2S` or `HS-`) and get exactly that species. The automatic,
  abundance-weighted acid–base treatment is designed in `SPEC.md` §6 but not
  built — see Future work item 1, which is the most important gap.

---

## Future work

Ordered by logical dependency first, then by effort within each tier. Earlier
items unblock later ones.

### Tier 1 — small, and builds directly on what exists

1. **pH speciation layer** *(the most important gap)*. Resolve a typed species
   onto its acid–base family, compute each member's fractional abundance from
   p$K_a$ values derived from the backend at the working temperature, and use
   the abundance-weighted group free energy. Matters most for sulfide, whose
   p$K_a$ of ~7.0 sits exactly at physiological pH, where picking a dominant
   species is a coin flip. Designed in `SPEC.md` §6.
2. **Supplemental formation-energy table** for species pyGCC lacks, chiefly
   glucose and pyruvate. Every entry carries a provenance string and a
   `verified` flag, and unverified values warn on use. Designed in `SPEC.md` §2.1.
3. **Curated reaction library** as YAML: denitrification step by step, DNRA,
   anammox, comammox, sulfur disproportionation, the four methanogenesis
   routes, AOM, acetogenesis, syntrophic propionate and butyrate oxidation,
   photoferrotrophy. Doubles as a much larger regression corpus.
4. **Jupytext-paired teaching notebooks** in `notebooks/`, one per figure type.
5. **Affinity ladder**: sorted bar chart of $\Delta G$ per electron for many
   metabolisms under one measured condition. Nearly free given `sweep`.
6. **Provenance export**: per-result record of pyGCC version, database file
   hash, and per-species source, dumpable as BibTeX.

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
