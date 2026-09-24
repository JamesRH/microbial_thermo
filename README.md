# microbial_thermo

Thermodynamic calculations and figures for microbial physiology and
biogeochemistry. Given a redox couple or a whole reaction, it balances the
chemistry, computes $E^\circ$, $E^{\circ\prime}$ and $\Delta G$ from a real
thermodynamic database, and draws the three figures a course needs: a stacked
half-reaction diagram, a redox tower, and an interactive free-energy explorer.

Built as a teaching aid. Every result can be asked to show its work.

**Code not reviewed, not production ready scientific software**

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
| `02_redox_tower_and_energy` | the tower and its scale bar, the interactive tower, hydrogen and temperature sweeps, the syntrophy window, the interactive explorer |
| `03_environmental_affinity` | pH speciation, the curated library, and the affinity ladder for a real porewater |
| `04_when_the_library_refuses_to_guess` | the three places the library asks instead of guessing: contested names, underdetermined equations, and averaged oxidation states |
| `photoAs`, `photoFe`, `photoNO2` | phototrophy: an uphill CO₂-fixation reaction per donor, then the photons that pay for it |
| `Scratchbook` | a worked problem end to end, starting from nothing but an equation string |

The numbered notebooks are committed with their outputs and are **executed by
the test suite** (`tests/test_notebooks.py`), so a library change that breaks
one is caught rather than discovered in class. They run from the `.py` pair,
never from the committed `.ipynb`, so the check does not depend on stored
outputs.

To regenerate the committed outputs after changing a notebook:

```bash
jupytext --sync --execute notebooks/02_redox_tower_and_energy.py
```

One caveat on notebook 02: the interactive tower is an ipywidgets view, which
needs a live kernel. Its widget state is embedded, so JupyterLab and nbviewer
render it, but **GitHub's notebook viewer shows that cell blank**. The
`plot_interactive_tower(..., save_html=...)` export in the same section is the
route that works with no kernel behind it.

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

**What the string may contain:** species names only, separated by a spaced `+`.
No electrons — they are inferred. No stoichiometric coefficients — the point is
that it is unbalanced. No H⁺ or H₂O unless you want them; they are supplied. A
bare `+` with no space after it is a charge (`Fe+3`, `VO2+`), which is why the
spacing matters.

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
leaves underdetermined.

`balance_equation` on its own has no couples to lean on, so where the nullspace
has more than one dimension it says so and hands back the basis in readable
form rather than picking for you:

```python
balance_equation("HS- + O2 -> SO4-2 + Sulfur(s) + H2O + H+")
# AmbiguousReactionError: 2 independent balanced solutions ... any combination
# of these balances: (0) 4·HS- , 5·O2(aq) , -2·SO4-- , -2·Sulfur(s) , -2·H2O;
# (1) 1·HS- , 2·O2(aq) , -1·SO4-- , -1·H+
```

Supply the missing information with `fix`, one entry per degree of freedom —
they set both the ratio between solutions and the overall scale. Moles are
positive; the sign comes from which side you wrote the species on:

```python
balance_equation(equation, fix={"HS-": 2, "Sulfur(s)": 1})
# HS- -2, O2(aq) -5/2, SO4-- 1, Sulfur(s) 1, H2O 1
```

Different constraints give different chemistry, which is the point — you are
choosing a reaction, not discovering one. Where the equation does *not* resolve to exactly one
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

Pass `per_atom=True` to `plot_half_reactions` and a species with a SMILES shows
each atom's own state instead of the mean — acetate as `-3, +3` rather than `0`,
since its methyl and carboxyl carbons really are four units apart and the mean
describes neither. Repeats collapse: butyrate reads `-3, -2(x2), +3`.

### Which phase a bare name means

`H2`, `O2`, `N2`, `CO2` and `CH4` are each claimed by a dissolved and a gaseous
entry. The **dissolved** form wins, declared in `species.yaml` rather than
falling out of file ordering, because that is the form a cell sees. It matters:
H₂(aq) and H₂(g) are 91 mV apart at pH 7.

```python
mt.resolve("H2")               # H2(aq)
mt.resolve("H2", phase="g")    # H2(g)
mt.default_registry().ambiguities()          # every contested name
mt.default_registry().undeclared_ambiguities # must stay empty; a test enforces it
```

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

### Where the numbers come from

Four sources, consulted in this order. **A later source can only ever add a
species; it never changes a number an earlier one already produced.**

| order | source | route | what it carries |
|---|---|---|---|
| 1 | `speq23.dat` | HKF, direct | 1597 aqueous species, gases, common minerals |
| 2 | `thermo.com.dat` | GWB log K + HKF basis | 2937 entries, chiefly minerals |
| 3 | `supcrtbl.dat` | Holland & Powell (HP11) | 266 species nothing else has |
| 4 | supplemental table | hand-entered | 4 species no database carries |

`speq23.dat` replaced pyGCC's default `speq21.dat` and the upgrade was free:
a strict superset, +3 species, and **not one formation energy of the 85
species this library exposes moved by as much as 1e-9 kJ/mol**. That is
asserted, not assumed. Pinning the file by name also lets the provenance
record say which one was read.

**`supcrtbl.dat` is a supplement, never a base.** SUPCRTBL (Zimmer *et al.*
2016) revised SUPCRT92's mineral end-members against Holland & Powell (2011),
so it is not a newer edition of the same numbers — it disagrees with
`speq23`/GWB by real amounts (hematite 1.7, goethite 2.4, magnetite 3.4
kJ/mol). Consulting it last means those disagreements cannot silently move a
published figure. It also holds only 444 species and lacks 1458 that
`speq23` has, so it could not be a base even if we wanted it.

It needs a different equation of state, and this is the trap worth knowing: a
SUPCRT92 mineral record carries Maier–Kelley coefficients in **calories**, a
SUPCRTBL one carries Holland & Powell parameters in **kJ** — pyrite is
`-38293.0` against `-160.16`. pyGCC evaluates both but has to be told which,
and it converts internally, returning calories either way. Converting again
lands you 4.184× out on a number that still looks perfectly plausible.

What it buys: **scorodite** (FeAsO₄·2H₂O, −1287 kJ/mol), the phase that
controls arsenic solubility in oxidised sediments, along with arsenopyrite
variants, barium arsenates and 260-odd other minerals.

### Species pyGCC does not have

Hydroxylamine, glucose and pyruvate are in none of pyGCC's databases, so they
come from a small hand-entered table instead. The biomass placeholder
`<CH2O>` lives there too, though it is a different kind of thing again — see
below.

**These values are of a different kind from everything else here.** They are
literature figures typed in by hand, not computed from an equation of state and
not cross-checked against a second source the way the mineral data was. Each
carries a `verified` flag, an unverified one **warns every time it is used**, and
any figure resting on one is footnoted in red. **Two of the four are now traced.** Glucose and pyruvate come from the OBIGT
database of CHNOSZ, and through it from Amend & Plyasunov (2001) and Canovas &
Shock (2016) respectively. OBIGT is SUPCRT-lineage, so the standard state
matches the rest of this library. Tracing them also gained each an enthalpy,
which is what makes a van 't Hoff correction possible away from 25 °C. Glucose
moved 4.7 kJ/mol in the process (−917.2 → −912.5), which is a useful reminder
of what an untraced "commonly cited" number is worth.

Hydroxylamine is still unverified, and its provenance now records where it was
looked for and not found, so nobody repeats the search. `Biomass(aq)` is
unverified permanently and by design.

**`Biomass(aq)` is not a compound at all.** `<CH2O>` is a stand-in for cell
carbon at oxidation state zero, and its value is this library's own glucose
divided by six, so the two are at least mutually consistent. It gets two
things wrong on purpose: real biomass is nearer CH₁.₈O₀.₅N₀.₂ — slightly more
reduced, so a little costlier to make — and it contains nitrogen this term
ignores entirely, so the cost of assimilating N is missing. Treat any
autotrophy yield built on it as a statement about sign and scale. The example
script below takes a `--biomass-dgf` override precisely so you can check
whether a conclusion survives the choice.

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

### Vanadium

Aqueous ions from the HKF database, oxides and the sulfate ion pair through the
log K route — no hand-entered values. `V+2`, `V+3`, `VO+2` (vanadyl),
`VO2+` (dioxovanadium), `VO4-3` (vanadate), `V2O4`, `V3O5`, `VOSO4(aq)`.

Two naming traps:

- **`VO2` is the solid**, vanadium(IV) oxide, stored as `V2O4` — the doubled
  formula of the same compound. **`VO2+` is not that**: it is the aqueous
  dioxovanadium(V) cation, one oxidation state higher. One character apart,
  V(IV) against V(V).
- **`VOSO4(aq)` is an ion pair**, not a redox form, and it is tabulated at 25 °C
  only. Pairing it against `VO2+` fails conservation because it carries sulfur;
  use `VO+2` for redox work. Likewise `V2O4` against `VO+2` is a dissolution,
  not a couple — both are V(IV) — and the library says so rather than inventing
  a potential.

The two reduction couples come out at +0.173 V (V(V)/V(IV)) and −0.486 V
(V(IV)/V(III)) at pH 7, so the simplest vanadate reduction is just

```python
mt.Reaction.from_equation("H2 + VO2+ -> VO+2")
# 2 H+ + H2(aq) + 2 VO2+ -> 2 VO++ + 2 H2O
```

To keep VOSO₄ explicitly, `from_equation` cannot infer it — vanadium's state is
unreadable in that formula — so name the couple yourself, carrying sulfate along
on the same side so sulfur conserves:

```python
from microbial_thermo.reaction import Couple, Reaction

Reaction.from_couples(
    donor=Couple.make("H2(g)", "H+"),
    acceptor=Couple.make(["VOSO4(aq)"], ["VO2+", "SO4-2"], key_element="V"),
    conditions=conditions, normalize_to="integer",
)
# 2 H+ + H2(g) + 2 VO2+ + 2 SO4-- -> 2 VOSO4(aq) + 2 H2O,  ΔG°′ = -141.6 kJ/mol
```

VOSO₄ then goes unannotated on the half-reaction diagram, since its oxidation
state cannot be assigned, but everything else is labelled normally.

### Arsenic

`As(III)` and `As(V)`, both from pyGCC — the arsenate ladder as direct HKF
data, As(III) by the log K route. No hand-entered values.

```python
mt.Reaction.from_couples(donor=("H2(g)", "H+"), acceptor=("As(III)", "As(V)"))
mt.Reaction.from_couples(donor=("As(III)", "As(V)"), acceptor=("H2O", "O2(aq)"))
```

Available: `As(OH)3(aq)`, `H2AsO3-`, `HAsO2(aq)`, `AsO2-` for As(III);
`H3AsO4(aq)`, `H2AsO4-`, `HAsO4--`, `AsO4---` for As(V). A bare `As(III)` means
`As(OH)3(aq)` and `As(V)` means `HAsO4--`.

Arsenic is worth having because **the couple is run in both directions**.
As(V)/As(III) sits at **+0.013 V** at pH 7 — almost exactly in the middle of
the tower, above the sulfur and carbon donors and below the nitrogen and oxygen
acceptors. So arsenate pays as a respiratory acceptor against hydrogen
(−82 kJ/mol per electron pair) *and* arsenite pays as a lithotrophic donor
against oxygen (−163 kJ/mol), and both organisms exist. Five metabolisms are
catalogued, three reducing and two oxidising.

Two things to know before quoting a number:

- **As(III) has two representations** in SUPCRT-lineage data, differing by one
  water: `As(OH)3` / `H2AsO3-` and `HAsO2` / `AsO2-`. Both are registered and
  they agree to 0.21 kJ/mol over 0–100 °C. That agreement is worth something,
  because the two arrive by *different routes* — `As(OH)3(aq)` through the GWB
  log K path, `HAsO2(aq)` as direct HKF. `HAsO2(aq)` is the one species where
  the two databases disagree at all (0.26 kJ/mol), which is why the bare name
  resolves to `As(OH)3(aq)`.
- **E°′ depends on which arsenate you write.** The three arsenate forms differ
  in proton count, so at pH 7 the couple reads +0.160 V against `H3AsO4`,
  +0.020 V against `H2AsO4-`, and +0.013 V against `HAsO4--`. Published E°′
  values for As(V)/As(III) disagree largely for this reason. The standard-state
  form is the unambiguous anchor: H₃AsO₄ + 2 H⁺ + 2 e⁻ → H₃AsO₃ + H₂O comes out
  at **+0.574 V** against a published +0.560 V, agreeing to 14 mV.

Arsenate's second p*K*a is **6.76**, right on physiological pH, so at pH 7 it is
a 36:64 mixture of `H2AsO4-` and `HAsO4--` and neither form alone is "arsenate".
Use the family for anything quantitative:

```python
mt.Conditions(pH=7.0, total_concentrations={"arsenate": 1e-6})
```

Arsenite, by contrast, is 99.4% the neutral `As(OH)3` at pH 7 — its p*K*a is
9.2 — which is the reason As(III) is the mobile, membrane-permeant and more
toxic state. It only ionises in genuinely alkaline water, which is where the
Mono Lake organisms live.

### Phototrophy: paying for an uphill reaction with light

A phototroph's defining trick is running a reaction that does not pay.
Photoferrotrophy, photoarsenotrophy and nitrite-driven photoautotrophy all fix
CO₂ with a donor far too weak to do it in the dark.

```python
from microbial_thermo.photons import photon_energy, photons_required
from microbial_thermo.figures import plot_light_budget
from microbial_thermo.reaction import Couple, Reaction

fixation = Reaction.from_couples(
    donor=Couple.make("As(III)", "As(V)"),
    acceptor=Couple.make("Biomass", "CO2"),   # (reduced, oxidized), as always
    conditions=conditions,
)
fixation.delta_G_standard_prime          # +80.3 kJ/mol — uphill
photon_energy("P870")                    # 137.5 kJ/mol per mole of photons
photons_required(fixation.delta_G_standard_prime)   # 0.73
plot_light_budget(fixation)
```

The reaction direction comes from which couple is the **acceptor**, not from
the order within a couple. `acceptor=Couple.make("Biomass", "CO2")` reduces
CO₂, so the reaction is written in the CO₂-fixation direction. Writing the
couple backwards still gives the right energy — it only swaps the
`reduced`/`oxidized` labels and the oxidation-state annotations.

At 870 nm (bacteriochlorophyll *a*, the pigment of purple bacteria), the three
donors rank:

| donor | ΔG°′ (kJ/mol per 2 e⁻) | photons needed |
|---|---|---|
| Fe(II) → ferrihydrite | +51.2 | 0.52 |
| As(III) → As(V) | +80.3 | 0.73 |
| NO₂⁻ → NO₃⁻ | +156.4 | 1.28 |

So one photon per electron pair covers iron and arsenite but **not** nitrite.
Two things worth knowing before quoting any of that:

- **`photon_energy` is an upper bound.** It returns $N_A hc/\lambda$; radiation
  carries entropy and a reaction centre captures only part of the excitation.
  Pass `efficiency=` to ask a less generous question.
- **`photons_required` is a floor, not a quantum requirement.** Real anoxygenic
  phototrophs run cyclic electron flow and reverse electron transport and spend
  several photons per electron. This answers only "how many photons' worth of
  energy is the reaction short by".

Photoferrotrophy against ferrihydrite is the interesting one: it crosses zero
near **pH 9 in the dark**, so above that the light is buying rate rather than
thermodynamics.

The `notebooks/photoAs`, `photoFe` and `photoNO2` notebooks work each donor
through the same steps.

### An arbitrary answer, when you want one

`balance_equation` refuses an underdetermined equation by default, because
conservation genuinely does not pick an answer. `choose="minimal"` overrides
that: it solves an integer program for the smallest whole-number coefficients
and **warns that it chose**.

```python
balance_equation("HS- + O2 -> SO4-2 + Sulfur(s) + H2O + H+", choose="minimal")
# UserWarning: ...has 2 independent balanced solutions and choose='minimal'
#              picked one of them by minimising the coefficients. That is an
#              arithmetic preference, not a chemical one...
# 5 HS- + 7 O2 -> 3 SO4-2 + 2 S + 2 H2O + H+
```

Use it to explore, never for a number you intend to quote — `fix=` is how you
say which reaction you actually mean, and it gives a different, equally valid
answer.

This is the same integer program `chempy.balance_stoichiometry(...,
underdetermined=None)` solves, through the same free CBC solver, and the two
agree on this equation. It is solved here rather than delegated because
chempy parses its own formula strings and this library's backend names are not
formulas — `Acetate`, `Methane(aq)`, `SO4--` — so the round trip would fail or
silently re-parse a species into something else.

**Solver note.** CBC is COIN-OR, open source. `pulp` from PyPI vendors a
binary; conda-forge's does not but installs `cbc` on PATH, which is why
chempy's own underdetermined mode fails on a conda install with the misleading
message *"check permissions on cbc"*. This library asks for whichever is
available, so it works either way.

### Reference state and ionic strength

```python
mt.Conditions.biochemical()                    # pH 7, 25 °C, I = 0.25 M, b-dot
mt.Conditions.biochemical(ionic_strength=0.7)  # seawater
reaction.activity_correction()                 # what the coefficients are worth
```

`Conditions.biochemical()` is the biochemical reference state: pH 7, 25 °C and
**I = 0.25 M**, the ionic strength the transformed biochemical tables are
quoted at and close to physiological.

**It will not reconcile a number here with an Alberty-convention table**, and
it is worth being precise about why. An ionic strength scales activity
*coefficients*. ΔG°′ is defined at unit activity, so there is nothing for a
coefficient to multiply and the shift is not small but exactly **zero** —
measured, not argued. What remains is the *convention* difference, and no
choice of ionic strength touches it.

Where ionic strength does bite is a reaction with stated concentrations.
`activity_correction()` reports exactly how much:

| ionic strength | ΔG (kJ/mol) | activity correction |
|---|---|---|
| 0 (ideal) | −10.97 | 0 |
| 0.10 | −10.49 | +0.48 |
| 0.25 | −10.32 | +0.65 |
| 0.70 (seawater) | −10.12 | +0.85 |

Under a kJ/mol from fresh water to seawater for this reaction — a real
correction, but not the reason two sources disagree by ten.

### Provenance

```python
record = reaction.provenance()
record.to_text()      # human-readable block to paste under a figure
record.to_dict()      # JSON-ready
record.to_bibtex()    # a citable @misc entry
record.unverified     # species resting on hand-entered values
```

Records the four things a computed number rests on: software versions, the
database files **with their hashes**, a per-species source line, and the
standard-state conventions applied. The hash matters more than it looks —
pyGCC ships several SUPCRT and GWB files, they are revised between releases,
and nothing in a computed value otherwise records which one was read.

```
Species
  As(OH)3(aq) <- thermo.com.dat (log K route)
  Biomass(aq) <- hand-entered supplemental table (A PLACEHOLDER…)  [UNVERIFIED]
  CO2(aq)     <- pyGCC speq21.dat (pyGCC default)
  H2O         <- IAPWS-95 via pygcc.iapws95
```

### The curated metabolism library

Thirty-seven named metabolisms, so you need not remember which couples to pair.
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

#### The syntrophy window

Two partners on one axis, with the overlap shaded where each is exergonic.

```python
from microbial_thermo.figures import plot_syntrophy_window
from microbial_thermo.library import reaction

producer = reaction("syntrophic_propionate_oxidation", sediment)
consumer = reaction("hydrogenotrophic_methanogenesis", sediment)

figure, window = plot_syntrophy_window(producer, consumer)
print(window.low, window.high, window.decades)   # 2.3e-06 2.1e-04 1.95
print(window.best_shared)                        # (2.2e-05, -5.5) bar, kJ/mol
```

The window edges are interpolated in log pressure off whichever curve is the
binding constraint at that edge, so the reported width is a property of the
chemistry rather than of how finely the axis was sampled — a 13-point and a
25-point grid agree to six decimal places.

`best_shared` gives the pressure at which the *worse-off* partner does best,
and what it gets there. For propionate against hydrogenotrophic methanogenesis
at pH 7 that is −5.5 kJ/mol: at any hydrogen pressure, one partner or the other
is doing at least that badly, which is well inside the energy-quantum band.
That is the figure's point.

Both partners must be normalised to the same electron count; sharing a y axis
between a per-6-electron and a per-8-electron free energy would look fine and
mean nothing, so the function refuses it.

#### The interactive tower

The static tower draws one set of conditions, which invites the reading that
the tower is a fixed table. It is not — couples move at different rates and
change places.

```python
from microbial_thermo.tower import tower_grid
from microbial_thermo.figures import interactive_tower, plot_interactive_tower

grid = tower_grid()                          # a few seconds per temperature
interactive_tower(grid)                      # ipywidgets, live in a notebook
plot_interactive_tower(grid, save_html="tower")   # standalone HTML, no kernel
```

Three sliders — pH, temperature, and the oxidised:reduced activity ratio — over
a precomputed grid. Nothing calls the backend while a slider moves. The grid
samples pH finely and temperature coarsely, because that asymmetry is real:
formation energies are cached per temperature, so a new pH costs about 10 ms
and a new temperature about 2.5 s. The activity ratio is not stored at all; it
is a closed-form Nernst shift, $(RT/nF)\ln(\text{ratio})$, applied on lookup, so
it stays continuous rather than quantised to grid steps.

The crossing worth showing students: at pH 7 and 25 °C oxygen sits above iron,
and raising the ferric:ferrous ratio past about $10^{1.5}$ puts iron on top,
because a one-electron couple moves twice as fast per decade as a two-electron
one.

Plotly's native sliders cannot express three *independent* dimensions — a
slider's steps have no access to the other sliders' positions — so the exported
page builds its own sliders and redoes the Nernst shift in JavaScript. Both
implementations agree to five decimal places, which is itself a cross-check.

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
python -m unittest discover -s tests     # 489 tests, ~180 s
MT_SKIP_NOTEBOOKS=1 python -m unittest discover -s tests   # skip the slow notebook runs
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

**Numbering is stable, and now actually is.** Completed items keep their
original number and move to the table above rather than being renumbered,
because `NOTES.md` and commit messages refer to them by number.

These were ordered-list items until 24 September 2026, which quietly broke that
promise: Markdown renumbers an ordered list from its first item regardless of
the digits written, so a file saying 24, 25, 26, 27 rendered as 18, 19, 20, 21
and nobody reading GitHub saw the numbers this repo cites. They are plain
bullets with literal `#N` labels now, which render the same everywhere.

**Next up:** item **20, compatible datasets** — high priority, direction
decided, and every other dataset item waits on it. `RESEARCH.md` carries the
measurements behind that decision and the survey of what else is out there.

### Completed

| # | Item | Where it lives |
|---|---|---|
| 3 | Syntrophy window figure | `figures/syntrophy.py`, `tests/test_syntrophy_figure.py`, notebook 02 |
| 5 | Underdetermined full-equation balancing | `balance_equation(..., fix=)`, `tests/test_ambiguity.py`, notebook 04 |
| 6 | Execute the notebooks in the test suite | `tests/test_notebooks.py` |
| 7 | Phase-ambiguous species names | `canonical:` in `species.yaml`, `registry.ambiguities()`, notebook 04 |
| 10 | Interactive redox tower | `tower.tower_grid`, `figures/interactive_tower.py`, notebook 02 |
| 11 | Per-atom oxidation states on the half-reaction figure | `oxidation.format_per_atom_states`, `per_atom=True`, notebook 04 |
| 1 | Ionic-strength reference state | `Conditions.biochemical()`, `Reaction.activity_correction()`, `tests/test_reference_state.py` |
| 4 | Provenance export | `microbial_thermo/provenance.py`, `Reaction.provenance()`, `tests/test_provenance.py` |
| 18 | `resolve()` covering every route | `backends/pygcc_backend.py`, `tests/test_provenance.py` |
| 19 | Reproducible HTML exports | `EXPLORER_DIV_ID` in `figures/explorer.py`, `tests/test_provenance.py` |
| 2 (part) | Glucose and pyruvate traced and verified | `data/supplemental_gibbs.yaml`, `tests/test_provenance.py` |
| 21 | Minimal-integer balancing as an opt-in | `balance_equation(..., choose="minimal")`, `tests/test_ambiguity.py` |

### Tier 1 — small, and builds directly on what exists

> **A standing note on the transformed (Legendre) standard state.** It was
> evaluated and deliberately declined; the reasoning is kept here so it does
> not get relitigated.
>
> - There are two conventions, not one. The *microbial bioenergetics*
>   literature — Thauer, Jungermann & Decker; Amend & Shock; LaRowe & Amend —
>   writes reactions with explicit species and explicit protons. **That is
>   already what this library computes**, and our −152.2 kJ/mol for
>   hydrogenotrophic sulfate reduction sits on the ≈−152 kJ/mol tabulated
>   there. The *biochemical* convention (Alberty, IUBMB, eQuilibrator) is the
>   one we do not match.
> - The two agree more than expected. For a reaction whose every reactant is a
>   single species, hydrogen conservation forces the transformed result to
>   equal our species-level ΔG°′ *exactly*; the transform is pure rebookkeeping
>   there. They diverge only through pseudoisomer grouping, bounded by the
>   mixing entropy $RT\ln(\text{populated forms})$ — at most ≈1.7 kJ/mol, for a
>   diprotic reactant sitting exactly on its p$K_a$.
> - **Ionic strength does not close the gap.** This was the original plan for
>   item 1, and measuring it killed the idea: an ionic strength scales activity
>   *coefficients*, and ΔG°′ is defined at unit activity, so the shift there is
>   not small but exactly **zero**. `Conditions.biochemical()` now exists and is
>   the right reference state for concentration-dependent work, but it does not
>   reconcile anything with an Alberty table.
> - The cost is not small: a second balancing path (dropping both the hydrogen
>   and charge rows from the conservation matrix), a transformed analogue of
>   the two-path cross-check, group identities throughout the reaction layer,
>   and figures that lose their H⁺ terms — which for a biogeochemistry course
>   is a regression, since proton stoichiometry is precisely what teaches why
>   pH moves the energetics.
>
> Revisit only to interoperate with eQuilibrator values or to do
> metabolic-pathway thermodynamics, where everything upstream is in Alberty's
> convention and mixing conventions silently corrupts a pathway sum.

- **#2 · Verify the supplemental values** — *partly done*. Glucose and pyruvate
   are now traced and verified (see the Completed table). **Hydroxylamine is
   still outstanding**: it is in none of the fifteen databases bundled with
   pyGCC and not in CHNOSZ's OBIGT either, in either its inorganic or organic
   tables. The likely primary source is the NBS tables — Wagman *et al.*
   (1982), *J. Phys. Chem. Ref. Data* **11**, Suppl. 2 — which NIST has since
   digitised; that should settle it. `Biomass(aq)` is on the unverified list
   **permanently** and correctly: ⟨CH₂O⟩ is a modelling placeholder, not a
   compound.

- **#20 · Compatible datasets** — **HIGH PRIORITY**, and the direction is decided
    (24 September 2026; the measurements behind it are in `RESEARCH.md`).
    **Layer, do not replace.** Keep the pyGCC/HKF backend, in three
    independently testable steps:

    1. `speq21.dat` → **`speq23.dat`**. A strict superset: 1594 → 1597 species,
       nothing dropped, adding epsomite, hexahydrite and kieserite. Free.
    2. Add **`supcrtbl.dat`** as a second mineral source — 444 species, 308 of
       them not in `speq21`, including arsenopyrite, scorodite, arsenolite and
       amorphous ferric arsenate. It cannot be a base (it lacks 1458 species we
       use) but it is the phase data the arsenic work actually needs. SUPCRTBL
       revised mineral end-members against Holland & Powell (2011), so
       disagreements with `speq21` are **expected** — declare them in
       `KNOWN_DISAGREEMENTS` with magnitudes, do not silence them.
    3. Add **OBIGT** (CHNOSZ) as the organics source. SUPCRT-lineage, so the
       standard state already matches, and every entry carries a citation.

    Nothing to gain on the GWB side: `thermo.com.dat` is already the richest
    file pyGCC ships. Small robustness fix that falls out of this work:
    `PygccBackend(database=…)` should reject a GWB file with a clear message
    rather than failing inside pyGCC's float parser.
- **#22 · External data sources, imported and provenance-tagged** *(moderate)*.
    Extend the mechanism behind the supplemental table from hand-entered
    one-offs to whole external sources, so a species we lack can come from
    OBIGT or elsewhere carrying its origin and a verification flag. This is
    the deliberate second half of the "borrow algorithms, not data" rule —
    allowed, but only where we have nothing, and never silently.

### Tier 2 — moderate, mostly new figures over existing machinery

- **#8 · Environmental gradient profiles**: read a CSV of depth, T, pH and
   concentrations and plot the affinity of many metabolisms against depth — the
   figure that makes redox zonation fall out of thermodynamics.
- **#9 · Two-dimensional contours** over pairs of variables ($p\mathrm{H_2}$ × pH,
   T × pH) with the $\Delta G = 0$ and energy-quantum contours drawn.
- **#12 · Problem generator and grader**: randomised conditions with worked
    solutions, built on show-your-work.

### Tier 3 — larger, or needing data the current backend lacks

- **#13 · Eh–pH (Pourbaix) diagrams** with water stability lines and the couples
    overlaid. **Borrow the geometry, not the data**: pymatgen has a mature
    `PourbaixDiagram`, but its pipeline wants a Materials Project API key and
    builds entries from DFT solid energies plus experimental ion energies —
    mixing those with our HKF aqueous species would put two provenances inside
    one diagram. `PourbaixEntry` can be built by hand, so read the algorithm
    and feed it our own ΔGf. CHNOSZ does exactly this in R.
- **#14 · Wider mineral support**: sulfides beyond pyrite, carbonates, and clays,
    all of which the GWB route already reaches — they need only registry
    entries and validation. **See `RESEARCH.md`** — the cheapest first step is
    `supcrtbl.dat`, which pyGCC already ships and this library can already
    load, and which carries arsenic minerals (arsenopyrite, scorodite,
    amorphous ferric arsenate) that `speq21.dat` lacks.
- **#15 · Uncertainty propagation** through formation-energy uncertainties, with a
    tornado plot showing which variable dominates. Note that SUPCRT-lineage
    databases mostly do not carry uncertainties, so this likely needs
    user-supplied values.
- **#16 · PHREEQC export**, emitting a reaction set as PHREEQC input so students can
    move from hand calculation to full speciation modelling.
- **#17 · Pressure beyond near-surface**, opening up the hydrothermal and deep
    subsurface range pyGCC is actually built for.
- **#24 · Latimer diagrams** for every element with more than two redox forms in
    the registry — the condensed chain of couples with E°′ on each arrow.
    **Nothing in Python draws these.** Cheap: it is the couples we already
    compute, laid out in oxidation-state order. Two prerequisites, both found
    by checking rather than assumed: elemental forms (`As`, `S`, `Fe`…) are in
    the backend but **not in the registry**, so `Couple.make("As", …)` fails
    today; and the diagram must use the species that actually dominates at the
    working pH, which the speciation layer can already decide.
- **#25 · Frost–Ebsworth diagrams** — volt-equivalent *N*·*E*° against oxidation
    state, where the slope between two points is the couple potential and
    convexity shows disproportionation. Also absent from Python. Feasibility
    is **confirmed**: computed for arsenic from our own data it gives +2.427 V
    for HAsO₄²⁻ at pH 0 against −0.472 V at pH 7, and recomputes cleanly at
    2 °C and 60 °C. One caution: the volt-equivalent convention for **negative**
    oxidation states needs deriving properly — a quick pass produced an AsH₃
    number that does not look right.
- **#26 · Interactive Latimer, Frost and Pourbaix** on the pattern item 10
    established: precompute the expensive axes on a grid, apply the cheap ones
    in closed form, ship ipywidgets for notebooks and standalone HTML with
    hand-built sliders for everyone else. Latimer and Frost are *cheaper* than
    the tower; Pourbaix is the expensive one and the best candidate for a
    precomputed grid. **This is the point of building them rather than
    borrowing**: a diagram that recomputes at the working pH and temperature
    is something no textbook version can do.
- **#27 · Cross-check oxidation states against an independent method** *(small)*.
    `oxidation.py` partitions bonds by electronegativity over SMILES, which is
    the right approach for organics and is where RDKit already earns its keep.
    For **inorganic solids** it is weakest, because it has no structural
    information. pymatgen's `Composition.oxi_state_guesses()` (ICSD statistics)
    and `BVAnalyzer` (bond valence) are more robust there. Use them as a
    **test-suite cross-check**, not a runtime dependency — assert that our
    assignment for arsenopyrite, pyrite and the manganese oxides agrees with an
    independent method. Same shape as the existing cross-database audit, and it
    would catch a class of error nothing currently catches.

---

## Examples

Runnable scripts that go beyond a single call. Each is library-first — import
the functions from a notebook, or run the file for a CLI.

### Arsenite-driven carbon fixation

`examples/arsenite_carbon_fixation.py` builds the **anabolic** half of
chemolithoautotrophy — As(III) → As(V) driving CO₂ into biomass — and plots it
against pH beside the catabolism that pays for it.

```bash
python examples/arsenite_carbon_fixation.py --save figures/as_fixation
python examples/arsenite_carbon_fixation.py --biomass-dgf -130 --no-catabolic
```

```python
from arsenite_carbon_fixation import fixation_vs_ph, plot_fixation_vs_ph

result = fixation_vs_ph(temperature_c=25.0)
figure, _ = plot_fixation_vs_ph(result=result)
```

The reaction is **endergonic** — about **+80 kJ/mol** per electron pair at
pH 7 — and that is the point: fixing carbon costs energy, and an autotroph has
to earn it back from catabolism. Three things the figure makes visible:

- The cost falls by **11.4 kJ/mol per pH unit**, which is exactly
  2 × RT·ln10. Two protons leave per electron pair, so the slope is
  predictable before you plot it — the script prints the observed slope beside
  the prediction as a self-check.
- The catabolism moves at the *same* rate and the same direction, so the net
  improves twice as fast, about 23 kJ/mol per pH unit. Alkaline water is
  thermodynamically kinder to these organisms, which is one reason the best
  studied arsenite oxidisers come out of soda lakes.
- At pH 4 a one-to-one budget nets only −14 kJ/mol — **inside** the biological
  energy quantum. The margin is real, not a rounding detail.

`--biomass-dgf` re-runs with a different placeholder energy. Across any
plausible value the sign does not change, which is the honest way to use a
placeholder: show that the conclusion does not rest on it.

```
$ python examples/arsenite_carbon_fixation.py --help
Usage: arsenite_carbon_fixation.py [OPTIONS]

  Plot arsenite-driven CO2 fixation against pH.

Options:
  --ph-low FLOAT                Lowest pH.  [default: 4.0]
  --ph-high FLOAT               Highest pH.  [default: 10.0]
  --points INTEGER              Grid points.  [default: 25]
  --temperature FLOAT           Temperature in °C.  [default: 25.0]
  --biomass-dgf FLOAT           Override ΔGf of the ⟨CH2O⟩ placeholder,
                                kJ/mol, for a sensitivity run.
  --catabolic / --no-catabolic  Also draw the arsenite/O2 catabolism that pays
                                for the fixation.  [default: catabolic]
  --save TEXT                   Write SVG and PNG to this path stem.
  --show / --no-show            Open an interactive window.
  --version                     Show the library and dependency versions, then
                                exit.
  --help                        Show this message and exit.
```

## Files

| path | purpose |
|---|---|
| `NOTES.md` | operational knowledge for anyone picking this up: pyGCC's quirks, the conventions that must not drift, and the traps |
| `RESEARCH.md` | survey of other libraries and databases — balancing, lookups, gas solubility, speciation, and where to get more minerals and organics |
| `SPEC.md` | design specification, including verified backend findings |
| `PROMPT.md` | the original project brief |
| `AGENTS.md` | repository conventions for AI agents |
| `microbial_thermo/` | the library |
| `tests/` | unittest suite |
| `notebooks/` | jupytext-paired teaching notebooks, committed with outputs and executed by the suite |
| `examples/` | runnable analysis scripts, library-first with a CLI |
