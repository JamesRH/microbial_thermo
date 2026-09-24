# Where to go for more chemistry

A survey of the libraries and databases that could replace or extend what this
project builds by hand, prompted by *Future work* item 14 — wider mineral
support — and the broader question of whether balancing, chemical lookups, gas
solubility and pH speciation are better solved by someone else's code.

**Two kinds of claim appear below and they are marked differently.** Things I
ran against this environment are marked **[checked]** and the command is
implied by the numbers. Everything else is from reading and is cited; treat it
as a lead, not a result.

---

## The short version

1. **OBIGT, the database inside CHNOSZ, is the single most useful thing here**
   and it is already paying off — item 2's glucose and pyruvate were closed
   with it. SUPCRT-lineage, so the standard state matches with no conversion,
   and every entry carries a citation key. **[checked]**
2. **`speq23.dat` is a free upgrade over `speq21.dat`** — a strict superset,
   +3 species, nothing dropped. **[checked]**
3. **`supcrtbl.dat` ships with pyGCC and loads today**, adding 308 species
   `speq21` lacks including the arsenic minerals. A supplement, not a base.
   **[checked]**
4. **A free solver for chempy already exists on this machine.** The failure is
   packaging, not licensing — one line of wiring fixes it. **[checked]**
5. **Nothing here should replace our balancing or speciation.** The libraries
   are more capable and would make the numbers less defensible.
6. **Latimer and Frost diagrams do not exist in Python.** They are cheap to
   build from what we already compute, and ours could do something textbook
   versions cannot: recompute at the working pH and temperature.

**Decided** (24 September 2026): layer `supcrtbl` and OBIGT on top of
`speq23`, keeping the pyGCC/HKF backend. Borrow algorithms from other
projects, not data — with imported-and-tagged external data left as a separate
future item.

---

## Equation balancing

**Keep what you have.** The sympy nullspace approach here already handles the
hard case — underdetermined systems, exposed through `fix=` — which is more
than the alternatives do.

- **chempy** `balance_stoichiometry` is the obvious library and this project
  already depends on it. It was dropped from the balancing path for a specific
  reason recorded in `NOTES.md`: underdetermined systems need an external CBC
  solver. Nothing found suggests that has changed.
- **scipy.linalg.null_space** computes a nullspace by SVD, in floating point.
  That is the wrong trade here: this library's stoichiometry is exact
  `Fraction` throughout, and normalising to an electron pair has to give exact
  quarters. Sympy's exact nullspace is the right tool and the slower one only
  in a sense that does not matter at this size.

The one thing worth borrowing is framing, not code: **pymatgen** has a
`Reaction` class and, more interestingly for us, **Pourbaix diagram**
generation — which is *Future work item 13*. Worth reading before writing that
from scratch.

## Chemical lookups

- **chemicals** (part of Caleb Bell's ChEDL) carries data for over 20,000
  chemicals as local databanks with on-demand loading, and is designed to be
  depended on by other packages. This is the right place to go for CAS numbers,
  formulae, molar masses, and physical properties — the identity and lookup
  layer `species.yaml` currently hand-maintains.
- **thermo** (same author) adds temperature- and pressure-dependent properties
  and phase equilibria.

These are chemical-engineering oriented: excellent for identity and bulk
properties, **not** a source of aqueous standard-state formation energies on
the SHE convention. They would supplement `species.yaml`, not `speq21.dat`.

## Gas solubility

- **Sander's compilation of Henry's law constants**, version 5.0.0 (2023), is
  the reference work: thousands of species, with temperature dependence given
  as *d ln H / d(1/T)* evaluated at 298.15 K. It is published with electronic
  supplementary data, so it is machine-readable even though no maintained
  Python package wraps it.
- **pyEQL** exposes solution-level properties and gets its speciation from
  PHREEQC.

For this library the gap is narrow. Gases are already handled through the
Maier–Kelley route with partial pressures, which is correct for the
`H2(g)`/`CO2(g)` work the course does. Sander's data would matter if you wanted
**dissolved** concentrations from measured headspace, which is a real teaching
use — students measure headspace and need aqueous activity.

## pH speciation

- **PHREEQC** is the reference implementation, with ion-association, Pitzer and
  SIT models.
- **phreeqpython** wraps the VIPhreeqc extension; **pyEQL** now talks directly
  to the USGS IPHREEQC libraries rather than going through phreeqpython.

**Recommendation: do not adopt these for the course.** The acid–base layer here
computes p*K*a values *from the backend's formation energies at the working
temperature*, so a calculation at 60 °C uses 60 °C p*K*a values. PHREEQC would
bring a second, independent thermodynamic dataset into the same calculation,
and the two would disagree in ways that are invisible in the output — precisely
the failure `tests/test_database_consistency.py` exists to prevent. The current
approach is more consistent, and for a teaching tool that matters more than
breadth.

Where PHREEQC *would* earn its place is **item 16, PHREEQC export** — emitting a
reaction set for students to carry into full speciation modelling. Generating
input is safe; importing results is not.

## More minerals and compounds, including organics — item 14

This is where the real opportunity is, and it is closer than expected.

### Already on disk **[checked]**

pyGCC bundles fifteen databases. The library uses `speq21.dat` (HKF) and
`thermo.com.dat` (GWB log K). Of the rest:

pyGCC has **two** database doors and they are not interchangeable — which is
the source of an error I reported wrongly in the first draft of this file:

- `dbaccess=` takes **SUPCRT** sequential files: `speq21.dat`, `speq23.dat`,
  `supcrtbl.dat`, `berman.dat`.
- `sourcedb=` + `sourceformat=` takes **GWB / EQ3-6 / PHREEQC** files:
  `thermo.com.dat`, `thermo.2021.dat`, `thermo_latest.tdat`, `data0.dat`,
  `phreeqc.dat`. This library has its own reader for the GWB ones
  (`backends/gwb.py`).

> **Correction.** I first reported `thermo.2021.dat` as failing to parse, with
> `ValueError: could not convert string to float: '*'`, and marked it
> **[checked]**. That was my mistake, not a pyGCC bug: I passed a GWB file
> through the SUPCRT reader. Through the right door it loads fine — 1023
> entries. pyGCC 1.5.3 is the current release and there is nothing to report
> upstream. The repository is on **Bitbucket**, not GitHub
> (`bitbucket.org/Tutolo-RTG/pygcc`), which is why searching GitHub found
> nothing.

**SUPCRT side**, through `dbaccess=` **[checked]**:

| file | species | vs `speq21.dat` |
|---|---|---|
| `speq21.dat` *(current)* | 1594 | — |
| `speq23.dat` | 1597 | **+3, nothing dropped** — a strict superset |
| `speq21_dimer.dat` | 1595 | +1 |
| `supcrtbl.dat` | 444 | +308 new, **−1458 missing** — a supplement, not a base |

`speq23.dat` adds epsomite, hexahydrite and kieserite — magnesium sulfate
hydrates. A free upgrade: strictly more, nothing lost.

`supcrtbl.dat` is Holland & Powell mineral data. It cannot be the base (it
lacks 1458 of the species we use) but it adds 308, including the arsenic set.

**GWB side**, through this library's own reader **[checked]**:

| file | entries | arsenic |
|---|---|---|
| `thermo.com.dat` *(current)* | 2937 | 29 |
| `thermo_latest.tdat` | 2904 | 29 |
| `thermo.com.tdat` | 2904 | 29 |
| `thermo.2021.dat` | 1023 | 0 |
| `data0.dat`, `phreeqc.dat`, `thermo_cemdata_mar.tdat` | 0 | — (reader does not understand these formats) |

**`thermo.com.dat` is already the richest GWB file we have.** There is nothing
to gain by switching GWB source — the gain is all on the SUPCRT side.

One small robustness fix falls out: `PygccBackend(database=…)` should reject a
GWB file with a clear message instead of letting it reach `float('*')`.

`supcrtbl.dat` is SUPCRTBL (Zimmer *et al.* 2016), which revised SUPCRT92 and
explicitly **added As-acid, As-metal aqueous species and As-bearing minerals**.
Loading it through `PygccBackend(database=…)` exposes `Ferric-As(am)`,
`Barium-As`, `Barium-H-As` and `As2O5(s)`; the raw file also carries
**arsenopyrite (FeAsS)**, **scorodite (FeAsO₄·2H₂O)** and **arsenolite**. For
the arsenic work just added, amorphous ferric arsenate is the phase that
actually controls arsenic mobility in oxidised sediments — this is not an
academic addition.

At 444 species it is a **supplement, not a replacement**. The natural design is
a second mineral source alongside the GWB route, with the cross-database
agreement test extended to cover it.

### Organics — the part no database solves

- **OBIGT** (CHNOSZ) contains 231 organic molecules and groups from Richard &
  Helgeson (1998), plus **protein group-additivity parameters**. Its
  `organic_aq.csv` is a plain CSV with reference keys resolving against
  `refs.csv`. **[checked]** — this is how glucose and pyruvate were traced.
- **AqOrg** (WORM portal, ASU) estimates thermodynamic properties and **HKF
  parameters for aqueous organic molecules by second-order group additivity**.
  This is the only approach here that extends to compounds *no* database
  contains, which is the actual shape of the organics problem.
- **Reaktoro** ships SUPCRT and SUPCRTBL in YAML with every organic species
  carrying an `organic` tag, so they can be filtered in or out.

### The wider ecosystem

- **Reaktoro** — C++ with Python bindings; reads PHREEQC, SUPCRT, SUPCRTBL,
  NASA and ThermoFun data. The most capable general engine.
- **ThermoFun** — C++/Python, fetches from the ThermoHub database; already
  adopted by GEMS and Reaktoro.
- **pyCHNOSZ** — Python wrapper over CHNOSZ; needs R.
- **AqEquil** — Python interface to EQ3/6, with database building.

**CHNOSZ deserves a paragraph of its own.** It is the closest thing to a prior
art for this entire project: an R package for thermodynamic calculations in
aqueous geochemistry *and geobiochemistry*, in continuous development since
2008, with balanced-reaction generation from basis species, standard molal
properties, activity and predominance diagrams, and protein thermodynamics.
Anyone extending this library should read it — not to copy it, but because it
has already made most of the design decisions once.

---

## Redox diagrams: Latimer, Frost–Ebsworth, Pourbaix

**No Python library draws Latimer or Frost diagrams.** Searching turns up
teaching material and no code. CHNOSZ draws predominance and activity diagrams
in R; pymatgen draws Pourbaix diagrams but from a different kind of data (see
below). So this would be genuinely new in the Python ecosystem — and it is
well within reach, because all three are just presentations of ΔGf values the
library already computes.

### Feasibility, checked

A Frost diagram plots the **volt-equivalent** *N*·*E*°(X^N/X⁰) against
oxidation state *N*. That is −ΔG/F for reducing the species to the element, so
it needs nothing beyond what `delta_Gf` already gives. Computed for arsenic
from this library's own data **[checked]**:

| species | state | volt-eq, pH 0 | volt-eq, pH 7 |
|---|---|---|---|
| HAsO₄²⁻ | +5 | +2.427 | −0.472 |
| H₂AsO₄⁻ | +5 | +2.027 | −1.286 |
| As(OH)₃ | +3 | +0.745 | −0.498 |
| HAsO₂ | +3 | +0.743 | −0.499 |
| As | 0 | 0 | 0 |

Lower is more stable. At pH 0 elemental arsenic sits lowest; at pH 7 both
oxidised forms drop below it. The same numbers recompute at 2 °C and 60 °C
without any extra machinery, because the temperature dependence is in the
formation energies.

**Three things this turned up that would otherwise bite during implementation:**

1. **Elemental forms are not in the registry.** `As`, and the other elements a
   Frost diagram needs, are in the backend's mineral dictionary but have no
   `species.yaml` entry, so `Couple.make("As", …)` fails. Any diagram work
   starts by adding them.
2. **The diagram depends on which species you write**, exactly as E°′ does —
   H₂AsO₄⁻ and HAsO₄²⁻ differ by 0.8 V in volt-equivalent at pH 7. A *correct*
   diagram must use the species that actually dominates at the working pH.
   **This library already has that machinery** in the speciation layer, which
   is precisely what textbook diagrams fixed at pH 0 or 14 cannot do. It is
   the strongest argument for building this rather than borrowing it.
3. **Negative oxidation states need care** in the volt-equivalent convention.
   My quick pass produced a number for AsH₃ I do not trust; that sign
   convention needs deriving properly rather than pattern-matching.

### Pourbaix, and why pymatgen is the wrong borrow

pymatgen has a mature `PourbaixDiagram`, but its data pipeline is
`MPRester.get_pourbaix_entries()` — **Materials Project API key required**, and
the entries are **DFT solid energies plus experimental ion energies**. Pulling
those in would put DFT-derived solids and HKF-derived aqueous species inside
one diagram: the same "two inconsistent datasets in one calculation" problem
that argues against PHREEQC for speciation.

`PourbaixEntry` *can* be constructed by hand, so the honest borrow is the
**geometry**, not the data — take the predominance/convex-hull algorithm, feed
it our own ΔGf values. CHNOSZ does exactly this with its own data. Adding
pymatgen as a dependency to use perhaps a few hundred lines of geometry is a
poor trade; reading its implementation is a good one. (pymatgen is not
currently installed here.)

### Could they be interactive?

**Yes, and the machinery already exists.** The interactive tower (item 10)
established the pattern: precompute a grid over the expensive axes, apply the
cheap ones in closed form, ship it as ipywidgets for notebooks and as
standalone HTML with hand-rolled sliders for everything else. Latimer and
Frost diagrams are cheaper than the tower — a handful of species per element
rather than a dozen couples — so the same approach applies directly. A Pourbaix
diagram is a 2-D predominance field and is the expensive one, but it is also
the one where a precomputed grid is most natural.

## A free solver for chempy — solved **[checked]**

`balance_stoichiometry(..., underdetermined=None)` fails here with
`PulpSolverError: PULP_CBC_CMD: Not Available (check permissions on cbc)`.
That message is misleading. The cause is packaging, not licensing or
permissions:

- chempy asks pulp for `PULP_CBC_CMD`, the CBC binary **vendored inside** the
  pip wheel.
- conda-forge's `pulp` ships **no** vendored binary — zero candidates on disk.
- But a system `cbc` **is** present at `…/envs/microbial-thermo/bin/cbc`, and
  pulp exposes it as `COIN_CMD`, which it lists as the one available solver.

Pointing one name at the other is enough:

```python
import pulp
pulp.PULP_CBC_CMD = pulp.COIN_CMD          # drive the system cbc
balance_stoichiometry({"HS-", "O2"}, {"SO4-2", "S", "H2O", "H+"},
                      underdetermined=None)
# -> 5 HS- + 7 O2  ->  H+ + 2 H2O + 2 S + 3 SO4-2      (verified balanced)
```

CBC is COIN-OR, open source under the EPL. So the answer to "is there a free
solver" is **yes, and it is already installed** — the wiring is the problem.

**But adopting chempy's behaviour would be a regression.** That answer is *a*
valid member of the solution family, chosen by an integer-minimisation
objective — not by chemistry, and returned without saying a choice was made.
This library refuses and asks, which is the whole point of item 5. The right
shape is an **opt-in** mode that returns a minimal-integer member and says
plainly that it picked one.

## Could these libraries shrink our codebase?

Honestly: **not much, and the places they could are mostly places we should
not.** Taking the pieces in turn.

| our code | could a library replace it? | verdict |
|---|---|---|
| balancing (`balance.py`) | chempy + CBC | **No.** Loses exact `Fraction`s and refuse-and-ask. Worth adding as an opt-in mode. |
| speciation (`speciation.py`) | PHREEQC via phreeqpython/pyEQL | **No.** pKa here comes from the backend at the working temperature; importing a second dataset breaks the single-provenance guarantee. |
| species registry (`species.yaml`) | `chemicals` (ChEDL, 20 000 compounds) | **Partly.** Good for identity, CAS, formulae, molar masses. Not a source of SHE-convention formation energies. Would reduce hand-maintenance. |
| formation energies (`backends/`) | Reaktoro, ThermoFun | Possible but a rewrite; see the dataset question. |
| oxidation states (`oxidation.py`) | RDKit *(already used)*, pymatgen | **Worth extending** — see below. |
| figures | CHNOSZ (R), pymatgen | Read them; do not depend on them. |

The uncomfortable conclusion is that the code most tempting to delete is the
code that encodes the project's judgement — exactness, refusal to guess,
single provenance. The libraries are bigger and more capable and would make
the numbers *less* defensible here, not more.

**Where CHNOSZ genuinely should change our approach** is not code but design.
It computes properties from basis species and generates balanced reactions
from them — a cleaner formulation than our couple-centric one, and the reason
it handles predominance diagrams so naturally. Worth reading before item 13.

### Would they give more robust formal redox states?

**For organics, no — we already use the better tool.** RDKit is installed and
`oxidation.py` already partitions bonds by electronegativity over SMILES,
which is the standard approach and gives per-atom states.

**For minerals, yes.** pymatgen has `Composition.oxi_state_guesses()`, which
enumerates plausible assignments using ICSD oxidation-state statistics, and
`BVAnalyzer`, which does bond-valence analysis on structures. Both are more
robust than electronegativity partitioning for inorganic solids — where our
approach is weakest, since it has no structural information. The natural use
is a **cross-check in the test suite** rather than a runtime dependency:
assert that our assignment for arsenopyrite, pyrite, the manganese oxides and
so on agrees with an independent method. That is the same shape as the
existing cross-database consistency test, and it would catch exactly the
class of error nothing currently catches.

---

## What I would do, in order

The dataset direction was **decided on 24 September 2026**: layer on top of
the existing backend rather than replace it. These are now *Future work* items
#20, #21, #22, #24, #25, #26 and #27 in `README.md`.

1. **Compatible datasets** — `speq21` → `speq23` (free), then `supcrtbl.dat`
   as a second mineral source, then OBIGT for organics. Extend
   `test_database_consistency.py` at each step. SUPCRTBL revised mineral
   end-members against Holland & Powell (2011), so disagreements with
   `speq21` are *expected* and must be declared, not silenced.
2. **Fix the chempy solver wiring** and expose it as an opt-in balancing mode
   that says it picked one answer out of many.
3. **Elemental species into the registry** — a prerequisite for any Frost or
   Latimer diagram, and cheap.
4. **The diagram family**: Latimer, then Frost–Ebsworth, then Pourbaix,
   sharing the interactive machinery the tower already established.
5. **Cross-check oxidation states against pymatgen** in the test suite, where
   our electronegativity method is weakest — inorganic solids.
6. **Read CHNOSZ before item 13.** Its basis-species formulation is cleaner
   than our couple-centric one and is why predominance diagrams fall out of it
   naturally.
7. **Leave balancing and speciation alone otherwise.** Both are better as they
   are, for reasons about consistency rather than capability.
8. **Consider AqOrg** if organics become central — the only route that
   generalises to compounds no database contains.

---

## Sources

- [pyEQL](https://github.com/KingsburyLab/pyEQL) · [engines](https://pyeql.readthedocs.io/en/v0.9.0/engines.html) · [phreeqpython](https://github.com/Vitens/phreeqpython)
- [Reaktoro thermodynamic databases](https://reaktoro.org/v1/thermodynamic-databases.html) · [loading databases](https://reaktoro.org/tutorials/basics/loading-databases.html)
- [ThermoFun (JOSS)](https://www.theoj.org/joss-papers/joss.04624/10.21105.joss.04624.pdf)
- [CHNOSZ](https://chnosz.net/) · [OBIGT vignette](https://chnosz.net/vignettes/OBIGT.html) · [customising the database](https://cran.r-project.org/web/packages/CHNOSZ/vignettes/custom_data.html)
- [pyCHNOSZ](https://github.com/worm-portal/pyCHNOSZ) · [AqEquil](https://github.com/worm-portal/AqEquil) · [WORM Portal](https://worm-portal.asu.edu/docs/)
- [SUPCRTBL (Zimmer et al. 2016)](https://www.sciencedirect.com/science/article/abs/pii/S0098300416300371) · [SupPhreeqc](https://www.sciencedirect.com/science/article/abs/pii/S0098300420305501)
- [Sander, Henry's law constants v5.0.0 (2023)](https://acp.copernicus.org/articles/23/10901/2023/acp-23-10901-2023.pdf)
- [chemicals (ChEDL)](https://chemicals.readthedocs.io/) · [thermo (ChEDL)](https://thermo.readthedocs.io/index.html)
- [chempy](https://github.com/bjodah/chempy) · [scipy.linalg.null_space](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.null_space.html) · [pymatgen](https://pymatgen.org/)
- [Amend (2019), demystifying microbial reaction energetics](https://enviromicro-journals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.14778)
