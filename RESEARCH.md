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
   and it is already paying off — items 2's glucose and pyruvate were closed
   with it. It is SUPCRT-lineage, so its standard state matches this library
   with no conversion, and every entry carries a citation key. **[checked]**
2. **`supcrtbl.dat` ships with pyGCC already and this library can load it**, and
   it carries arsenic minerals `speq21.dat` lacks. That is the cheapest possible
   progress on item 14. **[checked]**
3. **AqOrg** estimates HKF parameters for aqueous organics by group additivity,
   which is the only route that generalises to compounds no database contains.
4. **Do not replace the balancing.** It is not the weak part.

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

| file | loads? | species | note |
|---|---|---|---|
| `supcrtbl.dat` | **yes** | 444 | adds arsenic minerals; warns on 3 duplicate entries |
| `speq23.dat` | yes | 1597 | newer SUPCRT revision |
| `thermo.2021.dat` | **no** | — | `ValueError: could not convert string to float: '*'` |

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

## What I would do, in order

1. **Add `supcrtbl.dat` as a supplementary source.** Smallest change, real
   payoff, and it lands on the arsenic work already done. Extend
   `test_database_consistency.py` to cover the overlap — SUPCRTBL revised
   mineral end-members against Holland & Powell (2011), so disagreements with
   `speq21.dat` are *expected* and should be declared, not silenced.
2. **Use OBIGT as the standing source for anything pyGCC lacks.** It is
   convention-matched and citation-carrying. This is now recorded in `NOTES.md`.
3. **Read CHNOSZ before item 13 (Pourbaix) or item 15 (uncertainty).** Also
   look at pymatgen's Pourbaix implementation.
4. **Consider AqOrg** if organics become central — it is the only route that
   generalises.
5. **Leave balancing and speciation alone.** Both are better as they are, for
   reasons that are about consistency rather than capability.

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
