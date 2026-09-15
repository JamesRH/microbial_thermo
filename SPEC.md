# Software Specification — `microbial_thermo`

A Python library for microbial physiological thermodynamic calculations, built for
teaching first and research calculations second.

**Status:** draft v0.1 · **Branch:** `ai/claude/thermo-core`

---

## 1. Purpose and scope

Given a possibly-unbalanced redox reaction, optionally with a pH, temperature, and
non-unit concentrations, the library:

1. Balances the reaction for atoms, charge, and electrons.
2. Splits it into an oxidative half reaction (electron donor) and a reductive half
   reaction (electron acceptor).
3. Computes $E^\circ$, $E^{\circ\prime}$, $\Delta G^\circ$, $\Delta G^{\circ\prime}$, and
   in-situ $\Delta G$ under user conditions.
4. Assigns formal oxidation states to the atoms that change state.
5. Renders three figures: a stacked half-reaction diagram, a redox tower, and an
   interactive free-energy explorer.

### In scope for v1.0

| | |
|---|---|
| Backend | pyGCC only (equilibrator-api explicitly rejected) |
| Temperature | 0.01–100 °C |
| Pressure | 1 bar / PSAT (near-surface) |
| Phases | Aqueous species and gases; minerals read-only where the database has them |
| Organics | Low coverage acceptable; formate, acetate, lactate, methanol, ethanol, propionate, butyrate from pyGCC, plus a small supplemental table |

### Explicitly deferred

Everything in the *Future work* section of the README.

---

## 2. Backend findings that constrain the design

These were verified against pyGCC 1.5.3 and drive several decisions below.

1. **No arbitrary-reaction API.** `calcRxnlogK` only returns $\log K$ for a species'
   *database-defined* formation reaction. Arbitrary redox stoichiometry must be
   assembled by us as $\Delta G_\mathrm{rxn} = \sum_i \nu_i \Delta G_{f,i}$ using
   `species_eos.supcrtaq` per species. **This is our core numerical path.**
2. **Legacy SUPCRT/GWB naming.** `Fe++`, `SO4--`, `S2--`, `Methane(aq)`, `Sulfur`,
   `Acetate`, `Lactic_acid(aq)`. No `CH4(aq)`, no `Fe+2`. An alias layer is mandatory.
3. **No free-electron species.** pyGCC handles redox through an O₂(aq)/H⁺ basis.
   We compute $E$ from $\Delta G = -nFE$ ourselves.
4. **Activity coefficients are available.** `water_dielec(T, P, Dielec_method=...)`
   exposes `Ah`, `Bh`, `bdot` as functions of T and P. We use these for an extended
   Debye–Hückel (B-dot) per-species activity coefficient rather than writing the
   solvent model from scratch.
5. **Coverage gaps.** Glucose and pyruvate are absent from every bundled database.
   N₂O(aq) is absent. The default direct-access database (`speq21.dat`, 1594 species)
   *does* contain ionized organic forms — `Acetate`, `Formate(aq)`, `Lactate(aq)`,
   `Propanoate(aq)`, `Butanoate(aq)`, `Sulfur(s)` — so most fermentation intermediates
   are available after all.
6. **Water is not a database species.** `H2O` is absent from `dbaccessdic` because HKF
   parameterizes solutes, not the solvent. Water's Gibbs energy comes from
   `pygcc.iapws95(T=TC, P=P).G` (T in **Celsius**), which returns −56678 cal/mol at
   25 °C — within 9 cal/mol of the SUPCRT convention value and internally consistent
   with the rest of pyGCC. Since nearly every half reaction contains water, this is a
   required special case in the backend, not an optional one.
7. **Pressure has a cliff at 100 °C.** At exactly 100 °C, P = 1 bar is below the
   saturation pressure, so water would be steam and `supcrtaq` returns `NaN`. Policy:
   use P = 1 bar for T < 99 °C and pyGCC's saturation option (`P='T'`) at or above
   99 °C. The difference below 99 °C is <0.5 cal/mol, so the switch is numerically
   seamless. Passing a literal 1.013 also fails on the boundary — only the `'T'`
   sentinel is reliable there.
8. **T = 0 °C crashes** the IAPWS-95 water EOS. Floor is 0.01 °C.
9. **Speed: ~28 ms per species per call**, and `supcrtaq` returns `ndarray` of shape
   (1,), not a float. Passing a **vector** of T is *slower* than looping (507 ms for 10
   points vs 284 ms), so there is no vectorization win — use scalar calls behind a
   cache, and precompute grids for interactive widgets.

### 2.1 Consequence: the supplemental data table

Missing species are supplied by `data/supplemental_gibbs.yaml`, a small table of
literature $\Delta G_f^\circ$ values. Every entry carries a `provenance` string and a
`verified` flag. **Unverified entries raise a `UserWarning` at use** and are marked on
any figure that depends on them. The library must never silently substitute a value of
unknown origin — for a teaching tool that is worse than failing.

Supplemental entries are 25 °C only and do not carry HKF parameters, so they cannot be
extrapolated in T. Requesting a supplemented species at T ≠ 25 °C raises unless the
caller passes `allow_isothermal_extrapolation=True`, which applies a van 't Hoff
correction if $\Delta H_f$ is known and otherwise refuses.

---

## 3. Architecture

```
microbial_thermo/
  __init__.py            configure(), get_backend(), __version__
  units.py               pint registry, F, R, physical constants
  formula.py             formula string -> element counts + charge
  oxidation.py           per-atom oxidation states, NOSC
  species.py             Species dataclass, alias resolution, supplemental data
  speciation.py          acid-base families, pH-weighted transform, pKa from backend
  balance.py             full-reaction balancing, half-reaction splitting, normalization
  reaction.py            Reaction / HalfReaction; dG, E, Q; the two-path cross-check
  energetics.py          per-electron normalization, ATP yield, biological energy quantum
  showwork.py            LaTeX derivation emitter
  typeset.py             equation -> positioned token model (renderer-agnostic)
  backends/
    base.py              ThermoBackend ABC
    pygcc_backend.py     pyGCC implementation + caching
  figures/
    halfreaction.py      stacked aligned half-reaction diagram (matplotlib)
    tower.py             redox tower, static (matplotlib) + interactive (plotly)
    explorer.py          gapminder-style dG explorer (plotly)
    style.py             shared palette, ATP/BEQ axis helper
  data/
    aliases.yaml         friendly name <-> pyGCC species string
    acid_base.yaml       acid-base families and members
    supplemental_gibbs.yaml
    reactions.yaml       curated microbial metabolisms
  cli.py                 click entry points
```

Dependency rule: `figures/` may import from the core; the core must never import
`figures/` or plotting libraries. This keeps the numerical library importable in a
headless context and keeps figure work testable in isolation.

---

## 4. Core numerical model

### 4.1 The single source of truth

```
dGf(species, T, P)  ->  pint Quantity [kJ/mol]
```

Everything else derives from this. Backend returns cal/mol; the boundary converts once.

### 4.2 Reaction free energy

$$\Delta G^\circ(T,P) = \sum_i \nu_i \,\Delta G_{f,i}^\circ(T,P)$$

$$\Delta G = \Delta G^\circ + RT \ln Q, \qquad Q = \prod_i a_i^{\nu_i}$$

Activities $a_i = \gamma_i m_i$, with $\gamma_i$ from B-dot using pyGCC's `Ah`, `Bh`,
`bdot`. Water activity is 1 in dilute solution unless the caller overrides. Gases use
fugacity, computed from the gas/aqueous $\Delta G_f$ difference.

### 4.3 Half-reaction potentials

Written as reductions per IUPAC convention, then

$$E = \frac{-\Delta G_\mathrm{half}}{nF}$$

- $E^\circ$: all species at unit activity, including $[\mathrm{H}^+] = 1$ M.
- $E^{\circ\prime}$: same but $[\mathrm{H}^+] = 10^{-7}$, i.e. pH 7 (or the caller's pH).

Overall reaction: $\Delta E = E_\mathrm{acceptor} - E_\mathrm{donor}$ and
$\Delta G = -nF\Delta E$.

### 4.4 The two-path cross-check

Every `Reaction` computes $\Delta G$ **twice** — once by summing $\Delta G_f$ over the
full balanced reaction, once as $-nF\Delta E$ from the two half reactions — and asserts
agreement within a tolerance (default 0.01 kJ/mol). Disagreement raises
`ThermodynamicConsistencyError`. This catches sign-convention errors, electron
miscounts, and bad half-reaction splits automatically. It is the single most valuable
correctness feature in the library and runs on every calculation, not just in tests.

---

## 5. Balancing

`chempy.balance_stoichiometry` over an element + charge conservation matrix, exact
rational arithmetic via sympy. Free electrons are a first-class species (`e-`, charge −1).

**Underdetermined systems:** `chempy` defaults to returning a *symbolic* parametric
solution. We pass `underdetermined=None` to get smallest positive integers, but we
first check the nullspace dimension. If > 1 the reaction is genuinely ambiguous
(disproportionation being the common case) — we raise `AmbiguousReactionError` carrying
the basis, rather than silently picking one. This satisfies the zero-assumption policy.

**Half-reaction convention:** balance O with H₂O, H with H⁺, charge with e⁻.

### 5.1 Normalization

Default: **2 electrons**. Overridable via `normalize_to=`:

| Value | Meaning |
|---|---|
| `"electron_pair"` | default, n = 2 |
| `"electron"` | n = 1 |
| `"donor"` | 1 mol electron donor |
| `"acceptor"` | 1 mol electron acceptor |
| `"H2"` or any formula | 1 mol of that component |

Coefficients become `Fraction` after normalization; the typesetter renders them as
fractions, not decimals.

---

## 6. pH speciation (answers user question 1)

The caller never defines an acid-base group by hand. Writing `H2S` or `sulfide` both
resolve to the sulfide family. At the requested pH and T:

1. Look up the family in `acid_base.yaml`.
2. Get each $pK_a$ at T from the backend (verified path: $\Delta G_f$ difference).
3. Compute fractional abundance $f_i$ of each member.
4. Default `speciation="weighted"`: use the abundance-weighted group free energy —
   the correct treatment when pH ≈ p$K_a$.
5. `speciation="dominant"`: pick the majority species; a warning fires when the pH is
   within 1 unit of a p$K_a$, where this approximation is worst.
6. `speciation="explicit"`: use exactly the species the user typed, no transform.

The figure prints the dominant form with its percentage, and footnotes the reaction when
pH is within 1 unit of a p$K_a$. Sulfide at pH 7 is the motivating case — it is nearly
50/50 H₂S/HS⁻, so "pick the dominant one" is a coin flip there.

---

## 7. Oxidation states

RDKit's `CalcOxidationNumbers` fails conservation on neutral ethanol (verified — sums to
6, should be 0), so we implement our own:

- Per-atom, from a SMILES via RDKit's bond graph. Each bond's electrons are assigned to
  the more electronegative partner (Pauling scale); homonuclear bonds split evenly.
- **Validation invariant:** per-atom states must sum to the molecular charge. Enforced
  in code, not just in tests.
- Formula-only fallback when no SMILES is available: the mean state of the element,
  which is what the figure displays anyway when a species has several atoms of the
  element.
- **NOSC** for organics (LaRowe & Van Cappellen 2011):
  $$\mathrm{NOSC} = 4 - \frac{4a + b - 3c - 2d + 5e - 2f - Z}{a}$$
  for $\mathrm{C}_a\mathrm{H}_b\mathrm{N}_c\mathrm{O}_d\mathrm{P}_e\mathrm{S}_f^{\,Z}$.
  Test fixtures: CO₂ → +4, CH₄ → −4, acetate → 0, glucose → 0.

---

## 8. Figures

### 8.1 Half-reaction diagram (matplotlib → SVG + PNG)

Two equations stacked and aligned at their arrows, oxidative on top (donor),
reductive below (acceptor), with a connecting arrow annotated with the number of
electrons transferred.

- Far left: `Oxidation` / `Reduction` labels.
- Above the top equation and below the bottom equation: formal oxidation state of the
  element that changes state, centered on its species (mean if several atoms).
- Far right: $E^{\circ\prime}$ and $E^\circ$ per half reaction.

**Implementation:** `typeset.py` builds a renderer-agnostic token model — a list of
`(text, role, species_ref, x, y, width)` — measured with matplotlib's Agg renderer
(always Agg, even when exporting SVG, since the SVG backend's text metrics can differ).
Arrow alignment is achieved by measuring both equations' left-hand sides and padding the
shorter. Because layout is a data structure rather than hand-tuned coordinates, it is
unit-testable without rendering.

Mathtext is sufficient for chemical sub/superscripts; `usetex` is not required and would
add a system LaTeX dependency.

### 8.2 Redox tower

Static (matplotlib, SVG + PNG) and interactive (plotly). Both half reactions marked on
an $E^{\circ\prime}$ axis by default. Title reports overall $\Delta G^{\circ\prime}$.

Interactive sliders: pH, temperature, concentrations, gas partial pressures.

**ATP / energy-quantum scale (answers user question 2):** on the tower the y-axis is $E$,
where the ΔG mapping depends on $n$. So the ATP scale is rendered in the legend area as
an annotated ruler labeled for this reaction's own $n$, not as a twin axis. On the
explorer, where the y-axis *is* ΔG, it becomes a true twin right-hand axis, plus a shaded
band for the ~−20 kJ/mol biological energy quantum.

### 8.3 Free-energy explorer (plotly)

ΔG on y; x-axis variable selected from a dropdown (pH, temperature, H₂, any constituent),
gapminder style, via plotly `updatemenus`. SVG export via the modebar — **client-side, so
no kaleido and therefore no Chrome dependency**.

**Responsiveness:** pyGCC at ~30–80 ms/species/call cannot run per slider frame. We
precompute a grid over the swept variables and interpolate with scipy. The grid is baked
into the traces before export, so `fig.write_html()` produces a self-contained file with
working widgets and no live kernel — a file you can hand to students.

---

## 9. Units

`pint` at the API boundary: public function signatures, returned quantities, and figure
axis labels. Internal numeric hot paths (scipy interpolation, optimization) strip to bare
floats, since scipy does not accept `Quantity` objects. One shared registry in `units.py`
— multiple registries cause cryptic failures.

---

## 10. Show-your-work mode

`reaction.show_work()` emits the full derivation as LaTeX suitable for a notebook cell:
balancing, half-reaction split, each $\Delta G_f$ with its source, the $Q$ term by term,
each activity correction, and the Nernst step. Requested to follow immediately after core
functionality lands.

---

## 11. CLI

`click`, with logic decoupled from argument parsing so everything is importable into a
notebook. Every command supports `--version`, which prints the library version plus
pygcc, numpy, pandas, scipy, matplotlib, and plotly versions for traceability.

---

## 12. Provenance

Every computed result carries a `provenance` record: pyGCC version, database file name
and hash, per-species data source, and whether any supplemental or unverified value was
used. Dumpable as text or BibTeX.

---

## 13. Testing

Standard library `unittest`. Ground-truth fixtures:

- NOSC values for CO₂, CH₄, acetate, glucose (LaRowe & Van Cappellen table).
- p$K_a$ of H₂S at 25 °C = 6.99 (verified reproducible from pyGCC to 4 s.f.).
- Balancing: a set of hand-balanced microbial metabolisms.
- The two-path cross-check acts as a property test across the whole curated
  reaction library.
- Oxidation-state conservation: per-atom sum equals molecular charge, over a set of
  organics.

Textbook $E^{\circ\prime}$ values are used as *sanity ranges*, not exact assertions —
published tables differ in standard state and database vintage, and asserting equality
against them would encode someone else's conventions as truth.

---

## 14. Open questions

1. Does "edt" in the original prompt mean EDTA? Metal–EDTA complexes come from a
   different database lineage and would be a meaningful lift.
2. Which glucose $\Delta G_f^\circ$ source do you want to standardize on, given pyGCC
   has none?
