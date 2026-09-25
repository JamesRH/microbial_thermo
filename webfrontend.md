# Software Specification — `microbial_thermo` web frontend

A single-page Streamlit application exposing three teaching tools over the
existing `microbial_thermo` library, with no new thermodynamics and (mostly)
no new backend code — this is a UI layer over calls the library already makes.

**Status:** design only, nothing implemented. **Depends on:** `microbial_thermo`
0.1.0, as documented in the repo's [README.md](README.md).

Decisions below were confirmed with James on 2026-09-24; see §8 for the
questions that produced them and the reasoning kept for the record.

---

## 1. Purpose and scope

Three tools, one page, one shared conditions panel:

| # | Tool | Wraps |
|---|---|---|
| 1 | Compound $\Delta G_f^\circ$ lookup, with name suggestions | `mt.resolve`, `registry.suggest`, `backend.delta_Gf` |
| 2 | Half-reaction $E^\circ$ / $E^{\circ\prime}$ lookup | `mt.half_reaction` |
| 3 | Reaction balancer to an electron pair, two half-reaction strings, each with its potential, and a reveal button for $\Delta G^{\circ\prime}$ | `mt.Reaction.from_equation`, `Reaction.from_couples` |

Out of scope for v1: figures (tower, ladder, Pourbaix, Latimer, Frost), the
curated metabolism library, sweeps, provenance export, and multi-user
accounts. All of that stays notebook/CLI-only; §7 lists it as future work if
wanted.

---

## 2. Framework choice

**Streamlit**, single `app.py`, no separate frontend build step.

### Why, over the alternatives considered

- **Gradio** — comparable simplicity, but its component model is tuned for
  single-function ML demos (inputs → one output block). This app is three
  linked tools sharing one conditions panel and needs ordinary widget layout
  (columns, expanders, a sidebar) more than a fixed input/output contract.
  Streamlit's layout primitives fit that better.
- **Flask/FastAPI + HTML** — full control, but means writing and maintaining
  HTML/CSS/JS and a WSGI/ASGI deployment (gunicorn/uvicorn, a reverse proxy,
  static assets) for a three-tool internal tool. That is real ongoing surface
  area with no payoff here: nothing about this app needs a custom REST API,
  a JS framework, or non-form interactions.
- **Panel** (HoloViz) — a reasonable second choice, and the natural one if
  this app later needs to embed the existing Plotly/matplotlib figures
  interactively server-side. Noted in §7 as the thing to switch to if figures
  come into scope, not a reason to start there now.

### Why Streamlit fits this codebase specifically

- The library is already pure Python with no web code anywhere in it —
  Streamlit adds zero new languages or build tooling, matching AGENTS.md's
  "library-first" and KISS directives.
- `microbial_thermo` is a `pip install -e .` package already; the Streamlit
  app imports it exactly like a notebook does (`import microbial_thermo as mt`).
- Deployment path matches the "local now, hostable later" requirement
  directly (§4): `streamlit run app.py` locally, the same repo pushed to
  Streamlit Community Cloud for a shared URL later, no code changes between
  the two.
- Streamlit reruns top-to-bottom on every interaction, which is the wrong
  model for something CPU-heavy done from scratch each time — the caching
  strategy in §5 exists specifically to make that not matter here, since
  every backend call this app makes is a pure function of its arguments.

---

## 3. Dependency and environment impact

Add to `environment.yaml` and `pyproject.toml`'s `[project.optional-dependencies]`
as a new `web` extra — **not** a core dependency, since nothing else in the
library needs it:

```yaml
# environment.yaml, appended to the pip: block is unnecessary — streamlit is
# on conda-forge:
  - streamlit>=1.38
```

```toml
[project.optional-dependencies]
web = ["streamlit>=1.38"]
```

Per AGENTS.md §4, this is a `mamba install` addition, not a `pip install` —
unlike pyGCC, Streamlit is on conda-forge. No other new dependency: the app
imports `microbial_thermo`, `pandas` (already a dependency, used for the
species picker's searchable table), and Streamlit itself.

`mthermo --version`-style traceability is preserved by putting a small
"About / versions" expander in the app's sidebar that calls the library's
existing `mt.versions()` and adds `streamlit` to the printed dict — no new
version-tracking mechanism invented.

---

## 4. Deployment

**Local (day one):**

```bash
source setup.sh          # unchanged
mamba install -c conda-forge streamlit   # or: pip install -e ".[web]"
streamlit run app.py
```

Opens `http://localhost:8501` in a browser. This is the whole deployment
story for single-user, offline use, and it should stay the primary documented
path — it is the one that needs no accounts, no network exposure, and no
ongoing hosting cost.

**Shared, later, without rearchitecting:** two options, ordered by effort:

1. **Streamlit Community Cloud** (free tier). Point it at the GitHub repo and
   `app.py`; it builds from `environment.yaml`/`pyproject.toml` and gives a
   public `*.streamlit.app` URL. Zero server administration, but the app is
   then reachable by anyone with the link — fine for a course tool with no
   sensitive data (nothing here touches student records; it's chemistry), but
   worth a conscious decision, not a default. **Cold-start note:** pyGCC's
   database load happens on first backend use per process, so the first
   request after an idle period (free-tier apps sleep) will be several
   seconds slow before the resource cache in §5 takes over — acceptable for a
   teaching tool, worth a one-line note in the UI ("first load may take a
   moment").
2. **A department/university-run server**: `streamlit run app.py --server.port
   XXXX --server.address 0.0.0.0` behind the institution's existing reverse
   proxy (Apache/Nginx) with TLS terminated there, run under `systemd` or
   inside a container for restart-on-crash. This is the option if student
   data privacy policy rules out an external host; it costs a short runbook,
   not a rewrite.

Nothing in the app code should differ between these three; the only knobs are
how the process is started and what sits in front of it.

---

## 5. Architecture and caching

One process, one Streamlit script, three tabs (`st.tabs`) sharing one
sidebar. No database, no session persistence beyond Streamlit's own
per-browser-tab session state (used only for UI state like "which
suggestion did they click" and "is the reveal open").

### The backend is a shared, cached singleton

`mt.get_backend()` lazily builds a `PygccBackend` and caches gibbs/solvent
values in plain instance dicts (`microbial_thermo/backends/pygcc_backend.py`).
Streamlit serves concurrent users as separate threads inside **one process**
by default, so:

- Wrap backend construction in `st.cache_resource` so every session reuses
  the same `PygccBackend` instance (and its internal caches) instead of
  reloading pyGCC's database files per session:

  ```python
  @st.cache_resource
  def backend():
      return mt.get_backend()
  ```

- The backend's internal caches are plain `dict`s with no lock. Under
  CPython's GIL a `dict.get`/`dict[key] = value` pair cannot corrupt memory,
  but two threads racing on the same uncached `(species, T)` key can both
  miss the cache and both pay the compute cost once — a correctness non-issue
  (both threads compute the same pure value) and a minor, rare performance
  cost, not worth adding a lock for at this scale (a classroom, not a
  production multi-tenant service). Documented here so nobody "fixes" it
  into a deadlock later.

### Cache every pure lookup at the Streamlit layer too

Every value this app displays is a deterministic function of its inputs
(species names, phase, temperature, pH, concentrations, ...). Wrap each in
`st.cache_data` keyed on those arguments, e.g.:

```python
@st.cache_data(show_spinner=False)
def lookup_dgf(name: str, temperature_c: float, pressure_bar: float | None):
    return backend().delta_Gf(name, temperature_c, pressure_bar)

@st.cache_data(show_spinner=False)
def lookup_half_reaction(reduced: str, oxidized: str, conditions_key: tuple):
    ...
```

This matters concretely: the README records a new pH costing ~10 ms and a new
temperature ~2.5 s against this backend (from the tower's precompute note).
Without caching, moving a temperature slider naively would re-trigger a
multi-second recompute on every rerun; with `st.cache_data`, only genuinely
new argument combinations pay that cost, and slider drag-backs are instant.

`Conditions` objects are not hashable by default in a way `st.cache_data`
can key on cleanly (it's a mutable dataclass with dict fields) — pass the
*scalar arguments* into the cached function and construct `Conditions` inside
it, rather than passing a `Conditions` instance in. This also keeps the cache
key legible for debugging.

### Error surfacing, not exception surfacing

The library's exception hierarchy (`microbial_thermo/exceptions.py`) is
already designed for this — messages are written to explain, not just to
name a broken invariant. The app should catch by type and render, never show
a raw traceback:

| Exception | UI treatment |
|---|---|
| `SpeciesNotFoundError` | Red inline message with the exact wording, plus the `.suggestions` list rendered as clickable buttons that re-run the lookup with that name |
| `AmbiguousReactionError` | The message plus `.basis` rendered as a table of the independent solutions (mirrors the README's example table); a short explanation that the equation needs a `fix=`-equivalent constraint, with the couple-picker fallback offered (see §6.3) |
| `OutOfRangeError` | "This temperature/pressure is outside pyGCC's supported range (0.01–100 °C, near-surface pressure)" |
| `ThermodynamicConsistencyError` | Should never occur for a working reaction — if it does, show it verbatim and flag it as a bug to report, since the README treats this as an internal cross-check failure, not a user-input problem |
| `MissingDataError`, `BalancingError` | Show the message verbatim; these already read as prose |

---

## 6. Tool specifications

### 6.1 Compound $\Delta G_f^\circ$ lookup

**Inputs:** a text field for the species name, plus temperature (°C) and
optional pressure (bar) — *not* the full conditions panel, because formation
energy is a standard-state property of one species at a given $T$/$P$; pH,
concentrations and ionic strength have no effect on it and would be
misleading to show as if they did.

**Behavior:**

1. On every keystroke (debounced) or on submit, call
   `mt.default_registry().resolve(name)`.
2. On success: show the resolved `Species` — its canonical backend name,
   formula, phase, and display label — then `backend.delta_Gf(species.backend,
   temperature_c, pressure_bar)` in kJ/mol (and the raw calorie value in an
   expander, since that's the backend's native unit per `NOTES.md`).
3. If the species is one of the hand-entered, unverified table entries
   (hydroxylamine, `Biomass(aq)`; glucose and pyruvate are now verified per
   the README), show the same red "UNVERIFIED — hand-entered value" footnote
   convention the figures already use, sourced from `SpeciesRecord.verified`
   / `Provenance`. This is the one place the app *must* reproduce an existing
   library convention rather than inventing its own — the whole point of that
   flag is that it follows the number everywhere it's shown.
4. On `SpeciesNotFoundError`: render the message and suggestion buttons
   (§5). Clicking a suggestion re-runs the lookup with that name.
5. A small results table below the input lists the last N (e.g. 10) lookups
   this session, for comparing several compounds side by side — session
   state only, cleared on refresh, not persisted.

**Ambiguous names** (`H2`, `O2`, `N2`, `CO2`, `CH4`): if `resolve(name)`
succeeds but the name is in `registry.ambiguities()`, show a small phase
toggle (aqueous/gas) defaulting to the dissolved form per the library's own
convention, with a one-line note why (documented in the README: "that is the
form a cell sees").

### 6.2 Half-reaction $E^\circ$ / $E^{\circ\prime}$

**Inputs:** two species text fields (reduced, oxidized) — either side may
accept a comma-separated list for a multi-species couple, matching
`Couple.make`'s support for e.g. `Couple.make("Propanoate(aq)", ["Acetate",
"HCO3-"])`. Below that, the shared **Conditions panel** (§6.4), since this
tool's whole point is showing how $E$ moves off $E^\circ$ and $E^{\circ\prime}$.

**Output**, once both species resolve:

- The balanced half-reaction string as text (the same rendering
  `HalfReactionResult.__str__`/the CLI's `couple` command already produces —
  reuse it, don't reformat it).
- $n$ electrons and the oxidation-state change (`mt.mean_oxidation_state`) for
  the key element.
- Three potentials, always all three, clearly labeled and never conflated —
  this directly reflects the questionnaire answer that chose full `Conditions`
  exposure, and the README's own table is the right way to present it:

  | property | shown as |
  |---|---|
  | `.E_standard` | **E°** — pH 0, unit activity |
  | `.E_standard_prime` | **E°′** — conditions' pH, unit activity otherwise |
  | `.E` | **E** — full conditions (concentrations, partial pressures, ionic strength) |

  E°′ is the headline number (larger font / first column), since that is what
  the user asked for by name; E° and E sit alongside it rather than behind a
  click, because the panel's whole pedagogical value is watching E diverge
  from E°′ as conditions move away from standard.

**Errors:** `SpeciesNotFoundError` per §5 for either field independently (so
picking a bad acceptor doesn't clear a valid donor field).

### 6.3 Reaction balancer, two half reactions, reveal $\Delta G^{\circ\prime}$

**Primary input:** one free-text equation field, e.g. `NO3- + H2 -> NH2OH`,
calling `mt.Reaction.from_equation(equation, conditions=..., normalize_to="electron_pair")` —
`normalize_to="electron_pair"` is the library's own default and is exactly
"balance to 2 e⁻", so no override is needed for the common case. Below the
equation field, the same shared **Conditions panel** (§6.4).

**Fallback when the equation doesn't resolve cleanly** — this is the one
place where the free-text choice (over couple pickers) has a real cost, and
the UI needs to absorb it rather than dead-end the user:

- `AmbiguousReactionError` (disproportionation, or an underdetermined
  equation): render the solution basis as a table (mirroring the README's
  `balance_equation` example), and reveal a small "or specify the couples
  directly" expander with two `(reduced, oxidized)` input pairs that calls
  `Reaction.from_couples` instead. This is the same escape hatch the README
  itself recommends ("Name the couples with `from_couples` in that case") —
  the app should not invent a different one.
- If `infer_couples` cannot find exactly one donor and one acceptor at all,
  same fallback, different message.

**Output, once balanced**, in two clearly separated cards (one per half
reaction, donor left / acceptor right, matching the tower and half-reaction
figure's existing left-right convention):

- The half-reaction string (`reaction.donor_half`, `reaction.acceptor_half`),
  labeled oxidation/reduction.
- $E^{\circ\prime}$ for that half reaction, shown immediately (this is *not*
  behind the reveal — the user asked for "the two half reaction strings and
  the E°′ values for each" as the up-front output).
- $n$ electrons transferred (should read 2, given the electron-pair default)
  and the balanced overall equation as a single line underneath both cards.

**The reveal button**, directly under the two cards:

- A single `st.button("Reveal ΔG°′")` (or `st.toggle`, so it can be
  re-hidden) that, on click, shows `reaction.delta_G_standard_prime` in
  kJ/mol.
- Since full `Conditions` exposure was chosen, the same reveal also shows
  `reaction.delta_G` (in-situ, under whatever concentrations/pH/ionic
  strength the sidebar is set to) directly beneath ΔG°′, labeled distinctly,
  so a student comparing standard-state to real conditions doesn't have to
  hunt for the second number. ΔG°′ stays the first, larger line since it's
  what was explicitly asked for.
- This is a genuine answer-reveal pattern (consistent with the "problem
  generator and grader" already on the README's future-work list, item #12)
  — worth building as a real toggle component now, since #12 would want the
  identical interaction later.

**Unverified-species footnote:** if `reaction.provenance().unverified` is
non-empty, show the same red footnote convention as the figures, listing
which species it rests on (reusing `Provenance`, not a new mechanism).

### 6.4 Shared conditions panel (sidebar)

One `st.sidebar` form, live across all three tabs (tool 1 ignores everything
but temperature/pressure per §6.1), mapping directly onto `Conditions`'
fields — this is the "full `Conditions()` exposure" the questionnaire chose,
not a curated subset:

| Widget | `Conditions` field | Notes |
|---|---|---|
| Number input, °C | `temperature_c` | Default 25.0; range-checked against `OutOfRangeError` client-side (0.01–100) before calling the backend, so the error message is instant rather than round-tripped |
| Number input, pH | `pH` | Default 7.0 |
| Select: ideal / bdot / unit | `activity_model` | Default `"ideal"`, matching the dataclass default |
| Number input, ionic strength | `ionic_strength` | Only shown/enabled when `activity_model == "bdot"`, since `bdot` requires it and the dataclass raises otherwise — disable the incompatible combination in the UI instead of surfacing that `ValueError` |
| Repeatable key/value rows | `concentrations` | "species name → molality"; validated against `resolve()` on blur, same suggestion UX as §6.1 |
| Repeatable key/value rows | `total_concentrations` | "family name → total molal concentration"; a small help text lists the twelve supported families from the README (sulfide, carbonate, ammonia, phosphate, acetate, lactate, formate, propanoate, butanoate, sulfite, nitrite, sulfate) |
| Repeatable key/value rows | `partial_pressures` | "gas species → bar" |
| Number input | `pressure_bar` | Optional; blank means pyGCC's own default/saturation handling |
| Button | — | "Reset to standard state" — clears everything back to `Conditions()` defaults, since a panel this rich needs an obvious way back |

A one-line caption under the panel states the biochemical-convention caveat
verbatim from the README ("ionic strength does not reconcile this with an
Alberty-convention table; it moves ΔG at stated concentrations but not
ΔG°′, which is defined at unit activity"), so a student who cranks ionic
strength expecting ΔG°′ to move gets told why it didn't, in place, rather
than filing that as a bug.

---

## 7. Explicitly out of scope for v1 (future work)

Kept separate from §1 so a later pass has a ready list rather than a blank
page:

- **The three static figures and the interactive explorer** — these are
  matplotlib/Plotly objects the library already builds; Streamlit can render
  a matplotlib figure with `st.pyplot` and a Plotly one with
  `st.plotly_chart` with essentially no new code, so this is a natural v1.1,
  not a redesign. If it lands, re-open the Panel-vs-Streamlit question from
  §2 — Panel's native Bokeh/HoloViews integration may be worth it once
  figures are the majority of the app's surface, but not before.
- **`reaction.show_work()`** — an `st.expander("Show the full derivation")`
  rendering `.to_text()` would fit this app's teaching intent well and is a
  small addition once §6.3 exists, but the user asked specifically for the
  two half reactions, their potentials, and a ΔG°′ reveal — adding the full
  eight-step derivation now would be answering a question that wasn't asked.
- **The curated metabolism library** (`mt.library`) as a "pick a named
  metabolism" dropdown, feeding tool 3 instead of typing an equation.
- **`balance_equation(..., choose="minimal")`** as a one-click "just pick one"
  option in the ambiguous-equation fallback (§6.3), instead of only offering
  the couple-picker escape hatch. Deferred because the README is explicit
  that `choose="minimal"` is "an arithmetic preference, not a chemical one"
  and "never for a number you intend to quote" — exposing it in a teaching
  tool needs its own warning-and-confirmation UX, not a checkbox.
- **Provenance export** (`reaction.provenance().to_text()` /
  `.to_bibtex()`) as a download button next to the reveal — cheap to add,
  deferred only because it wasn't asked for.

### 7.1 Standalone interactive pages, as an alternative to rebuilding figures

*Added 2026-09-24, after the interactive element diagrams landed in the
library. This is a note for whoever picks up the figures bullet above.*

The library already exports **self-contained interactive HTML**: one file,
Plotly bundled inside it, no server, no Python, no network. Opening it in a
browser gives working sliders. That changes the calculus of the first bullet
in §7, because the choice is no longer "rebuild the figure's controls in
Streamlit or go without" — the third option is to serve a file the library
already knows how to write.

**What exists today**, all committed and regenerated by running the notebook
that owns them:

| Page | Written by | Notebook |
|---|---|---|
| `notebooks/tower_interactive.html` | `plot_interactive_tower(grid, save_html="tower_interactive")` | `02_redox_tower_and_energy` |
| `notebooks/explorer_methanogenesis.html` | `plot_energy_explorer(reaction, save_html="explorer_methanogenesis")` | `02_redox_tower_and_energy` |
| `notebooks/arsenic_interactive.html` | `plot_interactive_element("As", save_html="arsenic_interactive")` | `10_arsenic` |

**The recommended naming convention for the rest**, should the app want a
per-element page: `notebooks/<element-name>_interactive.html`, spelled out
rather than by symbol, to match `arsenic_interactive.html`.

```python
from microbial_thermo.figures import plot_interactive_element

for element, name in [
    ("C", "carbon"), ("N", "nitrogen"), ("S", "sulfur"), ("Fe", "iron"),
    ("Mn", "manganese"), ("As", "arsenic"), ("Se", "selenium"),
    ("Cr", "chromium"), ("U", "uranium"), ("Cu", "copper"),
]:
    plot_interactive_element(element, save_html=f"notebooks/{name}_interactive")
```

Each page carries Latimer, Frost and Eh–pH for that element, with sliders for
temperature, pH and dissolved activity, plus `show` and `label` dropdowns for
the Frost panel. What travels in the file is the *decomposition* — base
energy, proton and electron coefficients per species per temperature — and
the page rebuilds all three diagrams in JavaScript, convex hull included. The
test suite runs that script in node and checks it against the library, so the
page and the library agree to 1e-9.

**How the app would use them.** Either as plain links from the relevant tool,
or embedded with `st.components.v1.html(path.read_text(), height=700)`. No
new caching, no `st.pyplot`, no figure rebuilt per rerun — and it sidesteps
re-opening the Panel-vs-Streamlit question in §2 entirely, because the
interactivity is not Streamlit's problem.

**Three caveats, none fatal:**

1. **Size.** Each page is about 4.3 MB, because Plotly is bundled. Ten
   element pages is ~43 MB in the repo. `write_html` takes
   `include_plotlyjs="cdn"`, which drops each file to roughly 100 KB at the
   cost of needing the network — but the library currently hardcodes `True`
   in all three export paths, so exposing that is a one-line change per
   function and should be done before generating a full set.
2. **Fixed at export time.** The species list and the temperature grid are
   baked in when the file is written. A viewer can move the sliders the page
   was given and nothing else; there is no way to ask it about a species the
   exporter did not include. A Streamlit-rendered figure could take live
   input, so if the app ever wants "diagram for whatever the user just typed",
   this approach cannot do it and the §7 bullet stands.
3. **Not a substitute for the three v1 tools.** These are figures. The lookups
   and the balancer in §6 are text tools and remain the app's actual job.


---

## 8. Open questions resolved, and why

Recorded so these don't get relitigated (per the standing convention in
README's Future work §1):

1. **Deployment target** — "local by default, easy to host later." Resolved
   by choosing a framework (Streamlit) whose local and hosted deployment
   paths are the same code, and documenting both in §4 rather than picking
   one and bolting the other on later.
2. **Framework** — Streamlit, confirmed. Reasoning in §2.
3. **Conditions scope** — full `Conditions()` exposure, not fixed
   standard-state-only. This is the more expensive UI (§6.4 vs. a two-field
   version) but was chosen deliberately over the simpler "just E°′/ΔG°′"
   option; §6.2 and §6.3 both therefore show all three potentials / both
   free energies rather than hiding the standard-state ones behind a toggle.
4. **Reaction input** — free-text equation (`from_equation`), not couple
   pickers. This is the choice with real UX cost (ambiguity handling, §6.3's
   fallback), taken because it reads closer to how a student writes a
   reaction on paper. The couple-picker path is kept as a fallback, not
   dropped, exactly because free-text entry cannot express every reaction
   the library can (disproportionation with an unclear split, multi-product
   couples, VOSO₄-style cases the README documents as needing explicit
   `Couple.make`).

No open questions remain that block writing `app.py`. The one thing worth a
short conversation before implementation starts is §7's list — whether any
of those should actually be pulled into v1 rather than deferred.
