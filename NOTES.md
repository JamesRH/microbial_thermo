# Notes for whoever works on this next

`AGENTS.md` sets the rules; `SPEC.md` describes the design and §15 records where
the implementation departed from it; `README.md` is the user-facing guide. This
file is the accumulated operational knowledge — the things that cost time to
discover and are not obvious from reading the code.

---

## Getting going

```bash
source setup.sh          # builds the env if missing, activates, installs, registers the kernel
```

It must be **sourced**. `pygcc` comes from PyPI (not conda-forge) and is the one
sanctioned `pip` in the project.

```bash
python -m unittest discover -s tests     # full suite, 465 tests, ~165 s
MT_SKIP_NOTEBOOKS=1 python -m unittest discover -s tests   # ~130 s
python tests/test_balance.py             # one file, seconds — use this while iterating
ruff format microbial_thermo tests && ruff check microbial_thermo tests
```

The suite now executes the numbered teaching notebooks, which is the slowest
part (~25 s) and the part most likely to catch a rename you did not expect.

---

## What pyGCC will and will not do

These were all established by probing the installed package, not by reading its
docs, and several contradict what the docs imply.

**There is no arbitrary-reaction API.** `calcRxnlogK` only evaluates a species'
*database-defined* formation reaction. Every reaction energy in this library is
assembled as `sum(nu_i * dGf_i)` from per-species `supcrtaq` calls. That is the
core numerical path and there is no shortcut around it.

**Water is not a species.** `H2O` is absent from `dbaccessdic` because HKF
parameterises solutes, not the solvent. It comes from `iapws95(T=celsius, P=P).G`.

**Water needed a datum correction.** IAPWS-95 gives −237.140 kJ/mol; the rest of
the data assumes the SUPCRT convention of −237.18. A reference-state difference,
not an error, but it left a systematic 0.041 kJ/mol error *per mole of water*,
and water appears with coefficients of four or more in common half reactions.
`SUPCRT_WATER_GIBBS_CAL_25C` shifts it onto the right datum, preserving the
temperature dependence. `water_convention="iapws"` opts out.

**Gases and minerals are not HKF solutes.** They carry Maier–Kelley coefficients
and must go through `heatcap`, not `supcrtaq`. The discriminator is entry length:
13 parameters means aqueous HKF, anything else means `heatcap`.

**There is a cliff at exactly 100 °C.** At 1 bar water would be steam and
`supcrtaq` returns `NaN`. Only pyGCC's `P='T'` saturation sentinel works there —
a literal 1.013 also fails. Pressure follows the saturation curve at and above
99 °C. 0 °C exactly crashes IAPWS-95; the floor is 0.01.

**Vectorising over temperature is slower than looping** (507 ms for 10 points
against 284 ms). There is no vectorisation win. Use scalar calls behind the
cache, and precompute grids for anything interactive.

**Species names are legacy SUPCRT/GWB**: `Fe++`, `SO4--`, `S2--`, `Methane(aq)`,
`Sulfur(s)`, `Acetate`. Not `CH4(aq)`, not `Fe+2`. That is what `species.yaml`
exists to hide.

**`resolve()` only sees `species_dict`, not minerals.** `delta_Gf` handles all
three routes (HKF, GWB minerals, supplemental) but `backend.resolve()` checks
only the HKF dictionary, so it raises for Goethite, Pyrolusite and
`As(OH)3(aq)` while those species work perfectly well in reactions. The error
message then suggests the name you just passed, which looks like a bug and is
really this. Nothing downstream depends on `resolve` for those, so it has been
left alone -- but do not use `resolve` as an existence check.

**Arsenic has two As(III) representations differing by one water.**
`As(OH)3` / `H2AsO3-` (arsenous) and `HAsO2` / `AsO2-` (metarsenous) are the
same chemistry, and they agree to 0.21 kJ/mol over 0-100 C. They arrive by
*different routes* -- `As(OH)3(aq)` through GWB log K, `HAsO2(aq)` as direct
HKF -- so that agreement is real evidence. `HAsO2(aq)` is the single arsenic
species the two databases disagree on (0.26 kJ/mol, declared in
`KNOWN_DISAGREEMENTS`), which is why a bare `As(III)` resolves to
`As(OH)3(aq)`.

**Manganese oxides are in no direct-access database.** They come from GWB log K
values combined with HKF basis species — see `backends/gwb.py`. That mixing of
sources is justified in `tests/test_database_consistency.py`, which should be run
if the databases ever change. Three known inter-database disagreements are
declared there with their magnitudes.

---

## Conventions that must not drift

Breaking any of these silently produces plausible wrong numbers.

- **Couples are always `(reduced, oxidized)`.** Everywhere, without exception.
- **Half reactions are stored as reductions**, whatever direction they display
  in. That is what makes potentials comparable.
- **ΔE = E_acceptor − E_donor**, and ΔG = −nFΔE.
- **ΔGf(H⁺) = 0 and ΔGf(e⁻) = 0**, which puts everything on the SHE.
- **Stoichiometry is `Fraction`, never float.** Normalising to an electron pair
  must give exact quarters.
- **Figures normalise to an electron pair by default**, whatever the reaction was
  built with.
- **The two-path cross-check runs on every `Reaction` construction**, not only in
  tests. Do not make it optional to "speed things up"; redox sign errors are the
  characteristic failure here and this is what catches them.

---

## Traps

**`ruff format` reflows code between edits.** A string-replacement patch written
against what you last read will silently match nothing. Always confirm a patch
applied — several of mine reported success while changing nothing.

**Appending test classes after `if __name__ == "__main__"`** leaves them
unreachable when the file is run directly. `unittest discover` still finds them,
so it hides for a long time. One file had 28 tests and exposed 0 that way.

**`pkill -f <pattern>`** matches its own command line and will kill the shell
running it. It has already killed one session mid-command and lost an unwritten
file.

**A bare `+` is a charge, not a separator.** `Fe+3` is one species; only a `+`
with whitespace on both sides separates two. A `+` followed by a letter is a
missing-space error.

**Negative `n_electrons` means the couple was written backwards.** The potential
survives it — E = −ΔG/(nF) flips numerator and denominator together — so the
number looks right while the half reaction formats as an oxidation. That sign is
the only tell.

**A bare `H2` resolves to `H2(aq)`, not the gas.** The gas is 91 mV away at
pH 7. This is now *declared* rather than accidental: `canonical: true` in
`species.yaml` picks the winner for each contested name, and
`tests/test_ambiguity.py` fails if a new aqueous/gas pair is added without
declaring one. It used to depend on file ordering.

**Couple order does not set reaction direction.** `Couple.make(a, b)` is
always `(reduced, oxidized)`, but writing one backwards does *not* flip the
reaction: direction comes from which couple is passed as `donor=` and which as
`acceptor=`. A backwards couple gives identical n, E and dG, and only swaps
the `reduced`/`oxidized` labels and the oxidation-state annotations. So the
symptom is a mislabelled figure or a wrong `couple.reduced.backend`, never a
wrong energy. Confirmed by building both ways.

**`Couple` has no `.backend`.** The backend name lives on the `Species`, at
`couple.oxidized.backend` / `couple.reduced.backend`. For a multi-species side
use `couple.oxidized_side`, a tuple of `(Species, coefficient)` pairs -- note
a simple couple's side has length 1, so indexing `[1]` raises.

**An E°' is meaningless without naming the species.** The arsenate case makes
this unmissable: As(V)/As(III) at pH 7 reads +0.160 V written against H3AsO4,
+0.020 V against H2AsO4-, and +0.013 V against HAsO4--, purely because the
three differ in proton count. When checking a couple against a published
number, check the *standard-state* form first -- that convention is
unambiguous -- and only then argue about pH 7.

**Measure matplotlib text extents with the Agg renderer**, even when exporting
SVG; the SVG backend's metrics differ. And size equations by the *laid-out
extent*, not the sum of token widths — the latter omits inter-token gaps and lets
long equations overrun.

**`R` and `FARADAY` have magnitude 1.** They are defined as
`1 * ureg.molar_gas_constant` and `1 * ureg.faraday_constant`, so `R.magnitude`
is 1, not 8.314. Convert first: `R.to("J/(mol*K)")`. Likewise a bare float
temperature will not cancel against R — it needs `kelvin * ureg.kelvin`.

**nbclient kernels hang intermittently, and it is not the library's fault.**
Seen three times: a kernel starts, goes idle mid-notebook, and the client waits
forever on unchanged notebooks that pass on the next run. Once was a laptop
suspend; twice was not, so do not write it off as sleep the way I first did.
`tests/test_notebooks.py` now retries once on `CellTimeoutError` or
`DeadKernelError` and caps cell time at 240 s -- a `CellExecutionError` is a
real notebook bug and is never retried. When a run seems slow, check `ps` for a
live `ipykernel_launcher` at ~1% CPU before assuming the suite got heavier; the
giveaway is wall time far above the usual ~165 s with no new output.

**Plotly's native sliders cannot express independent dimensions.** A slider's
steps cannot read the other sliders' positions, so two sliders cannot select a
cell of a 2-D grid. The interactive tower's exported HTML builds its own
sliders in JavaScript instead. If you add another interactive figure with more
than one control, do not spend time trying to make `updatemenus` do it.

**The user runs JupyterLab against `notebooks/`.** Do not kill it. Their edits
arrive as working-tree changes; read a diff before assuming it is yours, and
never commit `-Copy1` duplicates (gitignored).

---

## Verify before asserting

The most useful discipline on this project, learned the hard way. Three claims I
wrote into prose or commit messages turned out to be wrong, each of which would
have read as plausible:

- proton counts in the notebook's pH commentary (said four; were two and ⅔),
- an acceptor-strength ordering that the printed numbers contradicted, and which
  was confounded anyway because the rows used different donors and pH,
- an equality check reported as failing when I had rounded the Faraday constant
  in the checking script.

Run the number before writing it down. For anything appearing in a notebook or a
figure, compute it independently rather than reusing the value under test.

Tests that parse rendered output back out — `tests/test_showwork.py` — exist for
this reason: they catch the narrative drifting from the arithmetic, which
recomputing the formula inside the test would not.

---

## Where the planned work is recorded

`README.md`, *Future work*, in three tiers by effort and dependency, with a
*Completed* table above them. **Item numbers are stable** — completed items keep
their number rather than being renumbered, because this file and the commit
messages refer to them by number. Items 3, 5, 6, 7, 10 and 11 are done.

Two of the remaining entries deserve reading before touching anything nearby:

- **Item 1** records the argument *against* the Legendre-transformed standard
  state, with the reasoning, so it does not get relitigated. The short version:
  this library already matches the microbial bioenergetics convention, the two
  agree exactly for single-species reactants, and ionic strength is the larger
  discrepancy.
- **Item 2** was the one real correctness debt. Glucose and pyruvate are now
  traced through CHNOSZ's OBIGT to Amend & Plyasunov (2001) and Canovas &
  Shock (2016); glucose moved 4.7 kJ/mol in the process, which is what an
  untraced "commonly cited" value is worth. **Hydroxylamine is still
  outstanding** and is in neither pyGCC nor OBIGT; its provenance records
  where it was looked for. `Biomass(aq)` is unverified permanently by design.
  **OBIGT is the place to look first** for anything pyGCC lacks: it is
  SUPCRT-lineage, so the standard state already matches, and every entry
  carries a citation key.
