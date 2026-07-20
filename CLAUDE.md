# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A proof-of-concept to replace DIYABC's C++ coalescent simulator
(`particuleC.cpp::dosimulpart`) with `msprime`. Originally scoped to the
simplest possible case (scenario 1 of the `human` dataset: 4
populations, merges only, no split/admixture), the pipeline now covers
the condensed SNP loci description format completely: all heritage
types (`<A>/<H>/<X>/<Y>/<M>`), split/admixture events and multi-scenario
weighted draws, and MAF filtering — validated against real DIYABC
output on `human`, `toy_example5`, and `toy_example3` (which exercises
split/admixture, scenario 3). The one remaining gap is the *detailed*
loci description format (one locus named per line, used by
sequences-mut datasets, as opposed to the condensed format) — not
implemented, see `loci_parser.py`. The goal is to demonstrate that a
`header.txt` → `msprime.Demography` → coalescent+mutation → summary
statistics pipeline built in Python produces a `reftable.bin`
structurally and statistically equivalent to the real DIYABC's.

`reference/human/` holds the ground truth produced by the real DIYABC
binary — never modify these files. `notes/exploration.md` is the running
research journal of everything reverse-engineered from the DIYABC C++
source (`particuleC.cpp`, `history.cpp`, `sumstat.cpp`, `data.cpp`,
`header.cpp`); read it before touching parsing or statistics code, it
records *why* things are done a certain way, including several
non-obvious findings (see Architecture notes below).

## Environment

Use the `diyabc_msprime` conda environment (Python 3.11, has msprime
1.4.2, tskit, numpy, scipy):

```bash
conda activate diyabc_msprime
```

The system Python (3.13, no conda env) does NOT have msprime installed —
always activate `diyabc_msprime` before running anything in this repo.

## Common commands

```bash
# Run the test suite
pytest tests/ -v

# Run a single test
pytest tests/test_scenario1_human.py::test_scenario1_events -v

# Lint / format (ruff, config in pyproject.toml)
ruff check .
ruff format .

# Pre-commit hooks (ruff check --fix + ruff format)
pre-commit run --all-files
```

Some tests (`test_compute_summary_statistics_scenario1`,
`test_run_reftable_simulation_scenario1`, and `test_write_reftable_bin`'s
`general_binary_path` arg) are skipped unless `DIYABC_GENERAL_PATH` is
set to a compiled DIYABC `general` binary — these compare against the
real C++ implementation and are optional for pure-Python pipeline work.

Regenerate `notes/tree.md`, `notes/commits.md`, `notes/api.md`, and
`notes/report.md` (a rollup of all three) with:

```bash
python3 tools/generate_report.py
```

## Architecture

The pipeline (`bridge/`) turns a DIYABC `header.txt` + `.snp` file into
a `reftable.bin`, entirely in Python, without any subprocess call to the
DIYABC binary. Data flows through these stages, each stage a separate
module with no cross-cutting logic:

1. **`scenario_parser.py`** — `header.txt` text → `Scenario` objects.
   Splits on `scenario N ... (...)` blocks and parses event lines
   (`sample`, `merge`, `varNe`, `split` for admixture — see `SplitEvent`
   in `scenario_types.py`) using the vocabulary from
   `history.cpp::ScenarioC::read_events`. Any vocabulary still
   unimplemented is skipped with a warning rather than failing the whole
   parse — only `NotImplementedError` is swallowed here, any other
   exception propagates.

2. **`prior_parser.py`** — extracts `Prior` and `OrderConstraint` objects
   from the `historical parameters priors` section.
   `is_constant_prior` replicates the exact near-degenerate-bounds rule
   from `readReftable.R` / `abcranger/readreftable.cpp` (constant if
   `(max-min)/max <= 1e-6`, never constant when `max == 0.0`) — DIYABC
   excludes these from `reftable.bin` columns.

3. **`loci_parser.py`** — parses the `loci description` section (condensed
   format only, not the detailed one-locus-per-line format used by
   sequences-mut). Handles both the single-heritage-type line
   (`"5000 <A> G1 from 1"`) and the multi-type line (`"70 <A> 10 <X> 10
   <M> 10 <Y> G1 from 1"`) — `LociDescription.total_loci` is a
   `dict[heritage_type, count]`, not a plain int, precisely to represent
   the multi-type case. `rewrite_loci_count` (used to test with a
   reduced loci count without editing header.txt by hand) is currently
   still single-type only and raises `NotImplementedError` on a
   multi-type line — not yet updated to match `parse_loci_description`.

4. **`parameter_sampling.py`** — draws one value per prior via
   rejection sampling until all `OrderConstraint`s are satisfied.

5. **`demography_builder.py`** — `Scenario` + drawn values →
   `msprime.Demography`. Populations are named `"pop1".."popN"` by their
   1-indexed position in `header.txt`. `get_parameter_names_used_by_scenario`
   determines which prior names a *specific* scenario actually references
   (important: `human/header.txt` declares 21 priors globally, but
   scenario 1 only uses 16 — the other 5 belong to scenarios 2-6; this
   subset, not the full prior list, must become the `reftable.bin`
   parameter columns, per a bug found via `readReftable.R` failing with
   "indice hors limites").

6. **`observed_data.py`** — maps population index (1,2,3,4 from
   `header.txt`, which never names populations) to real population name
   in the `.snp` file (e.g. `{1: "ASW", 2: "YRI", ...}`). The mapping is
   implicit: population *i* in the scenario = the *i*-th population by
   first-appearance order in the `.snp` file — there is no cross-reference
   in the DIYABC C++ source, this was verified by exhaustive code search.
   Never replace the order-preserving dict/Counter here with anything
   that could reorder keys (e.g. alphabetical sort) — it would silently
   break this mapping. Also parses two other things from the `.snp`
   header line: the sex-ratio (`parse_sex_ratio`, needed for `<X>/<Y>/<M>`
   ploidy/coalescence-coefficient) and the MAF threshold (`parse_maf_ratio`,
   `<MAF=N%>` or `<MAF=hudson>`, consumed by `with_maf_filter`/
   `with_maf_filter_shared_ancestry` in `ancestry_simulation.py`).

7. **`ancestry_simulation.py`** — simulates one independent tree per SNP
   locus (`simulate_independent_loci`, no recombination/linkage between
   loci) and mutates each with the Hudson algorithm
   (`simulate_snp_genotypes`): exactly one mutation per locus, placed on
   an edge chosen with probability proportional to branch length, fully
   vectorized over tskit's edge tables (not a per-node `branch_length()`
   loop). This guarantees every locus is globally polymorphic by
   construction — see `notes/exploration.md` for the DIYABC doc citation
   (section 2.4.3) and empirical validation. MAF filtering (`<MAF=N%>`,
   as opposed to `<MAF=hudson>`) is implemented (`with_maf_filter`/
   `with_maf_filter_shared_ancestry`, reject-and-resimulate a locus below
   the threshold, mirroring `ParticleC::mafreached`) — not needed for
   `human` (declares `<MAF=hudson>`, the no-op fast path) but used by
   some `toy_example3`/`toy_example5` test datasets.

8. **`snp_writer.py`** — writes simulated genotypes out in DIYABC `.snp`
   format (only used for the deprecated subprocess-based path, see below).

9. **`summary_statistics.py`** — pure-Python/numpy reimplementation of
   all 130 SNP summary statistics from `sumstat.cpp` (ML1-3, HW, HB,
   FST1-4, NEI, AML, F3, F4), validated column-by-column against the real
   `general` binary's output. This is the **current** path — no subprocess,
   no intermediate `.snp` file, ~2.4s/particle on 5000 loci single-threaded
   (was ~343s via subprocess before batch-size and vectorization fixes,
   see below; the intervening "~7s/particle" figure was stale — see
   "Known performance gap" below for the current, measured number and
   where that time actually goes).
   Formula provenance for each function is documented via exact
   `sumstat.cpp` function names in each docstring — check there before
   changing a formula.

10. **`pipeline.py`** — orchestrates 1-9 with no logic of its own.
   `compute_summary_statistics` is the main high-level entry point (its
   docstring still mentions delegating to the C++ binary; that's stale —
   it now calls `summary_statistics.compute_all_statistics` directly).

11. **`reftable_loop.py`** — runs `nrec` independent particles in
    parallel (`ProcessPoolExecutor`, one work directory per particle) and
    writes the DIYABC binary `reftable.bin` format (`write_reftable_bin`,
    verified against `reftable.cpp` / `readReftable.R` /
    `abcranger/readreftable.cpp`) as well as a human-readable
    `write_reftable_txt` (mirrors DIYABC's
    `first_records_of_the_reference_table_0.txt`, for direct diffing
    against real DIYABC output). Seeds are derived as `particle_index + 1`
    (msprime rejects `seed=0`). **Multi-scenario is supported**:
    `run_reftable_simulation` draws each particle's own scenario from a
    weighted list of candidates (`parameter_sampling.draw_scenario`,
    matching `particuleC.cpp::ParticleC::drawscenario`), and both writers
    handle rows that mix different `scenario_index` values —
    `write_reftable_bin` writes a variable-length record per row (only
    that row's own scenario's `nparam` columns, matching `reftable.cpp`'s
    own behavior, no NA-padding), `write_reftable_txt` uses a fixed union
    of parameter columns across all candidate scenarios. Validated
    end-to-end against real DIYABC multi-scenario reftables on
    `toy_example3` (exercises `SplitEvent`/admixture, scenario 3 —
    confirmed by the user 2026-07-20) — see also
    `tests/test_scenario1_human.py` for unit-level coverage of the
    weighted draw + `SplitEvent`/admixture translation.

### Two generations of architecture — mind the drift

An earlier version of this pipeline (`bridge/` history around
`compute_summary_statistics`) wrote simulated genotypes to a fake `.snp`
file and shelled out to the real DIYABC `general` binary
(`-R "ALL" -g 1`, `-g` is the *internal batch size*, not the loop count —
using `-g 50` caused a 29x slowdown by silently discarding 49/50
simulated particles) to compute statistics via `statobs_parser.py`. That
path is kept only for cross-validation (see the `DIYABC_GENERAL_PATH`
skip-marked tests) — it is no longer the default and is far slower.
`msprime_cpp/msprime_from_cpp.cpp` is an exploratory spike (embedding a
Python interpreter in a C++ process to call msprime directly, avoiding
subprocess overhead) — not wired into the pipeline, kept for reference
only.

### Known limitation (mostly resolved 2026-07-10 — see below)

`notes/exploration.md` (2026-07-03 entry) documents that an early
head-to-head statistical comparison (mean/median/stdev diffs +
two-sample Kolmogorov-Smirnov) between 1000-particle reftables from
real DIYABC vs. this pipeline found 124-126 of 125 variables differ
significantly.

**Root cause found 2026-07-10** (see `notes/exploration.md`'s
2026-07-10 entry and the persistent memory
`diyabc_header_trailer_line_bug`): the last line of header.txt/
headerRF.txt (`"scenario N1 N2 N3 ta ts ML1p_1 ..."`), which looks like
pure output-column documentation, is actually **re-read as input** by
`HeaderC::readHeaderAllStat` (header.cpp) to derive the historical
parameter count (`nparamhist = header_lastline.size() - 1 - nstat -
nparamut`) — counted from this line's token count, not from the real
`historical parameters priors` declarations. Test headers built by
copy-pasting this trailer line from a different scenario (e.g.
`toy_example5_scenario1`, used throughout the investigation, declared
"6 parameters" but only had 5 real priors, with two extra phantom
tokens `N4 r` left over) get a miscounted `nparamhist`, which corrupts
DIYABC's internal state and produces wildly wrong statistics (300% to
10000% off) for every population except the scenario's "hub" — exactly
the pattern chased throughout the investigation.

Rebuilding a test header with a correct trailer line (token count
matching the real declared+used priors) drops the significant-stat
count from ~48/55 to ~22/55 and single-population statistics (ML1p,
HWm) from 300-10000% off to <3%, not significant. **When hand-crafting
or copying a headerRF.txt/header.txt for a new test scenario, always
regenerate this trailer line to match that scenario's own priors —
never copy it from another scenario/template.**

A modest residual bias (~3-8%) remains on pairwise statistics
(FST2/NEI/F3/HB) even with a correct header — confirmed 2026-07-17 via
a from-scratch `jupyter nbconvert --execute` rerun of
`reference/toy_example5_1000loci/compare_reftables_te5_1000loci.ipynb`
(exact prior replay, 1000 particles × 650 loci) to **not be
statistically significant** (0/50 stats with KS p<0.05, min p=0.12) —
likely irreducible per-particle/per-locus noise at this loci count,
not a simulator bug. A much larger anomaly previously logged here
turned out to be unrelated to the simulators themselves — see
`notes/exploration.md`'s 2026-07-17 entry for detail. Confirmed by the
user 2026-07-20: the same paired-replay
validation also holds on `human` (5000 loci, the project's actual
reference dataset, scenario 1), plus `toy_example5` and `toy_example3`.

### Known performance gap (investigated 2026-07-20 — not a bug)

On `human_modif_scenario1_5000loci` (1000 particles, 5000 loci,
scenario 1, real DIYABC priors replayed via `replay_diyabc_priors.py`,
`max_workers=16`), this pipeline takes ~384s vs. ~137s for real DIYABC
on the same config (`-t 16`) — a ~2.8x gap. Investigated in full in
`notes/exploration.md`'s 2026-07-20 entry; summary:

- **Not the stats formulas or the simulator**: single-threaded
  per-particle cost is ~2.35s here vs. ~2.19s equivalent for DIYABC —
  nearly identical.
- **Not a `max_workers` misconfiguration**: this dev machine has 8
  physical cores / 16 logical threads (hyperthreading, of little help
  for this CPU-bound numeric workload) — `max_workers=16` still beats
  8 and 7 in absolute wall time, measured empirically.
- **The real cost**: ~93% of one particle's time is in the per-locus
  simulate+mutate loop, and within that, ~40% was Python/tskit overhead
  of materializing a full `TreeSequence`/`Tables` object (incl.
  re-decoding population metadata) for each of the 5000 independent
  loci — more than the actual msprime C coalescent engine itself
  (~30%).

**Partly fixed 2026-07-20**: `simulate_snp_genotypes`
(`ancestry_simulation.py`) now computes the population names + sample
IDs per population ONCE (from the first locus) instead of re-decoding
them for all 5000 — safe because every replicate of the same
`simulate_independent_loci`/`simulate_shared_ancestry_loci` call shares
the same `demography`/`samples`, only the coalescent topology differs
(verified empirically). Measured gain: ~17% single-threaded
(2.35s → 1.96s/particle), ~26% on the full 1000-particle parallel run
(384s → ~284s). Remaining gap vs. DIYABC (~284s vs. 137s) is the
incompressible part: the 8-physical-core ceiling and materializing a
`TreeSequence` per locus at all. `simulate_snp_genotypes` now accepts
an optional pre-computed `population_layout` (see `_population_layout`)
so `with_maf_filter`/`with_maf_filter_shared_ancestry`'s rejection loop
(maf>0 datasets) can also cache it across attempts — but the measured
gain there is much smaller (~2-3% on `toy_example3_scenario1`, vs. ~17%
on the `maf=0.0` fast path): each rejection-loop attempt already
simulates a single locus via its own `simulate_independent_loci`
call, so the eliminated per-locus metadata decode is a much smaller
share of that per-attempt cost than it was when one call handled 5000
loci at once.

**Also fixed 2026-07-20**: `build_samples_argument` was scanning the
`.snp` file twice internally (`population_index_to_name`, which itself
calls `count_samples_per_population`, followed by a second independent
call to `count_samples_per_population`) — now a single scan. And
`pipeline.py`'s `compute_summary_statistics`/`compute_summary_
statistics_from_values` no longer call `build_samples_argument` a
second time just to get `population_names`: they're read for free off
the keys of the already-simulated `genotypes_list`'s first locus
(`_population_names` helper), which `simulate_snp_genotypes` already
populates with the same names. Combined measured gain (all three fixes
together): ~1.84s/particle single-threaded (from 2.35s). Confirmed on
the full 1000-particle run (not a sub-sample extrapolation): **300s**
(from 384s) — **~22% overall**, 0 regressions (62/62 tests green
throughout). Remaining ~2.2x gap vs. DIYABC (~300s vs. 137s) is the
incompressible part described above.

Not blocking for this POC (goal: prove feasibility, already done) —
don't re-investigate the formulas or `max_workers` if this comes up
again, the cause is understood.

### `scripts/` — ad hoc investigation scripts

The many `.py` files under `scripts/` (`run_test.py`,
`generate_test_reftable.py`, `calibrate_reftable.py`, `validate_stats.py`,
`param_keepers.py`, `priors_keeper.py`, `benchmark_1000.py`, `profile_*.py`)
are ad hoc scratch/investigation scripts against `reference/human`, not
part of the `bridge/` package or its test suite — treat them as disposable
throwaways when reasoning about the actual architecture, but don't delete
them without checking with the user first since they double as informal
experiment logs.

They `import bridge` and use paths relative to the project root (e.g.
`"reference/human/header.txt"`), so run them as modules **from the repo
root**, not as bare scripts:

```bash
python3 -m scripts.run_test
```

Running `python3 scripts/run_test.py` directly fails with
`ModuleNotFoundError: No module named 'bridge'` — Python only puts the
script's own directory on `sys.path`, not the project root, when invoked
that way.
