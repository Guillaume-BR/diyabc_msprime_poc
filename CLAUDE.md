# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A proof-of-concept to replace DIYABC's C++ coalescent simulator
(`particuleC.cpp::dosimulpart`) with `msprime`, scoped to the simplest
possible case: scenario 1 of the `human` dataset (4 populations, merges
only, no split/admixture). The goal is to demonstrate that a
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
   (`sample`, `merge`, `varNe`) using the vocabulary from
   `history.cpp::ScenarioC::read_events`. Unimplemented vocabulary
   (`split`, needed for scenarios 2/3/5/6 of `human`) is skipped with a
   warning rather than failing the whole parse — only `NotImplementedError`
   is swallowed here, any other exception propagates.

2. **`prior_parser.py`** — extracts `Prior` and `OrderConstraint` objects
   from the `historical parameters priors` section.
   `is_constant_prior` replicates the exact near-degenerate-bounds rule
   from `readReftable.R` / `abcranger/readreftable.cpp` (constant if
   `(max-min)/max <= 1e-6`, never constant when `max == 0.0`) — DIYABC
   excludes these from `reftable.bin` columns.

3. **`parameter_sampling.py`** — draws one value per prior via
   rejection sampling until all `OrderConstraint`s are satisfied.

4. **`demography_builder.py`** — `Scenario` + drawn values →
   `msprime.Demography`. Populations are named `"pop1".."popN"` by their
   1-indexed position in `header.txt`. `get_parameter_names_used_by_scenario`
   determines which prior names a *specific* scenario actually references
   (important: `human/header.txt` declares 21 priors globally, but
   scenario 1 only uses 16 — the other 5 belong to scenarios 2-6; this
   subset, not the full prior list, must become the `reftable.bin`
   parameter columns, per a bug found via `readReftable.R` failing with
   "indice hors limites").

5. **`observed_data.py`** — maps population index (1,2,3,4 from
   `header.txt`, which never names populations) to real population name
   in the `.snp` file (e.g. `{1: "ASW", 2: "YRI", ...}`). The mapping is
   implicit: population *i* in the scenario = the *i*-th population by
   first-appearance order in the `.snp` file — there is no cross-reference
   in the DIYABC C++ source, this was verified by exhaustive code search.
   Never replace the order-preserving dict/Counter here with anything
   that could reorder keys (e.g. alphabetical sort) — it would silently
   break this mapping.

6. **`ancestry_simulation.py`** — simulates one independent tree per SNP
   locus (`simulate_independent_loci`, no recombination/linkage between
   loci) and mutates each with the Hudson algorithm
   (`simulate_snp_genotypes`): exactly one mutation per locus, placed on
   an edge chosen with probability proportional to branch length, fully
   vectorized over tskit's edge tables (not a per-node `branch_length()`
   loop). This guarantees every locus is globally polymorphic by
   construction — see `notes/exploration.md` for the DIYABC doc citation
   (section 2.4.3) and empirical validation. MAF filtering (`<MAF=N%>`,
   as opposed to `<MAF=hudson>`) is not implemented — not needed for the
   `human` dataset.

7. **`snp_writer.py`** — writes simulated genotypes out in DIYABC `.snp`
   format (only used for the deprecated subprocess-based path, see below).

8. **`summary_statistics.py`** — pure-Python/numpy reimplementation of
   all 130 SNP summary statistics from `sumstat.cpp` (ML1-3, HW, HB,
   FST1-4, NEI, AML, F3, F4), validated column-by-column against the real
   `general` binary's output. This is the **current** path — no subprocess,
   no intermediate `.snp` file, ~7s/particle on 5000 loci (was ~343s via
   subprocess before batch-size and vectorization fixes, see below).
   Formula provenance for each function is documented via exact
   `sumstat.cpp` function names in each docstring — check there before
   changing a formula.

9. **`pipeline.py`** — orchestrates 1-8 with no logic of its own.
   `compute_summary_statistics` is the main high-level entry point (its
   docstring still mentions delegating to the C++ binary; that's stale —
   it now calls `summary_statistics.compute_all_statistics` directly).

10. **`reftable_loop.py`** — runs `nrec` independent particles in
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
    of parameter columns across all candidate scenarios. Never validated
    end-to-end against a real DIYABC multi-scenario reftable yet, though
    (see `tests/test_scenario1_human.py` for the current unit-level
    coverage of weighted draw + `SplitEvent`/admixture translation).

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

A modest residual bias (~5-16%) remains on pairwise statistics
(FST2/NEI/F3/HB) even with a correct header — open, much smaller in
magnitude, not yet root-caused.

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
