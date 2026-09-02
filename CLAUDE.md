# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A proof-of-concept to replace DIYABC's C++ coalescent simulator
(`particuleC.cpp::dosimulpart`) with `msprime`. Originally scoped to the
simplest possible case (scenario 1 of the `human` dataset: 4
populations, merges only, no split/admixture), the pipeline now covers
the condensed SNP loci description format completely: all heritage
types (`<A>/<H>/<X>/<Y>/<M>`), split/admixture events and multi-scenario
weighted draws, MAF filtering, and PoolSeq (pooled read counts rather
than per-individual genotypes) — validated against real DIYABC output
on `human`, `toy_example5`, `toy_example3` (split/admixture, scenario
3), and `toy_example4` (PoolSeq). The *detailed* loci description
format (one locus named per line, used by MicroSat/sequences-mut
datasets, as opposed to the condensed format) is now parsed too — see
"MicroSat / sequences-mut header parsing" below. As of 2026-07-31, the
DNA-sequence side is fully wired end-to-end for one locus at a time
(see "DNA sequence substitution model" and "Mutation placement"
below): observed base frequencies, substitution model choice (`JK`/
`K2P`/`HKY`/`TN`), hierarchical `k1`/`k2`/`mus_rate` draws, per-site
rate heterogeneity (`mutsit`/invariant sites), and mutation placement
itself (via `msprime.sim_mutations`, not a hand-rolled port of
`particuleC.cpp`'s tree traversal) are all implemented and tested,
end-to-end, on `toy_example2_ms_dna`
(`bridge.ancestry_simulation.dna_mutation_simulation_per_locus`).
As of 2026-08-25, all 13 DNA-sequence-specific summary statistics
(`NHA`/`NSS`/`MPD`/`VPD`/`DTA`/`PSS`/`MNS`/`VNS` per-population,
`NH2`/`NS2`/`MP2`/`MPB`/`HST` per-pair) are implemented and tested too
— see "DNA sequence summary statistics" below. As of 2026-08-26, the
full DIYABC-replay pipeline (paired real-vs-msprime validation, same
architecture as the SNP side) is wired end-to-end for DNA sequences too
(`pipeline.compute_summary_statistics_dna`/`_from_values`,
`reftable_loop.replay_reftable_simulation_dna`) and cross-validated
against a real 1000-particle DIYABC reftable on `toy_example2_ms_dna` —
see "DIYABC-replay pipeline for DNA sequences" below. A 2026-08-27
follow-up investigation (50+50-loci stress dataset, then direct
`particuleC.cpp` source reading) identified a residual, narrow
statistical gap on mitochondrial (`<M>`) loci and initially attributed
it to a combinatorial admixture effect — not a port bug either way, but
a 2026-08-31 falsification test (no-admixture control scenario) showed
that attribution was incomplete and reopened the investigation. Fixed
2026-09-02: `<M>` loci weren't sharing a genealogy — see "RESOLVED
2026-09-02: G3 (`<M>`) variance deficit" below.
MicroSat itself (stepwise mutation model, `NAL`/`HET`-style summary
statistics) has no simulation-side code at all yet, only header
parsing. The goal is to demonstrate that a
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

# Lint / format (ruff, config in pyproject.toml)
ruff check .
ruff format .

# Pre-commit hooks (ruff check --fix + ruff format)
pre-commit run --all-files
```

Regenerate `notes/tree.md`, `notes/commits.md`, `notes/api.md`, and
`notes/report.md` (a rollup of all three) with:

```bash
python3 tools/generate_report.py
```

## Code style

Docstrings are written in **Google style** (`Args:`/`Returns:`/`Raises:`
sections, one-line summary first) — converted throughout `bridge/` on
2026-08-31/2026-09-01 (see `loci_parser.py` for a reference example).
Write any new or edited docstring the same way; don't reintroduce the
old free-prose style (no `Args:`/`Returns:` headers) even for small
functions.

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
   excludes these from `reftable.bin` columns. `_extract_historical_
   priors_section` falls back to an empty line as section terminator
   when `DRAW UNTIL` is absent (`toy_example1_ms` has no order
   constraints, hence no `DRAW UNTIL` line at all — see
   `header.cpp::readHeadersimHistParam`). Also parses the separate
   `group priors` section (`parse_group_priors` → `GroupPrior`, one per
   line, used by MicroSat/sequences-mut datasets) — see "MicroSat /
   sequences-mut header parsing" below for the non-obvious hierarchical
   dependency between some of these lines. Since 2026-07-30, also
   `get_parameter_used_by_model` (`GroupPrior` with `model=True` →
   `(k1_used, k2_used)` bool tuple, from the `name_model` token —
   `JK`→`(False,False)`, `K2P`/`HKY`→`(True,False)`, `TN`→`(True,True)`)
   — see "DNA sequence substitution model" below.

3. **`stats_group_parser.py`** — extracts the list of summary-statistic
   column names actually declared in the `group summary statistics`
   section (`parse_requested_statistic_names`), used to filter
   `summary_statistics.compute_all_statistics`'s output down to what
   `header.txt` really expects (otherwise `reftable.bin` ends up with
   extra/missing columns vs. the real DIYABC's, discovered on
   `toy_example5_modif`). Handles multiple `group Gx (N)` blocks
   (`_split_stats_blocks`) — the flattened result is order-preserved,
   not deduplicated across groups (deliberately: unclear yet whether two
   groups could legitimately declare the same column name for two
   different underlying statistics — collapsing them via `set`/
   `dict.fromkeys` would silently drop one).

4. **`loci_parser.py`** — parses the `loci description` section: the
   condensed format (single- or multi-heritage-type line, e.g.
   `"5000 <A> G1 from 1"` or `"70 <A> 10 <X> 10 <M> 10 <Y> G1 from 1"`,
   → one `LociDescription` with `loci_counts_by_heritage: dict[heritage_type, count]`)
   **and**, since 2026-07-28, the detailed one-locus-per-line format used
   by MicroSat/sequences-mut datasets (e.g. `"Locus_M_A_1_ <A> [M] G1 2
   40"` for a microsatellite locus, `"Locus_S_A_11_ <A> [S] G2 100"` for
   a DNA sequence locus → `list[LociDescriptionDetailed]`, one per line).

5. **`parameter_sampling.py`** — draws one value per prior via
   rejection sampling until all `OrderConstraint`s are satisfied.
   `_draw_one_value` supports `UN`/`LU`/`NO`/`LN`/`GA` (the full set DIYABC
   uses for historical parameters), each reproducing the exact
   `PriorC::drawfromprior`/`MwcGen::g*` formulas from `history.cpp`/
   `randomgenerator.cpp` — including the truncation-by-rejection loop for
   `NO`/`LN`/`GA` (redraw until inside `[min, max]`, never an unbounded
   draw) and the `GA` reparameterization (`gammavariate(shape=sdshape,
   scale=mean/sdshape)`, not `(mean, sdshape)` directly). Since
   2026-07-28, also draws `GroupPrior` values (`_draw_one_group_value`/
   `draw_group_parameter_values`, MicroSat/sequences-mut) — see
   "MicroSat / sequences-mut header parsing" below. Since 2026-07-30,
   also `sampling_group_local_param` (the second, per-locus tier of the
   `k1`/`k2`/`mus_rate` draw hierarchy, on top of the per-group
   `k1moy`/`k2moy`/`musmoy` from `draw_group_parameter_values`) and
   `sample_site_rates` (`mutsit`, per-site rate heterogeneity) — see
   "DNA sequence substitution model" below.

6. **`demography_builder.py`** — `Scenario` + drawn values →
   `msprime.Demography`. Populations are named `"pop1".."popN"` by their
   1-indexed position in `header.txt`. `get_parameter_names_used_by_scenario`
   determines which prior names a *specific* scenario actually references
   (important: `human/header.txt` declares 21 priors globally, but
   scenario 1 only uses 16 — the other 5 belong to scenarios 2-6; this
   subset, not the full prior list, must become the `reftable.bin`
   parameter columns, per a bug found via `readReftable.R` failing with
   "indice hors limites").

7. **`observed_data.py`** — maps population index (1,2,3,4 from
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
   `with_maf_filter_shared_ancestry` in `ancestry_simulation.py`). Since
   2026-07-30, also parses DNA sequence content from `.mss` files:
   `observed_sequences` (per-locus list of haplotype strings, correctly
   splitting diploid `<A>` loci into 2 entries vs. haploid `<M>` loci
   into 1, matching `do_sequence`'s `n=1`/`n=2` split — see "DNA sequence
   substitution model" below) and `base_frequency_by_locus` (`pi_A/C/G/T`
   per locus, pooled across all individuals/populations, matching
   `DataC::do_sequence`). Since 2026-07-31, also `observed_count_population`
   (individuals per population from a `.mss` file, `.mss`'s own genepop
   block-separator format — NOT reusable via `count_samples_per_population`,
   see "DNA sequence substitution model" below).

8. **`ancestry_simulation.py`** — simulates one independent tree per SNP
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
   some `toy_example3`/`toy_example5` test datasets. Since 2026-07-31,
   also the full DNA-sequence path: `build_transition_matrix` (`matQ`),
   `count_loci_per_group`, `build_group_local_param_per_locus`
   (`k1`/`k2`/`mus_rate` per locus), `build_matrix_per_locus` (`matQ`
   per locus), `build_rate_map`/`build_rate_map_per_locus` (`mutsit` →
   `msprime.RateMap`), and `simulate_dna_mutations`/
   `dna_mutation_simulation_per_locus` (the actual mutation placement,
   per locus, via `msprime.sim_mutations` — see "DNA sequence
   substitution model" and "Mutation placement" below). Unlike the SNP
   path, DNA sequence loci call `msprime.sim_ancestry` directly rather
   than `simulate_independent_loci` (which hardcodes `sequence_length=1`,
   wrong once `dnalength` varies per locus).

9. **`snp_writer.py`** — writes simulated genotypes out in DIYABC `.snp`
   format (only used for the deprecated subprocess-based path, see below).

10. **`summary_statistics.py`** — pure-Python/numpy reimplementation of
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

11. **`pipeline.py`** — orchestrates 1-10 with no logic of its own.
   `compute_summary_statistics` is the main high-level entry point (its
   docstring still mentions delegating to the C++ binary; that's stale —
   it now calls `summary_statistics.compute_all_statistics` directly).

12. **`reftable_loop.py`** — runs `nrec` independent particles in
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

### MicroSat / sequences-mut header parsing (started 2026-07-28 — header parsing only, simulation not started)

New reference datasets `reference/toy_example1_ms/` (microsat, no order
constraints) and `reference/toy_example2_ms_dna/` (microsat + DNA
sequences mixed) are used to extend header parsing to the *detailed*
loci format, not yet touching the simulation side. Don't assume
"MicroSat is supported" beyond parsing — `demography_builder.py`,
`ancestry_simulation.py`, and `summary_statistics.py` are all still
SNP-only.

- **`loci_parser.py`**: `parse_loci_description` now also handles the
  detailed one-locus-per-line format (`LociDescriptionDetailed`, one per
  line, fields `motif_size`/`motif_range` for `[M]` or `dnalength` for
  `[S]`). `dnalength` (the trailing number on `[S]` lines) is captured
  but is informational only — confirmed in `header.cpp:401-402` that
  DIYABC itself never reads it at this point (`dnalength` is already set
  from elsewhere earlier in the header; the comment there literally says
  `//inutile variable déjà renseignée`).

- **`prior_parser.py`**: `_extract_historical_priors_section` now falls
  back to an empty line as terminator when `DRAW UNTIL` is absent
  (`toy_example1_ms` has 0 order constraints, hence no `DRAW UNTIL` at
  all). New `parse_group_priors` → `GroupPrior` parses the separate
  `group priors` section (mutation-model hyperpriors, e.g. `MEANMU
  UN[1e-4,1e-3,5e-4,2]` / `GAMMU GA[1e-5,1e-2,Mean_u,2]` for `[M]`,
  `MODEL K2P 10 2.00` for `[S]`). **Non-obvious hierarchical dependency**:
  a `GAMxxx` line's declared "mean" (`Mean_u`, `Mean_P`, `Mean_u_SNI`...)
  is not a real bound — `readprior` (`history.cpp:89`) parses it via
  `atof` (silently 0.0 on non-numeric text), and it is unconditionally
  overwritten at draw time (`particuleC.cpp:819`,
  `priormutloc.mean = grouplist[gr].mutmoy`) with the value just drawn
  from the group's own `MEANxxx` prior. **Confirmed positional, not
  name-based**: `header.cpp::readHeadersimGroupPrior` reads a fixed
  sequence of `getline` calls per group and never inspects the name
  token (`ss1[0]`) at all — the text `MEANMU`/`GAMMU` is purely cosmetic
  in the file. `parse_group_priors`'s `else` branch (anything that isn't
  a `group` line or a recognized prior line is assumed to be a `MODEL`
  line) has no positive validation that it's really `MODEL ...` — a
  genuinely malformed line would be silently mis-parsed rather than
  raising a clear error (known gap, not yet fixed).

- **`stats_group_parser.py`**: `parse_requested_statistic_names` now
  handles multiple `group Gx (N)` blocks via `_split_stats_blocks`
  (shared block-splitting strategy, adapted for `group priors`' `group
  Gx [M]`/`[S]` header line, which has no count in parentheses).

- **`parameter_sampling.py`**: `_draw_one_value` gained `NO`/`LN`/`GA`
  (previously `UN`-only) — see item 5 above. `_draw_one_group_value`/
  `draw_group_parameter_values` draw `GroupPrior` values, reusing
  `_draw_one_value` by wrapping a `GroupPrior` into a throwaway `Prior`
  with a fictional `category="G"` (never matches `("N", "T")`, so no
  spurious rounding). `draw_group_parameter_values` resolves the
  `MEANxxx`→`GAMxxx` hierarchy positionally (previous drawn value in the
  group's list, matching DIYABC's own positional reading — see above),
  and offsets its seed by `_GROUP_PRIOR_SEED_OFFSET` (10,000,000) to
  never correlate with `draw_parameter_values` when both are called with
  the same base seed — same bug class as the scenario/prior correlation
  fixed 2026-07-16.

**Not yet done**: `rewrite_loci_count` still condensed-single-type only;
dedup-across-groups in `stats_group_parser.py` deliberately deferred (no
evidence yet whether two groups can legitimately share a column name);
DNA sequence *substitution model construction* (base frequencies,
model choice, `k1`/`k2` draw, `matQ`) is now done — see "DNA sequence
substitution model" below — but DNA sequence *mutation placement* along
the tree, and the entire MicroSat simulation side (demography rescaling
semantics, stepwise mutation model, MicroSat-specific summary
statistics catalog e.g. `NAL`/`HET` from `statdefs.cpp`), have not
started at all.

### DNA sequence substitution model (started 2026-07-29, mutation placement done 2026-07-31)

Picked up right after MicroSat/sequences-mut header parsing was declared
complete (chosen over MicroSat itself as "simpler to tackle first").
Covers the full path from `header.txt` + `.mss` to an actual mutated
`tskit.TreeSequence` per DNA sequence locus, validated end-to-end on
`toy_example2_ms_dna` (10 DNA sequence loci, model `K2P`): substitution
matrix construction, per-site rate heterogeneity, and — as of
2026-07-31 — mutation placement itself, achieved by handing the whole
problem to `msprime.sim_mutations` rather than reimplementing
`particuleC.cpp`'s tree traversal (see "Mutation placement" below for
why this is a faithful substitute, not an approximation).

`GroupPrior.model_bounds` (the single tuple field) was later split into
two separate named fields, `p_fixe: float | None` and `gams: float |
None` — same values, clearer call sites (`gp.p_fixe`/`gp.gams` instead
of `model_bounds[0]`/`model_bounds[1]`).

- **Observed base frequencies (`pi_A/C/G/T`) come from the `.mss` file,
  not from `header.txt`.** Verified against `~/Documents/Github/diyabc/
  src-JMC-C++`: `DataC::do_sequence` (`data.cpp:1444-1568`) computes them
  empirically per locus, pooling ALL populations/individuals together
  (never per-population) and excluding missing (`-`/`N`) positions from
  both numerator and denominator. Every simulated particle then reuses
  this same fixed value unchanged (`particleset.cpp:94`). A different
  header field (`readHeadersimLoci`, `header.cpp:1585`) does carry
  literal `pi` values in the header text, but that's for a
  no-observed-data launch mode this project doesn't use — irrelevant
  here. `observed_data.py`'s `observed_sequences`/`base_frequency_by_locus`
  reproduce `do_sequence`'s logic directly on the `.mss` file.

- **`p_fixe`/`gams` (now direct `GroupPrior` fields, see above) are
  `(proportion of invariant sites, gamma shape for per-site rate
  heterogeneity)`** — completely unrelated to `k1`/`k2`. `k1`/`k2` come
  from **separate** `GroupPrior` entries, `MEANK1`/`GAMK1`/`MEANK2`/
  `GAMK2`, present for every `[S]` group regardless of which model is
  chosen (even `K2P`, which only uses `k1`, still declares `GAMK2` —
  `header.cpp:640-670` reads a fixed sequence of lines per group, never
  conditional on the model). A first draft of `build_transition_matrix`
  mixed these two up (used the old `model_bounds` tuple as if it were
  `(k1, k2)`) — caught before commit.

- **`k1`/`k2`/`mus_rate` are all drawn hierarchically, two tiers each,
  mirroring MicroSat's `mutmoy`/`mutloc`**: `particuleC.cpp:840-867`
  draws `musmoy`/`k1moy`/`k2moy` once per particle per group (from
  `MEANMU`/`MEANK1`/`MEANK2`), overwrites `GAMMU`/`GAMK1`/`GAMK2`'s
  declared mean with that value, then EITHER draws an independent
  per-locus value (if the `GAMxxx`'s `sdshape > 0.001` AND the group's
  `nloc > 1`) OR reuses the group mean unchanged for every locus. `mus_rate`
  and `k1` both have the `nloc > 1` check; `k2` genuinely doesn't (an
  asymmetry in DIYABC's own source, not a bug in this port — replicated
  as-is via `check_nloc`). The per-group tier reuses `parameter_sampling.
  draw_group_parameter_values` unchanged — its existing positional
  `MEANxxx`→`GAMxxx` mean-substitution mechanism (already built for
  MicroSat) generalizes for free to any number of `MEANxxx`/`GAMxxx`
  pairs in the same group, no code changes needed. The per-locus tier
  is `parameter_sampling.sampling_group_local_param` (renamed from
  `sampling_kappa_per_locus` once it became clear `mus_rate` needed the
  exact same mechanism — nothing kappa-specific in its implementation,
  it just takes a `GroupPrior`/`k_moy`/`n_loci`/`check_nloc`/loci
  list/`rng`).

- **`get_parameter_used_by_model`** (`prior_parser.py`) maps a group's
  `name_model` to which of `k1`/`k2` are actually active
  (`JK`→neither, `K2P`/`HKY`→`k1` only, `TN`→both) — mirrors
  `header.cpp:680-694`'s `nparamvar` bookkeeping (which of `k1`/`k2`
  become `reftable.bin` columns). **DIYABC's own header token for
  Jukes-Cantor is `"JK"`, not `"JC"`/`"JC69"`** (the common
  phylogenetics-literature name) — got this wrong twice while writing
  this function and `build_transition_matrix` before catching it; see
  [[feedback_control_flow_chaining_bugs]] memory for the pattern.

- **`build_transition_matrix`** (`ancestry_simulation.py`) builds the
  4×4 matrix per `comp_matQ`'s four branches (`JK`/`K2P`/`HKY`/`TN`),
  base order `A/C/G/T` = index `0/1/2/3` (matching `data.cpp`'s `char
  base[] = "ACGT"`). Row-sum normalization needs `axis=1,
  keepdims=True` — omitting `keepdims` silently divides along the wrong
  numpy axis (a `(4,4)` array divided by a bare `(4,)` row-sum vector
  broadcasts across columns, not rows) and produces a matrix whose rows
  don't actually sum to 1, without raising any error.

- **`build_group_local_param_per_locus`** (`ancestry_simulation.py`,
  renamed from `build_kappas_per_locus` once it grew to also return
  `mus_rate`) orchestrates, for every `[S]` group in a `header.txt`:
  `draw_group_parameter_values` called ONCE for the whole `group_priors`
  dict with the plain particle `seed` (it applies its own internal
  offset, `parameter_sampling._GROUP_PRIOR_SEED_OFFSET`) to get
  `k1moy`/`k2moy`/`musmoy`, then `sampling_group_local_param` per
  parameter with its own `random.Random(seed + _KAPPA1_SEED_OFFSET)` /
  `_KAPPA2_SEED_OFFSET` / `_MUS_RATE_SEED_OFFSET` (offsets defined in
  `ancestry_simulation.py`, mirroring `parameter_sampling.py`'s — kept
  deliberately separate from `draw_group_parameter_values`'s own seeding
  so none of these draws correlate). Returns `{locus_name: (k1, k2,
  mus_rate)}` — a fixed-arity triplet for every locus regardless of
  which kappas the group's model actually uses (`0.0` filled in for the
  unused ones). Two bugs caught before commit: an early draft called
  `draw_group_parameter_values` per-parameter with an already-offset
  seed instead of once with the plain seed (wrong argument shape *and*
  wrong seed semantics), and `mus_rate` was only included in the
  triplet inside the `TN` branch — silently absent from the return
  value for any group using `K2P`/`HKY`/`JK` (i.e. absent for the
  entire `toy_example2_ms_dna` dataset, since it's `K2P`-only).

- **`build_matrix_per_locus`** (`ancestry_simulation.py`) is the full
  `header.txt` + `.mss` + `seed` → `{locus_name: matQ}` pipeline, tying
  together `build_group_local_param_per_locus`, `base_frequency_by_locus`,
  and `build_transition_matrix`. Validated end-to-end on
  `toy_example2_ms_dna` (10/10 sequence loci, every row-stochastic,
  reproducible across repeated calls with the same seed).

- **`parameter_sampling.sample_site_rates`** draws `mutsit` (per-site
  relative mutation rate, matching `header.cpp:707-738`): a
  `Gamma(shape=gams, mean=1)` draw per site, then the first `dnalength -
  nsv` sites forced to `0.0` (invariant), then the whole array normalized
  to sum to 1. **`gams == 0` is a valid, non-error input**: DIYABC's own
  `MwcGen::ggamma3` (`randomgenerator.cpp:199-204`) returns `mean` (`1.0`
  here) directly when `shape == 0.0` rather than dividing by it — an
  early draft raised `ZeroDivisionError` instead, which is wrong
  behavior, not just an unhandled edge case.

  **Non-obvious DIYABC bug, reproduced deliberately, not "fixed"**:
  `header.cpp:727-738` draws `dnalength - nsv` distinct random site
  indices into a `sitefix` array (with a proper duplicate-rejection
  loop) specifically to pick which sites are invariant — but the line
  that actually zeroes out `mutsit` uses the loop counter `i`, never
  `sitefix[i]`. The carefully-drawn random indices are computed and then
  never read. In practice, DIYABC's "invariant sites" are always the
  first `dnalength - nsv` sites in sequence order, not a random subset,
  despite the code's evident intent. Given this project's goal (bit-for-
  bit statistical fidelity to real DIYABC output, not an idealized
  reimplementation), `sample_site_rates` reproduces this exact behavior
  — the "intended" random-selection version is kept commented out in
  the function body, with a note that it's reserved in case a future
  decision (with the user's academic advisor) is to fix rather than
  replicate this bug.

- **`observed_data.observed_count_population`** counts individuals per
  population in a `.mss` file — the DNA-sequence/MicroSat equivalent of
  `count_samples_per_population`, but NOT a variant of it: `.mss` is
  genepop-format (`POP` as a block *separator* between populations, see
  `observed_sequences`), structurally different from `.snp`'s "IND SEX
  POP"/"POOL" column format, and the caller always already knows which
  file format it has (no need for a `.snp`-vs-`.mss` content-sniffing
  dispatcher analogous to `detect_snp_file_type`, which exists only
  because `.snp` itself can be either IND or POOL). Returns `{"pop1":
  N1, "pop2": N2, ...}`, same population-index-by-first-appearance
  convention and same msprime-facing naming as `build_samples_argument`.

- **`build_rate_map`/`build_rate_map_per_locus`** (`ancestry_simulation.py`)
  turn `mutsit` (relative, sums to 1) into an absolute-rate
  `msprime.RateMap` per locus: `rate[site] = mus_rate × dnalength ×
  mutsit[site]`, one breakpoint per site (`position=[0..dnalength]`).
  Validated end-to-end (`build_rate_map_per_locus`) against real
  `p_fixe`/`gams`/`mus_rate` on `toy_example2_ms_dna`: correct count of
  zero-rate (invariant) sites, distinct rate patterns across different
  loci (confirming independent per-locus `mutsit` draws, not the same
  pattern replayed), reproducible with the same seed.

### Mutation placement (2026-07-31) — via msprime.sim_mutations, not a hand-rolled port of particuleC.cpp

`particuleC.cpp`'s own mechanism (`put_mutations`, `init_dnaseq`,
`mute`, `draw_nuc` — `particuleC.cpp:1535-1775`) is: for each locus,
draw a Poisson number of mutations per branch (rate ∝ branch length ×
`mus_rate × dnalength`), assign each mutation to a site via a cumulative
walk over `mutsit`, and apply it via a cumulative walk over `matQ`'s row
for the site's current base — root sequence drawn from `pi_A/C/G/T` via
`draw_nuc`. This is textbook continuous-time-Markov-chain sequence
evolution along a tree, nothing DIYABC-specific, and `msprime` already
implements exactly this via `msprime.MatrixMutationModel(alleles,
root_distribution, transition_matrix)` + `msprime.sim_mutations(ts,
rate=..., model=...)` — confirmed by direct exploration (throwaway
scripts, not committed): observed mutation count matches the Poisson
expectation, ancestral-state distribution matches `pi`, zero silent
(self-transition) mutations given `matQ`'s zero diagonal, and the
observed transition/transversion ratio for a `K2P` matrix with `k1=8`
matches the theoretical `k1/2` (two transversion targets per row vs.
one transition target) almost exactly (4.09 observed vs. 4.0
theoretical). `msprime.RateMap` (passed as `rate=`) correctly implements
per-site heterogeneity: zero mutations in a zero-rate band, correct
ratio between two differently-rated bands. **Deliberately did NOT use
msprime's built-in named models** (`JC69`/`HKY`/`F84`/`GTR`) — their
internal rate-scaling convention allows/counts silent self-transitions
(normalizes by `max(row sum)`, non-zero diagonal), different from
`comp_matQ`'s convention (normalizes each row to sum to 1, zero
diagonal, every event is a real substitution) — using them would have
given a different effective substitution rate than DIYABC's for the
same nominal `mutation_rate`. The generic `MatrixMutationModel` fed our
own already-built `matQ`/`pi` sidesteps this entirely.

- **`simulate_dna_mutations`** (`ancestry_simulation.py`) is the narrow
  wrapper: `tree_sequence` + `pi` + `matQ` + `rate_map` + `seed` →
  mutated `tskit.TreeSequence`. Bugs caught before commit: `root_distribution`
  omitted entirely from an early draft (raises `TypeError` — msprime
  requires it), `random_seed` passed to `MatrixMutationModel` (which
  doesn't accept one at all — it's a static model description, not
  something that draws anything) instead of to `sim_mutations` (which
  silently got an auto-generated, non-reproducible seed as a result),
  and a `route_distribution` typo, plus passing the raw `pi` dict where
  an ordered `[pi_A, pi_C, pi_G, pi_T]` list matching `alleles` order
  was required.

- **`dna_mutation_simulation_per_locus`** (`ancestry_simulation.py`) is
  the full per-locus assembly: for every `[S]` locus, calls
  `msprime.sim_ancestry` **directly** (NOT `simulate_independent_loci`,
  which hardcodes `sequence_length=1` for the SNP/Hudson case and can't
  represent a locus-specific `dnalength`), with `sequence_length=
  locus.dnalength`, then `simulate_dna_mutations` with that locus's
  `matQ`/`pi`/`RateMap`. Samples come from `observed_count_population`
  on the `.mss` file. Needs its own per-locus seed offsets for BOTH the
  genealogy (`_ANCESTRY_SEED_OFFSET`) and the mutation placement
  (`_MUTATION_SEED_OFFSET`) — an early draft reused the identical
  `seed + offset` for every locus in the loop (no `+ i`), which produced
  byte-identical tree topologies across all loci of a dataset sharing
  the same `dnalength` (confirmed empirically: `ts1.tables.edges ==
  ts2.tables.edges` for two different loci) — silently violating the
  independent-loci-per-particle assumption this whole pipeline relies
  on for summary statistics. Validated end-to-end on
  `toy_example2_ms_dna`: 10/10 sequence loci (not the 10 MicroSat loci
  in the same header), independent topologies and mutation patterns
  across loci, reproducible with the same particle seed.

  **Ploidy/demography bug fixed 2026-08-24**: until this date,
  `dna_mutation_simulation_per_locus` called `msprime.sim_ancestry` with
  the raw `<A>` `demography` and `ploidy=2` for EVERY `[S]` locus,
  regardless of that locus's own heritage type — so `toy_example2_ms_dna`'s
  G3 group (`<M>`, mitochondrial, meant to be haploid with a rescaled
  demography) was actually being simulated as if it were `<A>` (diploid,
  unrescaled). The SNP path already gets this right, per-locus-type, in
  `simulate_genotypes_for_locus_type`. New helper
  `dna_ancestry_parameters_for_heritage` (`ancestry_simulation.py`)
  replicates that same dispatch for DNA sequences: `"A"` → demography
  unchanged, `ploidy=2`; `"H"`/`"M"` → `rescale_demography(demography,
  coalescence_coefficient(heritage, sex_ratio) / 2)`, `ploidy=1`; `"X"`/`"Y"`
  → `NotImplementedError` (deliberately, not deferred-by-oversight: `.mss`
  is genepop-format and carries no per-individual sex column the way
  `.snp` does, so `build_sex_stratified_samples_argument`/
  `build_male_only_samples_argument` have no equivalent to call here —
  sex-stratified DNA sequence loci would need a real per-individual sex
  source that doesn't exist in this file format). `sex_ratio` is read via
  the existing `parse_sex_ratio(mss_file_path)` — works unchanged on
  `.mss` because the `<NM=xNF>` token it looks for lives on the file's
  first line in exactly the same format as `.snp`, confirmed by direct
  inspection of `toy_example2_ms_dna`'s `.mss` file. `dna_mutation_
  simulation_per_locus` now calls this dispatch once per locus (a single
  group can mix heritage types across sequence loci, e.g. G2=`<A>`/
  G3=`<M>` in the same header, so the dispatch can never be hoisted
  outside the per-locus loop). Regression tests: `test_dna_ancestry_
  parameters_for_heritage` (dispatch itself) and `test_dna_mutation_
  simulation_per_locus_ploidy_matches_heritage` (`ts.num_samples` for an
  `<A>` locus is exactly 2x an `<M>` locus's, for the same population —
  this would have passed silently before the fix, since both were
  ploidy=2, hence checked directly against the post-fix dispatch, not
  just re-run of the pre-existing test).

**Not yet done**: no real reference dataset with a `JK`- or `TN`-model
DNA sequence group to cross-validate `build_transition_matrix` against
real DIYABC output (`toy_example2_ms_dna` only exercises `K2P`) — those
two branches are covered by hand-computed synthetic test values only.
MicroSat itself still has no simulation-side code at all, only header
parsing.

### DNA sequence summary statistics (started 2026-08-24, mentor mode — user-driven, reviewed/debugged with the assistant)

Picked up right after the ploidy/demography fix above. Goal: reproduce
the 13 DNA-sequence-specific statistics from `sumstat.cpp` (`cal_nha1p`/
`2p`, `cal_nss1p`/`2p`, `cal_mpd1p`, `cal_vpd1p`, `cal_mpw2p`, `cal_mpb2p`,
`cal_dta1p`, `cal_pss1p`, `cal_mns1p`, `cal_vns1p`, `cal_fst2p` —
confirmed against `toy_example2_ms_dna/headerRF.txt`'s `group summary
statistics` section, which requests exactly these 13 under the names
`NHA`/`NSS`/`MPD`/`VPD`/`DTA`/`PSS`/`MNS`/`VNS` per-population and
`NH2`/`NS2`/`MP2`/`MPB`/`HST` per-pair) in `bridge/summary_statistics.py`,
building on the mutated `TreeSequence`s from `dna_mutation_simulation_
per_locus`. Distinct from the MicroSat-specific `NAL`/`HET`/`VAR`/`MGW`/
`FST`/`LIK`/`DAS`/`DM2` stats declared in the same header's `G1` group —
those need allele-size arithmetic, not tskit genotypes, and are out of
scope here.

- **`compute_population_layout`** (`ancestry_simulation.py`, renamed
  from the private `_population_layout` since `summary_statistics.py`
  now needs it too) hit a real name-shadowing bug when made public: four
  call sites inside `simulate_snp_genotypes`/`with_maf_filter_shared_
  ancestry` (which both also have a **parameter** named `population_layout`)
  did `population_layout = population_layout(ts)` — the local parameter
  shadowed the module-level function, so this became `None(ts)` whenever
  the parameter defaulted to `None`, breaking 24 tests across three test
  files. Fixed by renaming the function only (not the parameter, which is
  documented at length in multiple docstrings) — see
  `feedback_name_shadowing_pattern` project memory, same bug class
  recurring.

- **`_genotype_matrix_by_population`** (`summary_statistics.py`): one
  `TreeSequence` (one DNA sequence locus) → `{pop_name: matrix}`, matrix
  shape `(n_sites, n_samples_pop)` — tskit's native `genotype_matrix()`
  convention, sliced per population via `compute_population_layout`.
  `genotype_matrix()` called once per `TreeSequence`, not once per
  sample (an early draft rebuilt the whole matrix per sample). Tested on
  both an `<A>` and an `<M>` locus of `toy_example2_ms_dna`: correct
  `n_sites`/`n_samples` shapes, no sample lost/duplicated across
  populations, and the 2:1 sample-count ratio between `<A>`/`<M>`
  confirms it composes correctly with the 2026-08-24 ploidy fix above.

- **Two-tier pattern established for the per-population ("1p") stats**,
  mirroring `sumstat.cpp`'s own per-locus/per-group split (`cal_*pl` +
  `cal_*1p` with an `nl` denominator): a `_count_*(matrix)` brick
  (one locus, one population → a scalar) plus a `mean_*_per_group
  (tree_sequences, population_names)` aggregator (mean over a **single**
  header `group Gx`'s loci — never loci from two different groups mixed
  together, since each group computes its own independent stat).
  Aggregators pre-fill `{pop_name: 0.0 for pop_name in population_names}`
  before accumulating, matching the C++'s `res = 0.0` declared before its
  `if (nl > 0)` guard — so every expected population always has a value,
  even for an empty `tree_sequences` list, rather than a population
  silently missing from the result dict. Documented, not fixed: this
  assumes every population in `population_names` is present on every
  locus of the group (divides by `len(tree_sequences)`, not a real
  per-population `nl` count) — true on `toy_example2_ms_dna` (checked
  empirically), not guaranteed in general; violating it raises `KeyError`
  rather than silently excluding that locus, unlike the C++.
  - `NSS` (`_count_segregating_sites` + `mean_segregating_sites_per_group`,
    `cal_nsspl`/`cal_nss1p`): a site is segregating for a population if
    not all its samples share the same base — vectorized as
    `np.any(matrix != matrix[:, [0]], axis=1)`, correct for any number of
    distinct values per site since "differs from sample 0" is equivalent
    to "not all identical" (not just "exactly 2 alleles").
  - `NHA` (`_count_distinct_haplotypes` + `mean_distinct_haplotypes_per_group`,
    `cal_nha1p`): number of distinct haplotypes = distinct **columns** of
    the matrix (`np.unique(matrix, axis=1)`). Bug caught and fixed: an
    early draft returned `.shape[0]` (number of *sites*) instead of
    `.shape[1]` (number of distinct haplotypes) — both happened to be
    `3` on the first hand-picked test matrix, masking the bug until a
    second matrix with `n_sites != n_distinct_haplotypes` was tried.
    `np.unique` on a `(0, n_samples)` matrix (locus with zero variable
    sites) correctly returns exactly 1 unique column for free, matching
    `cal_nha1p`'s explicit `dnavar == 0` → 1-haplotype special case
    without needing an explicit branch.
  - Both `_count_segregating_sites`/`_count_distinct_haplotypes` raise
    `ValueError` on `matrix.shape[1] == 0` (population with zero samples
    on a locus) rather than crashing obscurely or returning a silently
    wrong count — deliberately checked against `shape[1]` (samples), not
    `matrix.size` (an earlier draft used `.size`, which also triggers
    incorrectly on the *valid* `n_sites == 0` case, a locus with no
    mutations at all — see `feedback_control_flow_chaining_bugs`-style
    edge-case conflation).
  - A `ValueError` was independently added and removed **twice** from
    `mean_segregating_sites_per_group`'s empty-list handling during this
    session — once genuinely misplaced (inside `if num_loci > 0` instead
    of `else`, so it fired on the *normal* case), once syntactically
    correct but reintroduced the very "raise instead of 0.0-default"
    design this whole two-tier pattern was built to avoid. Kept as
    `0.0`-default, confirmed explicitly with the user both times — if
    this `raise` reappears a third time, check for an editor/autosave
    restoring a stale buffer rather than assuming it's an intentional
    edit.

- **`_pairwise_hamming_distances`** (`summary_statistics.py`, brick for
  `MPD`/`VPD`/`cal_mpdpl`/`cal_vpd1p`): one matrix → the 1D vector of
  Hamming distances for all `C(n_samples, 2)` pairs, via
  `(matrix[:, :, None] != matrix[:, None, :]).sum(axis=0)` then
  `np.triu_indices(..., k=1)` to keep `i < j` pairs only (no double-count,
  no diagonal). Verified against a hand-computed matrix. Same
  `matrix.shape[1] == 0` → `ValueError` guard added for consistency with
  `_count_segregating_sites`/`_count_distinct_haplotypes` (an initial
  draft silently returned `[]` instead, which would have propagated into
  `nan` + a numpy `RuntimeWarning` from `mean()`/`var()` on an empty
  array rather than a clear error at the source).

- **`MPD`/`VPD`** (`mean_pairwise_differences_per_group`/`variance_
  pairwise_differences_per_group`, `cal_mpd1p`/`cal_vpd1p`) needed a
  different exclusion regime than `NSS`/`NHA`: instead of a flat
  `num_loci` denominator, a **per-population** `valid_loci_count` dict,
  because a locus only contributes if it has at least 1 pair (`MPD`,
  `nd > 0`) or at least 2 pairs (`VPD`, `nd > 1`) — `_pairwise_hamming_
  distances` on a 1-sample matrix returns an empty vector, and
  `.mean()`/`.var()` on that gives `nan`, which would otherwise silently
  poison the whole group's sum. Not a bug in practice on this project's
  datasets (20-40 samples/population always) but the C++'s own guard
  reproduced for fidelity. A `> 1` vs `> 0` threshold confusion on
  `VPD`'s **final division** guard was flagged as a possible bug and
  turned out to be a false alarm — dividing by 1 is a no-op, so the two
  thresholds are numerically indistinguishable in every case; still
  switched to `> 0` for readability/consistency with the rest of the
  file, not because the `> 1` version was wrong.

- **`DTA`** (`_tajima_d_per_locus` + `mean_tajima_d_per_group`,
  `cal_dta1pl`/`cal_dta1p`) is the classic Tajima's D neutrality
  statistic, built directly on top of `MPD` (π) and `NSS` (S) — no new
  per-site logic needed beyond `_tajima_constants(n_samples)` (the
  `a1`/`e1`/`e2` coefficients, pure functions of sample size). Caught
  before commit: a parenthesization bug, `(n+1) / ((n-1)/3.0)` instead
  of `(n+1)/(n-1)/3.0`, made `b1` exactly 9× too large (verified
  numerically). **Two distinct, easy-to-conflate exclusion cases**:
  `n_samples < 2` excludes the locus entirely (`_tajima_d_per_locus`
  returns `None`, matching `OKK = false`); `n_samples >= 2` but `S == 0`
  (no segregating sites → the formula's denominator is 0) still
  **includes** the locus in the group average with a contributed value
  of `0.0` — the C++ never resets `OKK` in that second case (`cal_
  dta1pl` lines 1579/1594-1598). An initial draft used the wrong input
  shape entirely (SNP-style `genotypes_per_locus: list[dict]`, treating
  "number of loci" as "number of samples") before being rewritten to
  match the `tree_sequences`-based two-tier pattern of every other stat
  here.

- **`PSS`** (`_private_segregating_sites_per_locus` + `mean_private_
  segregating_sites_per_group`, `cal_pss1p`) is the one per-population
  stat that needs **every** population's matrix at once, not just the
  target's — a site counts as "private" to population `i` only if it's
  segregating in `i` and fixed in **every other population of the whole
  dataset** (not just the target's group). Required factoring a reusable
  `_segregating_sites_mask(matrix)` boolean helper out of `_count_
  segregating_sites` (previously computed the count directly). Since all
  populations' matrices for one locus come from the same underlying
  `genotype_matrix()` (just column-sliced), row `i` means the same
  physical site for every population — booleans compare directly by
  position, no index-matching search needed (the C++ does need one,
  `ssa[sample][j] == ssa[sa][k]`, because its per-population variable-
  site index lists are separate dynamic arrays). `nl` increments
  unconditionally every locus in `cal_pss1p` (no `samplesize > 0` guard
  at all, unlike `NSS`/`NHA`) — simple `num_loci` denominator.

- **`MNS`/`VNS`** (`_minor_allele_counts_at_segregating_sites` +
  `mean_minor_allele_count_per_group`/`variance_minor_allele_count_per_
  group`, `afs`/`cal_mns1p`/`cal_vns1p`): at each segregating site,
  `min(counts of each distinct base actually present)` — reproduces the
  C++'s "sort 4 slots ascending, skip zeros" (`afs`) via `np.unique(site,
  return_counts=True).min()` when `len(values) > 1`, simpler because
  `np.unique` only ever returns actually-present values (no need to
  handle the zero-slots explicitly). **`VNS` is a BIASED variance
  (`ddof=0`, division by `n` not `n-1`)** — confirmed against `cal_
  vns1p`'s `v = (sx2 - sx*sx/a) / a`, deliberately different from `VPD`'s
  `ddof=1`, easy to get wrong by pattern-matching against `VPD`. Both
  stats use the flat `num_loci` denominator (`nl` increments
  unconditionally in both `cal_mns1p`/`cal_vns1p`, like `PSS`) — no
  per-population exclusion needed.

- **Pairwise ("2p") stats reuse the 1p bricks on pooled/cross
  matrices**, all following the same aggregator skeleton (`{pair_key:
  0.0}` pre-filled, `num_loci` denominator, `"{i+1}.{j+1}"` keys via a
  plain double loop — not `_half_arrangements`, see below):
  - **`NH2`** (`mean_distinct_haplotypes_per_group_pairwize`, `cal_
    nha2p`): `_count_distinct_haplotypes` on `np.hstack(matrix_a,
    matrix_b)` — the two populations' matrices are column-slices of the
    same `genotype_matrix()`, so concatenation along the sample axis
    needs no realignment.
  - **`NS2`** (`mean_segregating_sites_per_group_pairwize`, `cal_
    nss2p`): identical trick with `_count_segregating_sites` instead.
  - **`MP2`** ("mean pairwise **within**",
    `mean_pairwise_differences_per_group_pairwize`, `cal_mpw2p`):
    `_pairwise_hamming_distances` computed **separately** on each
    population's own matrix (never concatenated), then pooled as a
    ratio of sums (`(sum_di_a + sum_di_b) / (nd_a + nd_b)`) — NOT a
    simple average of the two populations' own `MPD` values; only
    equal to that average when both populations have the same sample
    size (as they happen to in `toy_example2_ms_dna`).
  - **`MPB`** ("mean pairwise **between**",
    `mean_pairwise_differences_between_per_group_pairwize`, `cal_
    mpb2p`): new brick `_pairwise_hamming_distances_between(matrix_a,
    matrix_b)` — full cross-product `(matrix_a[:,:,None] !=
    matrix_b[:,None,:]).sum(axis=0)`, shape `(n_a, n_b)`, **no
    triangle extraction needed** (unlike the "within" case) since every
    `(p in a, q in b)` pair is valid, never a self-comparison. An early
    draft's docstring claimed this returned a flattened 1D vector of
    length `n_a*n_b`; it actually returns the 2D `(n_a, n_b)` matrix —
    `.mean()` on it is still numerically correct either way (numpy
    averages all elements regardless of shape), but a downstream `len(...)
    > 0` guard was checking the wrong thing (`n_a`, not `n_a*n_b`) —
    harmless in practice only because the brick's own guard already
    rejects 0-sample inputs before that check is ever reached.
  - **`HST`** (`mean_hst_per_group_pairwize`, `cal_fst2p` — note the
    lowercase C++ name, easy to confuse with MicroSat's unrelated
    `cal_Fst2p`): `(Hb - Hw) / Hb` where `Hb` = `MPB`-per-locus, `Hw` =
    `MP2`-per-locus. **Breaks the "mean over loci" pattern used by every
    other DNA stat** — it's a *ratio of sums* accumulated separately per
    pair across the whole group (`num[pair]`, `den[pair]`, divided once
    at the very end), the same style as `_fst_wc` (SNP side, already in
    this file), not `sum_of_per_locus_values / num_loci`. Took two
    attempts to get right: a per-pair aggregator needs `num`/`den` as
    **dicts keyed by pair**, not shared scalars reset once per locus —
    with a shared scalar, N>2 populations silently mix different pairs'
    contributions and only the last pair visited by the inner loop ever
    gets its dict entry updated (every other pair stays stuck at its
    `0.0` default). Completely invisible on this project's only real
    DNA-sequence dataset (2 populations = 1 pair, so "the shared scalar"
    and "the only pair" are the same thing by coincidence) — caught only
    via a synthetic 3-population mock test
    (`unittest.mock.patch` on `_genotype_matrix_by_population`). See
    `feedback_pairwise_accumulator_bug` project memory; the same bug
    class could in principle recur in `NH2`/`NS2`/`MP2`/`MPB` if ever
    exercised on a 3+ population dataset, even though those four
    happened to be written correctly.
  - `_half_arrangements` (built for `AML`/`F3`/`F4`'s asymmetric HALF
    ordering, where element order matters) was briefly misapplied to
    generate `NH2`'s plain symmetric pairs — not numerically wrong for
    `r=2` (its HALF filter happens to keep only the ascending-index
    permutation), but semantically confusing and needlessly expensive;
    switched to the plain double loop `compute_HW_HB`/`compute_FST2`
    already use elsewhere in this file.

**Resolved 2026-08-25**: the column-collision question below WAS a real
issue and IS handled — see `compute_all_statistics_dna` in the same
file, which embeds the group index in every column name
(`stats_group_parser.parse_requested_statistic_names` was fixed
alongside it), exactly mirroring how the real DIYABC reftable itself
names these columns (`NSS_2_1` for group G2, `NSS_3_1` for group G3,
never a bare `NSS_1`) — verified byte-for-byte against real `diyabc`
output on `toy_example2_ms_dna`. `compute_all_statistics_dna(header_text,
tree_sequences_by_locus, population_names)` is the top-level DNA entry
point (mirrors `compute_all_statistics`), now wired into `pipeline.py`
(see "DIYABC-replay pipeline for DNA sequences" below) — the
`stats_group_parser.py` docstring's old "dedup deferred, unclear if
legitimate" framing is stale, ignore it.

**Not yet done**: MicroSat's own stats (`NAL`/`HET`/`VAR`/`MGW`/`FST`/
`LIK`/`DAS`/`DM2`) not started at all — different data shape (allele
sizes, not tskit genotypes), no simulation side either. See
`notes/resume_stat_dna.md` for a biology-first (not code-first)
explanation of what each of the 13 DNA stats measures.

### DIYABC-replay pipeline for DNA sequences (2026-08-26) — cross-validated against a real reftable

Mirrors the SNP-side replay architecture (see item 12 above,
`run_reftable_simulation` vs. `replay_reftable_simulation`) exactly: a
`_dna`/`_from_values` sibling of each SNP function, never modifying the
SNP originals, so this can't regress the already-validated SNP path.
Six pieces, all in `bridge/`:

1. `reftable_loop.group_prior_column_names(header_text)` — real
   reftable column names for group priors (`µseq_2`, `k1seq_2`,
   `µmic_1`, `pmic_1`...). Verified scenario-INDEPENDENT (`nparamut` is
   a constant across scenarios, unlike `nparam` for historical params)
   — do NOT reuse `_kept_param_names_by_scenario` here, it expects
   `Prior` objects with `.bounds` (not strings) and filters by
   scenario, a concept group priors don't have.
2. `reftable_loop.parse_real_reftable_params_with_group_priors(path,
   priors, scenarios, group_priors_names)` — a NEW sibling of
   `parse_real_reftable_params`, returns triplets `(scenario_index,
   historical_values, group_priors_values)` with the two value dicts
   kept SEPARATE (they feed different downstream stages: demography
   vs. mutation model). Like the SNP original, handles the real
   reftable's per-row ragged column width (see `parse_real_
   reftable_params`'s own docstring for why a naive whitespace-split
   parser silently misaligns columns here).
3. `ancestry_simulation._group_prior_values_from_columns(group_priors_
   values, group_priors)` — reshapes the flat real-column dict into the
   nested `{group: {"MEANMU":.., "MEANK1":..}}` shape `draw_group_
   parameter_values` already produces, so `build_group_local_param_
   per_locus`'s existing body (model branching, per-locus `sampling_
   group_local_param` calls) is reused byte-for-byte in `build_group_
   local_param_per_locus_from_values`. Only the group-level (tier 1)
   draw is replaced with the real value; the per-locus (tier 2)
   dispersion around that mean is NEVER replaced — real DIYABC doesn't
   record it in the reftable, so there's nothing to replay, it keeps
   drawing from `seed`. Same principle propagates through `build_
   matrix_per_locus_from_values`/`build_rate_map_per_locus_from_values`/
   `dna_mutation_simulation_per_locus_from_values`.
4. `pipeline.compute_summary_statistics_dna`/`_from_values`,
   `reftable_loop._run_single_particle_dna`/`_from_values`,
   `reftable_loop.replay_reftable_simulation_dna` — each a direct
   mirror of its SNP sibling, no `num_loci`/`observed_reads_per_locus`
   params (no PoolSeq/loci-truncation concept for DNA sequences).
5. `write_reftable_txt` needed ZERO new code — already fully generic
   (`ParticleResult` + `priors`/`scenarios`), reused as-is for DNA
   results.
6. `scripts/replay_diyabc_priors_dna.py` (and its 50-loci-dataset copy
   `scripts/replay_diyabc_priors_dna_50loci.py`) — the runnable
   end-to-end script, mirrors `scripts/replay_diyabc_priors.py`. Run as
   `python3 -m scripts.replay_diyabc_priors_dna` from the repo root
   (see the `scripts/` section below for why). A 1000-particle real
   reftable replay runs in under 2 minutes.

**Validation result** (`toy_example2_ms_dna`, scenario 1, 1000
particles, 5+5 DNA loci): historical params match the real reftable
EXACTLY (`rdiff_mean = 0`, KS `p = 1.0` on `N1`/`t1`/`ta`/`ra`/`t2` —
confirms the replay plumbing itself is correct). Of the 42 DNA stat
columns, 11 show KS `p<0.05` — see the next section for the full
investigation of that gap.

### RESOLVED 2026-09-02: G3 (`<M>`) variance deficit in DNA-sequence stats (investigated 2026-08-27, reopened 2026-08-31, fixed by the user 2026-09-02)

On the 5+5-loci validation above, 11/42 DNA stat columns show a
significant KS difference, concentrated almost entirely in `MNS`/`VNS`/
`DTA`/`VPD` (variance-of-pairwise-differences-derived stats) and worse
on group G3 (`<M>`, mitochondrial/haploid) than G2 (`<A>`, autosomal
diploid). Full investigation in `notes/exploration.md`'s 2026-08-27
entry (source citations, tables, code excerpts); summary:

- **Tested and REFUTED**: "not enough loci" (the precedent from the
  2026-07-17 SNP entry above, where 10→650 loci per type resolved a
  similar-looking gap). Built `reference/toy_example2_ms_dna_50loci/`
  (new directory, `toy_example2_ms_dna` never touched — 50+50 DNA loci
  instead of 5+5, duplicated under new sequential names continuing the
  original numbering to avoid collisions), reran the same real-DIYABC
  replay validation: **17/47 columns significant, not fewer** — the
  proportion did not shrink. The corrected picture (a first comparison
  attempt was itself buggy — a naive `pd.read_csv(sep=r'\s+')` on the
  real reftable's ragged rows silently misaligns columns; fixed by
  reusing `_kept_param_names_by_scenario`/`group_prior_column_names` to
  parse each row by its own scenario) sharpened rather than changed the
  finding: G2 matches DIYABC very well (std ratio sim/real 0.94-1.04
  across all 8 stats), G3 is where the gap concentrates, and it's a
  **variance deficit** (std ratio sim/real 0.65-0.94, down to **0.26 on
  `DTA`**), not a location shift — `DTA`'s huge-looking `rdiff%` (up to
  +139%) is a red herring since its true mean is ~0, tiny absolute
  diffs blow up as a percentage.
- **Ruled out via two throwaway diagnostics measuring cross-locus CV
  (coefficient of variation) within real replayed particles**: ancestry
  alone (`msprime.sim_ancestry`, no mutation, `total_branch_length` as
  proxy) gives CV(G3)/CV(G2) ≈ **1.05**; the full mutation pipeline
  (real replayed group-prior means, `num_mutations` as proxy) gives ≈
  **1.02**. Both near 1: our own G2 and G3 have essentially the SAME
  relative inter-locus variance at every stage — the per-locus
  dispersion machinery, the `coalescence_coefficient`/`rescale_
  demography`/`ploidy` dispatch, and the mutation submodel are all
  internally consistent. The gap is specifically that real DIYABC's G3
  carries MORE variance relative to its own mean than G2, while our
  simulation's ratio stays flat around 1.
- **Root cause traced directly in `particuleC.cpp`** (sibling repo
  `~/Documents/Github/diyabc`, see `reference_diyabc_cpp_repo` project
  memory) rather than continuing to theorize from the Python side:
  - `ParticleC::coal_pop`'s continuous-approximation waiting-time
    formula (`coeffcoal * N / (k*(k-1)) * (-log(ra))`, lines ~1340)
    implies `Ne_effective = coeffcoal*N/2` — confirms our own
    `rescale_demography(factor=coeffcoal/2)` is an exact match on the
    MEAN behavior.
  - `ParticleC::evalcriterium` (lines 1251-1275): DIYABC has a SECOND,
    discrete "generation per generation" Wright-Fisher coalescent mode
    we don't replicate at all. Checked empirically (temporary trace
    instrumentation in `coal_pop`, rebuilt via the repo's own
    `CMakeLists.txt`, reverted after) — on this dataset's typical `N1`
    (1000-10000), the discrete mode never triggers for either heritage
    type. Ruled out as the cause, but a real, un-replicated code path
    to keep in mind for a dataset with much smaller `N`.
  - `ParticleC::split_pop` (lines 1513-1524): the admixture event
    (`ta split`) draws an INDEPENDENT Bernoulli coin per surviving
    lineage to assign it to one of the two ancestral populations.
    Traced 4 real particles: `<A>` arrives at `ta` with 5-13 surviving
    lineages (fairly evenly split), `<M>` with only 1-6 (one side is 0
    lineages in 3 of 4 cases) — with this few lineages the admixture
    partition becomes near-binary (all-or-nothing), a qualitatively
    different, more discrete regime than `<A>`'s. This plausibly
    explains why the earlier `total_branch_length`-based CV diagnostic
    found nothing: branch length is continuous and averages out an
    all-or-nothing population-assignment effect, while `DTA`/pairwise-
    difference stats are specifically sensitive to that exact kind of
    population-structure pattern.
  - **Final check — does OUR OWN msprime simulation reproduce this same
    effect?** `msprime.sim_ancestry(..., record_migrations=True)` +
    inspecting `ts.tables.migrations` (`dest` population at `t=ta`) is
    the direct msprime-side equivalent of tracing `split_pop`. Across
    30 real particles (1332 G2 loci, 968 G3 loci): G3 has one side at
    exactly 0 lineages in **55.6%** of loci vs. **30.6%** for G2 — same
    direction, comparable magnitude to the C++ trace. **Confirms our
    own port DOES reproduce the near-binary effect** (`msprime.add_
    admixture` does the same per-lineage independent draw as
    `split_pop`) — this is NOT a missing/underproduced effect in the
    port.

**Conclusion (2026-08-27, SUPERSEDED — see 2026-08-31 update below): no
bug found in the port.** The residual variance gap on G3 was attributed
to a real, source-confirmed combinatorial effect (few surviving
lineages at the admixture event → near-binary partition) present and
correctly reproduced on both sides. Investigation closed 2026-08-27.
**This attribution turned out to be incomplete** — the admixture
mechanism is real and correctly ported, but is NOT the (sole) cause of
the gap, see below.

**Update 2026-08-31 — falsification test reopens the investigation**:
`toy_example2_ms_dna` has 2 candidate scenarios drawn at equal weight
(`[0.5]`/`[0.5]`) — scenario 1 (the one studied above: merge, THEN `ta
split` admixture, THEN remerge) and scenario 2 (just a merge at `t1`,
**no admixture event at all**). If the admixture near-binary-partition
mechanism were the (main) cause, the G3 variance deficit should shrink
or vanish for particles that drew scenario 2. Tested directly on the
already-replayed 1000-particle real reftable (no new DIYABC run
needed): split the 1000 particles by their own `scenario_index` (488
scenario 1, 512 scenario 2, 0 scenario mismatches between the real and
replayed files — confirms row-by-row pairing is intact), recomputed the
std-ratio sim/real separately per scenario. Result: the deficit is
**essentially unchanged** between the two (mean G3-G2 std-ratio gap:
-0.177 with admixture vs. **-0.165 without**), and on `VPD`/`VNS` it is
actually **slightly worse without admixture** (0.714/0.778 vs.
0.909/0.830 with) — checked consistent across all 13 individual stat
families, not just an aggregate average. **This refutes admixture as
the main driver**: a mechanism tied specifically to the `ta split` event
cannot explain a deficit that survives intact in a scenario with no
such event. The true cause is more likely something general to `<M>`
itself (its reduced `Nₑ` means fewer surviving lineages at ANY point in
its history, not just at a discrete admixture event — plausibly fewer
independent coalescence events overall to average/smooth statistics
over, everywhere in the tree, not only at one partition point) — **not
yet investigated**. Full methodology and per-stat breakdown in
`notes/exploration.md`'s 2026-08-31 update (appended to the 2026-08-27
entry).

**Investigation status (as of 2026-08-31, before the fix below): REOPENED.**
The previously documented admixture explanation was incomplete.

**Fixed 2026-09-02.** Cause: `<M>` DNA-sequence loci (mitochondrial,
non-recombining) were each drawing their own independent genealogy in
`dna_mutation_simulation_per_locus`/`_from_values`
(`bridge/ancestry_simulation.py`), instead of sharing one — the SNP
side already got this right (`simulate_shared_ancestry_loci`, citing
`particuleC.cpp:2422-2435` `GeneTreeY`/`GeneTreeM`), the DNA-sequence
side never did. Averaging over 5 independently-drawn genealogies
instead of 1 shared one artificially shrank the inter-locus noise on
group-mean stats (`DTA_3`/`VNS_3`/`MNS_3`/`VPD_3`) — exactly the
deficit chased since 2026-08-27. Fix: new constant
`_SHARED_M_ANCESTRY_SEED_OFFSET`; a `<M>` locus now draws its ancestry
with a FIXED seed (`seed + _SHARED_M_ANCESTRY_SEED_OFFSET`, no `+i`),
so every `<M>` locus in the dataset shares the identical genealogy
(topology/node-times are seed-deterministic regardless of
`sequence_length` when there's no recombination, verified empirically).
`<A>`/`<H>` unaffected. Validated on a full 1000-particle real-reftable
replay: KS-significant columns dropped from 11/42 to 2/42, and the 2
remaining (`MPD_2_2`, `VPD_2_2`) are G2 columns with the ratio the
*other* way — ordinary sampling noise, not a deficit. Full
investigation writeup in `notes/exploration.md`. This fix was written
by the user (mentor mode), not the assistant.

**One separate, real bug found and fixed along the way** (unrelated to
the above — ruled out as its cause): `ancestry_simulation.build_group_
local_param_per_locus`/its `_from_values` twin recreated `random.
Random(seed + _KAPPA1_SEED_OFFSET)` (and `_KAPPA2_`/`_MUS_RATE_SEED_
OFFSET`) fresh inside the `for group in nloc_per_group` loop, with an
offset that didn't depend on the group — so any two groups sharing a
model needing the same parameter (G2 and G3 here, both `K2P`, both
using kappa1) replayed the IDENTICAL draw sequence, just recentred on a
different `k_moy`. Invisible at 5+5 loci (both groups' `GAMK1` happen
to declare the same shape, so relative dispersion looked identical by
coincidence); only surfaced once this investigation needed a real
multi-group same-model dataset to stress it. Fixed by constructing each
parameter's `rng` ONCE before the group loop (same pattern as `build_
rate_map_per_locus`'s already-correct `_SITE_RATE_SEED_OFFSET` usage).
Fixing it changed exactly the 15 G3/pairwise "golden value" test
assertions in `tests/test_summary_statistics.py`/`test_pipeline.py`
(regenerated) and left every G2-only assertion byte-identical — the
exact fingerprint confirming the fix is correctly scoped (G2 is always
first in the loop, so its RNG's starting state never changed). See
`feedback_seed_reuse_pattern` project memory — same bug class, 5th
occurrence in this codebase.

**Open, unexplained thread (flagged 2026-08-27, not investigated)**: on
the 50+50-loci dataset, `DTA_2_2` (and to a lesser extent `DTA_2_1`) —
both **G2** columns — also show a real, visually-confirmed distribution
difference, but it does NOT fit the G3-variance-deficit story above:
their std ratio sim/real is close to 1 (like the rest of G2), the
difference is instead a small but real MEAN shift (real ≈0.018, sim
≈0.043 for `DTA_2_2`). Not yet investigated — a legitimate loose end,
distinct from the G3 mechanism this section documents.

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
