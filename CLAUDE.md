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
**Not done yet**: wiring any of this into `reftable_loop.py`/the
summary-statistics side of the pipeline (no DNA-sequence summary
statistics implemented, no genotype extraction from the resulting
`TreeSequence`). MicroSat itself (stepwise mutation model, `NAL`/
`HET`-style summary statistics) has no simulation-side code at all
yet, only header parsing. The goal is to demonstrate that a
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
   → one `LociDescription` with `total_loci: dict[heritage_type, count]`)
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

**Not yet done**: `MPD`/`VPD` aggregation (mean/variance of
`_pairwise_hamming_distances`'s output — `VPD` additionally needs the
C++'s `nd > 1` per-locus exclusion, not just the `num_loci` denominator
the other stats use, since a locus contributes to `VPD`'s average only
if it has at least 2 valid pairs) plus `DTA`/`PSS`/`MNS`/`VNS` (the
remaining per-population stats) and all 5 pairwise stats (`NH2`/`NS2`/
`MP2`/`MPB`/`HST`) are unwritten. No test yet for the `MPD`/`VPD`
aggregation layer (only the `_pairwise_hamming_distances` brick is
tested). Also unresolved: `toy_example2_ms_dna/headerRF.txt` declares
`NSS`/`NHA`/etc. under the **same** column names (`NSS_1`, `NSS_2`,
...) in both `group G2` (`<A>`) and `group G3` (`<M>`) — the exact
column-name-collision-across-groups scenario `stats_group_parser.py`'s
docstring already flagged as an open question (dedup deferred, unclear
if legitimate) — still undecided how the eventual `compute_all_
statistics`-equivalent entry point for DNA sequences should key/merge
per-group results into `reftable.bin` columns when this happens. Not
wired into `reftable_loop.py` at all yet.

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
