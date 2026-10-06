# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

**This file holds rules and verdicts only.** The full narrative of every
chantier — investigations, falsified hypotheses, bugs caught in review,
source citations — lives in `notes/exploration.md`. Read that file before
touching parsing or statistics code; it records *why* things are done a
certain way. Do not grow this file back into a journal.

## What this is

A proof-of-concept replacing DIYABC's C++ coalescent simulator
(`particuleC.cpp::dosimulpart`) with `msprime`. The pipeline (`bridge/`)
turns a DIYABC `header.txt` + observed data file into a `reftable.bin`
entirely in Python — no subprocess call to the DIYABC binary — and must be
structurally and statistically equivalent to the real DIYABC's output.

| Data family | State | Validated against |
|---|---|---|
| SNP IndSeq | complete, 130 stats | `human` (5000 loci), `toy_example3`, `toy_example5` |
| SNP PoolSeq | complete | `toy_example4` |
| Heritage `<A>/<H>/<X>/<Y>/<M>` (SNP) | complete | `human`, `toy_example5` |
| MicroSat | complete, 11 stats + `AML` | `toy_example2_ms_dna_K2P`, `toy_example1_ms_modified` |
| DNA sequences | complete, 13 stats; `JK`/`K2P`/`TN` all covered | `toy_example2_ms_dna_{JK,K2P,TN}` (5 loci/group, 2 runs each), `toy_example2_ms_dna_50loci_K2P` (50/group, 2 runs) — the **G2 `MPD`/`VPD` residual is a DIYABC bug, located to the line** (its `MPD` denominator drops loci the simulation left monomorphic), see Open work |
| DNA sequences `<X>`/`<Y>` | implemented | synthetic fixture only — the real binary SIGSEGVs on this case |
| Serial/temporal sampling (MicroSat) | complete | `toy_example1_ms` (1 population, 4 sampling times) |
| Serial/temporal sampling (SNP IndSeq) | complete | `human_seriel` (1 population, 4 sampling times, 130 stats, 0/130 significant) |
| Serial/temporal sampling (SNP PoolSeq) | complete | `toy_example4_seriel` (100 loci, `<MRC=5>`, 4/133 vs ~6.7 expected) |
| Mixed MicroSat + DNA in one header | complete | `toy_example2_ms_dna_TN` (73 columns, 2 replays × 2 scenarios) |

**Validation status (06/10): every dataset in the "Validated against" column above is
validated with the final code**
— the one in which all derived seeds go through `_locus_seed` and the serial
IndSeq replay is restored. Replayed on 05-06/10 and read with a KS and a sign test per
cell: `human` and `human_seriel` (0 KS flag in every cell, sign-test median
+0.05 % / −0.05 %), `toy_example1_ms_seriel` and `_modified`, `toy_example2_ms_dna_50loci_{JK,K2P,TN}`
(the `TN` one is the mixed MicroSat + ADN header, 73 columns, 0 flag in 4 cells),
`toy_example3` (100 loci: 9 flags in 6 cells, none shared between scenarios; 500 loci: 0),
`toy_example4`, and `toy_example4_seriel` at 1000 loci (3 flags in 4 cells over two
independent DIYABC runs, no column in common, `ML3p_2.3.4` absent). Any column flagged
once and not again is noise; the flags that were not were traced to their cause (see
Closed investigations: `DTA_2_2`, `ML3p_2.3.4`, the serial IndSeq replay).

Validation means a *paired* comparison: the real DIYABC priors are replayed
particle-by-particle through our pipeline (`replay_reftable_simulation_snp*`,
`scripts/replay_diyabc_priors*.py`) and the two reftables compared
column-by-column with a two-sample Kolmogorov-Smirnov test. That test is
**conservative** on a paired replay: both sides share the same prior draws, so
their statistics are correlated and the KS understates the gap. Judge a
residual flag by whether it **persists across two independent replays**, not by
its p-value alone — a real effect keeps its order of magnitude, noise moves to
other columns (see `notes/exploration.md`, 24/09, for a worked example).
**Persistence only counts if the msprime-side noise of the two replays is
independent — before 05/10 it was not.** Particle seeds are `particle_index + 1`,
identical in every replay, and the per-locus seeds were `seed + OFFSET + i`, so
neighbouring particles reused the same random streams. A residual could then
"persist" across replays, models and scenarios by being one realization of *our*
noise replayed each time: `DTA_2_2` was KS-flagged in **11 cells of 12** (3 models
× 2 runs × 2 scenarios; the 12th at p = 0.16) and was exactly that. See "`DTA_2_2` residual (05/10)"
under Closed investigations.

**A per-column KS is also structurally blind to a small systematic shift.** A
uniform -1% offset across every column stays well inside each column's
inter-particle variance and never shows up, while a **sign test on
`rdiff_mean`** exposes it immediately. Run both: the KS for per-column
divergences, the sign test for a global bias. Several datasets declared
validated on the KS alone carry such a bias; it turned out to be the
loci-count residual — see "Systematic negative bias" under Closed
investigations.

**Report the sign test's `k/n`, p AND median — and trust none of the three on
one replay.** The p swings wildly with how many columns sit near zero while the
amplitude barely moves: on `toy_example2_ms_dna_K2P`, scenario 1 gives median
+1.569% at p = 0.644 and scenario 2 median +1.514% at p = 9.4e-04 — 0.05 point
apart in amplitude, three orders of magnitude apart in p. That made the median
look like the stable estimator. **It is not, and the two scenarios of one run
are not two replicates** — they share one DIYABC reftable and one set of
replayed draws, so their agreement measures nothing. A fresh DIYABC run moved
`JK` scenario 1 from +2.361% to −0.537%, sign included. **An independent replay
means a new DIYABC reftable, never another scenario of the same one** — missing
that distinction produced two successive wrong verdicts on the DNA shift (see
Open work). Per-run noise there is ±2 points, so aggregate by run before
concluding anything: 6 run-level means, all positive, is the signal; 12 raw
medians spanning −0.76% to +4.69% is not. The p is also
**anti-conservative**, since the columns are far from independent (the six
`FST2m` share populations, `HWm`/`HWv` come from one computation), so the
effective number of units is well below `n` — the opposite bias to the KS's.
Exclude the prior columns from `n`: on a paired replay they are the same draws
on both sides (KS p = 1 exactly), and they only add coin flips. That
prior check is itself worth keeping as a two-line sanity test — if a prior
column departs from p = 1, the replay is not paired and nothing downstream
means anything.
**Since 01/10 the mutation-parameter columns (`µmic_N`, `pmic_N`, `snimic_N`,
`µseq_N`, `k1seq_N`, `k2seq_N`) belong to that same family** — we now emit them,
and in a replay they are DIYABC's own injected values, so they come out at
`rdiff_mean = 0` and `ks_stat = 0` exactly. Exclude them from `n` too, and keep
them as the **second** pairing sanity test: a departure from zero there means the
mutation-parameter replay chain is broken, which the historical-prior check
cannot see (see the `nparamut` entry under Closed investigations).
**That second test is only as good as the text format it is read through, and it
was blind until 02/10.** `write_reftable_txt` and `rewrite_real_reftable_txt` both
wrote every column with `%12.6f`, whose precision is **absolute**, so it destroys
small values — and mutation rates are small. Measured on
`toy_example2_ms_dna_50loci_K2P`, `µseq_2 ~ UN[1e-8,1e-6]`, 1000 particles:
`%12.6f` left **2 distinct values, 508 of them crushed to exactly 0**, so the
pairing test compared 0 to 0 on half the particles and would have passed with our
injected µ off by a factor 50. `snimic_1` was wiped out entirely. The fix is
**per column family**: parameters (historical and mutational) in `%.6g` — 6
*significant* digits, scientific notation when needed, which is what DIYABC
writes itself (`8.749e-07`) — and statistics in `%12.8f`, where values of 0.02 to
12 are exact and the columns stay aligned. `%12.8f` on µ is **not** enough: it
leaves 100 distinct values and up to 19% relative error on the smallest, i.e.
precisely the low-`θ` particles where the G2 residual lives. Do not reunify the
two formats. `_split_mutation_and_stat_names` is what lets
`rewrite_real_reftable_txt` tell the families apart — the mutation columns cannot
be counted from `priors`/`scenarios`, which carry no `group priors` section, so it
reads the name: a statistic always starts with an uppercase ASCII letter
(`NAL_1_1`, `MPD_2_1`, `ML3p_1.2.3`), a mutation parameter never (`µseq_2`,
`pmic_1`, `k1seq_2`). It **raises** rather than shifting formats silently.
The `.bin` was never affected (`struct.pack` float32), so the structural
deliverable was always sound — only the diagnostic path lost information.
**Compare against DIYABC's text at DIYABC's precision: round our statistic
columns to 6 decimals before the KS.** DIYABC's text export writes 6 decimals, so a
statistic of order 1e-5 (`F3v`, `F4v`, `NEIv`, `FST2v`) takes a few dozen distinct
values on its side (`F3v_1.2.3`: 37 over 487 particles, against 395 for ours at 8
decimals). The KS then measures the mass of one quantization step (0.05 to 0.13) and
flags the column while the means agree. Measured: `human_seriel` scenario 2, 9 flags
(`F3v`/`F4v`) → 0 after rounding; the same artefact explains most flags on
small-amplitude columns elsewhere. The sign test is unaffected (it reads means).
**And when a format changes, restart the Jupyter kernel before regenerating
anything**: a long-lived kernel keeps the pre-edit module in memory, which
produced a `first_records_clean` file *newer* than the replay yet written by the
old code, and a pairing test that failed for a purely cosmetic reason.

**A `rdiff` gradient along a variable that moves the mean is not an effect until
you have looked at it in ABSOLUTE difference.** This is the quantitative form of
the near-zero-mean rule, and it cost this project its largest open chantier. On
`toy_example2_ms_dna_50loci_K2P` G2, `MPD`'s relative error runs −69.5% → −2.3%
across strata of the drawn `N1` — a factor **30**, read for five weeks as an
`N`-dependent genealogy effect. The absolute gap over those same strata is
+0.051 / +0.073 / +0.077 / +0.060 / +0.032: flat, CV 28%. The whole "gradient"
was one near-constant offset divided by a mean that varies 20x. So whenever a
residual is stratified by anything, **report the absolute difference next to the
`rdiff`, and read the sign and the trend off the absolute one.**

`reference/` holds ground truth produced by the real DIYABC binary —
**never modify these files.**

## Environment

Use the `diyabc_msprime` conda environment (Python 3.11, has msprime 1.4.2,
tskit, numpy, scipy):

```bash
conda activate diyabc_msprime
```

The system Python (3.13, no conda env) does NOT have msprime installed —
always activate `diyabc_msprime` before running anything in this repo.

The DIYABC C++ source is a sibling repo at `~/Documents/Github/diyabc`
(`src-JMC-C++/`) and is the ground truth for every formula. Don't nest it
inside this project.

## Common commands

```bash
# Run the test suite
pytest tests/ -v

# Lint / format (ruff, config in pyproject.toml)
ruff check .
ruff format .

# Pre-commit hooks (ruff check --fix + ruff format)
pre-commit run --all-files

# Regenerate notes/tree.md, notes/commits.md, notes/api.md, notes/report.md
python3 tools/generate_report.py
```

## Code style

Docstrings are **Google style** (`Args:`/`Returns:`/`Raises:`, one-line
summary first) — see `loci_parser.py` for a reference example. Write any new
or edited docstring the same way; don't reintroduce the old free-prose style.

### Vocabulary: population vs. sample

DIYABC's own distinction, and the codebase now follows it. **A sample is not a
population** — the discriminating test, which settles every case mechanically:
*on a serial dataset with 1 population and 4 samples, does this variable hold 1
value or 4?* Four → sample; one → population.

- **Population** keeps `pop`: `scenario_parser.py`, `demography_builder.py`,
  `event.pop`, the msprime API (`add_population_split`, the `population=`
  kwarg), and the string literals `"pop1".."popN"`, which are msprime's own
  population names — **never rename those**, the demography resolves them by
  name.
- **Sample** for everything that comes out of the `.snp`/`.mss` file and
  everything indexed by a statistics column: `count_individuals_per_sample`,
  `sample_index_to_name`, `counts_per_sample`, `sexes_per_sample`.
- **Spelling: `_per_sample`, singular.** A plural `_per_samples` coexisted for a
  while and is one letter away from the singular — unreadable. Don't
  reintroduce it.
- A **file-format** name stays the format's: the `.snp` column called `POP`
  stays `POP` (`pop_column_index`) even though it carries a sample identifier.
  Where format vocabulary and model vocabulary diverge, both words belong in the
  same function.
- `compute_population_layout` and `compute_sample_layout` **both keep their
  name**: the former really does slice by `ts.samples(population=...)`. Don't
  merge them — their equality on non-serial data is what proves the substitution
  neutral. Same verdict for a local built from a **tskit** attribute
  (`ts.node(...).population`, `ts.populations()`): that is a real population, so
  `pop_of_ind` / `samples_per_population` in the tests are correct as they stand.
  In `test_summary_statistics.py` the latter is even a useful tell — the test
  indexes a population-keyed dict by sample name, which only works because the
  dataset is non-serial. Renaming it would erase that clue.

**Done**: `observed_data.py`, `ancestry_simulation.py`, `reftable_loop.py`,
`pipeline.py`, `header_dataclasses.py`, `summary_statistics.py` (~1250
occurrences, 9 private functions) and their tests — commits `85471a0`,
`e5119b6` (a `docformatter` pass), `5245722`.

**The rename is provably output-neutral**, and this is why: every reftable
column name is a literal prefix followed by an **integer** (`f"ML1p_{i + 1}"`,
`f"FST2m_{key}"` with `key = f"{i + 1}.{j + 1}"`). No Python identifier ever
reaches an output string, so no column can move. The golden tests, keyed by
column name, confirm it empirically.

**Renaming prose is not renaming identifiers**, and it bit twice here. A
`population` → `échantillon` substitution in French breaks agreement (gender
changes: 69 cases of `de échantillon`, `la échantillon`, `une échantillon`,
plus non-adjacent ones like `toutes les échantillons … sont présentes`). Worse,
it **merges two notions under one word**: `compute_MPD`'s docstring became "un
locus où un échantillon a moins de 2 échantillons", where the second word means
msprime sample *nodes*; and `_genotype_matrix_by_sample`'s said "un échantillon
n'est plus sa propre échantillon", destroying the very sentence that explained
why `layout=` exists. An identifier cannot be ambiguous with itself, a sentence
can. Reach for a third word when tskit's own `num_samples` is meant: gene
copies, sample nodes, sequences.

## Architecture

Each stage is a separate module with no cross-cutting logic:

1. **`scenario_parser.py`** — `header.txt` → `Scenario` objects. Event
   vocabulary from `history.cpp::ScenarioC::read_events`: `sample`, `merge`,
   `varNe`, `split` (admixture). Only `NotImplementedError` is swallowed
   (block skipped with a warning); any other exception propagates.

2. **`prior_parser.py`** — `Prior` / `OrderConstraint` from the `historical
   parameters priors` section. `is_constant_prior` replicates the exact
   near-degenerate-bounds rule from `readReftable.R` / `abcranger`: constant
   if `(max-min)/max <= 1e-6`, **never** constant when `max == 0.0` — DIYABC
   excludes these from `reftable.bin` columns.
   `_extract_historical_priors_section` falls back to an empty line as
   terminator when `DRAW UNTIL` is absent (a header with 0 order constraints
   has no such line at all). `parse_group_priors` → `GroupPrior` handles the
   separate `group priors` section; `get_parameter_used_by_model` maps a
   group's `name_model` to which of `k1`/`k2` are active.

3. **`stats_group_parser.py`** — the summary-statistic column names actually
   declared in `group summary statistics`, used to filter computed stats down
   to what the header expects. Handles multiple `group Gx (N)` blocks; the
   flattened result is order-preserved and **deliberately not deduplicated**
   across groups.

4. **`loci_parser.py`** — the `loci description` section, both formats: the
   condensed one (`"70 <A> 10 <X> 10 <M> 10 <Y> G1 from 1"` → one
   `LociDescription` with `loci_counts_by_heritage`) and the detailed
   one-locus-per-line one (`"Locus_M_A_1_ <A> [M] G1 2 40"` →
   `list[LociDescriptionDetailed]`). The trailing `dnalength` on `[S]` lines
   is informational only — `header.cpp:401-402` never reads it there.

5. **`parameter_sampling.py`** — one value per prior by rejection sampling
   until all `OrderConstraint`s hold. `_draw_one_value` supports
   `UN`/`LU`/`NO`/`LN`/`GA`, each reproducing `PriorC::drawfromprior` /
   `MwcGen::g*` exactly — including truncation by rejection for `NO`/`LN`/`GA`
   (redraw until inside `[min,max]`, never an unbounded draw) and the `GA`
   reparameterization (`gammavariate(shape=sdshape, scale=mean/sdshape)`).
   Also `draw_group_parameter_values` (group priors),
   `sampling_group_local_param` (the per-locus tier) and `sample_site_rates`
   (`mutsit`).

6. **`demography_builder.py`** — `Scenario` + drawn values →
   `msprime.Demography`. Populations are named `"pop1".."popN"` by 1-indexed
   position in `header.txt`. `get_parameter_names_used_by_scenario` gives the
   prior names a *specific* scenario references — this subset, not the full
   prior list, must become the `reftable.bin` parameter columns (`human`
   declares 21 priors, scenario 1 uses 16).

7. **`observed_data.py`** — maps population index to the real name in the
   `.snp` file. The mapping is implicit: population *i* = the *i*-th
   population **by first-appearance order** in the file. **Never replace the
   order-preserving dict/Counter with anything that could reorder keys** (an
   alphabetical sort would silently break it). Also `parse_sex_ratio`,
   `parse_maf_ratio`, `parse_mrc_ratio`, and the `.mss` readers
   (`observed_sequences`, `observed_microsatellites`,
   `base_frequency_by_locus`, `observed_count_sample`,
   `individual_sexes_from_locus_genotype`). `.mss` is genepop-format (`POP` as
   a block separator), structurally different from `.snp`'s column format —
   the two families never share a reader.

8. **`ancestry_simulation.py`** — the simulation core (~2900 lines). SNP: one
   independent tree per locus, mutated with the Hudson algorithm — exactly one
   mutation per locus, placed on an edge chosen proportionally to branch
   length, vectorized over tskit's edge tables. This guarantees every locus is
   globally polymorphic by construction (DIYABC doc §2.4.3). MAF filtering
   (`with_maf_filter`) rejects and resimulates below threshold, mirroring
   `ParticleC::mafreached`; `<MAF=hudson>` is the no-op fast path. Also holds
   the whole DNA-sequence and MicroSat mutation machinery (see "Domain
   knowledge" below).

9. **`snp_writer.py`** — `.snp` output, only used by the deprecated
   subprocess path.

10. **`summary_statistics.py`** — pure numpy reimplementation of all stats:
    130 SNP (ML1-3, HW, HB, FST1-4, NEI, AML, F3, F4), 13 DNA-sequence, 11
    MicroSat + `AML`. Each function's docstring names the exact `sumstat.cpp`
    function it transcribes — **check there before changing a formula.**

11. **`pipeline.py`** — orchestrates 1-10, no logic of its own.
    `compute_summary_statistics_snp` (+ `_dna` / `_microsat` / `_from_values`
    variants) is the high-level entry point.

12. **`reftable_loop.py`** — runs `nrec` particles in parallel
    (`ProcessPoolExecutor`), writes the binary `reftable.bin` and the
    human-readable `write_reftable_txt`. Seeds are `particle_index + 1`
    (msprime rejects `seed=0`). Multi-scenario is supported: each particle
    draws its own scenario from the weighted list
    (`parameter_sampling.draw_scenario`, matching `ParticleC::drawscenario`),
    and `write_reftable_bin` writes a variable-length record per row (only
    that row's own scenario's `nparam` columns, no NA-padding), matching
    `reftable.cpp`. Both writers emit the mutation-parameter columns between the
    historical parameters and the statistics; `nparam[i]` in the binary header
    counts **both** families (`nparamhist + nparamut`, `reftable.cpp:64`) and
    there is no separate header field for `nparamut` — see Closed
    investigations.

### Two generations of architecture — mind the drift

An earlier version wrote simulated genotypes to a fake `.snp` file and shelled
out to the real DIYABC `general` binary to compute statistics. That path is
kept only for cross-validation (the `DIYABC_GENERAL_PATH` skip-marked tests)
— it is no longer the default and is far slower. If you see `-g`, note it is
the *internal batch size*, not a loop count (`-g 50` silently discards 49/50
simulated particles).

### Signatures: mind the ReplayContext refactor (2026-09-16/18)

`SnpReplayContext` / `MicrosatReplayContext` / `DnaReplayContext`
(`header_dataclasses.py`) are built **once per run**, not once per particle,
and hold everything previously re-read from disk. Consequence: most per-locus
builders took `(header_text, mss_file_path, ...)` and now take
`(context, ...)` — e.g. `build_matrix_dna_per_locus`,
`build_matrix_microsat_per_locus`, `microsat_mutation_simulation_per_locus`,
`dna_mutation_simulation_per_locus`. Not all moved:
`build_local_param_dna_per_locus`, `build_local_param_microsat_per_locus`
and `build_rate_map_per_locus` still take `(header_text, seed)` since they
read nothing off disk. **Any signature quoted in an older note may be stale —
check the real one before calling it from a notebook or script.**

Field-selection rule: a value belongs in the context if and only if computing
it involves an actual disk read, *not* merely because it varies per
locus/particle.

## Hard rules

- **`header.txt`/`headerRF.txt` trailer line.** The last line
  (`"scenario N1 N2 ... ML1p_1 ..."`) looks like output-column documentation
  but is **re-read as input** by `HeaderC::readHeaderAllStat` to derive
  `nparamhist = tokens - 1 - nstat - nparamut`. A trailer copied from another
  scenario miscounts it and corrupts DIYABC's internal state, producing stats
  300%–10000% off for every population except the scenario's hub. **When
  hand-crafting or editing a header, always regenerate this trailer line to
  match that scenario's own priors — never copy it from another scenario.**
  This has bitten the project **five** times, the last two on 30/09 when
  `toy_example2_ms_dna_50loci_TN` and `_50loci_JK` got `K2P`'s trailer: `k1seq`
  tokens under `MODEL JK` (two too many) and no `k2seq` under `MODEL TN` (two too
  few). Both reftables were unusable and had to be regenerated; the `TN` one had
  already produced an analysis that was retracted.
  **`reftable_loop.check_header_trailer_line` now guards this**, called from all
  seven run/replay entry points: it compares the trailer's mutation-parameter
  tokens against `group_prior_column_names`, which applies
  `get_parameter_used_by_model` (`JK` → no `k`, `K2P`/`HKY` → `k1`, `TN` → `k1`
  and `k2`). It is a no-op on a header with no `group priors` section, i.e. every
  SNP dataset. Note **DIYABC writes the `mus_rate` prefix two ways** — `µ`
  (U+00B5) on the `te2` family, plain ASCII `mu` on `toy_example1_ms` — and the
  reftable always follows its own header, so the guard normalizes both.
  `group_prior_column_names` hardcodes `µ`: on a `mu`-style header it returns
  column names that match nothing, which is latent (see Open work).
- **`DRAW UNTIL` exists if and only if `nconditions > 0`.** The
  `historical parameters priors (N,C)` header gives `C` order constraints;
  `readHeaderHistParam` (`header.cpp:274`) consumes the `DRAW UNTIL` line
  **inside** `if (this->nconditions > 0)`. With `C = 0`, that `getline` never
  runs, the stray line stays in the stream, and **every subsequent `getline`
  shifts by one** — `readHeadersimLoci` then reads `DRAW UNTIL` where it
  expects `loci description (N)` and fails with an unrelated-looking message.
  Reference: `human` (4 constraints) and `toy_example1_ms_modified` (3) have
  the line; `toy_example1_ms` (0) does not. Same trap family as the trailer
  line: a line that looks decorative is in fact a **positional token** in a
  sequential `getline` read, and the format tolerates no offset.
- **Check a real reftable against its header by CONTENT, never by timestamp.**
  `reftable_loop.check_real_reftable_matches_header` compares the reftable's own
  first line — its column names, written by `bintotxt` from the header's trailer —
  against that trailer, and is called from all three replay entry points. It
  catches both a header edited *after* generation and a reftable produced from a
  corrupt header, which no date can tell apart.
  **The old timestamp rule was wrong and is retired**: it required the reftable to
  be newer than `headerRF.txt`, but DIYABC rewrites `headerRF.txt` *after* the
  reftable in the same run — measured at 6 s on `toy_example2_ms_dna_50loci_TN` —
  so a correct run always fails it, and it distinguishes a 6-second-old pair from
  a five-week-old one not at all. It caught neither of the two corruptions of
  30/09; content did. Treat timestamps as a hint, never as a check.
- **`reference/` is read-only ground truth.** Never modify it.
- **Run scripts as modules from the repo root**: `python3 -m scripts.run_test`.
  `python3 scripts/run_test.py` fails with `ModuleNotFoundError: No module
  named 'bridge'` — Python only puts the script's own directory on `sys.path`.
- **A test that writes a modified `header.txt` to `tmp_path` must build a
  fresh `*ReplayContext`** from that text. Since the ReplayContext refactor,
  functions never re-read the header from disk, so "write a new file, then
  call the function on that directory" silently uses the original.
- **Never read a raw DIYABC reftable with a generic whitespace parser.**
  `first_records_of_the_reference_table_0.txt` is **ragged**: `bintotxt` omits
  the parameters a row's own scenario does not use, so those rows carry fewer
  fields than the header line. `pandas.read_csv(sep=r"\s+")` left-aligns them,
  which shifts **every** later column — mutation parameters and statistics
  included — on those rows only. Measured on `toy_example2_ms_dna_K2P`: header
  71 fields, scenario-1 rows 71, scenario-2 rows 68 (no `ta`/`ra`/`t2`), so the
  column labelled `ta` actually held a microsat mutation rate and `k1seq_2` held
  a statistic. Always read the rectangular rewrite produced by
  `rewrite_real_reftable_txt` (`first_records_clean.txt`, unused parameters
  written as `nan`, which is what our own `write_reftable_txt` emits too). Same
  trap family as the trailer line and `DRAW UNTIL`: a line whose *length* is
  load-bearing. It cost a wrong conclusion on 30/09 — and note that every
  historical validation of this project was run with `scenario = 1` hardcoded in
  the comparison notebook, i.e. on the only rows the raggedness never touched.

- **Never derive a seed by addition.** `seed + OFFSET + i` gives locus `i` of
  particle `s` the same seed as locus `i - 1` of particle `s + 1`, because particle
  seeds are consecutive: every random stream was shared by up to `nloci` neighbouring
  particles, which were therefore not independent. **Every derived seed in `bridge/`
  now goes through `configuration._locus_seed(seed, offset, index)`**
  (`np.random.SeedSequence([seed, offset, index])`, folded into `[1, 2**32 - 1]`
  because msprime rejects 0; `index` is the locus, batch, attempt or group, and 0 for a
  draw made once per particle). Converted on 05/10: the per-locus tree and mutation
  seeds (DNA, MicroSat), `with_mrc_filter` and `with_maf_filter` (batch and attempt
  seeds), `LIK`/`AML` (`_LIKELIHOOD_SEED_OFFSET`, `_GROUP_STAT_SEED_OFFSET`), and the
  single-offset draws (`random.Random(seed + OFFSET)` for µ/k/site rates/MicroSat,
  shared `<M>`/`<Y>` ancestry, scenario draw, group priors, binomial reads,
  `seed_for_type`) — the latter never overlapped below 1M particles, converted for
  uniformity and so that no new site reintroduces the pattern. **What stays raw**: the
  particle seed itself (`random.Random(seed)` for the historical priors,
  `random_seed=seed` for the SNP trees), unique per particle and not derived. Adding a
  new offset constant means a new entry in `configuration.py`, never a sum. **Every
  value gated by a seed moved with this change**: stored comparisons and reftables
  made before 05/10 are no longer comparable to a fresh run. Every validated dataset
  was re-replayed after the full conversion (see the validation status under the
  table at the top).

## Closed investigations — do not reopen

Each line is a verdict; the full investigation is in `notes/exploration.md`
under the date given.

- **Residual statistical bias (17/07)** — not a simulator bug. Shrinks and
  vanishes as loci count grows; 0/130 significant on `human` at 5000 loci.
- **Performance gap vs. DIYABC (20/07)** — understood, not a bug. ~300s vs
  ~137s on `human` 1000 particles × 5000 loci. Per-particle single-threaded
  cost is nearly identical to DIYABC's; the gap is the 8-physical-core ceiling
  plus materializing a `TreeSequence` per locus. **Don't re-investigate the
  stat formulas or `max_workers`.**
- **G3 `<M>` variance deficit in DNA stats (02/09)** — cause: `<M>` loci were
  each drawing an independent genealogy instead of sharing one. Fixed with
  `_SHARED_M_ANCESTRY_SEED_OFFSET`; KS-significant columns 11/42 → 2/42. The
  earlier "combinatorial admixture effect" attribution was wrong and was
  falsified by a no-admixture control scenario.
- **`FST` microsat gap, ~2x too low (21/09)** — cause: `np.bool_ + np.bool_`
  is a **logical OR**, not an integer sum, so `nA = sum((p[0]==al)+(p[1]==al))`
  undercounted every homozygote. Catastrophic on `<M>` (haploid duplication
  makes every individual a "homozygote" → negative FST). Fixed at both layers
  (`int()` in `_length_by_sample_and_individuals` *and* in
  `_compute_ni_nA_AA_for_one_sample`); all 11 MicroSat stats now match
  real DIYABC within noise. **Not** the genealogy, **not** admixture,
  **not** the formula — all six re-verifications of `cal_Fst2p` were correct.
- **SNI as the cause of the `FST` gap** — falsified twice, from both
  directions: forcing `SNI≈0` in a real DIYABC dataset didn't move the gap,
  and implementing SNI in our own simulation didn't either.
- **Admixture as the cause of either gap** — falsified for both the G3 `<M>`
  deficit and `FST`, by per-scenario splits on no-admixture control scenarios.
- **"DIYABC switches to a simplified substitution model below some threshold"
  (21/09)** — **no such thing exists in the C++**, checked exhaustively for
  SNP, MicroSat and DNA. `mutmod` is written only when the header token is
  read; `comp_matQ` branches on `mutmod` alone. The only genuine second mode
  anywhere is the discrete generation-by-generation coalescent selected by
  `evalcriterium` (`particuleC.cpp:1251-1275`) — a coalescent-side switch, not
  replicated. **Correction (25/09): it IS reachable.** The criterion is
  `ra = nLineages / N`, and the continuous approximation is kept only when
  `ra < 0.5` on segments longer than 100 generations (more permissive below:
  `ra < 0.033*nGen + 1.7` for `nGen <= 100`). A serial dataset concentrating a
  large sample in ONE population reaches it easily — `toy_example4_seriel` has
  400 individuals (800 gene copies) with `Npresent <= 1000`, so the 200→500
  segment is in discrete mode for **every** particle. The earlier "never
  triggered" claim was true of the datasets then in the repo, where samples
  were small relative to `N`.
  **The 30/09 verdict "it explains a real discrepancy" is RETIRED — falsified
  02/10 by direct measurement, and the chantier it opened is closed with it.**
  What the discrete branch actually is (`particuleC.cpp::coal_pop`, 1398-1500): a
  **haploid Wright-Fisher on `Ne = (int)(0.5*coeffcoal*N + 0.5)` gene copies**
  (`coeffcoal = 4` for `<A>` at `sexratio = 0.5`, so `Ne = 2N`). Each lineage
  draws a uniform parent in `1..Ne` via `rand1`, and every lineage sharing a
  parent merges into a single node at that generation — so several lineages can
  merge at once, producing polytomies. It is **the same process as the continuous
  branch, not a second model**: the continuous wait `coeffcoal*N/(nl(nl-1))` is
  exactly `2*Ne/(nl(nl-1))`, its Kingman limit.
  Measured with a faithful Python port of that loop against
  `msprime.StandardCoalescent` and `msprime.DiscreteTimeWrightFisher`, on the
  mutation-free shape statistic `r = <mean pairwise distance> / <total branch
  length>` — the π/`S` of the low-mutation limit, worth `1/a_n` under Kingman —
  with n = 80 gene copies and 2000 replicates:

  | | `r` ≈ π/`S` | `L` ≈ `S` | π ≈ `r·L` |
  |---|---|---|---|
  | `ra = 3.2` | 0.969 | 1.039 | 1.007 |
  | `ra = 0.8` (the real switch point of G2's priors) | 0.994 | 1.000 | 0.995 |

  **Nul**, and DTWF tracks both Kingman and DIYABC for `Ne` from 50 to 3000. The
  reason is structural, so a different setup will not give a different answer: the
  multiple mergers happen at **near-zero branch length**, so π and `S` absorb them
  proportionally and the ratio is preserved. The port was validated from the other
  side — at `ra = 3.2` it drops 80 lineages to 41 in a single generation with
  merges of 4, and at `ra = 0.05` it degenerates to Kingman one merge at a time —
  so the polytomy mechanism really is there, and really does cost nothing.
  Consequences: **not replicating the discrete coalescent costs nothing
  measurable**; the `N1` gradient formerly tabulated here is not its signature
  (see the G2 entry under Open work for what it is); and the two obstacles that
  blocked the chantier — DTWF requires `ploidy = 2`, and the criterion is dynamic
  since `nLineages` is unknown before simulating — are moot.
  The criterion itself, kept for reference (`particuleC.cpp:1250-1273`), is
  evaluated per segment and per locus and returns 1 when the continuous
  approximation is acceptable:

  ```cpp
  if (seqlist[iseq].t1 < 0) return 1;              // last, infinite segment: always continuous
  nGen = seqlist[iseq].t1 - seqlist[iseq].t0;
  ra   = (double)nLineages / (double)seqlist[iseq].N;
  if      (nGen <=  30) OK = ra < (0.0031*nGen*nGen - 0.053*nGen + 0.7197);
  else if (nGen <= 100) OK = ra < (0.033*nGen + 1.7);
  else                  OK = ra < 0.5;
  ```

  Note `ra` is built on **`N`**, not on `Ne`, and that the `nGen <= 30` branch is
  quadratic, non-monotone and discontinuous with the next one. Only **bounded**
  segments can go discrete, so a single-population constant-size history never
  does.
  Don't reopen the "simplified substitution model" rumor; if a precise source
  turns up, check it against those line numbers first.
- **`LIK` pseudo-count `nal` too large (24/09)** — `_compute_LIK_for_one_locus`
  counted every allele appearing in `variant.alleles`, i.e. every state ever
  produced by a mutation, **including states no sample carries any more**
  (overwritten by a later mutation on the same lineage). `cal_lik2p` counts an
  allele only when its frequency summed over all samples is non-zero
  (`frt > 0.000001`), so our `b = 1/nal` was up to **2x too small**. Measured on
  `toy_example1_ms` G2: `nal` 13/6/14/12/18 against the C++'s 7/3/11/10/9.
  Invisible on every 2-sample dataset; surfaced only on a group with a 5x
  higher mutation rate (`MEANMU UN[5e-4,5e-3]` vs G1's `UN[1e-4,1e-3]`), where
  many more states are created then lost. The function's own docstring
  asserted the old behaviour was "exactly the C++'s `nal` since this project
  only has 2 populations" — that claim was wrong in a second way too: the
  `count_i ∪ count_j` union it described is a **no-op**, since
  `_length_by_sample` gives every population the same key list (taken from
  `variant.alleles`). Fixed by summing counts across all populations of
  `length_by_pop` and keeping only non-zero ones. Golden `LIK` values
  regenerated.
- **Systematic negative bias across datasets (25/09)** — a small shift of
  msprime below DIYABC, invisible to the per-column KS and visible only through
  a **sign test on `rdiff_mean`**. Measured on nine datasets: −2.65%
  (`toy_example5`, 70 loci) down to +0.06% (`human_seriel`, 5000 loci).
  **Verdict: it is the residual bias closed on 17/07, and it shrinks with loci
  count** — now measured single-variable on `toy_example4_seriel`, everything
  else held constant: 100 loci → sign test p = 2.6e-03, median −0.495%, KS
  4/133; 1000 loci → p = 0.185 (not significant), median −0.076%, KS 1/133.
  A factor 6.5 for 10x the loci, and doubly telling since more loci make the KS
  *more* sensitive, not less. Do not treat a sign-test bias on a small-loci
  dataset as a port defect; check the loci count first.
- **PoolSeq replay simulated twice the gene copies (25/09)** —
  `compute_summary_statistics_snp_from_values` built `counts_by_sample` **before**
  the IND/POOL branch, with the IndSeq formula
  (`count_per_samples.values()`, no `// 2`), and its PoolSeq branch used it as is.
  On a POOL file `count_per_samples` already counts **haploid gene copies**, so the
  replay simulated 400 copies per sample instead of 200 (1600 nodes / 800
  individuals instead of 800 / 400). Halving the sampling noise on every allele
  frequency lowered every differentiation statistic: a tight additive offset of
  −0.0025 on all six `FST2m` pairs, plus NEI, ML2p/ML3p, always in the same
  direction. **No guard could catch it**: `sum(counts) == ts.num_individuals`
  and `sum(layout) == ts.num_samples` both held — the inconsistency existed only
  against the observed file, which nothing re-reads at that point. Fixed by
  making the formula family-dependent. After the fix on `toy_example4_seriel`:
  `FST2m` offset −0.00254 → −0.00052 with the window now straddling zero,
  KS 10/133 → 4/133, sign-test median −1.06% → −0.495%. Found because the user
  asked why one branch computed what the other received from elsewhere — an
  asymmetry, not a failing test.
- **Text-reftable column order (21/09)** — `reftable.cpp::bintotxt` walks the
  header's **trailer line** (`entetehist`) and looks each name up *by name*,
  so a text reftable's column order is trailer order, **not** prior-declaration
  order. `_historical_columns_order` implements this; it raises with a
  symmetric difference on a trailer typo rather than silently shifting every
  value. The `.bin` format is different (scenario's own `histparam` order,
  constants excluded) and was left untouched.

- **Mutation parameters (`nparamut`) absent from our reftables (01/10)** —
  implemented and validated; was the long-standing Open work item. The values
  were drawn inside the worker, used for the mutation model and discarded.
  `ParticleResult.group_priors_values` now carries the **flat**
  `{column_name: float}` dict: the two SNP runners pass `{}` (a SNP header has no
  `group priors` section, and `group_prior_column_names` **raises** there rather
  than returning `[]` — `check_header_trailer_line` catches that `ValueError`
  deliberately), the four DNA/MicroSat ones fill it.
  Chain: `build_{group,microsat}_local_param_per_locus` return
  `(params_per_locus, values)`; the MicroSat path carries them up through
  `build_matrix_microsat_per_locus` (now a **triplet**), the DNA path through
  `build_rate_map_per_locus` — DNA has **two** branches and
  `build_matrix_dna_per_locus` is the one that drops the value.
  `_group_prior_columns` returns `(column, group, prior)` triplets and
  `group_prior_column_names` is now only its projection on the first element:
  **one authority decides which columns exist**, so writer and reader can no
  longer disagree on the count — which was the whole point, since
  `parse_real_reftable_params_with_group_priors` reads by position.
  **The `.bin` format carries no column names and no separate `nparamut` field**:
  `nparam[i]` absorbs both families (`reftable.cpp:64` checks
  `nparam[i] == nparamvar + nparamut`, and the write loop at `reftable.cpp:199`
  has a single counter). Writing an extra header int shifts the whole file — it
  was done once and caught by `test_write_reftable_bin_multi_scenario`.
  Validated by **strict identity** (`rdiff_mean = 0`, `ks_stat = 0`) on a replay
  of `toy_example1_ms` (6 columns, pure MicroSat) and
  `toy_example2_ms_dna_TN` (9 columns — the richest case: `k1` *and* `k2`, plus
  the inter-type ordering). The *drawing* path cannot be validated by identity and
  rests on a unit test only.
  Two traps met on the way, both worth knowing: `group_priors` is the **parsed
  declaration** (lists of `GroupPrior`) while `values`/`group_priors_values` are
  the **drawn values** (nested `{group: {prior: float}}`) — both live at once
  inside `build_rate_map_per_locus`; and `draw_group_parameter_values` does **not**
  filter by type, so a DNA run also carries a mixed dataset's MicroSat columns,
  with identical values (same seed). That is correct, not a bug.

- **A mixed MicroSat + DNA header produced only half its statistics (02/10)** —
  implemented and validated. The two stat paths are **structurally disjoint**:
  `compute_all_statistics_dna` skips every `ms_or_seq != "S"` locus
  (`summary_statistics.py:3994`), `compute_all_statistics_microsat` every `!= "M"`
  (`:4096`), and no runner called both. Measured on `toy_example2_ms_dna_TN`: the
  header declares **58** statistics (16 on G1 microsat, 42 on G2+G3 DNA), a
  MicroSat run produced 16, a DNA run 42, **never 58**. So the `.bin` was not
  structurally equivalent on this dataset family, and every past validation of
  `toy_example2_ms_dna_*` covered a subset — true separately, false in
  conjunction. `stats_filter="HEADER"` was simply unusable there (it raises on the
  42 missing).
  **Five properties made the merge possible, all measured before writing code**:
  the two families' column names never collide; 16 + 42 = exactly the 58 the
  header asks for; and at equal seed both paths draw **identical** historical and
  group parameters. Above all, the four simulation loops enumerate the **full**
  `list_loci` and skip with `continue` (`ancestry_simulation.py:2130, 2544, 2973,
  3244`), so `i` is the locus index in the header — a `[M]` and an `[S]` locus can
  never share it, hence **no seed collision between families** and values
  identical to what the separate runs produced.
  Shape: `compute_summary_statistics_mixed{,_from_values}` run steps 1-2-3
  (demography, scenario, sample sets) **once**, then 4-5-6 per family, then
  `_filter_statistics` **once** on the merged dict — that final filter is what
  guarantees header order and turns any missing statistic into a `ValueError`.
  The single draw is not a performance matter (the duplicate cost 0.12 ms against
  160 ms per particle, 0.075 %) but a **correctness** one: both halves of a
  reftable row now describe the same demographic history by construction, not by
  reproducible coincidence.
  Four helpers, drawing and replay kept as **siblings** per the project rule — a
  `from_values: bool` flag was tried and hid three divergences at once (the
  function called, its argument list, and its return arity); only the first had
  been wired.
  Validated by replay, 2 independent DIYABC runs × 2 scenarios: **73 columns**
  (1 + 5 historical + 9 mutation + 58 statistics), historical and mutation columns
  at `rdiff = 0` (the pairing proof), and 10/232 KS flags (4.3 % against 5 %
  expected). **9 of the 10 are in G2**, all in `MPD`/`VPD`/`DTA`/`MNS`, and the
  only two passing the persistence criterion (`VPD_2_1`, `MPD_2_1`, flagged in
  both replays) belong to that family — i.e. the documented G2 residual, see Open
  work. The lone non-G2 flag, `MGW_1_1`, appears in 1 cell of 4: noise. G1 being
  clean at 15/16 is the point, since the MicroSat half is the one whose call path
  changed.
  **`values` fell out of the replay chain four times** while building this (the
  drawing helper used instead of the replay one, `values` overwritten by an
  unpack, then absent from the signature altogether). It is the one datum the
  replay variant *receives* and the drawing variant *produces*, so copying the
  drawing shape drops it every time. Only three of the four failures were loud.

- **`DTA_2_2` residual (05/10)** — **correlated noise, not a formula or simulator
  defect.** KS-flagged in 11 cells of 12 on `toy_example2_ms_dna_50loci_{JK,K2P,TN}`
  (p ~ 1e-4 to 1e-2; the 12th at 0.16), `DTA_2_1` in 2 of 12 at the margin, ours
  above DIYABC's mean in 12 of 12. Cause: the
  additive per-locus seeds above, replayed identically in every file (seeds
  1..1000). Fixed by `_locus_seed`; **validated: 3 models × 2 replays, every
  statistic without a KS flag** (as reported by the user — the per-`θ`-quintile
  ratio and the `MPD <= 0.75*NSS` violation rate were not re-measured separately).
  How it was established, by elimination — keep the order, it is the method:
  1. `cal_dta1pl` and `_tajima_d_per_locus` are identical line by line (constants,
     `dnavar < 1` → 0 and counted, `n < 2` → excluded). The `.mss` has no missing
     data, and both samples have n = 20. The `MPD` denominator bug cannot reach
     `DTA`, which counts monomorphic loci with D = 0.
  2. The two samples are **exchangeable** in scenario 2 (same `N1`, symmetric merge),
     so `DTA_2_2 − DTA_2_1` must be 0 on average. In DIYABC it is (−0.011..+0.004,
     mixed signs); in ours it was **+0.011..+0.021 in 6 files of 6**, while
     `NSS`/`MPD` stayed symmetric. 2600 particles: +0.0132, t = 4.2.
  3. tskit's `Tajimas_D` gives the same asymmetry (+0.0131) → our D computation is
     exonerated. Tree only, no mutation, 300 000 replicates at fixed parameters:
     branch π 20038.4 vs 20037.1 (t = −0.12) → genealogy symmetric. Mutation at
     fixed parameters, 3 rates × 100 000 replicates: `D(2−1)` t = −0.60 / +1.46 /
     −0.50 → mutation symmetric. So the defect was not in any component tested alone.
  4. Seeds spaced by 1000 (no overlap), 2000 particles: `DTA_2_2 − DTA_2_1` =
     **−0.0068** (t = −2.15) — the sign flipped, and both means fall inside DIYABC's
     range. Windows 1-600 / 601-2600 / spaced gave +0.011 / +0.013 / −0.007: not a
     stable effect, and the standard error assumed independence that did not exist.
  **Traps met on the way.** (a) Swapping the two samples' layouts is **vacuous**:
  every statistic is a function of a sample's node set, so the swap exchanges the
  columns exactly and flips the sign whatever the cause. (b) A large `rdiff_mean` on
  `DTA` proves nothing: means 0.005–0.04 against a standard deviation of ~0.12 (the
  near-zero-mean rule again). (c) `_length_by_sample_and_individuals` duplicates a
  haploid's single allele into `(x, x)` by design, so its counts exceed
  `_length_by_sample`'s by one; compare against `_genotypes_by_sample_and_individuals`,
  which keeps true ploidy. (d) An unexplained `TSK_ERR_BAD_OFFSET` from
  `sim_mutations` occurred once in 2000 particles and never again on the same seeds;
  not understood, not reproduced.
  Consequences: the golden values moved twice on 05/10 (33 tests with `_locus_seed`
  on the per-locus sites, then 38 when the conversion was extended to every derived
  seed); each hand-derived expectation was recomputed from the new allele counts, not
  pasted. **The `te2` fixture no longer exercises the `num_sites == 0` guard**: with
  the final seeds G2 and G3 both have their 5 loci polymorphic (it had 1, then 2, loci
  without mutation under the earlier seeds, hence the 5/4 and 5/3 factors that the
  `MPD`/`VPD` test comments used to quote). The guard now has its own seed-independent
  test, `test_mean_pairwise_differences_excludes_loci_without_mutation`, which fails
  when the guard is removed (checked). Any stored comparison or reftable produced
  before 05/10 used the old seeds and is no longer comparable to a fresh run.

- **Serial IndSeq replay broken from 25/09 to 06/10 (`KeyError: Population with name
  'pop2' not found`).** `simulate_particle_genotypes_from_values`, the replay sibling
  of `simulate_particle_genotypes`, no longer built its `SampleSet` list and per-sample
  counts from the scenario: lost at commit `392d193` (25/09 17:25, "symétrisation
  INDSEQ POOLSEQ"), forty minutes after `851ba78` where the IND branch still had them.
  On a serial dataset msprime then received the `pop1..popN` names of the `.snp`'s POP
  blocks against a one-population demography. Non-serial datasets were unaffected
  (`sample_sets=None` falls back to the names), which is why nothing failed. Restored by
  copying the drawing sibling's block; guarded by
  `test_simulate_particle_genotypes_from_values_serial_indseq` (synthetic 1-population,
  4-date context, both MAF paths; fails without the block, checked). **The
  `human_seriel` validation predated this loss, so it described a code path that was
  broken for ten days; it was rerun on 06/10 on the real dataset (5000 loci, 2
  scenarios: 0 KS flag after the 6-decimal rounding, sign-test median +0.05 % /
  −0.05 %).** Lesson: a refactor that
  touches one half of a drawing/replay pair must be checked against the other half, and
  a feature validated by replay needs a test that replays it.

## Domain knowledge

### SNP

Genotypes are simulated per locus and mutated with exactly one Hudson
mutation. `<M>` and `<Y>` loci share a single genealogy across all loci of
their type (`simulate_shared_ancestry_loci`, `particuleC.cpp:2422-2435`
`GeneTreeM`/`GeneTreeY`). Non-`<A>` heritage types get a rescaled demography
(`rescale_demography(coalescence_coefficient(heritage, sex_ratio) / 2)`) and
`ploidy=1`. PoolSeq reads go through `with_mrc_filter` against a global pool.

### MicroSat

- **Header parsing.** The `MEANxxx` → `GAMxxx` hierarchy is **positional, not
  name-based**: `header.cpp::readHeadersimGroupPrior` reads a fixed sequence of
  `getline` calls per group and never inspects the name token, so the text
  `MEANMU`/`GAMMU` is purely cosmetic. A `GAMxxx` line's declared "mean" is not
  a real bound — it is unconditionally overwritten at draw time with the value
  just drawn from the group's own `MEANxxx` prior (`particuleC.cpp:819`).
- **GSM mutation model.** At each mutation event a step of `d` repeat units is
  drawn geometrically from `Pgeom` (`Pgeom=0` is the SMM special case), clamped
  to `[kmin, kmax]`. Ancestral state = midpoint of `[kmin, kmax]`.
  `msprime.TPM` with `p→ε` reproduces this channel exactly, with
  `m = 1 - Pgeom` (both `Pgeom=0` and `Pgeom=1` are valid DIYABC inputs that
  map to msprime's forbidden `m=1`/`m=0` literals, hence the epsilon clamp).
  **The TPM grid must be anchored on `root`, not on `kmin`** — the ancestral
  state does not necessarily fall on a `motif_size`-multiple offset from
  `kmin`, and GSM steps happen from the current state. Accepted consequence:
  TPM's bounds can differ from DIYABC's literal `kmin`/`kmax` by up to
  `motif_size - 1` bp.
- **SNI channel.** A single combined Poisson process at `mut_rate + sni_rate`;
  each event is classified by a Bernoulli draw
  `p_sni = sni_rate / (sni_rate + mut_rate)` into a GSM step or a ±1 bp SNI
  step. Because SNI shifts the allele into *any* residue class modulo
  `motif_size`, the transition matrix needs a **dense** grid (one row per
  integer in `[kmin, kmax]`), unlike GSM alone. Performance: build exactly
  `motif_size` shared `msprime.TPM` matrices (one per residue class — the
  state count is constant within a class but differs across classes), not one
  per dense state; the naive version cost ~34x, this one ~1.5x.
  Validating `sni_rate=0` against the GSM-only matrix is **not** a same-shape
  comparison — extract the dense rows/columns matching
  `build_transition_matrix_microsat(...).alleles` first.
- **Statistics.** Two distinct data representations are required, not one:
  `_length_by_sample_and_individuals` (for `FST`) returns a 2-tuple per
  individual, duplicating a haploid's single allele into `(taille, taille)` —
  safe because `cal_Fst2p` itself treats a haploid copy as a doubled
  diploid-like pair, making the ANOVA ploidy-invariant. `_genotypes_by_pop_and_
  individuals` (for `LIK`) keeps **true** ploidy, because `cal_lik2p` applies a
  genuinely different formula per ploidy. Every other stat uses the flat
  `_length_by_sample` `(taille, compte)` table.
  `LIK` is the only **asymmetric** stat — it iterates all ordered pairs
  `i != j` (the header declares `LIK 1.2 2.1` explicitly).
  Aggregation falls into three families, easy to conflate: mean of per-locus
  values with a per-population or per-pair valid-loci denominator; a **ratio of
  sums** with no per-locus averaging (`MGW`, `DAS`, `FST`); and `DM2`'s
  stateful variant. Picking the wrong one produces a plausible but wrong
  number with no exception.
  A locus with `num_sites == 0` (zero mutations drawn) is a real, fully
  monomorphic locus, not an error — tskit represents it as no site at all, so
  all three helpers special-case it with an arbitrary placeholder value.
- **`AML` (`cal_Aml3p`).** Bisection over a mixture proportion `a`, driven by
  the **sign** of the finite-difference slope at both bounds, with three
  boundary cases (flat → uniform random `a`; always negative → `a=0`; always
  positive → `a=1`). No pseudo-count, unlike `LIK` — a zero-frequency term is
  skipped. Frequencies must be hoisted out of the bisection loop
  (`_prepare_loci_for_admixture`): ~20 evaluations per triplet, measured 15.2s
  → 0.77s per particle. Seeds must vary per triplet **and** per group.

### DNA sequences

- **Base frequencies come from the `.mss` file, not the header.**
  `DataC::do_sequence` computes `pi_A/C/G/T` empirically per locus, pooling all
  populations and excluding missing (`-`/`N`) positions from both numerator and
  denominator; every particle then reuses that fixed value. The literal `pi`
  values in the header belong to a no-observed-data launch mode this project
  doesn't use.
- **`p_fixe`/`gams` are `(proportion of invariant sites, gamma shape)`** —
  completely unrelated to `k1`/`k2`, which come from separate `MEANK1`/`GAMK1`/
  `MEANK2`/`GAMK2` group priors present for *every* `[S]` group regardless of
  model.
- **`k1`/`k2`/`mus_rate` are drawn hierarchically, two tiers.** The per-locus
  tier is used only if the `GAMxxx`'s `sdshape > 0.001` **and** the group's
  `nloc > 1`; `mus_rate` and `k1` have the `nloc > 1` check, `k2` genuinely
  does not — an asymmetry in DIYABC's own source, replicated via `check_nloc`.
- **`gams == 0` is a valid, non-error input** to `sample_site_rates`: `MwcGen::ggamma3` (`randomgenerator.cpp:199-204`) returns `mean` (`1.0`) directly when `shape == 0.0` rather than dividing by it. Raising `ZeroDivisionError` there would be wrong behaviour, not an unhandled edge case.
- **`dnatrue` selects two different sequence representations, and the reftable
  path uses the one nobody checks.** With `dnatrue = true` (how `statobsRF.txt` is
  computed, `header.cpp:1930`) a node's sequence is `dnalength` real bases. With
  **`dnatrue = false`, which is the reftable/`-R` path**, `init_dnaseq` builds an
  *artificial sequence made only of variable sites*: length `k` = the number of
  distinct physical sites hit by a mutation, `tabsit` mapping physical site →
  compressed index, and `sitmut[i] = tabsit[sitmut2[i]]`. The representation is
  faithful (several mutations at one physical site share an index, so they
  overwrite along the tree as they should) and gives identical π and `S` — **but
  it behaves differently when a locus has no mutation at all**: `dna` stays empty,
  and the empty string is `SEQMISSING`. Consequence to keep in mind generally:
  **a verification run on observed data does not validate the simulated path**,
  since the two take different branches of `init_dnaseq` and of `cal_numvar`. That
  is how the G2 `MPD` bug survived a 6-decimal check (see Open work).
- **DIYABC's Jukes-Cantor token is `JK`**, not `JC`/`JC69`. Got this wrong
  twice.
- `build_transition_matrix_dna`'s row-sum normalization needs
  `axis=1, keepdims=True` — omitting `keepdims` divides along the wrong numpy
  axis and yields rows that don't sum to 1, with no error raised.
- **Mutation placement is delegated to `msprime.sim_mutations`** with a generic
  `MatrixMutationModel` fed our own `matQ`/`pi` — **deliberately not** msprime's
  named models (`JC69`/`HKY`/`F84`/`GTR`), whose rate-scaling convention allows
  and counts silent self-transitions (non-zero diagonal), unlike `comp_matQ`'s
  (each row sums to 1, zero diagonal, every event a real substitution).
- DNA loci call `msprime.sim_ancestry` **directly**, not
  `simulate_independent_loci` (which hardcodes `sequence_length=1`, wrong once
  `dnalength` varies per locus). A MicroSat locus, conversely, *is*
  `sequence_length=1` — a single site with a ~39-state alphabet, not many sites.
- `<M>` **and** `<Y>` loci share one genealogy across all loci of their type,
  each via its **own dedicated** seed offset (`_SHARED_M_ANCESTRY_SEED_OFFSET` /
  `_SHARED_Y_ANCESTRY_SEED_OFFSET`) — sharing one offset would make every `<Y>`
  locus share the `<M>` tree too.
- **`<X>`/`<Y>` sex inference**: `.mss` has no SEX column; DIYABC infers sex
  from the genotype's own ploidy *at that locus* — everyone starts female and
  flips to male when an `<X>` locus shows a haploid genotype (unconditionally,
  even for the missing-data placeholder) or a `<Y>` locus shows a present
  genotype. The missing-sequence token is `<[]>` (empty between brackets), so
  the regex must be `\S*`, not `\S+`. Sample sets must therefore be dispatched
  **inside** the per-locus loop for `<X>`/`<Y>`.
- **Statistics.** Two-tier pattern: a stateless per-locus brick plus an
  aggregator over one header group's loci (never mixing groups). Three
  denominator regimes coexist: flat `num_loci` (`PSS`, `MNS`, `VNS`),
  per-population valid-loci (`MPD`, `VPD`), per-pair valid-loci. `VNS` is a
  **biased** variance (`ddof=0`), deliberately different from `VPD`'s `ddof=1`
  — easy to get wrong by pattern-matching. `HST` and `MP2` are **ratios of
  sums** accumulated per pair, not means of per-locus ratios; their
  accumulators must be **dicts keyed by pair**, never shared scalars (a shared
  scalar is invisible on 2-population data and wrong on 3+).
- **Column naming**: whether a stat's column gets a `_<group>_` segment
  (`multi_group`) depends on the number of distinct groups in the **whole**
  header, any type — not the number of groups of the stat's own type.

### Serial/temporal sampling

A scenario may sample the **same** population at several dates
(`reference/toy_example1_ms`: `0 sample 1`, `50 sample 1`, `200 sample 1`,
`500 sample 1`, 4 `POP` blocks of 20 individuals, statistics declared on
indices 1..4). The underlying gap was never "temporal sampling" — it is that
**a sample is not a population**, and serial data is simply the first case
where the two differ:

- `history.cpp::read_events` counts `nsamp` separately from `npop` (`nn0`) and
  assigns `event[i].sample = ++nsamp` in file order; nothing ties them.
- `data.cpp` speaks only of samples (`nsample`, `samplesize[ech]`,
  `ssize[locustype][sa]`). Statistics indices are **sample** indices.
- `particuleC.cpp:1185/1213/1528` reads a sample's size from the observed data
  and tags each node `gt.nodes[i].sample = sa + 1`; a SAMPLE event activates
  exactly the nodes carrying its own index. The data block `sa` therefore maps
  to the `sa+1`-th `sample` line **positionally**, with no name cross-reference
  — and `buildSuperScen` (`header.cpp:1057`) only takes maxima, so DIYABC
  guarantees nothing about this ordering *across* scenarios.

**THE INVARIANT the whole feature rests on** — node order = `SampleSet` order =
`sample`-event order = `POP`-block order = `samples_default` key order =
statistics column indices. Verified empirically: msprime allocates sample node
IDs contiguously in the order of the `SampleSet` list, individuals too, and
this holds under heterogeneous ploidy.

Two bricks and a wiring:

- **`build_sample_sets_from_scenario(scenario, values, )`**
  (`ancestry_simulation.py`) — one `msprime.SampleSet(n, population, time)` per
  `sample` event, in order. ``'s **keys are ignored**: only the
  order of its values matters, the k-th count going to the k-th event. Indexing
  it by `f"pop{event.pop}"` returns the same count for every serial sample —
  that bug was written, caught by a deliberately unequal-count test, and is now
  locked by it. A sample time may be a parameter name, not a literal
  (`particuleC.cpp:599-605`), so this must be called **after** the prior draw;
  `build_demography` already evaluates those expressions (and discards them),
  which is why no constant-prior failure mode is introduced here. Guards on
  `len(sample_events) != len(counts)` — DIYABC does not, and a silent
  misalignment is the worst possible outcome.
- **`compute_sample_layout(ts, )`** — the sample-aware twin of
  `compute_population_layout`, same return shape so the two are
  interchangeable. Slices `ts.individuals()` by cumulative counts and
  concatenates their `.nodes`; slicing `ts.samples()` by `count × ploidy`
  instead would break on `<X>`, whose sample is described by **two**
  `SampleSet`s (females `ploidy=2`, males `ploidy=1`). On a non-serial dataset
  it returns exactly what `compute_population_layout` returns — that equality
  is a test, and it is what makes the substitution provably neutral on every
  already-validated dataset.
- **Wiring**: the 4 layout-consuming helpers of `summary_statistics.py` take an
  optional `layout=`; the 27 stat functions that call them take a parallel
  `layouts=` list (`zip(..., strict=True)` everywhere); `compute_all_statistics_
  dna`/`_microsat` take `layouts_by_locus`, a dict **keyed by locus name** —
  never a flat list, because the dispatch filters and reorders loci per group
  and two divergent comprehensions would misalign silently. `pipeline.py`
  builds both the `SampleSet` list and the layouts, per particle, from the
  scenario it re-parses plus `counts_by_sample_for_locus` (`samples_default` for
  `<A>/<H>/<M>/<X>`, male counts for `<Y>`).

**The SNP paths are wired too** (25/09), IndSeq and PoolSeq. PoolSeq follows
the same shape with two specifics: `context.count_per_samples` gives the **haploid**
pool size (gene copies, not individuals), hence `poolseq_counts_by_sample`
and its `// 2` — a single definition used by both `pipeline` and
`simulate_poolseq_reads_with_mrc_filter`, because the two candidate dicts have
the SAME total and `_check_layout_matches` cannot tell them apart. And
`pool_sizes`, consumed by `compute_all_statistics_poolseq` for the read-bias
correction, stays in **haploid** units — the two scales coexist deliberately.
Both MRC entry points (the `mrc <= 0` fast path and the rejection loop) forward
the counts; `observed_reads_per_locus[locus_index : locus_index + 1]` pins the
observed coverage to its locus, so a rejected attempt never consumes one.

The SNP shape differs from MicroSat/DNA in one way worth
knowing: the 130 SNP statistics take `genotypes_per_locus`, a list of
`{name: [genotype]}` dicts — the layout is consumed **once**, inside
`simulate_snp_genotypes`, to build those dicts. There is therefore **no
`layouts=` threading through the stat functions** as on the MicroSat/DNA side;
everything happens in `ancestry_simulation.py`. What travels down is
`` (an ordered `{name: individual count}` dict), not a
precomputed layout: the layout depends on **ploidy**, which varies per heritage
type (`<A>` 2, `<H>`/`<M>` 1), so a single layout built upstream would be wrong
for every type but one. Each MAF loop derives its own from the counts, where
the `TreeSequence` is born and the ploidy is known. Both MAF loops have a
`maf == 0.0` fast path that must forward the counts too — `<MAF=hudson>`
datasets take *only* that path.

**Never derive `` from `samples`.** They are different objects:
`samples` goes to msprime (`dict[str, int]` *or* `list[SampleSet]`),
`` feeds `compute_sample_layout`. A `SampleSet` list carries no
names, and an `<X>` locus describes **one** sample with **two** `SampleSet`s.

**Deliberately not done, now guarded**: `<X>`/`<Y>` loci override the
sample sets inside the per-locus loop with the sex-stratified builders, which
are **not** serial-aware — on both the MicroSat/DNA and SNP sides. Those four
builders rebuild their msprime keys from the **position of the observed block**
(`f"pop{i}"`, `i` = i-th `POP` block; SNP side via `sample_index_to_name`,
MicroSat/DNA side inside `individual_sexes_from_locus_genotype`) instead of the
population its `sample` event names. On a serial scenario `i` runs past the
number of sampled populations.

Why that had to raise rather than let msprime complain — the two outcomes are
wildly unequal and **the dangerous one is reachable**. Measured on
`toy_example2_ms_dna_K2P`, whose scenario 1 declares 5 populations for 2 samples and
whose `build_demography` creates `pop1..pop5` all `initially_active=True`:
`samples={"pop1": 5, "pop2": 5, "pop3": 5}` is accepted with **no error at all**
(15 nodes instead of 10), the samples silently attached to an
ancestral, never-sampled population. Only once the invented index exceeds the
maximum does msprime raise a loud `KeyError` (`pop6`).

`reftable_loop.raise_if_serial_with_sex_linked_loci` refuses the combination at
run setup — once per run, not per particle, so the failure precedes particle 0
rather than hitting particle 743 — and raises as soon as **any** scenario of the
header is serial, without waiting to see which one `draw_scenario` picks. The
predicate is `scenario_parser.is_serial_scenario`, i.e.
`len(sample_events) > len({e.pop for e in sample_events})`: "the same population
is sampled twice", which is the definition of serial. **Do not rephrase it as a
comparison with `npop`** — a scenario may declare never-sampled populations
(`toy_example2_ms_dna_K2P`: 5 for 2), and that comparison would misclassify it.
No reference dataset triggers the guard: the serial ones are `<A>`/`<M>`, the
`<X>`/`<Y>` ones are not serial — which is what makes adding it provably neutral
on everything already validated.

Related hazard, measured:
a layout whose total length disagrees with `ts.num_samples` produces silently
wrong genotypes rather than an error, because `simulate_snp_genotypes` does a
membership test (`s in derived_samples`), never an index — a non-existent node
id simply yields `0`.

### Replay (`_from_values`) chain

Every `_from_values` function is a **sibling** of its drawing counterpart,
never a modification of it — the originals are validated against real
reftables and must not be touched. Only the **group-level** (tier 1) draw is
replaced with the real DIYABC value; the per-locus (tier 2) dispersion around
that mean is never replaced, because DIYABC doesn't record it in the reftable.
`group_prior_column_names` is scenario-**independent** (`nparamut` is constant
across scenarios, unlike `nparam`) — do not reuse
`_kept_param_names_by_scenario` for it.

## Deliberately reproduced C++ bugs

These are **not** to be "fixed" — the goal is fidelity to real DIYABC output,
not an idealized reimplementation. All three are judged unintentional, kept
because DIYABC's own output depends on them.

- **`MPD`/`VPD` denominator** (`particuleC.cpp:1777`, `sumstat.cpp::cal_mpd1p` /
  `cal_vpd1p`) — added 02/10. DIYABC's `init_dnaseq` leaves a locus's sequence
  as the empty string when the simulation drew no mutation and `dnatrue` is
  false (the reftable path); the empty string is its missing-data marker
  `SEQMISSING`; so `cal_mpd1p` (`if (nd > 0)`) and `cal_vpd1p` (`if (nd > 1)`)
  drop the locus from their denominator while `cal_nss1p`, whose `OK` comes from
  the *observed* `samplesize()`, keeps it. Reproduced by the
  `if ts.num_sites == 0: continue` guard in `compute_MPD`/`compute_VPD`, which is
  **commented as a deliberate bug with the lines to delete**. Full investigation
  under "G2 `MPD`/`VPD` residual" in Open work.

- **`mutsit` / `sitefix`** (`header.cpp:727-738`): DIYABC draws random distinct
  site indices to choose the invariant sites, then zeroes `mutsit` using the
  **loop counter** instead of `sitefix[i]`. In practice the invariant sites are
  always the *first* `dnalength - nsv` sites. The "intended" version is kept
  commented out in `sample_site_rates`.
- **`DM2` / `cal_dmu2p`**: `moy[]` is recomputed only when both populations
  have samples at a locus, but the line consuming it sits **outside** that
  guard — so on a locus where one population is absent, the stale `moy[]` from
  the previous locus is reused, divided by the current locus's `motif_size`,
  while the denominator `nl` does not count that locus. Reproduced by threading
  an explicit `previous_moy` between locus calls, which makes `DM2` the one
  stat where the order of `tree_sequences` matters.

## Design decisions

- **Sibling duplication over parameterization.** The SNP → DNA → MicroSat
  function families are deliberately kept as near-identical siblings rather
  than factored into one parameterized function. Reaffirmed twice with the
  user; the rule is "never modify an already-validated original". Don't
  refactor these without asking.
- **No dedup across stat groups** in `stats_group_parser.py` — it is still
  unclear whether two groups could legitimately declare the same column name
  for two different underlying statistics, and collapsing via `set` would
  silently drop one.
- **Review-checklist patterns** (numpy bool addition, seed reuse, pairwise
  accumulators, duplicate `def` names, name shadowing, signature refactors)
  live in the persistent memory files, not here.

## Open work

- **Serial sampling is not wired for `<X>`/`<Y>` — now guarded, not fixed.**
  The per-locus dispatch still overrides the serial `SampleSet`s, on both the
  SNP and MicroSat/DNA sides. Since 28/09 the combination is refused up front by
  `reftable_loop.raise_if_serial_with_sex_linked_loci`, called from all seven
  run/replay entry points; the predicate is
  `scenario_parser.is_serial_scenario`. Implementing it for real was
  deliberately deferred: **no dataset would validate it**, and this project
  never declares a path correct without a paired replay against the real
  DIYABC. See "Serial/temporal sampling" under Domain knowledge.
- **`ML3p_2.3.4` on `toy_example4_seriel`: RESOLVED 05/10 — an artefact of
  overlapping seeds, not a real residual (verdict at the end of this entry).**
  Significant in **three** replays (p = 0.0001 / 0.0013 / 0.0002) with a stable
  −10 to −14% rdiff, and it survives a 10x increase in loci count that wipes out
  every other flagged column. The persistence criterion was met at the time (the
  05/10 verdict below shows why that was not enough). `ML2p_2.4` is intermediate: significant in
  both 100-loci replays, not at 1000, amplitude decaying — the loci-count
  residual, with a consistently negative sign worth noting. Left unexplored by
  the user's call (cost vs. stakes) until 05/10.
  **Localised 30/09: it is the oldest SAMPLE, not a population, and not
  admixture.** Splitting the whole ML family by index on run 0 scenario 1 gives a
  perfect separation — every column carrying index **4** diverges, no other one
  does:

  | clean | rdiff | KS p | | affected | rdiff | KS p |
  |---|---|---|---|---|---|---|
  | `ML1p_1` | −0.30% | 1.00 | | `ML1p_4` | −3.58% | 0.107 |
  | `ML1p_2` | +0.19% | 0.97 | | `ML2p_1.4` | −10.48% | 0.067 |
  | `ML1p_3` | +0.49% | 1.00 | | `ML2p_2.4` | −9.27% | 0.009 |
  | `ML2p_1.2` | +0.08% | 1.00 | | `ML2p_3.4` | −3.28% | 0.246 |
  | `ML2p_1.3` | +0.09% | 1.00 | | `ML3p_1.2.4` | −13.68% | 0.011 |
  | `ML2p_2.3` | +0.20% | 0.98 | | `ML3p_1.3.4` | −22.06% | 0.029 |
  | `ML3p_1.2.3` | +0.28% | 1.00 | | `ML3p_2.3.4` | −14.32% | 0.0001 |

  7/7 affected, all negative, against 0.96–1.00 KS p for the other five. **This
  is not a small-mean artifact**: the KS is scale-invariant, so the separation
  cannot come from the index-4 columns having smaller means (0.008–0.083 vs
  0.24–0.44). The earlier hypothesis — a link to `toy_example4`'s observed data
  being built under scenario 3 with pop4 admixed — **is falsified**: this header
  declares **one** population sampled at four dates (`0/50/200/500 sample 1`), so
  there is no population 4. Index 4 is the **oldest sample**, at t = 500.
  The number of populations in the statistic is irrelevant: `ML3p_1.2.3` is clean
  at +0.28% while `ML1p_4` is not. So this is a **serial-sampling** lead, not an
  admixture or ML-combinatorics one, and the discrete coalescent is already
  falsified for it (see the `Npresent` stratification below).
  **The discrete coalescent is NOT the cause — falsified 30/09.** Stratifying the
  paired replay by the drawn `Npresent` (3 runs, including the 1000-loci one)
  gives the *opposite* gradient to the one that signature produces. `ML3p_2.3.4`
  relative error:

  | `Npresent` | <150 | 150–350 | 350–600 | 600–1000 |
  |---|---|---|---|---|
  | run 0 | **+8.1%** | −45.8% | −30.0% | −32.3% |
  | run 1 | **−1.7%** | −31.8% | −31.3% | −59.4% |
  | run 0, 1000 loci | **−3.6%** | −21.8% | −20.7% | −16.3% |

  Every particle here is already in discrete mode (800 gene copies with
  `Npresent <= 1000`, so `ra >= 0.8`), and within that regime a smaller `N`
  widens the discrete/continuous gap — yet **that is exactly where the two
  simulators agree best**. On `toy_example2_ms_dna_50loci_K2P` the `N1`
  gradient runs the other way (π error −72% at `N1 < 160` down to −3% above 4000)
  — and that one turned out not to be a coalescent effect at all, see the G2 entry
  under Open work. `FST2m_2.4` stays within ±4% in every stratum, so the pairing
  and the stratification are sound.
  Read the `<150` column as saturation, not accuracy: tiny `N` collapses
  diversity and drives a likelihood statistic to a boundary value on both sides.
  The falsification was later confirmed by a second route: `Npresent` is the
  wrong variable, since sample 4 (t = 500) lives in the **`Npast`** regime
  whenever `tbn < 500`. Stratifying by `Npast` *conditional* on `tbn < 500`, the
  residual is **worst where `Npast` is largest** (`ML3p_2.3.4` −70.0% / −56.8% /
  −32.6% above 20000, against +8.0% / −61.3% / −12.3% in 1600–6000) — whereas
  discrete mode needs `ra = 800/N >= 0.5`, i.e. `Npast < 1600`, a stratum holding
  only 6 particles. The residual lives where **both** simulators are continuous.

  **The failing configuration, specified** (3 runs, scenario 1). Splitting by
  `tbn`, the residual collapses once the size change is older than the oldest
  sample — `ML3p_2.3.4`:

  | `tbn` | run 0 | run 1 | run 0, 1000 loci |
  |---|---|---|---|
  | <200 | −47.1% | −51.9% | −20.0% |
  | 200–500 | −41.8% | −57.5% | −24.6% |
  | **>500** | **−2.2%** | **−8.9%** | **−7.3%** |

  So it takes all three of: **(a)** `tbn < 500`, i.e. the oldest sample drawn
  beyond the `VarNe`, inside `Npast`; **(b)** `Npast` large — with
  `Npresent <= 1000` and `Npast ~ UN[10,50000]` that is a 20–50x expansion into
  the past, i.e. a forward-time bottleneck; **(c)** the ML family, since
  `FST2m_2.4` stays within ±4% in every stratum, so sample 4's allele frequencies
  are not globally wrong. `ML1p_1` (the present-day sample) is clean in every
  stratum of every split, at |rdiff| < 1.7% and KS p > 0.85.
  In one sentence: **a sample taken from an ancestral population much larger than
  the recent one.** First steps if anyone picks this up: re-read `cal_ml3p` in
  `sumstat.cpp`, and compare a single particle at fixed `tbn ≈ 100`,
  `Npast ≈ 40000`, `Npresent ≈ 500` — the allele-frequency spectrum of sample 4
  is where the two sides must be made to disagree visibly.
  **Verdict 05/10: the residual does not survive decorrelated seeds.** The three
  replays that made it "persistent" all used seeds 1..N on the PoolSeq path, whose
  `with_mrc_filter` seeds were additive (`seed + batch_index * 20` for the tree
  batches, `seed + attempt + _MRC_REJECTION_SEED_OFFSET` for the reads): the same
  mechanism as `DTA_2_2`, see Closed investigations. Both sites now go through
  `_locus_seed` (`_MRC_BATCH_SEED_OFFSET` added to `configuration.py`). Replayed on
  **two independent DIYABC runs at 1000 loci**, two scenarios each: `ML3p_2.3.4`
  KS p = **0.996** (run 1, scenario 1) and 0.30 (scenario 2), no ML column among the
  flags of run 0, and no column carrying the oldest sample flagged in any of the four
  cells. Flag counts 0 / 6 / 0 / 2 of 120 (about 6 expected by chance), with **no
  column in common** between the two runs. Read the localisation and the `tbn` /
  `Npast` / `Npresent` stratifications above as the record of what overlapping
  seeds produced; they are **not** evidence of a coalescent or `cal_ml3p` defect, and
  the "first steps if anyone picks this up" are void.
- **G2 `MPD`/`VPD` residual — CAUSE FOUND, LOCATED TO THE LINE, AND REPRODUCED
  IN `compute_MPD`/`compute_VPD` (02/10). Replay validation done 05/10.**
  In one sentence: **`cal_nss1p`'s per-locus denominator is derived from the
  OBSERVED data, `cal_mpd1p`'s from the SIMULATED data**, and on a locus where
  the simulation drew no mutation the two disagree.

  **The chain, link by link.**
  1. `init_dnaseq` (`particuleC.cpp:1777`) starts with `string dna = "";` and
     ends `if (nmutot > 0) {…} else { if (dnatrue) {…fills dnalength bases…} }`.
     The reftable path runs with **`dnatrue = false`** (see the `dnatrue` entry
     under Domain knowledge), so on a locus with **no mutation** the `else`
     branch does nothing and `dna` stays **empty**.
  2. `cree_haplo` copies that empty string into all of `haplodna[sa][ind]`.
  3. **`#define SEQMISSING ""`** (`particuleC.hpp:17`) — the empty string *is*
     the missing-data marker. The whole locus becomes "missing".
  4. `cal_mpdpl` gates its pairs on `haplodna[…] != SEQMISSING`, so no pair
     passes and `ndd = 0`, hence `*nd = 0`.
  5. `cal_mpd1p` counts the locus only `if (nd > 0)` → **the locus vanishes from
     the denominator.**
  6. `cal_nss1p` instead counts it through `OK`, which `cal_nsspl` sets from
     `samplesize()` — and `samplesize()` reads `dataobs.ssize` plus the
     **observed** missing-haplotype list, never the simulated sequences. The
     locus **is** counted, contributing 0 sites.

  So `MPD` = (Σ π over polymorphic loci) / (**number of polymorphic loci**) while
  `NSS` = (Σ `S`) / (**all loci**). Measured in a debug run: **1168 of 18832
  simulated loci carry no mutation** (6.2%, and far more inside low-`θ`
  particles).

  **Why that produces exactly the observed signature.** As `θ → 0`,
  `MPD → θ/(a_79·θ) = 1/a_79 = 0.202`, a **constant** — which is why the gap
  looked like a fixed additive +0.06 independent of `θ`: `MPD` does not go to
  zero as it should, it plateaus. Confirmed: in the lowest `θ` quintile DIYABC's
  `MPD` mean is **0.1942** against the predicted 0.202. And correcting the
  denominator by `f = 1 − exp(−a_79·θ)` collapses `MPD` back onto the coalescent
  prediction (ratio to theory 2.434 → 1.010 / 1.057 / 1.032 / 0.996 across
  strata) and takes the **`MPD <= 0.75*NSS` violations from 4.6% to 0.1%**. That
  `f` over-corrects in the lowest quintile only (0.805), because the per-locus µ
  is Gamma-dispersed (CV 0.71) and a mixture of Poissons has more zeros than a
  single Poisson of the same mean — an imperfect estimator of `f`, not a flaw in
  the mechanism.

  **Scope: `MPD` and `VPD`, and nothing else.** Predicted from the code, then
  confirmed blind on the reftable. Ratio DIYABC/ours by `θ` quintile:

  | `θ` | `MPD` | `VPD` | `NSS` | `PSS` | `NHA` | `MNS` | `MP2` | `MPB` | `HST` |
  |---|---|---|---|---|---|---|---|---|---|
  | ≤ 0.20 | **1.951** | **1.907** | 0.973 | 0.980 | 0.992 | 0.952 | 0.951 | 0.955 | 0.989 |
  | 0.20–0.52 | **1.206** | **1.177** | 0.979 | 0.997 | 0.992 | 0.982 | 0.967 | 0.965 | 0.999 |
  | 1.0–1.75 | 1.013 | 1.016 | 0.989 | 1.004 | 0.997 | 0.976 | 0.978 | 0.983 | 1.007 |
  | 1.75–3.88 | 0.995 | 0.974 | 0.986 | 1.001 | 0.992 | 0.985 | 0.981 | 0.987 | 1.020 |

  Only the two whose denominator reads the simulated sequences move;
  the seven whose `OK` comes from `samplesize()` are flat at ~1 with no `θ`
  dependence. `cal_vpd1p`'s gate is **`nd > 1`**, not `nd > 0` — one notch
  stricter, so it also drops loci with a single usable pair.
  `cal_dta1p` is **not** affected: its early return on `dnavar < 1` leaves
  `OKK = true`, so monomorphic loci are counted with D = 0, as ours are.

  **Why every earlier verification missed it.** Transcribing
  `cal_nsspl`/`cal_mpdpl` and running them on the **observed** G2 sequences
  reproduces `statobsRF.txt` to 6 decimals. That check passed because
  `statobsRF` is computed with `dnatrue = true`, where a monomorphic locus still
  gets its full 100 bases and is never "missing". **A verification on observed
  data does not validate the simulated path** — the two run different branches
  of `init_dnaseq` and of `cal_numvar`.

  **The frequency spectrum of DIYABC's sequences is correct** — measured directly
  on the dumped simulated sequences of locus 10, 200 (particle, sample) blocks:
  23.3% singletons against 24.1% expected, mean minor-allele count 7.77 against
  6.40, and Kingman-consistent in **every** stratum of `S` including `S = 1`. So
  the earlier inference of "a spectrum shifted to k ≈ 13–27 instead of 5.5" was
  **wrong**: that was the denominator artefact read through π, not a genealogy
  difference. Don't look for one.

  **`DTA` — a retraction and a small open item.** `DTA` was briefly written up as
  "a second defect of inverse sign" on the strength of a DIYABC/ours ratio of
  0.12–0.67. **That was the near-zero-mean trap.** The absolute means are 0.0015
  to 0.0509 against a standard deviation of 0.07–0.14, so the mean sits 20–50x
  below the spread and the ratio is noise. `cal_dta1p` has the correct
  denominator. The KS flag on `DTA_2_2` that remained (4 cells of 4, +0.011 to
  +0.031) was **not** a distributional difference: it was correlated noise from
  additive per-locus seeds, resolved 05/10 — see "`DTA_2_2` residual (05/10)" under
  Closed investigations. The `S = 1` guard candidate was never needed (`_tajima_d_per_locus`
  already returns 0 when the denominator is not positive).

  **Reproduced (02/10)**, as the fidelity rule demands — it needed a denominator,
  not a simulator. `compute_MPD` and `compute_VPD` now `continue` past any locus
  with `ts.num_sites == 0`, so such a locus leaves their denominator while
  `compute_NSS` (which divides by `len(tree_sequences)`) still counts it. The
  guard is a two-line block at the top of each loop, **preceded by a comment that
  states the defect and the exact lines to delete to restore the correct
  behaviour** — the user asked for the previous (correct) logic to stay recoverable.
  `num_sites == 0` is the right test: it is equivalent to DIYABC's `nmutot == 0`,
  since the `MatrixMutationModel` we use has a zero diagonal, so every mutation
  is real and several at one site still count as one site.
  Effect on the golden tests: `test_mean_pairwise_differences_per_group` and
  `test_variance_pairwise_differences_per_group` changed on **G2 only**, by
  exactly **5/4 = 1.25** on both populations — that fixture has 5 G2 loci of which
  1 carries no mutation — while G3 (5 of 5 polymorphic) stayed bit-identical. That
  asymmetry is the unit-level proof that the guard keys on the absence of
  mutation and nothing else (**superseded 05/10**: that fixture changed with the
  seeds and no longer has a monomorphic locus; the unit-level proof is now
  `test_mean_pairwise_differences_excludes_loci_without_mutation`). The new golden values carry a comment saying not to
  restore the old ones without removing the guard.
  **Validated against the real DIYABC on 05/10** (3 models × 2 replays, after the
  seed fix, every statistic without a KS flag, as reported). The expectations
  below were the ones set at the time; the ratio per `θ` quintile and the
  `MPD <= 0.75*NSS` violation rate were not re-measured separately. Original text:
  rerun the two replays of `toy_example2_ms_dna_50loci_K2P` (restart the Jupyter
  kernel first) and expect:
  `MPD`/`VPD` ratio DIYABC/ours **flat at ~1 across `θ` quintiles** instead of
  1.95 → 0.995; `MPD_2_1`, `MPD_2_2`, `VPD_*` **dropping out of the KS flags**;
  our own `MPD <= 0.75*NSS` violations **rising to ~4.6%** to match DIYABC's (it is
  a property of the bug, so a faithful port must now reproduce it too);
  `DTA_2_2` **remaining flagged** (it did not: see above).
  If the ratio does not flatten, the guard is not equivalent to DIYABC's and the
  per-sample vs per-locus distinction in `cal_mpdpl` is the first thing to recheck.
  **The other DNA datasets change too** (`toy_example2_ms_dna_{JK,K2P,TN}` at 5
  loci/group, the mixed path): anything with monomorphic loci. Their stored
  comparisons predate this and are no longer comparable to a fresh run.

  **How the sequences were obtained**, because this cost a whole session: the
  dump lives in `cal_nss1p` under `debuglevel == 9` **and** `kloc == 10` (the
  first G2 `<A>` `[S]` locus). The invocation is
  `general -p ./ -r 2 -R ALL -a 9`. **`-R ALL` is the load-bearing part**:
  without it `randomforest` stays false, `readHeader` goes into
  `readHeaderEntete` instead of `readHeaderAllStat` and **hangs at 100% CPU**,
  which is what made three earlier attempts (`-r 1 -g 1`, `-r 5 -g 5`, `-r 1000`)
  look like a slow simulation when they had never left the header. `-g` defaults
  to 100, so `-r 2` still simulates 100 particles — 19 s, 173k log lines. Use
  `nice -n 19` and one process at a time.

  Falsified along the way, so don't retry them: the discrete coalescent (see the
  `evalcriterium` entry under Closed investigations); `put_mutations` is
  textbook — per-branch `poisson(length * mus_rate * dnalength)`, so the
  allocation of mutations across branches is unbiased; `coal_pop`'s continuous
  branch is correct, including the subtlety that the new parent node carries the
  target `pop` before both `draw_node` calls yet can never be drawn as its own
  child because it holds the highest index; `cal_numvar` clears `nuvar` between
  loci; `liberednavar` runs after the whole stat loop; `sample_site_rates` is a
  faithful port, and `nsv = 90` means the **first 10** sites are the invariant
  ones with 90 mutable; `dnalength` is 100 both declared and real in the `.mss`;
  the tier-2 per-locus µ draw is active and unbiased on both sides; and "`MPD`
  normalised by polymorphic loci" was the right shape but was dismissed in error
  because `PSS` exceeds 1 — `PSS` is the count of **private** segregating sites,
  not a proportion.

  **A second latent bug found in passing, not triggered here.** In `mute`'s
  non-ACGT error branch, `n` is reused as the print loop counter and ends at
  `dna.length()`, and the `exit(1)` below it is commented out — so `switch
  (dna[n])` reads `dna[size()]` (`'\0'`), no case matches, `dnb` keeps the
  original base, and `dna[n] = dnb` writes out of range. Effect if it fires: the
  mutation is silently lost plus an out-of-range write. Same family as the
  documented `mutsit`/`sitefix` bug — **a loop counter reused in place of the
  intended variable**. Ruled out empirically: zero "probleme" lines in
  `diyabc_run.log` (432 lines, the full 1000-particle run). Also latent:
  `init_dnaseq`'s inverse-CDF loop `while (s < ra) { k++; s += mutsit[k]; }` is
  unbounded, so a normalized `mutsit` summing to 1−ε lets `k` run past
  `dnalength`; probability ~1e-16 per mutation, no guard.

  **Practical rule.** Beyond `θ ≈ 1` the ratio is 1.013 then 0.995, so the defect
  bites only on low-diversity loci. On such a dataset **DIYABC's `MPD` and `VPD`
  are not a trustworthy reference**, and a disagreement there is not evidence
  against our port.

- **The DNA global shift is the loci-count residual, and it shrinks as
  documented.** Same dataset family, sign-test median aggregated per run (the
  only independent unit): at 5 loci per group, `K2P` gave +1.542 % and +0.071 %
  (scatter 1.47 points); at 50, +0.613 % and +0.309 % (scatter 0.30 points). A 5x
  drop in scatter and a halved magnitude for 10x the loci — the signature of
  "Systematic negative bias" under Closed investigations, positive-signed here.
  Twelve measurements at 5 loci ranged −0.76 % to +4.69 %, which is why two
  successive verdicts were written and retracted on that dataset before the loci
  count was raised. **Do not read a global bias off a 5-loci-per-group dataset.**
  The `JK` and `TN` branches of `build_transition_matrix_dna` are validated by the
  per-column agreement at 5 loci (14 flags / 504 comparisons = 2.8 % against
  25.2 expected), which never depended on the global shift.
- **`_genotypes_by_sample_and_individuals` not audited** for the `np.int64` leak
  that caused the `FST` bug. Harmless today (`LIK`'s formula does no boolean
  addition), but the next consumer of those tuples inherits the trap.
- **`group_prior_column_names` hardcodes the `µ` prefix** (U+00B5) while DIYABC
  also writes plain ASCII `mu` — `toy_example1_ms`'s trailer and reftable both use
  `mumic_1`. **The 30/09 verdict "inert, not latent" is RETIRED — since 01/10 the
  spelling is active in the text path.** It rested on "the names it returns are
  never looked up by name anywhere", and that stopped being true when
  `write_reftable_txt` began emitting these columns: its header line carries our
  `µ` spelling, and the comparison against a real reftable is **by column name**.
  Measured on `toy_example1_ms`: 4 of the 6 columns match by name, the two `µ`
  ones do not (`mumic_1`/`mumic_2` on DIYABC's side).
  **Normalize in the comparator, not in `bridge/`** — `_normalize_mu` already
  exists for it. The decisive reason: the `.bin`, which is the project's
  structural deliverable and what `abcranger` reads, contains **no column names at
  all**, so the spelling cannot affect it. Making `_group_prior_columns` follow the
  header's own spelling would be needed only to make our text file diffable as is.
  Still true, and still the reason not to "fix" the spelling hoping to fix a bug:
  `parse_real_reftable_params_with_group_priors` reads these columns **by
  position** (`tokens[1 + len(priors_param_names) + i]`) and uses
  `group_priors_names` only as the returned dict's keys. So the spelling cannot
  cause a wrong *read*.
  **The real exposure the positional read creates is different**: a wrong *count*
  from `group_prior_column_names` shifts every value silently, labelling them with
  keys that do not belong to them — which is exactly the 30/09 corruption
  (`MODEL JK` with a `K2P` trailer: 7 names expected, 5 columns present). Nothing
  downstream checks that the column read is the one intended. That is what
  `check_header_trailer_line` and `check_real_reftable_matches_header` now guard,
  upstream. Don't "fix" the `µ` spelling expecting to fix a bug; do reach for the
  guards if a group-prior value ever looks shifted.
- **`parse_group_priors`'s `else` branch** assumes any unrecognized line is a
  `MODEL` line, with no positive validation — a malformed line is silently
  mis-parsed rather than raising.

## `scripts/` — ad hoc investigation scripts

The `.py` files under `scripts/` are ad hoc scratch/investigation scripts, not
part of the `bridge/` package or its test suite — treat them as disposable
when reasoning about the architecture, but **don't delete them without asking
the user**, since they double as informal experiment logs.
