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
| DNA sequences | complete, 13 stats; `JK`/`K2P`/`TN` all covered | `toy_example2_ms_dna_{JK,K2P,TN}` (5 loci/group, 2 runs each), `toy_example2_ms_dna_50loci_K2P` (50/group, 2 runs) — a **G2 residual is open**, see Open work |
| DNA sequences `<X>`/`<Y>` | implemented | synthetic fixture only — the real binary SIGSEGVs on this case |
| Serial/temporal sampling (MicroSat) | complete | `toy_example1_ms` (1 population, 4 sampling times) |
| Serial/temporal sampling (SNP IndSeq) | complete | `human_seriel` (1 population, 4 sampling times, 130 stats, 0/130 significant) |
| Serial/temporal sampling (SNP PoolSeq) | complete | `toy_example4_seriel` (100 loci, `<MRC=5>`, 4/133 vs ~6.7 expected) |

Validation means a *paired* comparison: the real DIYABC priors are replayed
particle-by-particle through our pipeline (`replay_reftable_simulation*`,
`scripts/replay_diyabc_priors*.py`) and the two reftables compared
column-by-column with a two-sample Kolmogorov-Smirnov test. That test is
**conservative** on a paired replay: both sides share the same prior draws, so
their statistics are correlated and the KS understates the gap. Judge a
residual flag by whether it **persists across two independent replays**, not by
its p-value alone — a real effect keeps its order of magnitude, noise moves to
other columns (see `notes/exploration.md`, 24/09, for a worked example).

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
   `base_frequency_by_locus`, `observed_count_population`,
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
    `compute_summary_statistics` (+ `_dna` / `_microsat` / `_from_values`
    variants) is the high-level entry point.

12. **`reftable_loop.py`** — runs `nrec` particles in parallel
    (`ProcessPoolExecutor`), writes the binary `reftable.bin` and the
    human-readable `write_reftable_txt`. Seeds are `particle_index + 1`
    (msprime rejects `seed=0`). Multi-scenario is supported: each particle
    draws its own scenario from the weighted list
    (`parameter_sampling.draw_scenario`, matching `ParticleC::drawscenario`),
    and `write_reftable_bin` writes a variable-length record per row (only
    that row's own scenario's `nparam` columns, no NA-padding), matching
    `reftable.cpp`.

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
`(context, ...)` — e.g. `build_matrix_per_locus`,
`build_matrix_microsat_per_locus`, `microsat_mutation_simulation_per_locus`,
`dna_mutation_simulation_per_locus`. Not all moved:
`build_group_local_param_per_locus`, `build_microsat_local_param_per_locus`
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
  (`int()` in `_length_by_pop_and_individuals` *and* in
  `_compute_ni_nA_AA_for_one_population`); all 11 MicroSat stats now match
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
  **Established 30/09 — it does explain a real discrepancy, and this is now the
  one known systematic difference between the two simulators.** Not replicating
  that switch shows up as an `N`-dependent divergence of genealogy *shape*, which
  is visible on any statistic sensitive to the frequency spectrum. Measured on
  `toy_example2_ms_dna_50loci_K2P` G2 (`<A>` sequences, 80 sampled gene copies,
  `N1 ~ UN[10,10000]`, all five populations at `N1`), stratifying the paired
  replay by the drawn `N1` — `MPD` (π) relative error, 2 independent DIYABC runs
  × 2 scenarios:

  | `N1` | n | run 0 sc1 | run 0 sc2 | run 1 sc1 | run 1 sc2 |
  |---|---|---|---|---|---|
  | <160 | 6–9 | −71.8% | −66.8% | −61.9% | −71.1% |
  | 160–500 | 16–29 | −52.8% | −53.4% | −54.4% | −49.5% |
  | 500–1500 | 42–52 | −28.0% | −28.6% | −26.5% | −21.3% |
  | 1500–4000 | 109–140 | −9.9% | −9.0% | −12.1% | −8.7% |
  | >4000 | 291–317 | −3.2% | −2.0% | −1.9% | −2.4% |

  Monotone, and **confirmed on all three substitution models** — the cause being
  coalescent-side, `JK` and `TN` must show the same gradient, and they do.
  `toy_example2_ms_dna_50loci_{K2P,TN,JK}`, 2 DIYABC runs × 2 scenarios each,
  12 cells, `MPD_2_1` relative error:

  | `N1` | `K2P` | `TN` | `JK` |
  |---|---|---|---|
  | <160 | −71.8 / −66.8 / −61.9 / −71.1 | −82.1 / −78.6 / −63.2 / −60.3 | −56.8 / −84.0 / −33.9 / −54.4 |
  | 160–500 | −52.8 / −53.4 / −54.4 / −49.5 | −59.6 / −34.2 / −44.8 / −45.6 | −51.1 / −42.2 / −52.8 / −55.4 |
  | 500–1500 | −28.0 / −28.6 / −26.5 / −21.3 | −27.9 / −23.7 / −29.6 / −16.8 | −30.6 / −27.1 / −20.7 / −22.3 |
  | 1500–4000 | −9.9 / −9.0 / −12.1 / −8.7 | −10.7 / −10.0 / −10.5 / −13.2 | −2.5 / −9.6 / −13.1 / −11.5 |
  | >4000 | −3.2 / −2.0 / −1.9 / −2.4 | −3.7 / +0.1 / −2.5 / −3.7 | −1.0 / −1.6 / −2.8 / −3.0 |

  Overall `MPD` amplitude overlaps completely across models (`JK` −1.9 to −5.4%,
  `K2P` −4.1 to −5.1%, `TN` −2.1 to −6.6%), so there is **no model component on
  top of the coalescent effect**. The `<160` stratum scatters because it holds
  5–17 particles. Carry the divergence-of-sign argument through `MPD`, negative in
  all 12 cells, rather than through `MNS`/`NSS`, whose excess is
  scenario-dependent.
  **First attempts on `TN` and `JK` had to be discarded**: both reftables came
  from headers whose trailer line had been copied from the `K2P` variant, and a
  `TN` analysis run on one of them was retracted — see the trailer-line rule under
  Hard rules, now guarded by `check_header_trailer_line` and
  `check_real_reftable_matches_header`. **The decisive
  quantity is π/`S`**, which depends on branch-length *shape* alone, not on the
  mutation rate or the loci count: ours is flat at ~0.24 across every stratum
  (the Kingman value), DIYABC's climbs 0.258 → 0.274 → 0.351 → 0.521 → 1.344 as
  `N1` falls. So it is DIYABC's genealogy that changes with `N`, not ours —
  multiple lineages merging within one generation produce polytomies and shift
  length from terminal to internal branches, raising π while `S` barely moves.
  That is why our `S` runs +1 to +3% while our π runs −4 to −5%: **more sites at
  rarer frequencies.** (The `N1 < 160` stratum reports π/`S` > 1, which is
  structurally impossible for a per-site process — worth its own look, but it
  rests on 6–9 particles and is not needed for the conclusion.)
  **Not replicated — now an open chantier**, see "Discrete coalescent" under Open
  work for the exact criterion (three branches, not two) and the msprime
  constraints. Until then, judge any `N`-dependent residual against this gradient
  before looking elsewhere, and prefer datasets whose priors keep
  `nLineages / N` small. Likely relevant to the open `ML3p_2.3.4` residual on
  `toy_example4_seriel`, which has 800 gene copies with `Npresent <= 1000` —
  i.e. discrete mode for *every* particle.
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
  `_length_by_population` gives every population the same key list (taken from
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
  `compute_summary_statistics_from_values` built `counts_by_sample` **before**
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
  `build_microsat_transition_matrix(...).alleles` first.
- **Statistics.** Two distinct data representations are required, not one:
  `_length_by_pop_and_individuals` (for `FST`) returns a 2-tuple per
  individual, duplicating a haploid's single allele into `(taille, taille)` —
  safe because `cal_Fst2p` itself treats a haploid copy as a doubled
  diploid-like pair, making the ANOVA ploidy-invariant. `_genotypes_by_pop_and_
  individuals` (for `LIK`) keeps **true** ploidy, because `cal_lik2p` applies a
  genuinely different formula per ploidy. Every other stat uses the flat
  `_length_by_population` `(taille, compte)` table.
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
- **DIYABC's Jukes-Cantor token is `JK`**, not `JC`/`JC69`. Got this wrong
  twice.
- `build_transition_matrix`'s row-sum normalization needs
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
not an idealized reimplementation. Both are judged unintentional C scoping
slips, kept because DIYABC's own output depends on them.

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

- **Write the mutation parameters (`nparamut`) into our reftables.** DIYABC's
  reftable carries `µmic_1`, `pmic_1`, `snimic_1`, `µseq_2`, `k1seq_2`… ; ours
  carries none of them. Measured on `toy_example1_ms`: our replay has 156 columns
  against DIYABC's 162, the six missing ones being exactly
  `mumic_1 pmic_1 snimic_1 mumic_2 pmic_2 snimic_2`, nothing extra on our side.
  `write_reftable_bin`'s own docstring already flags it ("Ne gère PAS les
  paramètres de mutation (absents de human) -- à ajouter … si un dataset avec
  microsatellites/séquences est traité plus tard"); that "later" arrived when
  MicroSat and DNA were validated, and the note stayed.

  **Two distinct costs.** *(a) Structural*: the project's goal is a `reftable.bin`
  equivalent to DIYABC's, and `abcranger` reads it expecting
  `nparam[scenario] = nparamhist + nparamut` floats per record while ours writes
  `nparamhist` only — a reader trusting the header would shift every row. It has
  never bitten because validation goes through the **text** comparison of a
  replay, never through `abcranger` on our output. *(b) A missing check*: in a
  replay the historical priors are verifiable (KS p = 1 exactly, which is what
  proves the pairing is real), but the mutation parameters are **injected without
  ever being checked**. If `parse_real_reftable_params_with_group_priors` read a
  shifted column — the positional-read risk noted under the `µ` prefix entry — the
  simulation would use wrong µ/k1/k2 and nothing would show it; the prior check
  would still pass. Emitting these columns turns that risk from silent into
  visible.

  **Three layers, in order.** `ParticleResult` gains a
  `group_parameter_values: dict[str, float]` field — today the values are drawn
  inside the worker, used for the mutation model and discarded, so the information
  never reaches the writers. Then the six `_run_single_particle*` carry them up,
  keeping the sibling symmetry (the `_from_values` variants receive DIYABC's real
  values, the drawing variants draw their own; both must store the same thing).
  Then `write_reftable_txt` and `write_reftable_bin` emit them **last, after the
  demographic parameters**, as that docstring already states from
  `readReftable.R`, in `group_prior_column_names` order — which is
  scenario-independent, unlike `nparamhist`.
  This touches validated originals, so the necessity is worth stating: it is real
  for `write_reftable_bin` (structural equivalence) and arguable for the text one.
- **Replicate DIYABC's discrete generation-by-generation coalescent.** This is
  the one known systematic difference between the two simulators, measured and
  reproducible — see the `evalcriterium` entry under Closed investigations for
  the `N1`-stratified gradient (π error −72% to −2%) and the π/`S` argument.

  **The criterion, verbatim** (`particuleC.cpp::ParticleC::evalcriterium`, 1250-1273),
  evaluated **per segment and per locus**, with `nLineages` the lineage count on
  entering that segment. It returns 1 when the *continuous* approximation is
  acceptable, 0 when DIYABC switches to discrete:

  ```cpp
  if (seqlist[iseq].t1 < 0) return 1;              // last, infinite segment: always continuous
  nGen = seqlist[iseq].t1 - seqlist[iseq].t0;
  ra   = (double)nLineages / (double)seqlist[iseq].N;
  if      (nGen <=  30) OK = ra < (0.0031*nGen*nGen - 0.053*nGen + 0.7197);
  else if (nGen <= 100) OK = ra < (0.033*nGen + 1.7);
  else                  OK = ra < 0.5;
  ```

  Note the `nGen <= 30` branch is **quadratic** and non-monotone (0.720 at
  nGen = 0, minimum ≈0.49 around nGen ≈ 8.5, 1.92 at nGen = 30), and that it is
  discontinuous with the next branch (2.72 at nGen = 31). Transcribe it as is;
  it is a heuristic, not a formula to tidy.

  **What msprime gives us** (verified on 1.4.2): `msprime.DiscreteTimeWrightFisher(duration=...)`
  exists and mixes with `StandardCoalescent` through a time-ordered `model=[...]`
  list, so a per-segment schedule is expressible. Two obstacles, both real:
  - **DTWF requires `ploidy = 2`** (`LibraryError: The DTWF model only supports
    ploidy = 2`). Every non-`<A>` heritage type runs at `ploidy=1` here, so
    `<H>`/`<M>`/`<Y>` and `<X>` males have no direct route. Whether a halved
    rescaled size at `ploidy=2` is equivalent needs proving, not assuming.
  - **The criterion is dynamic**: `nLineages` is not known before simulating, and
    it differs per locus. DIYABC evaluates it once on entering each segment with
    the live count. An msprime `model=[...]` schedule is fixed up front, so it can
    only approximate — e.g. from the expected lineage count at each segment, or by
    simulating in stages and re-deciding between them.

  **Success criterion, already in place**: rerun the `N1`-stratified replay on
  `toy_example2_ms_dna_50loci_K2P` and require the `MPD` gradient to flatten and
  DIYABC's π/`S` climb (0.258 → 1.344 as `N1` falls) to be matched instead of
  staying flat at the Kingman ~0.24. That measurement is the test; it exists and
  is reproducible over 2 runs × 2 scenarios.

- **Serial sampling is not wired for `<X>`/`<Y>` — now guarded, not fixed.**
  The per-locus dispatch still overrides the serial `SampleSet`s, on both the
  SNP and MicroSat/DNA sides. Since 28/09 the combination is refused up front by
  `reftable_loop.raise_if_serial_with_sex_linked_loci`, called from all seven
  run/replay entry points; the predicate is
  `scenario_parser.is_serial_scenario`. Implementing it for real was
  deliberately deferred: **no dataset would validate it**, and this project
  never declares a path correct without a paired replay against the real
  DIYABC. See "Serial/temporal sampling" under Domain knowledge.
- **`ML3p_2.3.4` on `toy_example4_seriel`: a real residual, not investigated.**
  Significant in **three** replays (p = 0.0001 / 0.0013 / 0.0002) with a stable
  −10 to −14% rdiff, and it survives a 10x increase in loci count that wipes out
  every other flagged column. The project's persistence criterion is met, so
  this is **not** a false positive. `ML2p_2.4` is intermediate: significant in
  both 100-loci replays, not at 1000, amplitude decaying — the loci-count
  residual, with a consistently negative sign worth noting. Left unexplored by
  the user's call (cost vs. stakes).
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
  simulators agree best**. On `toy_example2_ms_dna_50loci_K2P` the discrete
  signature runs the other way (π error −72% at `N1 < 160` down to −3% above
  4000). `FST2m_2.4` stays within ±4% in every stratum, so the pairing and the
  stratification are sound.
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
- **G2 frequency-spectrum residual — CAUSE FOUND, not fixed.** On
  `toy_example2_ms_dna_50loci_K2P`, two independent DIYABC runs × 2 scenarios ×
  2 samples, no exception: `MNS`/`NSS` **+1 to +3%** (8/8 positive) while `MPD`
  (π) is **−3.7 to −5.1%** (8/8 negative). Opposite directions rule out a rate
  mismatch — a rate error moves `S` and π together, as in the G3 control. KS
  flags 17 of 168 comparisons (10.1% against 5% expected), **all in G2**, all in
  `MPD`/`VPD`/`DTA`, with `MPD_2_2` and `DTA_2_2` in all four cells.
  **The cause is DIYABC's discrete generation-by-generation coalescent**, which
  this project does not replicate — see the `evalcriterium` entry under Closed
  investigations for the `N1`-stratified measurement and the π/`S` argument.
  Falsified along the way, so don't retry them: `sample_site_rates` is a faithful
  port (`p_fixe` as a percentage, `ggamma3(1.0, gams)` ≡
  `gammavariate(gams, 1/gams)` with mean 1 and variance `1/gams`, the `sitefix`
  bug, zeros applied before normalizing over *all* sites); and
  `_segregating_sites_mask` matches `cal_nsspl` exactly (≥2 distinct non-missing
  states within the sample), as `_pairwise_hamming_distances` matches
  `cal_mpdpl`.
  Fixing it means a coalescent-side rewrite, so the practical rule is the one in
  the closed entry: check the `N1` gradient first, and prefer priors that keep
  `nLineages / N` small.
  G2 is where it surfaces because its mutation rate is **10x lower** than G3's
  (`UN[1e-8,1e-6]` vs `UN[1e-7,1e-5]`): few mutations per locus, so the
  statistics are dominated by genealogy shape instead of saturating.
  **`DTA` carries no usable amplitude** — real 0.0200 vs simulated 0.0459 is an
  absolute shift of +0.026 on a mean-of-per-locus-D that nearly cancels to zero,
  hence `rdiff` of +130%. Judge this family on `MPD`. **Never read a `rdiff` on a
  statistic whose mean sits near zero.** This item absorbs the former `DTA_2_2`
  entry.

- **The DNA global shift is the loci-count residual, and it shrinks as
  documented.** Same dataset family, sign-test median aggregated per run (the
  only independent unit): at 5 loci per group, `K2P` gave +1.542 % and +0.071 %
  (scatter 1.47 points); at 50, +0.613 % and +0.309 % (scatter 0.30 points). A 5x
  drop in scatter and a halved magnitude for 10x the loci — the signature of
  "Systematic negative bias" under Closed investigations, positive-signed here.
  Twelve measurements at 5 loci ranged −0.76 % to +4.69 %, which is why two
  successive verdicts were written and retracted on that dataset before the loci
  count was raised. **Do not read a global bias off a 5-loci-per-group dataset.**
  The `JK` and `TN` branches of `build_transition_matrix` are validated by the
  per-column agreement at 5 loci (14 flags / 504 comparisons = 2.8 % against
  25.2 expected), which never depended on the global shift.
- **`_genotypes_by_pop_and_individuals` not audited** for the `np.int64` leak
  that caused the `FST` bug. Harmless today (`LIK`'s formula does no boolean
  addition), but the next consumer of those tuples inherits the trap.
- **`group_prior_column_names` hardcodes the `µ` prefix** (U+00B5) while DIYABC
  also writes plain ASCII `mu` — `toy_example1_ms`'s trailer and reftable both use
  `mumic_1`. **Verified 30/09 to be inert, not latent**: the names it returns are
  never looked up by name anywhere.
  `parse_real_reftable_params_with_group_priors` reads these columns **by
  position** (`tokens[1 + len(priors_param_names) + i]`) and uses
  `group_priors_names` only as the returned dict's keys; and
  `write_reftable_txt` does not emit them at all — on `toy_example1_ms` our
  replay has 156 columns against DIYABC's 162, the six missing ones being exactly
  `mumic_1 pmic_1 snimic_1 mumic_2 pmic_2 snimic_2`, with nothing extra on our
  side. So the spelling cannot cause a wrong read, and the comparison notebooks
  drop those columns on both sides.
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
