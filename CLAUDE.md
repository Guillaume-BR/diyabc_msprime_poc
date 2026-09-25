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
| MicroSat | complete, 11 stats + `AML` | `toy_example2_ms_dna`, `toy_example1_ms_modified` |
| DNA sequences | complete, 13 stats | `toy_example2_ms_dna` (`K2P` only) |
| DNA sequences `<X>`/`<Y>` | implemented | synthetic fixture only — the real binary SIGSEGVs on this case |
| Serial/temporal sampling (MicroSat) | complete | `toy_example1_ms` (1 population, 4 sampling times) |
| Serial/temporal sampling (SNP IndSeq) | complete | `human_seriel` (1 population, 4 sampling times, 130 stats, 0/130 significant) |
| Serial/temporal sampling (SNP PoolSeq) | complete | `toy_example4_seriel` (100 loci, `<MRC=5>`, 10/133 and 7/133 over two replays vs ~6.7 expected) |

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
validated on the KS alone carry such a bias — see "Systematic negative bias"
under Open work.

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
  This has bitten the project three times.
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
- **Freshness check after editing a header**: `stat -c '%y %n'` on
  `headerRF.txt` vs `first_records_of_the_reference_table_0.txt` — the real
  reftable must always be *newer* than the header it was generated from.
- **`reference/` is read-only ground truth.** Never modify it.
- **Run scripts as modules from the repo root**: `python3 -m scripts.run_test`.
  `python3 scripts/run_test.py` fails with `ModuleNotFoundError: No module
  named 'bridge'` — Python only puts the script's own directory on `sys.path`.
- **A test that writes a modified `header.txt` to `tmp_path` must build a
  fresh `*ReplayContext`** from that text. Since the ReplayContext refactor,
  functions never re-read the header from disk, so "write a new file, then
  call the function on that directory" silently uses the original.

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
  were small relative to `N`. Whether this explains any observed discrepancy is
  NOT established — the one experiment run on it was confounded (see
  `notes/exploration.md`, 25/09).
  Don't reopen on the same rumor; if a precise source turns up, check it
  against those line numbers first.
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

- **`build_sample_sets_from_scenario(scenario, values, counts_by_samples)`**
  (`ancestry_simulation.py`) — one `msprime.SampleSet(n, population, time)` per
  `sample` event, in order. `counts_by_samples`'s **keys are ignored**: only the
  order of its values matters, the k-th count going to the k-th event. Indexing
  it by `f"pop{event.pop}"` returns the same count for every serial sample —
  that bug was written, caught by a deliberately unequal-count test, and is now
  locked by it. A sample time may be a parameter name, not a literal
  (`particuleC.cpp:599-605`), so this must be called **after** the prior draw;
  `build_demography` already evaluates those expressions (and discards them),
  which is why no constant-prior failure mode is introduced here. Guards on
  `len(sample_events) != len(counts)` — DIYABC does not, and a silent
  misalignment is the worst possible outcome.
- **`compute_sample_layout(ts, counts_by_samples)`** — the sample-aware twin of
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
the same shape with two specifics: `context.count_samples` gives the **haploid**
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
`counts_by_samples` (an ordered `{name: individual count}` dict), not a
precomputed layout: the layout depends on **ploidy**, which varies per heritage
type (`<A>` 2, `<H>`/`<M>` 1), so a single layout built upstream would be wrong
for every type but one. Each MAF loop derives its own from the counts, where
the `TreeSequence` is born and the ploidy is known. Both MAF loops have a
`maf == 0.0` fast path that must forward the counts too — `<MAF=hudson>`
datasets take *only* that path.

**Never derive `counts_by_samples` from `samples`.** They are different objects:
`samples` goes to msprime (`dict[str, int]` *or* `list[SampleSet]`),
`counts_by_samples` feeds `compute_sample_layout`. A `SampleSet` list carries no
names, and an `<X>` locus describes **one** sample with **two** `SampleSet`s.

**Deliberately not done, and it is silent**: `<X>`/`<Y>` loci override the
sample sets inside the per-locus loop with the sex-stratified builders, which
are **not** serial-aware — on both the MicroSat/DNA and SNP sides. On a serial
dataset carrying such loci the serial `SampleSet`s would be ignored for them.
No dataset triggers it today, and no guard raises. Related hazard, measured:
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

- **Serial sampling is not wired for `<X>`/`<Y>`, and nothing guards it.** The
  per-locus dispatch overrides the serial `SampleSet`s silently, on both the
  SNP and MicroSat/DNA sides. See "Serial/temporal sampling" under Domain
  knowledge.
- **Systematic ~1% negative bias on several datasets, project-wide.** Visible
  only through a **sign test on `rdiff_mean`**, never through the per-column
  KS. Measured medians: `toy_example5` -2.65% (70 loci, validated since July),
  `toy_example4_seriel` -1.06% / -0.79% over two replays (100 loci),
  `toy_example1_ms` -0.87%, `toy_example3` -0.74%; and neutral on
  `toy_example4` +0.22%, `toy_example4_MRC1` +0.00%, `human_seriel` +0.06%
  (5000 loci). Reproducible per column (`r = 0.86` between replays — but note
  both replays share the same real reftable, so a high `r` is expected as soon
  as ANY real per-column difference exists; it proves the difference is not
  sampling noise, not that it is large). **Not specific to serial sampling**:
  `human_seriel` is serial and clean. Leading hypothesis, matching the
  residual-bias verdict closed 17/07: the bias shrinks with loci count
  (`toy_example3` -0.74% → -0.24% at 500 loci), but `toy_example5` barely moves
  (-2.65% → -2.33% at 350), so it is not settled. Next experiment: rerun
  `toy_example4_seriel` at 500-1000 loci, everything else unchanged.
- **Rename `pop` → `samp` for everything that comes from the `.snp`/`.mss`**,
  and **only then**. Two distinct families of `f"pop{...}"` coexist, and only
  one is misnamed:
  - *real msprime populations* — `demography_builder.py` (7 sites) and
    `build_sample_sets_from_scenario`'s `population=f"pop{event.pop}"`. These
    are correct and must keep the name.
  - *samples wearing a population name* — `observed_count_population`,
    `individual_sexes_from_locus_genotype`, `build_samples_argument`. Their keys
    come from file-block order, i.e. the **sample** index.

  **Partly unblocked since 25/09.** The pipeline now builds `SampleSet`s on both
  the SNP and MicroSat/DNA paths, so the msprime-facing use of
  `build_samples_argument` survives only as the `sample_sets is None` fallback
  for direct callers (tests, `scripts/`). Renaming its keys still breaks those,
  so the fallback has to go — or be translated — first. Conversely the
  `counts_by_samples` / layout / `population_names` side is now free of msprime:
  those strings never reach `sim_ancestry`, so they can be renamed safely. That
  is also the bulk of the work. Scope, for planning: 329 `population_names`,
  86 `_by_population`, 40 `population_layout`, 226 literal `"popN"` in tests —
  one commit, never a half-rename leaving two vocabularies side by side.
- **`DTA_2_2` mean shift** on the 50+50-loci DNA dataset: a small but real
  mean difference (real ≈0.018, sim ≈0.043) on a **G2** column, with a std
  ratio near 1 — so it does not fit the (now resolved) G3 variance story. Never
  investigated.
- **No `JK`- or `TN`-model reference dataset** exists, so those two branches of
  `build_transition_matrix` are covered by hand-computed synthetic values only
  (`toy_example2_ms_dna` exercises `K2P` only).
- **`_genotypes_by_pop_and_individuals` not audited** for the `np.int64` leak
  that caused the `FST` bug. Harmless today (`LIK`'s formula does no boolean
  addition), but the next consumer of those tuples inherits the trap.
- **`parse_group_priors`'s `else` branch** assumes any unrecognized line is a
  `MODEL` line, with no positive validation — a malformed line is silently
  mis-parsed rather than raising.
- **`rewrite_loci_count`** still handles the condensed single-type format only.

## `scripts/` — ad hoc investigation scripts

The `.py` files under `scripts/` are ad hoc scratch/investigation scripts, not
part of the `bridge/` package or its test suite — treat them as disposable
when reasoning about the architecture, but **don't delete them without asking
the user**, since they double as informal experiment logs.
