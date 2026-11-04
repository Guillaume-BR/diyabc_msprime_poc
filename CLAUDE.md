# CLAUDE.md

Guidance for Claude Code in this repository. **Only what an agent cannot guess from the
code goes here.** Anything else lives elsewhere:

- `notes/results.md` — validation status, closed investigations (verdicts), benchmarks, open work.
- `notes/exploration.md` — the narrative of every chantier (read it before touching parsing or statistics code).
- Docstrings of `bridge/` — each statistic names the exact `sumstat.cpp` function it transcribes.

Do not grow this file back into a journal: a measurement or a p-value belongs in `results.md`.

## What this is

A proof-of-concept replacing DIYABC's C++ coalescent simulator (`particuleC.cpp::dosimulpart`)
with `msprime`. `bridge/` turns a DIYABC `header.txt` + observed data file into a `reftable.bin`
entirely in Python (no call to the DIYABC binary) and must be structurally and statistically
equivalent to the real DIYABC output. Data families covered: SNP IndSeq / PoolSeq (heritage
`<A>/<H>/<X>/<Y>/<M>`), MicroSat, DNA sequences, mixed MicroSat + DNA, serial sampling.
Per-dataset state: `notes/results.md` §1.

`reference/` holds ground truth produced by the real DIYABC binary — **never modify it**
(and don't remind the user about modified files there).

## Environment and commands

Use the `diyabc_msprime` conda environment (Python 3.11, msprime 1.4.2). The system Python has
no msprime. The DIYABC C++ source is a sibling repo, `~/Documents/Github/diyabc/src-JMC-C++/`,
and is the ground truth for every formula (don't nest it here).

```bash
conda activate diyabc_msprime
pytest tests/ -v
ruff check . ; ruff format .
python3 tools/generate_report.py     # regenerates notes/tree.md, commits.md, api.md, report.md
python3 -m scripts.run_test          # scripts run AS MODULES from the repo root
```

## Code style

- Docstrings are **Google style** (`Args:`/`Returns:`/`Raises:`, one-line summary first); see `loci_parser.py`.
- **Population vs. sample** (DIYABC's own distinction). Test: *on a serial dataset with 1 population and 4
  samples, does this variable hold 1 value or 4?* Four → sample, one → population.
  - `pop` stays in `scenario_parser.py`, `demography_builder.py`, `event.pop`, the msprime API, and the string
    literals `"pop1".."popN"` (msprime population names, resolved by name: **never rename**).
  - `sample` for everything coming from the `.snp`/`.mss` file or indexed by a statistics column. Suffix
    `_per_sample`, **singular**, never `_per_samples`.
  - File-format names stay the format's (`.snp` column `POP`, `pop_column_index`).
  - `compute_population_layout` and `compute_sample_layout` both keep their name; don't merge them (their equality
    on non-serial data proves the substitution neutral). A local built from a tskit attribute
    (`ts.node(...).population`) is a real population.
- **Renaming prose is not renaming identifiers**: a blind French `population` → `échantillon` breaks agreement and
  merges two notions (a sample vs msprime sample *nodes*). Use "gene copies", "sample nodes", "sequences" for the latter.
- Reftable column names are a literal prefix + an integer, so renaming Python identifiers can never move an output
  column.

## Architecture — only the traps

Pipeline: `scenario_parser` → `prior_parser` / `stats_group_parser` / `loci_parser` → `parameter_sampling` →
`demography_builder` → `observed_data` → `ancestry_simulation` → `summary_statistics` → `pipeline` →
`reftable_loop`. One module per stage, no cross-cutting logic. Non-obvious points:

- **`scenario_parser`**: only `NotImplementedError` is swallowed (block skipped with a warning); anything else propagates.
- **`prior_parser.is_constant_prior`**: constant iff `(max-min)/max <= 1e-6`, **never** constant when `max == 0.0`
  (DIYABC excludes constants from `reftable.bin` columns). `_extract_historical_priors_section` falls back to an empty
  line as terminator when `DRAW UNTIL` is absent.
- **`stats_group_parser`**: the flattened statistics list is order-preserved and **deliberately not deduplicated**
  across groups.
- **`loci_parser`**: two formats (condensed `"70 <A> 10 <X> ... G1 from 1"` and one-locus-per-line). The trailing
  `dnalength` on `[S]` lines is informational only (`header.cpp:401-402`).
- **`parameter_sampling._draw_one_value`**: `UN/LU/NO/LN/GA` reproduce `PriorC::drawfromprior`, including truncation by
  rejection for `NO/LN/GA` and the `GA` reparameterization (`gammavariate(shape=sdshape, scale=mean/sdshape)`).
- **`demography_builder`**: populations are `"pop1".."popN"` by 1-indexed position in `header.txt`. The reftable
  parameter columns are the priors a *specific* scenario references (`get_parameter_names_used_by_scenario`), not the
  full prior list.
- **`observed_data`**: population *i* = the *i*-th population **by first appearance** in the file. Never replace the
  order-preserving dict/Counter with anything that could reorder keys. `.mss` (genepop, `POP` separators) and `.snp`
  (columns) never share a reader.
- **`ancestry_simulation`**:
  - SNP: one independent tree per locus, exactly one Hudson mutation on an edge chosen proportionally to branch
    length. MAF filtering (`with_maf_filter`) mirrors `ParticleC::mafreached`; `<MAF=hudson>` is the no-op fast path.
  - DNA/MicroSat genealogies come from one place, `_simulate_ancestries_for_loci(list_loci, family, ...)`, grouped by
    `(heritage, sequence length)` (length = `dnalength` for DNA, **1 for MicroSat**). `<M>`/`<Y>`: one shared tree per
    group; `<A>`: one `sim_ancestry(num_replicates=len(group))`; `<X>`/`<Y>`/`<H>` own samples, one call per locus.
    The mutation loops only read that dict and keep their own seeds (`i` = position in the **full** `list_loci`).
  - SNP IndSeq trees come from `simulate_light_replicates` (msprime private API: `ancestry._parse_sim_ancestry`,
    `_msprime.Simulator.run`), bit-identical to `sim_ancestry(num_replicates=n)`. Same for the light `<A>` mutation
    path (`_light_a_trees`, `_mutate_tree`, `LightTables`). Everything is behind the guard
    `msprime.__version__ == _VERSION_MS_PRIME_LIGHT_REPLICATES` (`configuration.py`); the public-API path is the fallback.
  - MAF loop: one unlimited tree stream per particle (`_MAF_BATCH_SEED_OFFSET`), mutation seed of attempt `k` =
    `_locus_seed(seed, _MAF_REJECTION_SEED_OFFSET, k)`.
  - Never route a shared-tree locus (`<M>`, `<Y>`) through the light `<A>` path without copying its tables per
    locus: lightweight tables are consumed by the mutation call.
- **`summary_statistics`**: the DNA dispatch (`_DNA_PER_SAMPLE_STATS`, `_DNA_PAIRWISE_STATS`) calls every statistic with
  the same keywords. **A new shared input is a new keyword on all 13 statistics** (the ones that don't read it accept
  and ignore it, and say so in the docstring); each consumer keeps an `if x is None:` fallback so it stays callable alone.
- **`pipeline`**: orchestrates, no logic of its own. `compute_summary_statistics_*` (+ `_dna`/`_microsat`/`_mixed`/
  `_from_values`) are the entry points.
- **`reftable_loop`**: `ProcessPoolExecutor`, seeds = `particle_index + 1` (msprime rejects `seed=0`); multi-scenario
  (`draw_scenario` mirrors `ParticleC::drawscenario`); `write_reftable_bin` writes a variable-length record per row (only
  that row's scenario's `nparam`, no NA padding). `nparam[i]` counts **both** `nparamhist + nparamut` and there is no
  separate `nparamut` header field (`reftable.cpp:64`): adding a header int shifts the whole file.
- **`snp_writer.py`** is only used by the deprecated subprocess path (shelling out to the DIYABC `general` binary, kept
  for the `DIYABC_GENERAL_PATH` skip-marked cross-validation tests). If you see `-g`, it is the *internal batch size*,
  not a loop count (`-g 50` silently discards 49/50 simulated particles).
- **ReplayContext refactor (16-18/09)**: `SnpReplayContext` / `MicrosatReplayContext` / `DnaReplayContext` are built once
  per run and hold what was re-read from disk. Most per-locus builders take `(context, ...)`; `build_local_param_*_per_locus`
  and `build_rate_map_per_locus` still take `(header_text, seed)`. A value belongs in the context iff computing it involves
  a disk read. **A signature quoted in an old note may be stale: check the real one.**

## Hard rules

- **`header.txt`/`headerRF.txt` trailer line.** The last line (`"scenario N1 N2 ... ML1p_1 ..."`) looks like documentation
  but is **re-read as input** (`HeaderC::readHeaderAllStat`): `nparamhist = tokens - 1 - nstat - nparamut`. A trailer
  copied from another scenario/model corrupts DIYABC's state (statistics 300 %–10000 % off). When editing a header,
  **regenerate the trailer from that scenario's own priors, never copy it**. It has bitten five times.
  `reftable_loop.check_header_trailer_line` guards it (seven entry points; no-op without a `group priors` section).
  DIYABC writes the `mus_rate` prefix as `µ` (U+00B5) on `te2`, ASCII `mu` on `toy_example1_ms`: `_detect_mu` reads the
  spelling from the trailer, `_normalize_mu` (in `configuration.py`) normalizes both.
- **`DRAW UNTIL` exists iff `nconditions > 0`** (`header.cpp:274`): with 0 order constraints that line must be absent, or
  every following `getline` shifts by one. Same family as the trailer: a decorative-looking line that is positional.
- **Check a real reftable against its header by CONTENT, never by timestamp** (`check_real_reftable_matches_header`
  compares the reftable's first line with the trailer). DIYABC rewrites `headerRF.txt` *after* the reftable, so timestamps
  prove nothing.
- **Never read a raw DIYABC reftable with a generic whitespace parser.** `first_records_of_the_reference_table_0.txt` is
  **ragged** (rows omit the parameters their scenario doesn't use), so `pandas.read_csv(sep=r"\s+")` shifts every later
  column. Read the rectangular rewrite from `rewrite_real_reftable_txt` (`first_records_clean.txt`, unused = `nan`), and
  compare **all** rows, not a hardcoded `scenario = 1`.
- **Never derive a seed by addition.** `seed + OFFSET + i` makes neighbouring particles share random streams (consecutive
  particle seeds) and produces fake "persistent" flags. **Every derived seed goes through
  `configuration._locus_seed(seed, offset, index)`**; a new offset is a new constant in `configuration.py`, never a sum.
  What stays raw: the particle seed itself. A batch seed must be unique per (family, batch) (`_FAMILY_SEED_OFFSET`,
  `batch_A_counter`): a KS cannot see two consumers sharing trees, so any new derived seed needs a test that two
  independent consumers get different trees (`test_independent_family_trees`).
- **Replay (`seed_offset`)**: particle seeds are `index + 1` in every replay, so a second replay re-uses the same streams.
  The four replay chains take `seed_offset`, and each `scripts/replay_diyabc_priors_*.py` exposes it as `SEED_OFFSET`
  (**script default 200 000, not 0**). Re-run a suspicious flag with an offset before investigating it: a flag that moves is noise.
- **MAF datasets: a leftover `.snp.bin` / `.snpbin.txt` makes DIYABC skip the conversion, so no `maf.txt`, and it runs
  WITHOUT the MAF filter, silently** (`data.cpp:1083, 1137`). It also **forces the sex ratio to 0.5** (`data.cpp:1649`).
  Before a DIYABC run on a dataset with `maf > 0` or `NM != 1`, and after any edit of the `.snp` first line: delete
  `*.snp.bin` and `*.snpbin.txt`, then check `maf.txt` exists, the run prints `THIS->MAF=<value>` and `reftable.log` ends in `END|n|`.
- **`-R ALL` ignores the `group summary statistics` block and rewrites `headerRF.txt` at start-up.** It regenerates all
  statistics of the data type (2 samples: 58 columns; 4 samples: 224, incl. `AML`/`SML`), so editing statistics out of a
  header does nothing. It reads the trailer as input (`header.cpp:797`): a trailer with extra words makes it take statistic
  names for parameters and duplicate them. **Check `wc -w` of the trailer (1 + hist + mut + nstat) before AND after a run.**
  It also removes unused historical parameters from the priors section but leaves them in the trailer: don't declare unused
  priors. `-R ALL` without `-r` gives an empty reftable; `-r N` without `-g N` simulates 100; `-R` without `ALL` hangs.
  The replay scripts use `STATS_FILTER = "HEADER"` so a missing statistic raises; a replay not going through them must
  pass `"HEADER"` itself (`"ALL"` means "all *implemented*", no error).
- **An independent replay needs a NEW DIYABC reftable, initialised in its own step**: `./diyabc -p ./ -n "t:16;s:<seed>;f:"`
  only writes `RNG_state_0000.bin` and exits; `-R ALL ...` is a second command. Three msprime seed sets against ONE reftable
  are one replica; **the two scenarios of one run are not two replicates either**.
- **Reading a validation (the method)**:
  - Paired replay: real DIYABC priors replayed particle by particle (`replay_reftable_simulation_*`), reftables compared
    column by column. Run **both** a KS per column and a **sign test on `rdiff_mean`** (a uniform −1 % is invisible to the
    KS). Report `k/n`, p **and** median; trust none on one replay (p swings by 3 orders of magnitude at equal amplitude and is
    anti-conservative as columns are correlated; the KS is conservative on a paired replay).
  - Aggregate per independent run (new DIYABC reftable) before concluding.
  - Sanity tests, excluded from `n`: historical-prior columns and mutation-parameter columns (`µmic_N`, `pmic_N`, `snimic_N`,
    `µseq_N`, `k1seq_N`, `k2seq_N`) must come out at `rdiff = 0`, `ks = 0`/p = 1 exactly; otherwise the replay is not paired.
  - Round our statistic columns to **6 decimals** before the KS (DIYABC's text export has 6): small-amplitude columns
    (`F3v`, `F4v`, `NEIv`, `FST2v`) get flagged by quantization alone.
  - A `rdiff` stratified by a variable that moves the mean is not an effect until read in **absolute difference**.
  - A near-zero mean makes `rdiff` meaningless.
  - Small loci count (< ~50 per group) gives a small negative/positive global shift that shrinks with loci count: don't
    read a global bias off such a dataset (the "loci-count residual").
- **Reftable text formats are per column family, do not unify**: parameters `%.6g`, statistics `%12.8f`.
  `_split_mutation_and_stat_names` tells the families apart (a statistic starts with an uppercase ASCII letter, a mutation
  parameter never) and raises rather than shifting silently. After a format change, **restart the Jupyter kernel** before
  regenerating anything.
- **Text-reftable column order is the header trailer order** (`bintotxt` looks names up by name; `_historical_columns_order`),
  not prior-declaration order. The `.bin` uses the scenario's own `histparam` order, constants excluded.
- **The private msprime API behind the light replicates stays guarded and tested.** Keep the version guard,
  `test_compare_light_replicates_with_reference`, `test_light_replicates_resume_after_interrupted_runs` (forces
  `_LIGHT_EVENT_CHUNK` to 20 and 1), `test_mutate_light_tables_equal_public_path` and
  `test_light_a_replicates_falls_back_on_other_msprime_version`. A single `Simulator.run` returns 1 (`MAX_EVENTS`) on a long
  simulation and `finalise_tables` then writes a **truncated tree without any error**.
- **A test that writes a modified `header.txt` to `tmp_path` must build a fresh `*ReplayContext`** from that text:
  functions never re-read the header from disk.
- **A genepop `.mss` is accepted only if every line before the first `POP` is a locus** (`DataC::testfile`): a blank line
  counts as a locus, and every individual line with a different field count is refused.
- **A refactor touching one half of a drawing/replay pair must be checked against the other half** (the replay broke
  unnoticed twice with green tests). A naming change is checked against every reader of those names, not only the writer.
- **`parse_group_priors`**: a `MODEL` line must follow a prior of its group and have exactly 4 fields
  (`MODEL <name> <p_fixe> <gams>`), else a clear `ValueError`.
- **Comparison notebooks** (`reference/*/compare_reftables_*.ipynb`) choose columns **by name**, never by a hand-written
  `N_STATS`/`N_PRIORS` position constant (it silently cut trailing statistics). DNA-only replays of mixed headers must use
  `replay_diyabc_priors_mixed.py`.
- **Serial + `<X>`/`<Y>` is refused** up front (`raise_if_serial_with_sex_linked_loci`, predicate
  `scenario_parser.is_serial_scenario`: "the same population is sampled twice"). Don't rephrase it as a comparison with
  `npop` (a scenario may declare never-sampled populations). The `<X>`/`<Y>` sample-set builders rebuild msprime keys from the
  position of the observed block, so on a serial scenario they would silently attach samples to an ancestral population.
- **DIYABC SNP: a header with 2 SNP loci groups core-dumps DIYABC** (its own limit).

## Domain knowledge (C++ semantics the code can't tell you)

### SNP

- `<M>`/`<Y>` loci share one genealogy per type (`simulate_shared_ancestry_loci`, `particuleC.cpp:2422-2435`), each with its
  own dedicated seed offset (sharing one offset would make `<Y>` share the `<M>` tree).
- Non-`<A>` types: rescaled demography (`rescale_demography(coalescence_coefficient(heritage, sex_ratio) / 2)`), `ploidy=1`.
  `<A>` (ploidy 2, IndSeq and PoolSeq) is rescaled by `coalescence_coefficient("A", sex_ratio) / 4`, skipped at sex ratio 0.5.
  **We follow the real `<NM=...>`** (what a freshly converted dataset does in DIYABC); the factor is guarded by two pytest spies.
- PoolSeq reads go through `with_mrc_filter` against a global pool. On a POOL file `count_per_samples` counts **haploid gene
  copies**: `poolseq_counts_by_sample` (`// 2`) is the single definition used by `pipeline` and
  `simulate_poolseq_reads_with_mrc_filter`; `pool_sizes` for the read-bias correction stays in haploid units.
- Population layout: IndSeq computes it once, inside `simulate_snp_genotypes`, from the per-sample counts (the ploidy varies by
  heritage type, so a layout built upstream would be wrong). **No `layouts=` threading through the 130 statistics.** Both MAF
  loops have a `maf == 0.0` fast path that must forward the counts too (`<MAF=hudson>` takes only that path).
- **Never derive the per-sample counts from `samples`**: `samples` goes to msprime (dict or `SampleSet` list), the counts feed
  `compute_sample_layout`. An `<X>` locus describes one sample with two `SampleSet`s.
- A layout whose length disagrees with `ts.num_samples` yields silently wrong genotypes (membership test, not index).

### MicroSat

- **Header**: the `MEANxxx` → `GAMxxx` hierarchy is **positional, not name-based** (`header.cpp::readHeadersimGroupPrior`).
  A `GAMxxx` line's declared mean is overwritten at draw time with the value just drawn from the group's `MEANxxx`
  (`particuleC.cpp:819`).
- **GSM**: geometric step from `Pgeom` (`Pgeom=0` is SMM), clamped to `[kmin, kmax]`, ancestral state = midpoint.
  `msprime.TPM` with `m = 1 - Pgeom` reproduces it (epsilon clamp for `m=1`/`m=0`). **The TPM grid must be anchored on the
  root, not on `kmin`**.
- **SNI**: one Poisson process at `mut_rate + sni_rate`, each event a GSM step or ±1 bp SNI step with probability
  `sni_rate / (sni_rate + mut_rate)`; the matrix is **dense** (one row per integer in `[kmin, kmax]`); built by blocks per
  residue class (`n_alleles >= 2 * motif_size`, below it the old function raises too).
  `_gsm_matrix_numpy`: the vector of powers must be built **by numpy** (a Python list comprehension rounds differently).
- **Statistics**: `_length_by_sample_and_individuals` (for `FST`) duplicates a haploid allele into `(x, x)`;
  `_genotypes_by_sample_and_individuals` (for `LIK`) keeps true ploidy; everything else uses the flat `_length_by_sample`.
  `LIK` is the only asymmetric statistic (all ordered pairs). Aggregation families: mean of per-locus values, **ratio of sums**
  (`MGW`, `DAS`, `FST`), and `DM2`'s stateful variant. A locus with `num_sites == 0` is a real monomorphic locus (tskit has no
  site), special-cased in all helpers. **`np.bool_ + np.bool_` is a logical OR**, not a sum (the `FST` root cause): cast with `int()`.
- **`LIK` `nal`** counts only alleles whose total frequency over all samples is non-zero (`frt > 0.000001`).
- **`AML` (`cal_Aml3p`)**: bisection driven by the sign of the finite-difference slope at both bounds, three boundary cases
  (flat → uniform random `a`; always negative → `a=0`; always positive → `a=1`), no pseudo-count (zero-frequency term skipped),
  except a heterozygote with ONE zero frequency adds the log of the other (`pente_lik`, `sumstat.cpp:1252-1260`).
  Frequencies are hoisted out of the bisection (`_prepare_triplet_for_admixture`); seeds vary per triplet **and** per group.

### DNA sequences

- **Base frequencies come from the `.mss`** (`DataC::do_sequence`, empirical per locus, pooling populations, missing positions
  excluded); the `pi` literals in the header are for a no-observed-data mode.
- **`p_fixe`/`gams` = (proportion of invariant sites, gamma shape)**, unrelated to `k1`/`k2` (from `MEANK1`/`GAMK1`/
  `MEANK2`/`GAMK2` group priors, present for every `[S]` group whatever the model).
- **`k1`/`k2`/`mus_rate` are drawn in two tiers**; the per-locus tier needs `sdshape > 0.001` **and** the group's `nloc > 1`;
  `k2` genuinely has no `nloc > 1` check (an asymmetry of the C++, replicated via `check_nloc`).
- `gams == 0` is valid (`MwcGen::ggamma3` returns the mean): don't raise `ZeroDivisionError`.
- **`dnatrue` selects two sequence representations; the reftable path (`dnatrue = false`) builds an artificial sequence of
  variable sites only**, where a locus with no mutation keeps an empty string = `SEQMISSING`. **A verification on observed data
  does not validate the simulated path.**
- The Jukes-Cantor token is **`JK`**, not `JC`/`JC69`.
- `build_transition_matrix_dna` row-sum normalization needs `axis=1, keepdims=True`.
- **Mutations go through `msprime.sim_mutations` with a generic `MatrixMutationModel`** fed our own `matQ`/`pi`, **not** msprime's
  named models (their rate scaling allows silent self-transitions, unlike `comp_matQ`).
- DNA loci call `msprime.sim_ancestry` **directly**, not `simulate_independent_loci` (which hardcodes `sequence_length=1`).
  A MicroSat locus *is* `sequence_length=1` with a ~39-state alphabet.
- **`<X>`/`<Y>` sex inference**: no SEX column in `.mss`; everyone starts female, flips to male when an `<X>` locus shows a
  haploid genotype (even the missing placeholder) or a `<Y>` locus a present genotype. Missing token `<[]>` → regex `\S*`.
  Sample sets are dispatched **inside** the per-locus loop.
- **Statistics**: stateless per-locus brick + aggregator over one header group (never mixing groups). `VNS` is a **biased**
  variance (`ddof=0`), unlike `VPD` (`ddof=1`). `HST`/`MP2` are **ratios of sums** with accumulators **dicts keyed by pair**.
  **Column naming**: a `_<group>_` segment appears iff the **whole** header has more than one group, any type.
- **`SML`**: triplets from `_half_arrangements(n, 3)`; haplotype ids come from `_haplotype_ids_per_node` (one numbering per
  locus, equal id iff equal column of `genotype_matrix()`; compare equality only, never write an expected id in a test; sample
  nodes must be `0..n−1`); `cal_aml3p` is `AML` on haplotype ids, with the guard `sumstat.cpp:2100-2106`. `compute_all_statistics_dna`
  needs `seed` from the caller.
- **The discrete coalescent (`evalcriterium`, `coal_pop`) is not replicated, on purpose**: it is the same process as the
  continuous one (verdict in `results.md`). Don't reopen.

### Serial/temporal sampling

**A sample is not a population.** The scenario samples the same population at several dates; `read_events` counts `nsamp`
separately from `npop`; the data block `sa` maps to the `sa+1`-th `sample` line **positionally**, with no name cross-reference
(DIYABC guarantees nothing about this ordering across scenarios).

**THE INVARIANT**: node order = `SampleSet` order = `sample`-event order = `POP`-block order = `samples_default` key order =
statistics column indices.

- `build_sample_sets_from_scenario(scenario, values, ...)`: one `SampleSet(n, population, time)` per `sample` event, in order.
  The counts dict's **keys are ignored**, only value order matters (indexing by `f"pop{event.pop}"` gives every serial sample the
  same count). Call **after** the prior draw (a sample time may be a parameter name). Guards on `len(sample_events) != len(counts)`.
- `compute_sample_layout(ts, ...)`: slices `ts.individuals()` by cumulative counts (slicing `ts.samples()` by `count × ploidy`
  breaks on `<X>`). On non-serial data it equals `compute_population_layout`.
- Wiring: the stat functions of MicroSat/DNA take `layouts=` (zip `strict=True`); `compute_all_statistics_*` take
  `layouts_by_locus`, a dict **keyed by locus name**, never a flat list (dispatch filters and reorders loci).
  `pipeline.py` builds `SampleSet`s and layouts per particle.

### Replay (`_from_values`) chain

Every `_from_values` function is a **sibling** of its drawing counterpart, never a modification of it. Only the group-level
(tier 1) draw is replaced by DIYABC's value; the per-locus (tier 2) dispersion is never replaced (not recorded in the reftable).
`group_prior_column_names` is scenario-**independent** (don't reuse `_kept_param_names_by_scenario`); it **raises** on a SNP header
(`check_header_trailer_line` catches it deliberately). `ParticleResult.group_priors_values` is a flat `{column_name: float}` dict.
`values` is the one datum the replay variant *receives* and the drawing variant *produces*: copying the drawing shape drops it.

## Deliberately reproduced C++ bugs

Not to be "fixed": the goal is fidelity to real DIYABC output.

- **`MPD`/`VPD` denominator** (`particuleC.cpp:1777`, `cal_mpd1p`/`cal_vpd1p`): a locus with no mutation has an empty sequence =
  `SEQMISSING` and leaves the `MPD`/`VPD` denominator, while `cal_nss1p` keeps it. Reproduced by
  `if ts.num_sites == 0: continue` in `compute_MPD`/`compute_VPD` (commented with the lines to delete; guarded by
  `test_mean_pairwise_differences_excludes_loci_without_mutation`). Only `MPD`/`VPD` are affected, and beyond `θ ≈ 1` the defect
  vanishes.
- **`DM2` / `cal_dmu2p`**: `moy[]` is recomputed only when both populations have samples, but the consuming line sits outside
  that guard, so a stale `moy[]` from the previous locus is reused. Reproduced via an explicit `previous_moy` threaded between
  locus calls: `DM2` is the one statistic where the order of `tree_sequences` matters.
- **`mutsit`/`sitefix`** (`header.cpp:722-748`) — **the one exception, NOT reproduced (decision of the user, 03/11)**: DIYABC zeroes
  `mutsit` with the loop counter instead of `sitefix[i]`, so its invariant sites are always the first ones. `sample_site_rates`
  zeroes `rng.choice(dnalength, size=nb_sites_fixes, replace=False)` (the intended behaviour); the faithful line
  `mutsit[:nb_sites_fixes] = 0.0` is kept in a comment. **Keep it, do not "fix" it back.**

## Design decisions

- **Sibling duplication over parameterization.** SNP → DNA → MicroSat families are kept as near-identical siblings; the rule is
  "never modify an already-validated original". Don't refactor without asking. **One exception (user, 02/11): in `pipeline.py` the
  DNA and MicroSat paths share four helpers** (`_draw_common_data`, `_simulate_and_compute_statistics` and their `_from_values`
  counterparts), dispatched on `type_of_data` through `_MUTATION_BY_TYPE`, `_MUTATION_FROM_VALUES_BY_TYPE`, `_STATS_BY_TYPE`. Drawing
  and replay never share a helper (a `from_values` flag once hid three divergences), nor do the SNP/PoolSeq paths or
  `reftable_loop.py`. A new data family goes through the same dispatch dicts, not a copy.
- **No dedup across stat groups** in `stats_group_parser` (unclear whether two groups could legitimately declare the same column).
- A new function whose gain is under ~1 ms per particle is not worth adding.
- Every optimization touches code validated against DIYABC. An output-neutral one must stay **identical** at equal seed
  (`scripts/compare_particles_{dna,microsat,snp}.py`, then the golden tests); one that changes the random stream is validated by a
  paired replay read with KS, sign test **and a re-run with `seed_offset`**. Ask first (sibling rule).

## Measuring performance — traps

Benchmarks and verdicts: `notes/results.md` §3. Method:

- Profile **one particle in-process**, never the replay (`cProfile` sees only the main process). `cProfile` ranks, `perf_counter`
  decides: it misranks small-numpy-call code by x2-5. A synthetic micro-benchmark lies too: capture the real arguments with a spy.
  A micro-benchmark on objects already touched lies (`ts.tables` is cached): build fresh objects per repetition.
- Single-particle timings carry ±7 % noise: run 10-20 and take min/median. Judge an optimization on the **full replay at 16 workers,
  alternated against the reference, ≥ 5 times**. Measure `HEAD` (clean `git worktree add --detach`) and the modified tree
  alternately.
- SNP simulation is a **lazy generator**: time the call **and** the `list(...)` consuming it. MicroSat is eager.
- To time a statistic, patch the **dispatch dicts** (`_MICROSAT_*_STATS`), not the module attribute. A timer wrapper is not free,
  and nested wrappers double-count.
- **A drawing/replay pair must be optimized AND checked on both halves.** After sharing a computation, count calls with
  `pstats.Stats.print_callers`.
- Running DIYABC for a timing: copy the directory first (never run in `reference/`), pass `cwd=` the copy (`-p ./` is relative to the
  process), initialise the RNG with `-n "t:16"`, always `-R ALL`, check the log ends in `END|<nrec>|`. Before comparing, check the
  machine (`ps -eo pcpu,args --sort=-pcpu | head`) and the reference (`maf.txt`, `THIS->MAF=`); a stored `ref_*.json` is void once its
  dataset directory changed or a statistical change happened. Run one heavy computation at a time, with `nice` and a timeout; never
  `pkill -f` / `pkill python` (it killed the Jupyter kernel once).
- **Identity check**: `compare_particles_microsat save ref.json` on `HEAD` in a clean worktree, `compare ref.json` after; exact equality.
  When `reference/` is tracked, `ln -s` inside the worktree creates the link *inside* it: remove the checked-out directory first.
- **Tests of an output-neutral optimization**: (1) the oracle is the **old code copied into the tests**, imported from `bridge/` for
  the new one; (2) **at least two different loci**; (3) a cache needs a **call-count spy**; (4) run the **mutation check** on a *copy*
  of the module, never on the working tree. Compare with `==`, not `approx`.
- **Golden tests pin the seeded simulation, not the formula**: regenerating them proves nothing; a formula test is built on a hand-made
  `TreeSequence`. Test statistic bricks with real tskit dtypes.
- **Tooling**: tracked in this repository, `scripts/compare_particles_{microsat,dna,snp,snp_pool}.py` (drawing path only; **no identity
  script covers the replay paths**), `profile_particle_{microsat,dna}.py`, `time_aml_microsat.py`. The stage timer, the benchmark scripts
  and the `timings` plumbing live on branch **`bench/timings`** (worktree `../diyabc-msprime-poc_timings`). **Never merge it into
  `main`** (it makes `replay_reftable_simulation_*` return a 5-tuple and breaks the replay scripts); bring `main` into it with
  `git merge main`; run replays from `main`; don't recopy `bridge/` between directories.

## Closed — do not reopen (verdicts only; evidence in `notes/results.md` §2)

- Residual statistical bias vs DIYABC: it is the **loci-count residual** (shrinks with loci count); not a simulator bug.
- "Simplified substitution model below some threshold": **does not exist in the C++**.
- Discrete coalescent (`evalcriterium`) as the cause of any discrepancy: **falsified by direct measurement**; not replicated.
- SNI and admixture as the cause of the `FST`/G3 gaps: falsified. The `FST` gap was `np.bool_ + np.bool_`.
- G2 `MPD`/`VPD` residual: a **DIYABC bug**, reproduced (see above).
- `DTA_2_2`, `ML3p_2.3.4`, `AML_2_3.1.2`: correlated noise from overlapping seeds (fixed by `_locus_seed`) or seed noise.
- `LightRateMap`: tried and removed. Reusing one `Simulator` across batches: falsified (changes the stream).

## `scripts/`

The `.py` files under `scripts/` are ad hoc investigation scripts, not part of `bridge/` nor its tests. Treat them as disposable
when reasoning about the architecture, but **don't delete them without asking the user** (they double as experiment logs).
Open work and pending validations: `notes/results.md` §4.
