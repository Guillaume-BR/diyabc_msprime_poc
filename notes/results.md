# results.md

Results, verdicts and measurements moved out of `CLAUDE.md` on 04/11 (moved verbatim, not rewritten).
Narrative of each investigation: `notes/exploration.md`. Rules live in `CLAUDE.md`.

# 1. Validation status (was "What this is")

### What this is

A proof-of-concept replacing DIYABC's C++ coalescent simulator
(`particuleC.cpp::dosimulpart`) with `msprime`. The pipeline (`bridge/`)
turns a DIYABC `header.txt` + observed data file into a `reftable.bin`
entirely in Python — no subprocess call to the DIYABC binary — and must be
structurally and statistically equivalent to the real DIYABC's output.

| Data family | State | Validated against |
|---|---|---|
| SNP IndSeq | complete, 130 stats | `human` (5000 loci), `toy_example3`, `toy_example5` |
| SNP PoolSeq | complete | `toy_example4` |
| Heritage `<A>/<H>/<X>/<Y>/<M>` (SNP) | complete; real sex ratio followed, `<A>` rescaled (22/10) | `human`, `toy_example5`; per type and mixed on `toy_example5_with_maf` against DIYABC freshly converted (sex ratio 0.3): 0/50 flagged for `<A>`, `<X>`, `<Y>` and the mix |
| MicroSat | complete, 11 stats + `AML` | `toy_example2_ms_dna_K2P`, `toy_example1_ms_modified` |
| DNA sequences | complete, 13 stats + `SML` (triplets, only with ≥ 3 samples); `JK`/`K2P`/`TN` all covered | `toy_example2_ms_dna_{JK,K2P,TN}` (5 loci/group, 2 runs each), `toy_example2_ms_dna_50loci_K2P` (50/group, 2 runs) — the **G2 `MPD`/`VPD` residual is a DIYABC bug, located to the line** (its `MPD` denominator drops loci the simulation left monomorphic), see Open work |
| DNA sequences `<X>`/`<Y>` | implemented | synthetic fixture only — the real binary SIGSEGVs on this case |
| Serial/temporal sampling (MicroSat) | complete | `toy_example1_ms` (1 population, 4 sampling times) |
| Serial/temporal sampling (SNP IndSeq) | complete | `human_seriel` (1 population, 4 sampling times, 130 stats, 0/130 significant) |
| Serial/temporal sampling (SNP PoolSeq) | complete | `toy_example4_seriel` (100 loci, `<MRC=5>`, 4/133 vs ~6.7 expected) |
| Mixed MicroSat + DNA in one header | complete | `toy_example2_ms_dna_TN` (73 columns, 2 replays × 2 scenarios) |
| Serial MicroSat + DNA, 4 samples (`AML`, `SML`) | complete, 224 columns (24 `SML`) | `toy_example2_ms_dna_50loci_TN_seriel`: 3 independent DIYABC reftables, several `seed_offset`s, 12 cells: 4 `SML` flags out of 288 columns, none repeated; MicroSat G1 (10 loci) shows a same-sign deficit in all four cells of two reftables, both scenarios (median −0.75 to −1.91 % over its 76 statistics, 54-63/76 negative, 9 of the 12 KS flags), DNA groups −0.70 to +0.47 %: probably the loci-count residual, **not demonstrated, open** (02/11, 04/11) |
| Non-serial 4 populations with real admixture, MicroSat + DNA (`AML`, `SML`) | complete, 224 columns | `toy_example2_ms_dna_50loci_K2P_4pop`: ONE DIYABC reftable, 2 scenarios, 0 KS flag out of 233 statistics each, `SML` smallest p 0.07, sign test significant in scenario 2 at −0.27 % on the first replay; a second replay (04/11; whether it was a new DIYABC reftable was not recorded) gives 96/222, p = 0.051, median +0.122 %, KS clean: the −0.27 % was noise (03-04/11) |
| Admixture with k ≥ 2 ancestral populations (`split`, no DIYABC counterpart for k > 2) | implemented (06/11), rates are absolute probabilities, last ancestor = complement | k = 2: paired replay against the real reftables, **every value matches**; k = 3: no DIYABC reference exists, checked by a χ² of the lineages' destinations against the declared rates (`test_buiild_demography_with_three_ancestral_populations`, positive and permuted-rates twin). Dirichlet prior for the rates: implemented (12/11), see Open work |

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

**Validation status after the 13-14/10 changes (`<A>` replicates, shared `<M>` tree,
seed-independent replays).** The `<A>` batching changes every DNA and MicroSat
value at equal seed, so these were re-replayed against the real DIYABC:
`toy_example2_ms_dna_50loci_{JK,K2P,TN}` (two replays each, no KS flag; `TN` is the
mixed header), the non-serial MicroSat sets (no divergence reported) and
`toy_example1_ms_seriel`. On the latter, `AML_2_3.1.2` and `AML_2_1.3.4` were flagged
with the base seeds and both **vanished with `seed_offset`** (p = 0.005 → 0.79 for
the first); the sign test stays at −0.4 % to −1.3 %, unchanged by the batching (the
known loci-count residual, 50 loci here). SNP code is untouched by the batching.
Statistics optimizations of the same days are output-identical and were checked by
`scripts/compare_particles_dna.py` / `compare_particles_microsat.py`, not by replay.

**Validation status after the 19-21/10 changes (IndSeq light replicates, single tree
stream).** The no-filter branch (`maf == 0.0`) is **bit-identical** to the code before:
`compare_particles_snp` gives 0 differences on `human` (30 particles, 3900 values),
`human_seriel` (10, 1300) and `toy_example5` (30, 1530; it exercises `<X>`). The MAF branch
is identical up to the first batch of attempts and **changes the tree stream beyond it**
(one stream per particle instead of one fresh `Simulator` per batch of
`max(20, num_loci // 4)`); validated by a paired replay against a DIYABC reftable generated
**with the filter on** (`toy_example5_with_maf`, 1000 particles, 50 statistics, rectangular
read): 11 KS flags (p < 0.05) after the change, 12 before, both led by the `FST*m` family, same
sign-test median (−3.6 % / −3.5 %), i.e. neutral. The reference reftable first used for it was wrong, see Closed
investigations ("`maf.txt`"). **Final MAF code (light replicates, lean `run`, single stream, `_derived_mask` prefix, `<A>`
rescale) re-replayed on 22/10 against a regenerated `toy_example3_500loci` reftable** (500 `<A>`, `maf = 0.05`, `NM=1`,
converted in the same run: `maf.txt`, `.snp.bin` and `reftableRF.bin` share the minute; 1000 particles, six scenarios,
priors paired): 0 KS flag out of 130 statistics (smallest p = 0.20), sign test 63/130 (p = 0.79), median +0.01 %
by my count; the user reports 0 flag, 55/127 (p = 0.155), +0.070 % on 127 columns (the 3 not kept are the `F4v` columns, cut by the notebook's hand-written `N_STATS`, see
Closed investigations, 04/11). One reftable, no `seed_offset` replica. **`human_seriel`** (5000 `<A>`, one population at four dates,
`MAF=hudson`, regenerated 22/10 with `.snp.bin`, `maf.txt` and the reftable in the same minute; two scenarios, 1000
particles, priors paired) with the final no-filter code: 0 KS flag out of 130 (smallest p = 0.078), sign test 50/130
negative (p = 0.011), median `rdiff` +0.04 % by my count; the user reports 41/130 (p = 3.1e-5), +0.075 %. The sign test
is significant but the amplitude is under 0.1 % (median |`rdiff`| 0.075 %), the same order as the +0.05 % of the
06/10 validation: the loci-count residual, positive-signed here, not a defect.

**Validation status after the 22/10 change (`<A>` rescaled by the real sex ratio).** Decision: follow the
real `<NM=...>` (what a freshly converted dataset does in DIYABC). Neutral at sex ratio 0.5, so `human`,
`human_seriel`, `toy_example3_500loci` (IndSeq) and `toy_example4`, `_hudson`, `_seriel` (PoolSeq,
`scripts/compare_particles_snp_pool.py`) stay **bit-identical**; `toy_example5` changes (1482 of 1530 values).
Paired replay against DIYABC freshly converted (no `.snp.bin`, sex ratio 0.3), `toy_example5_with_maf`,
1000 particles: `<A>` −6.7 % (37/50 flagged) → **−0.03 % (0/50)**, `<X>` −0.02 % (0/50), `<Y>` −0.2 % (0/50),
`<M>` +1.5 % / +0.5 % / +0.1 % on three independent reftables (2, 0, 0 flagged), and the **mixed header**
(350 `<A>` + 50 each) −0.10 % (0/50, sign test 34/50, p = 0.015). One reftable per cell except `<M>`.
**PoolSeq** at a sex ratio ≠ 0.5: `toy_example4` with its `<NM=...>` edited, DIYABC freshly converted (no `.snp.bin`),
recorded 04/11 (sex ratio 0.2308, 1000 particles, 100 loci, one DIYABC reftable): scenario 1 65/127 negative (p = 0.86), median −0.076 %, KS clean; scenario 2 44/127 (p = 6.9e-4), median +0.636 %, 3/140 KS flags (`FST1m`, `FST1v`, `HBm`); a second replay gave 35/127 (p = 4.4e-7), median +0.986 % (the scenario was noted ambiguously, probably 2). **Open**: no baseline at sex ratio 0.5 on the same 100 loci, see Open work. The factor itself is
guarded by two pytest spies (`..._rescales_autosomal_demography_by_sex_ratio`, IndSeq and PoolSeq, mutation `/ 4` →
`/ 2` checked crosswise); the 0.5 case (sizes untouched) is not covered by them.

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


# 2. Closed investigations

Each line is a verdict; the full investigation is in `notes/exploration.md`
under the date given.

- **Residual statistical bias (17/07)** — not a simulator bug. Shrinks and
  vanishes as loci count grows; 0/130 significant on `human` at 5000 loci.
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

- **`<A>` genealogies batched with `num_replicates` (13/10) — statistically equivalent,
  not identical.** One `sim_ancestry(num_replicates=n)` per `(heritage, length)` group
  replaces `n` calls (about 1.25 ms of fixed msprime cost each, against 0.15 ms
  batched, measured on the real demography; the SNP path already used this
  technique, `simulate_independent_loci`). Replicates draw from one stream, so
  **every DNA and MicroSat statistic changed at equal seed** (44 golden values in 28
  tests regenerated by a script that rewrites each literal from the actual value,
  median relative move 12 %, 90th percentile 78 % on 5-10 loci); validated by replay, see the validation
  status at the top. Two things it must keep: the batch index (previous bullet under
  Hard rules) and the iteration: the generator is lazy, consumed with
  `zip(loci, replicates, strict=True)` inside the helper, so a locus count that
  disagrees with the replicate count raises instead of shifting trees. The `<M>`
  shared tree (one call per `(heritage, length)`) is the exact half of the change:
  bit-identical to the old per-locus calls (`compare_particles_*` 0 differences).
- **`AML_2_3.1.2` / `AML_2_1.3.4` on `toy_example1_ms_seriel` (14/10) — correlated
  noise, not a defect.** Flagged at p = 0.005 (run 0, scenario 2) and p = 0.001 (run 1,
  scenario 1) with the base seeds, **identically before and after the batching**
  because they sit in G2, the 5 `<M>` loci, which are bit-identical to the old code
  (86 of 162 columns of the replay identical). With `seed_offset` 100 000 / 200 000
  they disappear (p = 0.79 for the first) and other columns appear (`MGW_2_3`,
  `MGW_2_4`, `FST_1_1.2`, `DM2_2_3.4`, p between 0.007 and 0.047, none present with
  the base seeds), well below the ~7.6 flags per 152 columns expected at 5 %. G2 has
  5 loci sharing one tree: one random realization per particle.

- **IndSeq light replicates and single tree stream (19-21/10) — adopted.** `human` particle 1.52 s →
  1.27 s, `toy_example3_500loci` particle ~265 → ~136 ms, both validated as in the status above. The
  levers, with what each gave and why: `ts.tables` rebuilds a whole `TableCollection` on the **first**
  access of a fresh `TreeSequence` (`_draw_single_mutation_edge_child` 35 → 13 µs by reading
  `ts.edges_child` / `edges_parent` / `nodes_time`; −16 % of a particle); the tree objects
  (`finalise_tables`, `asdict`, `fromdict`, `load_tables`, `TreeSequence.__init__`) were ~45 % of a
  tree (125 → ~60 µs); calling `Simulator.run` directly skips `model._as_lowlevel()` (`asdict` +
  `deepcopy`, −7 µs per tree); `_derived_count` returns 1 for a sample node (45 % of draws);
  one stream instead of batches saves a `Simulator` creation (~1.8 ms, mostly C allocation) per batch,
  ~21 ms of 130. **Falsified**: reusing one `Simulator` across batches with
  `random_generator.seed = s; reset()` does **not** reproduce a fresh `Simulator`'s trees (~10 % of
  replicates differ, even on a never-used simulator: allocation state decides which lineages are
  drawn); `sim_ancestry` options (`discrete_genome`, `model`) change nothing; the batch size
  (125 vs optimum ~240) is worth ~2 % and moves the stream.
- **"50 of 57 columns KS-flagged" on `toy_example5_with_maf` (21/10) — the DIYABC reference had no
  MAF filter.** Its folder still held the `.snp.bin` and `.snpbin.txt` of an earlier run, so DIYABC did
  not convert the `.snp`, wrote no `maf.txt` and ran Hudson (statistics of a filter-less run: `HWm_1` 0.088 against 0.142-0.149 with the filter, DIYABC 3 s instead of 7 s).
  Same inputs, regenerated with `maf.txt`: our replay agrees. The "7.2x" timing was the same
  artefact (filter-on msprime against filter-off DIYABC); the real ratio is 3.3x.
- **"Persistent −3.5 % on `toy_example5_with_maf`" (21/10) — not established.** One DIYABC reftable
  read with three msprime seed sets; a second DIYABC reftable shows no divergent column. The per-type
  `<M>`/`<Y>` residual seen next was DIYABC's sex-ratio mode, not a defect of ours (Open work, 22/10).

- **Per-type `<M>`/`<Y>`/`<A>` residuals on `toy_example5_with_maf` (21-22/10) — DIYABC's sex-ratio mode plus
  one real defect (`<A>` not rescaled), fixed.** `<M>` −7 % / `<Y>` +9 % came from comparing our 0.3 with a
  DIYABC that had loaded a `.snp.bin` (forced 0.5): with `<NM=1NF>` on our side 6 reftables fall to 0-1 flagged
  columns of 50. The real defect was `<A>` at a sex ratio ≠ 0.5 (pair rate `2 / (16 s (1 − s) N)` in DIYABC
  against our `1 / (2 N)`), −6.7 % against a freshly converted DIYABC; fixed and validated (status above).
  Claiming "the coefficients equal DIYABC's" from our reading of the file and from the `sexratio=` line of the
  DIYABC log was the error: the value that counts is the one applied after `loadfromfile`.

- **DNA serial "x2 against x7" (02/11) — not a defect of the serial path: unequal loads, then DIYABC's cost.**
  `-R ALL` made DIYABC compute 224 statistics (24 `SML`) while we produced 200. At equal load: DIYABC 20.5 s,
  msprime 50.7 s (2.5x), 44.3 s after the `SML` optimisation (≈ 2.2x); DIYABC ×13 between 2 and 4 samples, us ×4
  (our floor is the ~100 `sim_mutations` calls per particle). DIYABC's own variance on this dataset is ±12 %.
  `SML` ≈ 28 % of DIYABC's time (`gdb` stub, 2 runs); `AML` nil. Falsified: discrete coalescent (`Npresent ~
  UN[5000,10000]` leaves DIYABC at 20 s), parasitic measurement. The 4-population non-serial control was never built.
- **`_log_likelihood_admixture` heterozygote with one zero frequency (23/10) — real infidelity, fixed.** The C++
  adds `log(fx)` or `log(fy)` for the non-zero one (`pente_lik`, shared by `cal_Aml3p` and `cal_aml3p`); we added
  0. 344 of 720 `AML` values change on `ms_modified`, nothing else. Replays on `toy_example1_ms_seriel`: the `AML`
  flags of 13/10 vanish (smallest `AML` p 0.001 → 0.056); only `AML_2_3.1.2` keeps a stable sign (+2.3 to +3.8 %,
  p 0.028-0.162 over three seed offsets) on ONE DIYABC reftable, absent on the other, and not flagged on a further
  replay (04/11): noise. Beware: the replay script's default `SEED_OFFSET` is 200 000, not 0.
- **Site rates drawn with numpy (02/11) — same law, new stream, validated by replay.** `sample_site_rates` takes a
  `np.random.Generator` (one `default_rng(_locus_seed(seed, _SITE_RATE_SEED_OFFSET, 0))` per particle) and
  `build_rate_map` works on arrays: `build_rate_map_per_locus` 12.2 → 6.0 ms, DNA particle ≈ −5 %. **Every DNA
  statistic changed at equal seed**: 15 golden tests were regenerated after an independent check with tskit
  (NSS, NHA, MPD, DTA, PSS, NS2, NH2, MPB), and anything produced before 02/11 (reftables, stored comparisons) is no
  longer comparable at the bit level. Replay, 4 datasets × 2 scenarios: KS flags against DIYABC unchanged, new-vs-old
  mean of the relative differences on the DNA columns +0.016 % (sd 0.077 % between cells, t = 0.6). Trap: tskit's
  `segregating_sites(mode="site")` counts each derived allele of a multi-allelic site (8 where `cal_nss` counts 7
  sites): do not use it to check `NSS`/`NS2`/`DTA`.
- **A genepop `.mss` is accepted only if every line before the first `POP` is a locus (03/11).** `DataC::testfile`
  (`data.cpp`) counts as loci ALL the lines between the title and the first line containing `POP` (case-insensitive,
  so a blank line counts), then refuses any individual line whose number of fields after the comma differs
  (`Unreckognized file format`, printed right after the dataset title, before the header is read). Seen on a hand-built
  4-population file: one blank line before `POP` (111 "loci" instead of 110, every individual refused) and 8 duplicated
  rows at the start of the first block. Check with `testfile`'s own rule before blaming the header.
- **`LightRateMap` (03/11) — tried and removed: no gain on the real replay.** Skipping the `msprime.RateMap` object
  (31 µs a piece) on the light mutation path saved 4 ms per particle in isolation, outputs byte-identical, but
  the 1000-particle replays, 5 alternations against `HEAD`, came out +0.52 s (mixed) and +0.22 s (DNA) slower, never
  faster. **Judge an optimisation on the full replay at 16 workers, alternated against the reference, ≥ 5 times**: a step
  timing and a single-process particle both missed what the replay showed.
- **`µ`/`mu` prefix of `group_prior_column_names` (03/11) — fixed.** `_group_prior_columns` wrote `µ` in hard, so a
  `mu`-style header (`toy_example1_ms_*`) got column names matching nothing in DIYABC's reftable. `_detect_mu` now reads
  the spelling from the mutation-parameter tokens of the header's trailer line (`mu` if one starts with `mu`, fallback
  `µ` when the trailer is missing, as in a fresh `header.txt` before the first `-R ALL`). Cosmetic: the `.bin` carries
  no column names and `parse_real_reftable_params_with_group_priors` reads by position, so a wrong spelling never
  caused a wrong *read*; the real exposure is a wrong *count*, guarded by `check_header_trailer_line` and
  `check_real_reftable_matches_header`. Guarded by `test_detect_mu` and
  `test_group_prior_column_names_follows_trailer_spelling` (expected names written in hard).
  **Regression caught 04/11 by the first replay of `toy_example1_ms_seriel` (`KeyError: 'µmic_1'`), not by the 635 tests:**
  `_group_prior_values_{dna,microsat}_from_columns` (`ancestry_simulation.py`) looked up `µseq_N`/`µmic_N` by a
  hand-written name in the dict the replay builds from `group_prior_column_names`, which now says `mumic_1`; every
  replay-path fixture used a `µ` header (the `te2` family). Fixed by moving `_MU_VARIANTS` and `_normalize_mu` to
  `configuration.py` (`reftable_loop` imports `ancestry_simulation`: the reverse would be circular) and normalizing the
  keys at the top of both readers; `test_group_prior_values_*_from_columns_mu_spelling` (`mu` header on MicroSat and
  DNA, `µ` unchanged), mutation-checked on a copy (the two `mu` tests go red). **Lesson: a naming change is checked
  against every reader of those names, not only the writer; see the replay-path hole in Open work.**
- **`parse_group_priors` took any unknown line for a `MODEL` line (03/11) — fixed.** The `else` branch ignored the
  first word and read `last_group_name` only from the prior branch, so a misspelt `MODEL`, a missing or non-numeric field
  gave a raw `IndexError`/`ValueError` or a silent misread, and a `MODEL` right under a `group` line was filed in the
  PREVIOUS group (no error, the worst case). Now `last_group_name` is set on the `group` line, a `MODEL` line must follow
  at least one prior of its group, have exactly 4 fields (`MODEL <name> <p_fixe> <gams>`, last two numeric) and start with
  `MODEL`, otherwise a clear `ValueError` quotes the line. Guarded by `test_parse_group_priors_model_error` (the four cases).
  Not checked: whether DIYABC tolerates extra words after `gams`, we refuse them.
- **`_genotypes_by_sample_and_individuals` and the `np.int64` leak (03/11) — audited, none.** The sizes go through
  `np.array(decoded[0])[decoded[1]]` and `.tolist()` before the tuples are built (monomorphic locus: `(0, 0)` literals),
  and `_decode_locus` uses `int(a)`, so every element is a Python `int` (since the 15-16/10 array reads).
  `test_genotypes_by_sample_and_individuals_matches_reference` asserts `type(x) is int` on every element; mutation check
  on a copy of the module (`.tolist()` → `list(...)`): 8 of the 16 polymorphic cases go red (the 8 monomorphic ones never
  reach `.tolist()`).
- **Comparison notebooks cut statistics silently (04/11) — fixed in all 23.** Each `reference/*/compare_reftables_*.ipynb`
  took its columns by POSITION with a hand-written `N_PRIORS`/`N_GROUP_PRIORS`/`N_STATS`, and "colonnes écartées : 0"
  compared with that constant, not with the header. When the constant was lower than the declared count the LAST columns
  vanished, with no error: `F4v` (3) on `toy_example3*`/`toy_example4*` (127 of 130), `F3v`/`F4m`/`F4v` (10) on
  `toy_example4_seriel` (120), `AML_2_*` (6) on `toy_example1_ms_seriel` (146 of 152), `SML_3_3.2.4`/`SML_3_4.2.3` on
  `TN_seriel` and `K2P_4pop` (222 of 224); `toy_example2_ms_dna_50loci_JK` also carried `N_GROUP_PRIORS = 9` (the `TN`
  value, JK has 5), which filed four DNA statistics under "group priors". Now the columns are chosen BY NAME
  (`parse_requested_statistic_names`, `group_prior_column_names` normalized by `_normalize_mu`, `[]` for a SNP header,
  priors = the rest); an `assert` fails when DIYABC's reftable lacks a declared statistic and a warning lists the ones
  the replay lacks. The earlier "0 flag on 127 columns" statements therefore exclude those columns, and the stored
  outputs and `comparaison_summary_*.csv` are stale until each notebook is rerun. Checked by executing every notebook on
  a copy of its folder (20 run to the end; `toy_example2_ms_dna_JK`, `_XY` and `toy_example5_1000loci` lack their replay
  file). **The `toy_example2_ms_dna_50loci_{JK,K2P}` and `ms_dna_{JK,K2P,XY}` replays are DNA-only (42 of 58 statistics,
  the 16 MicroSat G1 ones missing): they validate the DNA groups only; rerun them with
  `replay_diyabc_priors_mixed.py`.** `scripts/replay_diyabc_priors_{snp,ms,dna,mixed}.py` now use `STATS_FILTER =
  "HEADER"`; first run (04/11): `snp` and `mixed` clean, `ms` raised the `µ`/`mu` regression above, `dna` raised
  "statistiques non calculées" on the 16 MicroSat statistics of a mixed header (the guard doing its job: use `mixed`).
  `replay_diyabc_priors_dna_50loci.py` (obsolete, folder renamed) was deleted by the user.
- **`SML` (02/11) — implemented and validated.** `cal_aml3p` is `AML` with one haplotype per allele (`cal_freq`),
  so `_compute_AML_one_triplet` is reused unchanged on haplotype ids; the only addition is the guard of
  `sumstat.cpp:2100-2106` (statistic = 0.5 when no locus shows two distinct haplotypes among the two parents
  together). Details, tests, validation table: `notes/exploration.md`, 23/10-02/11.


# 3. Performance — measured state

## Benchmark of 04/11 (latest; the 22/10 table below is kept for history)

Même protocole que ci-dessous, une exécution par ligne. DIYABC : `benchmark_diyabc_results.csv` ; msprime :
`benchmark_results.csv` (hors écriture du reftable).

| Case | Directory | Loci | DIYABC | msprime | msprime / DIYABC | msprime sim / stats (sum over workers) |
|---|---|---|---|---|---|---|
| IndSeq Hudson | `human` | 5000 | 132.1 s | 185.7 s | 1.41x | 98.9 / 1.1 % |
| IndSeq MAF > 0 | `toy_example3_500loci` | 500 | 5.95 s | 23.5 s | 3.95x | 98.2 / 1.8 % |
| IndSeq serial | `human_seriel` | 5000 | 156.3 s | 173.6 s | 1.11x | 98.8 / 1.2 % |
| PoolSeq Hudson | `toy_example4_hudson` | 100 | 28.8 s | 20.6 s | 0.71x | 98.7 / 1.3 % |
| PoolSeq MRC > 0 | `toy_example4` | 100 | 96.0 s | 42.9 s | 0.45x | 99.4 / 0.6 % |
| PoolSeq serial | `toy_example4_seriel` | 1000 | 600.0 s | 238.5 s | 0.40x | 99.8 / 0.2 % |
| DNA only (42 stats) | `toy_example2_dna_50loci_TN` | 100 | 1.81 s | 11.6 s | 6.43x | 59.6 / 40.4 % |
| DNA only serial (148 stats) | `toy_example2_dna_50loci_TN_seriel` | 100 | 20.7 s | 37.3 s | 1.80x | 27.6 / 72.4 % |
| MicroSat | `toy_example1_ms_modified` | 50 | 2.59 s | 9.9 s | 3.84x | 43.6 / 56.4 % |
| MicroSat serial | `toy_example1_ms_seriel` | 50 | 2.67 s | 10.6 s | 3.96x | 40.9 / 59.1 % |

Les lignes DNA (jeux DNA seul, 100 loci, 42 et 148 statistiques) et la répartition sim / stats DNA et MicroSat ne sont pas
comparables avec celles du 22/10.

## Performance — measured state (**full 1000-particle benchmark of 22/10**, one run per row; `benchmark_results.csv` on `bench/timings`)

Verdicts only; each case below is meant to be worked in its **own session**.
Narrative and raw numbers: `notes/exploration.md`.

**Protocol.** Replay of the 1000 real DIYABC particles of each `reference/`
directory (`replay_reftable_simulation_*`, `max_workers=16`) against
`time ./diyabc -p ./ -R ALL -r 1000 -g 1000 -m -t 16`, same machine (8 physical
cores / 16 threads). Every DIYABC run was checked: `reftable.log` ends in
`END|1000|`. msprime wall time **excludes** the reftable write. **One run per
case: the variance is not measured.**

| Case | Directory | Loci | DIYABC | msprime | msprime / DIYABC | msprime time: sim / stats |
|---|---|---|---|---|---|---|
| IndSeq Hudson | `human` | 5000 | 116.5 s | **157.2 s** (238.4 s before 21/10; 291.6 s on 07/10) | **1.35x** | 98.8 / 1.2 % |
| IndSeq MAF > 0 | `toy_example3_500loci` | 500 | 5.3 s | **19.7 s** (39.7 s before 21/10) | **3.72x** (6.7x before) | 98.2 / 1.8 % |
| IndSeq MAF > 0, `<X>/<M>/<Y>` | `toy_example5_with_maf` | 500 (350 `<A>`, 50 each) | 7.0 s (`-t 16`; 3.0 s **without** the filter, see Hard rules) | 23.0 s and 23.7 s (26.8 s at 8 workers) | **3.3x** | not re-measured |
| IndSeq serial | `human_seriel` | 5000 | 146.2 s | **150.0 s** (228.7 s before 21/10) | **1.03x** (1.6x before) | 98.7 / 1.3 % |
| PoolSeq Hudson | `toy_example4_hudson` | 100 | 26.2 s | 19.2 s | 0.73x | 98.6 / 1.4 % |
| PoolSeq MRC > 0 | `toy_example4` | 100 | 86.7 s | 40.3 s | 0.47x | 99.4 / 0.6 % |
| PoolSeq serial | `toy_example4_seriel` | 1000 | 524.9 s | 219.8 s | 0.42x | 99.8 / 0.2 % |
| DNA (mixed MicroSat + DNA, 58 stats) | `toy_example2_ms_dna_50loci_TN` | 110 | 2.05 s (1.6 s on 22/10, same machine, other day) | **18.3 s** mixed replay on 02/11 (DNA-only replay: 13.86 s, 100 loci and 42 stats: **not the same load**; 12.3 s on 22/10 was DNA-only against DIYABC's full 110 loci / 58 stats, so the old 7.57x flattered us) | **≈ 8.9x** at equal load before the 02/11 light `<A>` mutations, **≈ 6.3x** after the three 02/11 changes (SML early return, light `<A>` mutations, numpy site rates): mixed replay **13.14 s** against DIYABC 2.1 s, measured by the user on 03/11; DNA-only replay 11.81 s = 5.6x (not the same load). Path: 8.9x (18.3 s) → 7.1x after the light mutations alone (14.5 s) → 6.3x. ±1 with DIYABC's ±12 % variance; the same replay varies by ±10 % within an hour on this machine | not re-measured |
| DNA serial (mixed, 224 stats incl. 24 `SML`) | `..._50loci_TN_seriel` | 110 | 20.5 s (16-20.5 s over 4 runs; the 1.6 s of 07/10 had 58 statistics, `-R ALL` gives 224 from 4 samples) | **44.3 s** (50.7 s before the 02/11 `SML` optimisation; 30.8 s was without `SML`: unequal load) | **≈ 2.2x** (2.4-3.1x before, DIYABC variance ±12 %) | not re-measured |
| MicroSat | `toy_example1_ms_modified` | 50 | 2.1 s | **9.7 s** (20.4 s on 12-14/10; 228.7 s before) | **4.59x** (8.9x before) | 49.9 / 50.1 % |
| MicroSat serial | `toy_example1_ms_seriel` | 50 | 2.2 s | **10.3 s** (23.6 s on 12-14/10; 256.7 s on 07/10) | **4.66x** (9.8x before) | 47.2 / 52.8 % |

The 14/10 morning benchmark (`benchmark_results.csv` on `bench/timings`) ran the code
with the `<A>` batching and the shared `<M>` tree, **before** the DNA statistics work of
the same day; the DNA 13.4 s is a later single run after it. The SNP gains are not
isolated (candidates: `record_provenance=False` on every `sim_ancestry` /
`sim_mutations`, and the sample layout read by array). One run per row.

- **Orchestration is not the problem.** In all ten cases
  `(sum sim + sum stats) / 16` lands within 0.5-1.5 % of the wall time, so
  pickling the context, worker start-up and queueing cost ~1 %. **Don't
  re-investigate `max_workers`** (16 beat 8 and 7 on `human`) **nor the SNP stat
  formulas** (at most 1.3 % of the time in every SNP case).
- Per-particle seconds measured at 16 workers are inflated by hyperthreading:
  trust shares and ratios, not the absolute cost of a particle (use a
  `max_workers=1` run for that).
- **IndSeq (21/10), verdict: done.** `human` 1.8x → **1.35x**, MAF > 0 6.7x → **3.72x**, serial 1.6x → **1.03x**
  (benchmark 22/10; levers in Closed investigations). What is left per particle is msprime's own: on `human` `Simulator.run` 43 %,
  `finalise_tables` + `reset` + `asdict` ~23 %, `_derived_mask` 11 % before its fix (the only part that is ours);
  on the MAF loop `Simulator.run` 26 %, table conversion 18 %, `_locus_seed` + `Random(seed)` ~10 %
  (they define the stream, do not touch). **Per core, DIYABC does a particle of `toy_example5_with_maf`
  in ~48 ms of CPU and we need ~166 ms (3.5x)**; one tree costs us ~110 µs of which the C
  coalescent is ~50, DIYABC's whole pipeline ~32. Going lower means leaving msprime. **Done 22/10**:
  `_derived_mask` returns a one-node mask for a sample node and scans only the edges whose parent is not
  younger than the mutated node (a prefix, since edges are sorted by parent time; `side="right"` is
  load-bearing: the node's own child edges have its time): 30.7 → 16.5 µs per locus (leaf case alone
  saves as much as the prefix), masks identical on 800 loci, `human` particle 1.18-1.23 → 1.09-1.12 s
  (−5 to −9 %), nothing measurable on the MAF loop (the mask only runs for the ~500 accepted loci).
  **Measured and dropped**: the same prefix in `_derived_count` — slower on 80 genes (3.7 → 4.1 µs),
  break-even near 160 genes (318 edges, one constant-size population), +17 % at 240 genes on the
  4-population `human` demography, where the function is not even called; ~2 ms per particle at best
  for a threshold that depends on the demography. **Not done**: reusing one `Simulator` across batches
  (−10 %, changes the stream, see the falsified reuse above). 8 physical cores: 16 workers gain only 8-17 % over 8, for DIYABC (7.6 → 7.0 s)
  as for us (26.8 → 23 s); **no parallelism problem**.
- **SNP IndSeq: the simulation is everything (07/10 state).** On `human`, ~93 % of a particle
  is in `simulate_snp_genotypes`, ~40 % of it the per-locus `TreeSequence`
  overhead (population-layout caching already done). **MAF > 0 is the outlier
  (8.3x against 2.2x)**: the rejection / resimulation path of `with_maf_filter`
  carries the extra cost (batching already done, x26.8 over the original).
- **PoolSeq: msprime is already faster** (0.45-0.79x). DIYABC's slowness on the
  serial case is plausibly its discrete coalescent (see the `evalcriterium`
  entry), **unverified**. Nothing to do.
- **MicroSat: `ms_modified` replay 228.7 s → 41.0 s (x5.6); gap to DIYABC 100x → 18x
  (12/10).** One run, variance not measured. KS unchanged against the real reftable;
  **the sign test on `rdiff_mean`, run 04/11 on `ms_modified`: median −0.713 % (106/150 negative, p = 4.4e-7) and −0.723 % (114/150, p = 1.2e-10), KS 0/161 in both scenarios (one reftable): the known loci-count residual at 50 loci, same sign and order as `toy_example1_ms_seriel`.** Each step is
  output-neutral at equal seed (`compare_particles_microsat`: 0 differences on 30
  particles x 152 values). One particle in-process without a profiler (seed 1, drawing
  path): 1.364 → 0.255 s. DNA (60-64x) is untouched: it has its own matrices and no AML.
  - **AML**: `_extract_mutual_data` computes `(length_by_sample, genotypes_by_sample)`
    once per group, not once per triplet (12 x 50 extractions → 50);
    `_prepare_triplet_for_admixture` only filters; `_vectorize_triplet` turns a
    triplet's frequencies into arrays (haploid / homozygous / heterozygous) so
    `_log_likelihood_admixture` evaluates each `a` on whole arrays. A heterozygote with
    ONE zero frequency contributes the log of the other frequency (`pente_lik`,
    `sumstat.cpp:1252-1260`; corrected 23/10, the note that said "0 as a whole genotype"
    was wrong, see Closed investigations).
  - **SNI matrices**: `_gsm_matrix_numpy` replaces `msprime.TPM(p=ε, m, lo=0, hi)`,
    whose `general_microsat_rate_matrix` is a pure-Python O(n²) loop (152 100
    `gamma_microsat` calls per particle). Tested against msprime on 161 `(hi, m)` pairs,
    max difference 1.1e-15. `Pgeom` is drawn per locus (50 distinct values over 50
    loci), so memoizing across loci is useless. Wired in BOTH siblings (drawing and
    `..._from_values`), see the traps below.
  - **Shared reads**: `compute_all_statistics_microsat` reads `_length_by_sample` and
    `_genotypes_by_sample_and_individuals` once per locus and passes them to every
    statistic as `length_by_sample_per_locus` / `genotypes_by_sample_per_locus` (11
    reads → 1, 2 → 1). A dispatch dict calls all its statistics alike, so the ones that
    do not read them (`FST` among others) accept the keyword and ignore it; their
    docstrings say so.
  - **Reads by array**: `_individuals_by_sample`,
    `_length_by_sample_and_individuals` and `_genotypes_by_sample_and_individuals` read
    `ts.individuals_nodes` (padded with -1) instead of `ts.individual(k)` per
    individual (16 150 calls per particle). `_individuals_by_sample` raises
    `ValueError` for an individual in no sample or in several. Kept on purpose: a
    monomorphic locus gives `(0, 0)` to EVERY individual, haploid ones included, in
    `_genotypes_by_sample_and_individuals`. Changed: elements are Python `int`, no
    longer `np.int64`. The old loops are the `_reference_*` oracles of
    `tests/test_summary_statistics.py`.
  - **Series 15-19/10: particle 133 → ~66 ms (x2.0), statistics 74 → ~30 ms,
    simulation ~60 → ~36 ms** (median in-process on `ms_modified`, `HEAD` vs modified
    measured alternately). Every step is **output-identical**
    (`compare_particles_microsat`, 30 particles x 152 values). Narrative, tables and the
    bugs caught: `notes/exploration.md`, sections 15/10, 16/10 and 19/10.

    | Change | Gain |
    |---|---|
    | `FST`: `_allele_counts_by_sample` (`ni`, `nA`, `AA` once per sample) | 19-26 → 9-10 ms |
    | `FST`: pairs derived from the shared genotypes (`_pairs_from_genotypes`), no second tree read | 8.4 → 3.3 ms |
    | `_length_by_sample`: list conversion out of the allele loop, then `np.bincount` | ~8 + 0.5 ms |
    | `_individuals_by_sample`: node → sample table instead of one `np.isin` per sample | 3-5 ms |
    | `AML`: empty genotype blocks skipped; frequency cache per parent; vectorization cache per (focal, parent), both created in `compute_AML_microsat` | 16.4 → ~9 ms |
    | `H2P`: dict of `pop_2` instead of a `next(...)` scan per allele | 2.7 → 1.4 ms |
    | `build_transition_matrix_microsat_with_sni_numpy`: dense matrix built **by blocks** (`dense[r::motif, r::motif] = local matrix of residue r`, tridiagonal SNI, mixing on the whole matrix) instead of 79 rows in Python | 33.5 → 8.1 ms (**particle -25 ms**) |
    | `_gsm_matrix_numpy`: one power per distinct distance instead of one per cell | -0.5 ms |
    | `LIK`: `nal` (`_count_alleles`) once per locus, `count_j`/`total_count_j` once per sample and locus | 8.2 → 5.9 ms |
    | `_decode_locus`: the variant of a locus decoded once, shared by `_length_by_sample` and `_genotypes_by_sample_and_individuals` | statistics -2.1 ms |

    **What each one must keep** (none is visible to a unit test on one locus):
    the order of `unique_alleles` in `FST` (built from dict keys in first-appearance
    order: it fixes the order of the float sums); `_pool_allele_counts_for_two_samples`
    iterates `pop_1` and ignores alleles only in `pop_2`; `_pairs_from_genotypes`
    gives `(x, x)` to a haploid and relies on the monomorphic locus being `(0, 0)` for
    everyone in **both** representations; `_individuals_by_sample` now raises on any
    node shared by two samples (stricter than the old count); `bincount` equals a
    count by size only because allele sizes of a locus are distinct (checked: 0
    duplicates over 1 499 loci), with `minlength=len(alleles)` or zero-count alleles
    vanish silently. The `AML` caches live **one call** of `compute_AML_microsat`
    (never module-level) and share arrays that nothing downstream mutates.
    **Not done on purpose**: evaluating `a` and `a + 0.001` together in
    `_log_likelihood_admixture` (changes the float summation order, and the bisection
    branches on the sign of a difference of close values).

    **Matrices, `LIK`, decode (19/10) — what they must keep.** The block construction
    is valid for `n_alleles >= 2 * motif_size` (below it the *old* function already
    raises: a residue with one state gives a 1x1 matrix and a division by zero); the
    SNI matrix is `0.5` on both off-diagonals plus the two corners `[0, 0]` and
    `[-1, -1]` (the folded border mass), and the mixing `gsm*(1-p) + sni*p` is
    elementwise, hence bit-identical. In `_gsm_matrix_numpy` the vector of powers must
    be built **by numpy** (`(1 - m) ** np.arange(...)`): a Python list comprehension
    rounds differently, 15 of 45 `(hi, m)` pairs. `LIK`: the `nal` rule is unchanged (an
    allele whose total over all samples is 0 is **not** counted: `test_count_alleles`
    pins it); each optional argument of `_compute_LIK_for_one_locus` has a `None`
    fallback. `_decode_locus` returns `(sizes per allele, allele index per node)` — the
    **indices** (`variant.genotypes`, read-only, stable after the decode loop: checked),
    not sizes per node; `_genotypes_by_sample_and_individuals` rebuilds
    `np.array(sizes)[indices]` itself; `None` for a monomorphic locus.
    **Float order facts measured on this machine** (for any future vectorization of
    `LIK`): `np.log10` of a scalar and of an array are bit-identical (36 000 real
    values), `math.log10` is not; `np.sum` differs from the sequential subtraction of
    the code in 51 % of 20-term segments, `np.cumsum(...)[-1]` in 0 %.
  - Replay path before the numpy matrix reached it (3 particles, `time_stages_microsat`):
    matrices 61 %, `sim_ancestry` 19 %, `sim_mutations` 4 %, all statistics together
    13 %. `ms_seriel` was not re-measured.
- **DNA: 104.1 s → 13.4 s (x7.8), one particle 426 → 90 ms in-process (14/10).** Every
  step output-identical at equal seed except the `<A>` batching (Closed investigations).
  Floor check: 1000 particles x 90 ms / 8 physical cores = 11.3 s, measured 13.4 s
  (x1.19: pool start-up, hyperthreading); **the benchmark scales with physical
  cores**, not 16. Steps, with the gain measured *in situ* (single particle, then
  `perf_counter`): provenance off (3-4 %, `cProfile` said 10 %); genotype matrix read
  once per locus (1290 → 100 calls); `_minor_allele_counts_at_segregating_sites`
  without the per-site `np.unique` (65 → 13 ms under `cProfile`); `_triu_indices`
  cached; `compute_sample_layout` from `ts.individuals_nodes` (111 → 3 µs per call);
  `rescale_demography` once per heritage type (1 copy instead of 50: **6 ms real**,
  56 ms under `cProfile`); shared `<M>` tree and batched `<A>` replicates
  (`sim_ancestry` 100 → 2 calls; one particle 0.294 → 0.123 s); `_count_distinct_haplotypes` as a
  base-K integer code (54 → 9 µs per call); `_pairwise_hamming_distances` as one
  indicator matrix product per state (22 → 14 µs) and computed once per
  (locus, sample) and shared (926 → 200 calls per particle, -11 ms).
  - `_count_distinct_haplotypes`: columns coded in base K = max state + 1 after a
    shift by the minimum (in `int64`: `int8` extremes overflow otherwise); **the
    overflow guard `K ** n_sites > 2**62` must be evaluated with Python ints**, numpy
    ints wrap silently; above it, falls back to `np.unique(axis=1)`. Zero sites → 1.
  - `_pairwise_hamming_distances`: K ≤ 8 states → sum over states of
    `(shifted == k).T @ (shifted == k)` in `float32` (exact), distance = sites − agreements,
    upper triangle by cached `_triu_indices`; K > 8 → the old broadcast. **Output is
    `int64`** (a float output flips every downstream type). Zero sites → zeros of length
    `n(n−1)/2`, not an empty vector.
  - **What remains (90 ms per particle):** `sim_mutations` 36 % (100 calls x 0.3 ms: ~0.18 ms
    fixed msprime cost plus ~0.115 ms of table conversions Python↔C (measured 02/11; the `RateMap` itself costs 0.6 µs to convert); cannot be
    batched, each locus has its own matrix, rate map and seed), statistics 34 % (13
    statistics, none above ~3 ms, plus ~9 ms of shared preparation), site rates 13 %
    (done 02/11 with `Generator.gamma`: 5.8 → 0.8 ms, the stream changed and was revalidated by replay, see
    the Closed investigations entry), `sim_ancestry` 12 %, matrices 2 %. **Verdict (02/11):
    statistics at the numpy call-overhead floor, see `notes/exploration.md` §11; DNA is done unless the architectural
    batching of loci is taken on.** The remaining ratio to DIYABC (~7-8x) is the one-library-call-per-locus
    model, not a local defect.

**Open leads, one session each**: (1) MicroSat, **timing work paused on 19/10**. A
particle is ~66 ms: statistics ~30 ms (`AML` 9.5, `LIK` 5.9, shared reads 3.8, `FST` 3.5,
`H2P` 1.5, eight small ones ~4.3, decode and rest 1.7), simulation ~36 ms
(`_simulate_ancestries_for_loci` 13.5, `sim_mutations` 12.7, matrices 7.5). **Remaining
overhead, by nature**: ~9 ms is fixed `msprime` cost per call (`sim_mutations` 130 us per
call whether the model has 2 or 79 alleles, x 50; `sim_ancestry` 1.45 ms per call to
build, x 2) — not reducible while loci are simulated one by one; ~17 ms is `msprime`'s
real work (by subtraction, not timed directly; the 65 us per tree measured is for 20
haploid samples in one population, the pipeline has 4 populations); **~38 ms (57 %) is
our own small-array numpy calls and Python loops** (0.4-0.9 us per numpy call on 20
elements; `AML`'s likelihood is 574 evaluations x ~12 us, i.e. call overhead, not
arithmetic). The only way left is **batching**, ~10-15 ms at best, each with an
arithmetic-order risk: the 24 `AML` triplets bisected in lock-step, `LIK` per locus
(**tried 20/10 and dropped**: an exact array version of `compute_LIK` (0 differences on 300 random
loci and 30 particles) saved only 0.8 ms in situ, 1.2 % of a particle, for three more functions, a
fallback branch and a dozen tests; rule kept: below ~1 ms per particle a new function is not worth it), the matrices of several loci in
one 3-D array. Small exact ones: ~1 ms each (the 8 small statistics, the per-individual
tuples of the shared reads). **Measured 22/10: 9.7 s against DIYABC 2.1 s, 4.59x** (the 8-10 s, ~4x estimate held). Per-statistic timer:
`scripts/bench_time_stats.py`; the simulation split (patching the outer functions
`build_matrix_microsat_per_locus`, `_simulate_ancestries_for_loci`,
`msprime.sim_mutations`) was not kept as a script. (2)-(3) IndSeq (MAF > 0, Hudson, serial): **closed on 21/10**, see the IndSeq verdict above. (4) PoolSeq: none. (5) DNA: see the verdict above.
**Every optimization touches code validated against DIYABC**: an output-neutral one
must stay *identical* at equal seed (`compare_particles_dna` / `_microsat`, then the
golden tests); one that changes the random stream (the `<A>` batching) is validated by
a paired replay read with the KS, the sign test **and a re-run with `seed_offset`**. The
sibling rule applies unless the user decides otherwise — ask first.

**Measuring, traps already met**
- Profile **one particle in-process**, never the replay: under a
  `ProcessPoolExecutor`, `cProfile` sees only the main process. It also inflates
  call-heavy code ~1.7x (`ts.individual()`, `gamma_microsat`): use it to rank,
  confirm with `perf_counter`. `python -m cProfile -m` needs Python >= 3.12;
  the environment is 3.11.
- SNP simulation is a **lazy generator**: `simulate_particle_genotypes_from_values`
  returns before simulating anything. Time the call **and** the `list(...)`
  that consumes it. (The MicroSat path is eager.)
- Wrapping a statistic to time it must patch the **dispatch dicts**
  (`_MICROSAT_*_STATS`), which hold their own references; patching the module
  attribute changes nothing.
- Running DIYABC for a timing: copy the directory first (never run in
  `reference/`), pass `cwd=` the copy (`-p ./` is relative to the **process**),
  initialise the RNG with `-n "t:16"` (`;`-separated `t:` / `c:` / `s:` / `f:`
  options, not a bare number), always `-R ALL`, and check the log ends in
  `END|<nrec>|`.
- **A drawing/replay pair must be optimized AND checked on both halves.** The numpy
  matrix was wired into `build_matrix_microsat_per_locus` only;
  `build_matrix_microsat_per_locus_from_values` (the replay path) kept `msprime.TPM`
  until 12/10, while every profile and identity check ran on the drawing path. Symptom:
  a single-particle gain far larger than the replay gain (x4 against x1.4). Time the
  replay path itself (`time_stages_microsat`).
- **A useless computation is invisible to the tests and to the identity check.** `N2P`
  briefly recomputed genotypes it never read (~9 % of a particle under `cProfile`,
  identical outputs).
  After sharing a computation, count the calls with `pstats.Stats.print_callers`:
  `_length_by_sample` and `_genotypes_by_sample_and_individuals` must each run once
  per locus.
- Identity check of an optimization: `compare_particles_microsat save ref.json` on the
  code BEFORE, `compare ref.json` after. Exact equality, no tolerance. **Take the
  reference on `HEAD` in a clean worktree** (`git worktree add --detach <tmp> HEAD`,
  symlink `reference` into it, `save` there, `compare` from the main tree) rather than
  with `git stash`, which also stashes unrelated work; **a reference older than a
  statistical change is void** (`ref_stats_dna.json` stopped being valid at the `<A>`
  batching).
- **`cProfile` ranks, `perf_counter` decides, and the synthetic micro-benchmark lies
  too.** Small numpy calls are inflated 2-5x: provenance (10 % predicted, 3-4 % real),
  `rescale_demography` (12 % predicted, 2 % real), `_pairwise_hamming_distances`
  (66 µs predicted by a miscounted loop, 22 µs real). A function benchmarked on a
  synthetic matrix of 11 sites x 20 columns (9 µs) cost 22 µs on the real matrices,
  which include 40- and 80-gene merges. Capture the real arguments with a spy (keep
  the original function in a variable *before* patching, or the spy times itself) and
  time the function on them; count calls with `print_callers`.
- **Single-particle timings carry ±7 % noise** (the same code gave 0.410 and 0.382 s):
  run 10-20 particles and take the minimum or the median; a gain under ~5 % is not
  detectable on one.
- **A timer wrapper is not free, and nested wrappers under one label double-count**
  (the tell is a part longer than its whole). Wrap only outer functions, and measure
  `HEAD` (clean worktree) and the modified tree **alternately**, twice, 20 particles: the
  machine drifts by ~1 ms between passes, a figure from another session is not a
  baseline. When `reference/` is tracked, `ln -s` into the worktree creates the link
  *inside* it: remove the checked-out directory first. `cProfile` misranks small-call
  code by x4-5 (`H2P` 13 ms against 2.8 real): rank with it, decide with `perf_counter`.
- **Tests of an output-neutral optimization** (all four were missed at least once in the
  15-16/10 series, each time found by `compare_particles_microsat`, not by the tests):
  (1) the oracle is the **old code copied into the tests**, and the test imports the
  function from `bridge/` — a test redefining the function locally checks its own copy;
  (2) **at least two different loci**, or a loop whose `append` fell out of the `for`
  keeps the last locus and passes; (3) a cache needs a **call-count spy**
  (`monkeypatch` on the module attribute: 4 parents, 12 pairs), because a cache that
  computes before it reads gives the right answer and saves nothing; (4) run the
  **mutation check** (break the thing, see the test go red) on a *copy* of the module,
  never on the working tree. Compare with `==`, not `approx`.
- **Golden tests pin the seeded simulation, not the formula.** They stayed green
  through all the output-neutral refactors (their use) and moved with the first statistical
  change; regenerating them proves nothing about correctness (the replay does). A
  formula test should be built on a hand-made `TreeSequence`.

- **A micro-benchmark on objects you already touched lies.** `ts.tables` is cached per
  `TreeSequence`: timing it on trees warmed by earlier repetitions gave 12 µs for a function that
  costs 35 µs on the fresh trees of the real loop. Build fresh objects for every repetition.
- **Before timing or comparing with DIYABC, check the machine and the reference.** `ps -eo pcpu,args
  --sort=-pcpu | head` (a Jupyter kernel at 130 % CPU contaminated a scaling test and the
  `load average` was 5.7 before it started); `maf.txt` and `THIS->MAF=` (Hard rules); a stored
  `ref_*.json` is void once its dataset directory changed (one recorded at 100 loci was compared to
  the 500-loci directory: 1496 differences out of 1530, unrelated to the code) or if it was
  rewritten by a mistaken rerun.

**Where the tooling lives.** Tracked in this repository:
`scripts/compare_particles_microsat.py` and `scripts/compare_particles_dna.py`
(`save` / `compare` a JSON of N particles at equal seed: the identity check, exact
equality, drawing path only; **the replay paths are covered by no identity script**),
`scripts/profile_particle_microsat.py` and `scripts/profile_particle_dna.py` (one
particle, `cProfile` plus a wall time without profiler) and
`scripts/time_aml_microsat.py`. The **stage timer** `scripts/time_stages_microsat.py`,
the benchmark scripts (`benchmark_1000_msprime.py`, `bench_1000_diyabc.py`),
`scripts/profile_particle.py` and the `timings` plumbing of `pipeline.py` /
`reftable_loop.py` live on the branch **`bench/timings`**, checked out as the git
worktree `../diyabc-msprime-poc_timings` (its `reference` is a symlink to this
repository's). **Never merge `bench/timings` into `main`**: the plumbing makes
`replay_reftable_simulation_*` return a 5-tuple, which breaks
`scripts/replay_diyabc_priors_*.py` (`AttributeError: 'list' object has no attribute
'group_priors_values'`). Bring `main` into it with `git merge main`; run replays from
`main`. Don't recopy `bridge/` between directories (it overwrites uncommitted work):
the worktree exists for that.


# 4. Open work

- **`split` with k ancestral populations (06/11, Dirichlet 12/11) — implemented end to end, k = 3 run through `reftable_loop` on human.**
  Syntax `t split D A1 … Ak r1 … r(k-1)` (`len(args) = 2k`, k ≥ 2; the 4-token form is k = 2). `SplitEvent.ancestral_pops: list[int]`,
  `admixture_rates: list[str]`. Rates are absolute probabilities, the k-th ancestor gets `max(0, 1 - Σr)`; a negative rate or
  `Σr > 1 + 1e-9` raises `ValueError` (which `parse_header_scenarios` does **not** swallow, unlike `NotImplementedError`).
  Joint prior of the rates: `ra,rb A DI[min,max,0.0,0.0]`, a symmetric Dirichlet(1,…,1) (`DirichletPrior`, sibling of `Prior`);
  `min`/`max` bound **every** share including the complement, mean/sdshape ignored. Drawn by k normalised Gamma(1,1) with rejection,
  `max_attempts = 1000` (`ValueError` if exceeded). `check_header_with_dirichlet_prior(priors, scenarios)` (called by the 4 `run_*` and
  4 `replay_*` entry points) requires every `split` that names a Dirichlet rate to name exactly that prior's names. Identity check at
  equal seed on a header without `DI`: nothing changed. Narrative: `notes/exploration.md`, 2026-11-06 and 2026-11-12.
  Still open: acceptance-probability check at parse time, non-symmetric α, global `min > max` test coherence.

- **DIYABC's sex-ratio mode (22/10) — we follow the real ratio; `toy_example5` is validated in that mode.**
  DIYABC forces `this->sexratio = 0.5` when it loads an existing `.snp.bin` (`data.cpp:1636-1642`) and keeps the
  real `<NM=...>` only on the run that converts the `.snp` (1644-1664); coefficients at 0.3: first run `A` 3.36,
  `X` 2.91, `Y` 0.6, `M` 1.4; later runs `A` 4, `X` 3, `Y` 1, `M` 1. `<A>` is rescaled by
  `coalescence_coefficient("A", s) / 4` (ploidy 2, skipped when `s == 0.5`) in `simulate_genotypes_for_locus_type`
  and `simulate_poolseq_reads_with_mrc_filter`. The 0.5 mode shows up against us as `<M>` −7 %, `<Y>` +9 % (and
  `<A>` exact); a reference to validate a dataset with `NM != 1` must be regenerated with no
  `.snp.bin`/`.snpbin.txt` in the folder (Hard rules). **Done 22/10**: `reference/toy_example5` and
  `reference/toy_example5_with_maf` regenerated that way (`.snp.bin`, `maf.txt` and `reftableRF.bin` carry the same
  minute, i.e. converted in the same run); replay of the 1000 particles, 50 statistics, priors paired, rounded to
  6 decimals: `toy_example5` 0 KS flag, sign 21/50, median +0.09 %; `toy_example5_with_maf` 0 flag, sign 32/50
  (p = 0.065), median −0.18 %. The −2.65 % of `toy_example5` (70 loci) filed under "Systematic negative bias" is
  **not reproduced** on the regenerated reference; it is consistent with the old reftable being in the
  first-run mode while our `<A>` was not rescaled (ours would sit low), but the old reftable's mode is not
  recorded, so that is a hypothesis, not a verdict.
- **Serial sampling is not wired for `<X>`/`<Y>` — now guarded, not fixed.**
  The per-locus dispatch still overrides the serial `SampleSet`s, on both the
  SNP and MicroSat/DNA sides. Since 28/09 the combination is refused up front by
  `reftable_loop.raise_if_serial_with_sex_linked_loci`, called from all seven
  run/replay entry points; the predicate is
  `scenario_parser.is_serial_scenario`. Implementing it for real was
  deliberately deferred: **no dataset would validate it**, and this project
  never declares a path correct without a paired replay against the real
  DIYABC. See "Serial/temporal sampling" under Domain knowledge.
- **`ML3p_2.3.4` on `toy_example4_seriel` (05/10) — RESOLVED: an artefact of overlapping seeds, not a residual.** It was
  KS-flagged in three replays (p = 0.0001 to 0.0013, −10 to −14 %, only the columns carrying the oldest sample, t = 500)
  because the PoolSeq `with_mrc_filter` seeds were additive (`seed + batch_index * 20`, `seed + attempt + offset`): the
  `DTA_2_2` mechanism, see Closed investigations. Both sites now go through `_locus_seed` (`_MRC_BATCH_SEED_OFFSET`).
  Two independent DIYABC runs at 1000 loci, two scenarios each: `ML3p_2.3.4` KS p = 0.996 and 0.30, 0 / 6 / 0 / 2 flags
  of 120 (about 6 expected), no column in common, none carrying the oldest sample. The localisation, the `tbn` /
  `Npast` / `Npresent` stratifications and the "first steps if anyone picks this up" are the record of what
  overlapping seeds produced, **not evidence** of a coalescent or `cal_ml3p` defect (the discrete coalescent was
  falsified on 30/09). `ML2p_2.4` (intermediate, decaying with loci count) is the loci-count residual. Full text:
  `notes/exploration.md`, 25/09, 05/10 and the 03/11 archive.
- **G2 `MPD`/`VPD` residual (02/10) — a DIYABC bug, located to the line, reproduced, validated by replay 05/10.**
  `cal_nss1p` takes its per-locus denominator from the OBSERVED data (`samplesize()`), `cal_mpd1p` (`nd > 0`) and
  `cal_vpd1p` (`nd > 1`) from the SIMULATED sequences. On the reftable path (`dnatrue = false`, `init_dnaseq`) a locus with
  no mutation keeps an empty sequence, which is `SEQMISSING` (`particuleC.hpp:17`): it leaves the `MPD`/`VPD`
  denominator and stays in `NSS`'s. As `θ → 0`, `MPD` plateaus at `1/a_n` instead of going to 0 (6.2 % of simulated loci
  carry no mutation). Scope: `MPD` and `VPD` only (`NSS`, `PSS`, `NHA`, `MNS`, `MP2`, `MPB`, `HST` flat at ~1 by `θ`
  quintile; `DTA` counts monomorphic loci with D = 0). Reproduced by `if ts.num_sites == 0: continue` in
  `compute_MPD`/`compute_VPD` (see "Deliberately reproduced C++ bugs"; guarded by
  `test_mean_pairwise_differences_excludes_loci_without_mutation`). **A verification on observed data does not validate
  the simulated path** (`statobsRF` is computed with `dnatrue = true`). **Practical rule: beyond `θ ≈ 1` the ratio is 1.013
  then 0.995, so the defect bites only on low-diversity loci; there DIYABC's `MPD`/`VPD` are not a trustworthy
  reference.** Falsified, do not retry: the discrete coalescent, `put_mutations`, `coal_pop`, `cal_numvar`,
  `sample_site_rates`, a shifted frequency spectrum. To dump DIYABC's simulated sequences: `general -p ./ -r 2 -R ALL -a 9`
  (`debuglevel == 9`, `kloc == 10` in `cal_nss1p`; **`-R ALL` is load-bearing**, without it the header read hangs at
  100 % CPU). Two latent C++ bugs found in passing, not reproduced and not triggered here (`mute`'s reused loop
  counter, `init_dnaseq`'s unbounded inverse-CDF loop). Full text and tables: `notes/exploration.md`, 02/10 and the
  03/11 archive.
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

- **MicroSat G1 deficit on `toy_example2_ms_dna_50loci_TN_seriel` (04/11) — open.** Four cells (reftables 0 and 1,
  both scenarios): over the 76 MicroSat G1 statistics (10 loci) `rdiff` median −1.21 / −1.91 / −1.34 / −0.75 %,
  54-63 of 76 negative, `LIK`/`MGW`/`N2P` −1 to −4 % in every cell, 9 of 12 KS flags in that group; the DNA groups
  stay within −0.70 / +0.47 %. Not scenario-2-specific. Compatible with the loci-count residual (−0.4 to −1.3 % at 50
  loci on `ms_seriel`, −0.72 % on `ms_modified`) but the decrease is not 1/loci (x5 fewer loci, x1.5 to 2.5 more bias),
  so **not demonstrated**. Caveat: the reftable-0 CSVs (23/10) predate the 02/11 regeneration of that reftable. A
  discriminating test needs a dataset with more MicroSat loci in the same group.
- **PoolSeq at sex ratio 0.2308, scenario 2 (04/11) — open.** Positive shift on two replays (+0.64 %, +0.99 %, 3 KS
  flags on the differentiation statistics `FST1m`, `FST1v`, `HBm` on the first), scenario 1 clean. No baseline: run the
  same 100-loci replay at sex ratio 0.5 (`toy_example4`) and a `seed_offset` replay before attributing it to the
  rescaling (the factor itself is guarded by the two pytest spies). Dataset: `toy_example4_sex_ratio_diff`.
- **No identity or smoke test covers the replay (`_from_values`) paths — open.** `compare_particles_*` checks the drawing
  path only. The replay broke unnoticed twice (serial IndSeq for ten days, the `µ`/`mu` readers on 04/11, both with green
  tests). Idea: a test replaying 2-3 particles per data family (IndSeq, PoolSeq, MicroSat, DNA, mixed) from the first rows
  of a real reftable, asserting that it runs and returns exactly the declared statistic names, on at least one `mu` and
  one `µ` header.
- **Notebooks to rerun (04/11).** All 23 `reference/*/compare_reftables_*.ipynb` (stored outputs and summary CSVs stale);
  the three without a replay file need one first; the DNA-only `ms_dna_50loci_{JK,K2P}` / `ms_dna_{JK,K2P,XY}` replays
  need the mixed script.

