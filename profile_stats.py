import cProfile, pstats, io
from bridge.pipeline import run_poc_for_directory
from bridge.ancestry_simulation import build_samples_argument
from bridge.summary_statistics import (
    compute_ML1,
    compute_ML2,
    compute_ML3,
    compute_HW_HB,
    compute_FST1,
    compute_FST2,
    compute_NEI,
    compute_AML,
    compute_F3_F4,
    compute_FST3_FST4_FSTG,
)

# Simuler une fois (5000 loci), hors profilage
genotypes_per_locus, _ = run_poc_for_directory("reference/human", 1, 5000, 1)
genotypes_list = list(genotypes_per_locus)
pop_names = list(
    build_samples_argument("reference/human/human_snp_all22chr_maf5.snp").keys()
)

funcs = {
    "ML1": compute_ML1,
    "ML2": compute_ML2,
    "ML3": compute_ML3,
    "HW_HB": compute_HW_HB,
    "FST1": compute_FST1,
    "FST2": compute_FST2,
    "NEI": compute_NEI,
    "AML": compute_AML,
    "F3_F4": compute_F3_F4,
    "FST3_FST4": compute_FST3_FST4_FSTG,
}

import time

print(f"{'Fonction':<12} {'Temps (s)':>10}")
print("-" * 24)
total = 0
for name, f in funcs.items():
    t0 = time.time()
    f(genotypes_list, pop_names)
    dt = time.time() - t0
    total += dt
    print(f"{name:<12} {dt:>10.3f}")
print("-" * 24)
print(f"{'TOTAL':<12} {total:>10.3f}")
