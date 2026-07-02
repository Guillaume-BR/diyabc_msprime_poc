import os, time
from bridge.reftable_loop import run_reftable_simulation
from bridge.reftable_loop import write_reftable_bin
from bridge.prior_parser import parse_priors
from bridge.scenario_parser import parse_header_scenarios
from pathlib import Path

header_text = open("reference/human/header.txt").read()
priors, _ = parse_priors(header_text)
scenarios = parse_header_scenarios(header_text)
scenario1 = next(s for s in scenarios if s.index == 1)

Path("./tmp/bench_particles").mkdir(parents = True, exist_ok=True)

t0 = time.time()
results = run_reftable_simulation(
    reference_directory="reference/human",
    scenario_index=1,
    num_loci=5000,
    nrec=1000,
    general_binary_path=None,
    base_work_directory="./tmp/bench_particles",
    stats_filter="ALL",
    max_workers=8,
)
t1 = time.time()

write_reftable_bin(results, priors, scenario1, "./tmp/bench_reftable.bin")
t2 = time.time()

print(f"Simulation 1000 particules : {t1-t0:.1f}s ({(t1-t0)/60:.1f} min)")
print(f"Écriture reftable.bin      : {t2-t1:.1f}s")
print(f"TOTAL                      : {t2-t0:.1f}s ({(t2-t0)/60:.1f} min)")
