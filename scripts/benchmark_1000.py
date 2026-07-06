"""
Pour mesurer le temps de calcul de la simulation de 1000 particules de 5000 loci chacune, et l'écriture du fichier reftable.bin.
"""

import time

from bridge.prior_parser import parse_priors
from bridge.reftable_loop import (
    run_reftable_simulation,
    write_reftable_bin,
    write_reftable_txt,
)
from bridge.scenario_parser import parse_header_scenarios

with open("reference/Exemple5/headerRF.txt") as f:
    header_text = f.read()
priors, _ = parse_priors(header_text)
scenarios = parse_header_scenarios(header_text)
scenario1 = next(s for s in scenarios if s.index == 1)

t0 = time.time()
results = run_reftable_simulation(
    reference_directory="reference/Exemple5",
    scenarios=[scenario1],
    num_loci=100,
    nrec=1000,
    stats_filter="ALL",
    max_workers=16,
)
t1 = time.time()

write_reftable_bin(results, priors, [scenario1], "./tmp/bench_reftable.bin")
t2 = time.time()

write_reftable_txt(results, priors, [scenario1], "./tmp/bench_reftable.txt")
t3 = time.time()


print(f"Simulation 1000 particules : {t1 - t0:.1f}s ({(t1 - t0) / 60:.1f} min)")
print(f"Écriture reftable.bin      : {t2 - t1:.1f}s")
print(f"Écriture reftable.txt      : {t3 - t2:.1f}s")
print(f"TOTAL                      : {t3 - t0:.1f}s ({(t3 - t0) / 60:.1f} min)")
