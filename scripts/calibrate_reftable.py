import time

from bridge.reftable_loop import run_reftable_simulation
from bridge.scenario_parser import parse_header_scenarios

with open("reference/human/header.txt") as f:
    scenarios = parse_header_scenarios(f.read())

start = time.time()
results = run_reftable_simulation(
    reference_directory="reference/human",
    scenarios=scenarios,
    num_loci=5000,  # le vrai nombre de loci de human, pas 10 -- important pour une estimation réaliste
    nrec=10,
    stats_filter="ALL",
)
elapsed = time.time() - start

print(f"Temps pour 10 particules (5000 loci chacune) : {elapsed:.1f} secondes")
print(f"Temps moyen par particule : {elapsed / 10:.1f} secondes")
print(f"Estimation pour 1000 particules : {elapsed / 10 * 1000 / 60:.1f} minutes")
