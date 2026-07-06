from bridge.prior_parser import parse_priors
from bridge.reftable_loop import run_reftable_simulation, write_reftable_bin
from bridge.scenario_parser import parse_header_scenarios

with open("reference/human/header.txt") as f:
    header_text = f.read()

priors, _ = parse_priors(header_text)
scenarios = parse_header_scenarios(header_text)

results = run_reftable_simulation(
    reference_directory="reference/human",
    scenarios=scenarios,
    num_loci=10,
    nrec=5,
    stats_filter="ALL",
)

write_reftable_bin(results, priors, scenarios, "./tmp/notre_reftable.bin")
print("reftable.bin écrit.")
