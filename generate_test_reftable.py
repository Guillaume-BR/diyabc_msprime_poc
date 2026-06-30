import os
from bridge.reftable_loop import run_reftable_simulation, write_reftable_bin
from bridge.prior_parser import parse_priors
from bridge.scenario_parser import parse_header_scenarios

GENERAL_BINARY_PATH = os.environ["DIYABC_GENERAL_PATH"]

header_text = open("reference/human/header.txt").read()

results = run_reftable_simulation(
    reference_directory="reference/human",
    scenario_index=1,
    num_loci=10,
    nrec=5,
    general_binary_path=GENERAL_BINARY_PATH,
    base_work_directory="./tmp/reftable_particles",
    stats_filter="ALL",
)

priors, _ = parse_priors(header_text)
scenarios = parse_header_scenarios(header_text)
scenario1 = next(s for s in scenarios if s.index == 1)

write_reftable_bin(results, priors, scenario1, "./tmp/notre_reftable.bin")
print("reftable.bin écrit.")
