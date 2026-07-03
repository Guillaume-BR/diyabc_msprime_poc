from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.scenario_parser import parse_header_scenarios

with open("reference/human/header.txt") as f:
    header_text = f.read()
priors, _ = parse_priors(header_text)

scenarios = parse_header_scenarios(header_text)
scenario1 = next(s for s in scenarios if s.index == 1)

used_names = get_parameter_names_used_by_scenario(scenario1)
kept = [p.name for p in priors if not is_constant_prior(p) and p.name in used_names]

print("Nombre de paramètres gardés (filtrage combiné):", len(kept))
print(kept)
