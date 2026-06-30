from bridge.prior_parser import parse_priors, is_constant_prior
from bridge.scenario_parser import parse_header_scenarios
from bridge.demography_builder import get_parameter_names_used_by_scenario

header_text = open("reference/human/header.txt").read()
priors, _ = parse_priors(header_text)

scenarios = parse_header_scenarios(header_text)
scenario1 = next(s for s in scenarios if s.index == 1)

used_names = get_parameter_names_used_by_scenario(scenario1)
kept = [p.name for p in priors if not is_constant_prior(p) and p.name in used_names]

print("Nombre de paramètres gardés (filtrage combiné):", len(kept))
print(kept)
