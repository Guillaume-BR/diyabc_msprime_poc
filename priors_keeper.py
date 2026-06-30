from bridge.prior_parser import parse_priors, is_constant_prior

header_text = open("reference/human/header.txt").read()
priors, _ = parse_priors(header_text)

kept = [p.name for p in priors if not is_constant_prior(p)]
print("Nombre de priors gardés par notre code:", len(kept))
print(kept)
