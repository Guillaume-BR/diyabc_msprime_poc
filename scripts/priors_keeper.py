from bridge.prior_parser import is_constant_prior, parse_priors

with open("reference/human/header.txt") as f:
    header_text = f.read()
priors, _ = parse_priors(header_text)

kept = [p.name for p in priors if not is_constant_prior(p)]
print("Nombre de priors gardés par notre code:", len(kept))
print(kept)
