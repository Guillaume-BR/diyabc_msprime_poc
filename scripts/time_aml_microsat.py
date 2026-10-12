"""
Mesure, sur UNE particule MicroSat en processus unique, comment se répartit
le temps d'AML entre la préparation des loci (`_prepare_loci_for_admixture`)
et la log-vraisemblance (`_log_likelihood_admixture`).

Question posée : le cache de `(length_by_sample, genotypes_by_sample)` suffit-il
(gain borné par la part de la préparation), ou faut-il vectoriser la
vraisemblance ?

Usage, depuis la racine du dépôt :
    python3 -m scripts.time_aml_microsat [dossier_reference] [nb_particules]

Pièges évités (voir CLAUDE.md, « Measuring, traps already met ») :
- une seule particule à la fois, dans CE processus (cProfile ne voit pas
  les workers d'un ProcessPoolExecutor) ;
- on enveloppe avec `perf_counter`, pas avec cProfile (qui gonfle ~1.7x le
  code riche en petits appels) ;
- les fonctions chronométrées sont appelées par leur NOM depuis
  `summary_statistics`, donc patcher l'attribut du module suffit ;
  le dict de dispatch `_MICROSAT_TRIPLET_STATS`, lui, garde sa propre
  référence : il est patché explicitement.
"""

import statistics
import sys
from collections import defaultdict
from pathlib import Path
from time import perf_counter

from bridge import pipeline
from bridge import summary_statistics as ss
from bridge.header_dataclasses import MicrosatReplayContext
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    allele_bounds_per_locus,
    observed_count_sample,
    observed_microsatellites,
    parse_sex_ratio,
)

REFERENCE_DIR = Path(
    sys.argv[1] if len(sys.argv) > 1 else "reference/toy_example1_ms_modified"
)
N_PARTICLES = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SCENARIO_INDEX = 1

elapsed: dict[str, float] = defaultdict(float)
calls: dict[str, int] = defaultdict(int)


def timed(label, fn):
    """Enveloppe `fn` : cumule durée et nombre d'appels sous `label`."""

    def wrapper(*args, **kwargs):
        t0 = perf_counter()
        try:
            return fn(*args, **kwargs)
        finally:
            elapsed[label] += perf_counter() - t0
            calls[label] += 1

    return wrapper


# Fonctions appelées par leur nom dans summary_statistics : patch du module.
ss._prepare_triplet_for_admixture = timed(
    "aml_prepare", ss._prepare_triplet_for_admixture
)
ss._log_likelihood_admixture = timed("aml_loglik", ss._log_likelihood_admixture)
# Dict de dispatch : sa propre référence, à remplacer à la main.
ss._MICROSAT_TRIPLET_STATS["AML"] = timed("aml_total", ss.compute_AML_microsat)
# Étape « statistiques » vue depuis pipeline (nom importé dans son espace).
pipeline.compute_all_statistics_microsat = timed(
    "stats_total", pipeline.compute_all_statistics_microsat
)

header_text = pipeline.read_header_text(REFERENCE_DIR)
mss_path = REFERENCE_DIR / header_text.splitlines()[0].strip()
list_loci = parse_loci_description(header_text)
microsat_observed = observed_microsatellites(mss_path, list_loci)
context = MicrosatReplayContext(
    header_text=header_text,
    mss_path=mss_path,
    list_loci=list_loci,
    microsat_observed=microsat_observed,
    bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
    samples_default=observed_count_sample(mss_path),
    sex_ratio=parse_sex_ratio(mss_path),
)

# Une particule d'échauffement (imports paresseux, caches), non comptée.
pipeline.compute_summary_statistics_microsat(context, SCENARIO_INDEX, seed=1000)
elapsed.clear()
calls.clear()

per_particle = []
for seed in range(1, N_PARTICLES + 1):
    before = dict(elapsed)
    t0 = perf_counter()
    pipeline.compute_summary_statistics_microsat(context, SCENARIO_INDEX, seed=seed)
    wall = perf_counter() - t0
    delta = {k: elapsed[k] - before.get(k, 0.0) for k in elapsed}
    delta["wall"] = wall
    per_particle.append(delta)

print(f"\n{REFERENCE_DIR} -- {N_PARTICLES} particule(s), 1 processus\n")
print(
    f"{'étape':<18}{'médiane (s)':>9}{'% du wall':>11}{'% AML':>10}{'appels/part.':>14}"
)
med = {k: statistics.median(p[k] for p in per_particle) for k in per_particle[0]}
rows = [
    ("wall", "particule entière"),
    ("stats_total", "toutes les stats"),
    ("aml_total", "AML"),
    ("aml_prepare", "  AML: préparation"),
    ("aml_loglik", "  AML: loglik"),
]
for key, label in rows:
    pct_wall = 100 * med[key] / med["wall"]
    pct_aml = (
        100 * med[key] / med["aml_total"] if key.startswith("aml") else float("nan")
    )
    n_calls = calls[key] / N_PARTICLES if key in calls else float("nan")
    print(
        f"{label:<18}{med[key]:>9.3f}{pct_wall:>11.1f}{pct_aml:>10.1f}{n_calls:>14.0f}"
    )

other = med["aml_total"] - med["aml_prepare"] - med["aml_loglik"]
print(f"\nAML hors préparation/loglik (dichotomie, boucle de triplets) : {other:.3f} s")
print(f"Simulation (wall - stats) : {med['wall'] - med['stats_total']:.3f} s")
print(
    "\nLecture : le gain maximal d'un cache est borné par la ligne "
    "'AML: préparation'.\n"
    "Imbrication : aml_prepare et aml_loglik sont INCLUS dans aml_total, "
    "qui est inclus dans stats_total."
)
