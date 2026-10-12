"""
Profil d'UNE particule MicroSat en processus unique : où va le temps, une fois
l'AML allégé (cache de `_extract_mutual_data`) ?

Deux lectures complémentaires, car chacune a un défaut connu :
- cProfile CLASSE les fonctions (qui coûte le plus) mais gonfle ~1.7x le code
  riche en petits appels (`ts.individual()`, `gamma_microsat`) ;
- `perf_counter` donne les durées réelles des grandes étapes, mais seulement
  celles qu'on a enveloppées.
Lire les parts du premier, confirmer avec le second.

Usage, depuis la racine du dépôt :
    python3 -m scripts.profile_particle_microsat [dossier_reference] [seed]

Pièges (voir CLAUDE.md, « Measuring, traps already met ») : une seule
particule, dans CE processus ; les statistiques sont dans des dicts de
dispatch qui gardent leurs propres références, d'où le chronométrage par
cProfile plutôt que par patch.
"""

import cProfile
import pstats
import sys
from pathlib import Path
from time import perf_counter

from bridge import pipeline
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
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
SCENARIO_INDEX = 1
TOP_N = 25

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

# Échauffement hors mesure (imports paresseux, caches).
pipeline.compute_summary_statistics_microsat(context, SCENARIO_INDEX, seed=1000)

# 1. Durée réelle, sans profileur.
t0 = perf_counter()
pipeline.compute_summary_statistics_microsat(context, SCENARIO_INDEX, seed=SEED)
wall_plain = perf_counter() - t0

# 2. Même particule sous cProfile.
profiler = cProfile.Profile()
profiler.enable()
pipeline.compute_summary_statistics_microsat(context, SCENARIO_INDEX, seed=SEED)
profiler.disable()

stats = pstats.Stats(profiler)
total_profiled = stats.total_tt

print(f"\n{REFERENCE_DIR}, seed {SEED}, scénario {SCENARIO_INDEX}")
print(
    f"Sans profileur : {wall_plain:.3f} s | sous cProfile : {total_profiled:.3f} s "
    f"(gonflement x{total_profiled / wall_plain:.2f})\n"
)

print(f"--- Top {TOP_N} par temps CUMULÉ (inclut les sous-appels) ---")
stats.sort_stats("cumulative").print_stats("bridge|msprime|tskit", TOP_N)

print(f"--- Top {TOP_N} par temps PROPRE (hors sous-appels) ---")
stats.sort_stats("tottime").print_stats(TOP_N)
