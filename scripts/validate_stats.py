# validate_stats.py
"""
Simule un petit reftable PoolSeq (toy_example4) et calcule nos 130
statistiques dessus -- test de fumée, pas une comparaison au binaire
DIYABC : `general` ne peut pas recalculer de stats sur nos reads
simulés tant qu'il n'existe pas d'équivalent PoolSeq de
bridge.snp_writer.write_snp_file (celui-ci n'écrit que le format IndSeq
IND SEX POP). Voir validate_observed_stats_poolseq.py pour une vraie
comparaison contre le binaire réel -- sur les données OBSERVÉES
(statobs.txt), pas simulées.
"""

from pathlib import Path

from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_poolseq_reads_with_mrc_filter,
)
from bridge.loci_parser import parse_loci_description
from bridge.pipeline import build_random_demography_for_scenario_index, read_header_text
from bridge.summary_statistics import compute_all_statistics_poolseq

REFERENCE_DIR = Path("reference/toy_example4")

header_text = read_header_text(REFERENCE_DIR)
demography, values = build_random_demography_for_scenario_index(
    header_text, scenario_index=1, seed=42
)

snp_filename = header_text.splitlines()[0].strip()
snp_path = f"{REFERENCE_DIR}/{snp_filename}"
num_loci = parse_loci_description(header_text).total_loci["A"]

reads_list = list(
    simulate_poolseq_reads_with_mrc_filter(demography, snp_path, num_loci, seed=42)
)
pool_sizes = build_samples_argument(snp_path)
population_names = list(pool_sizes.keys())

our_stats = compute_all_statistics_poolseq(reads_list, population_names, pool_sizes)

print(f"Paramètres tirés : {values}")
print(f"{len(reads_list)} loci simulés, {len(our_stats)} statistiques calculées.\n")
for key in sorted(our_stats):
    print(f"{key:<15} {our_stats[key]:>15.6f}")
