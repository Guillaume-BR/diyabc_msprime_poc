# validate_observed_stats_poolseq.py
"""
Valide nos formules de stats PoolSeq (bridge/summary_statistics.py) en
les comparant au vrai statobs.txt de DIYABC -- calculé sur les données
OBSERVÉES réelles, sans aucune simulation msprime. Isole la correction
des FORMULES de celle du simulateur (contrairement à validate_stats.py,
qui teste la simulation mais ne peut pas comparer au binaire réel faute
d'un writer .snp PoolSeq).

Nécessite un statobs.txt déjà généré par DIYABC pour REFERENCE_DIR (ex:
en lançant `./diyabc -p ./ -R ALL -r 1 -g 1` une fois dans ce dossier).

Résultat (22/07/2026, toy_example4) : 0/130 divergence >1% -- mais
SEULEMENT une fois observed_reads() corrigée pour purger les loci sous
le seuil MRC avant de prendre les num_loci premiers (voir
notes/exploration.md, entrée du 22/07/2026, et la mémoire persistante
project_poolseq_support_in_progress -- sans cette purge, 130/130
divergeaient).
"""

from bridge.ancestry_simulation import (
    _reindex_reads_by_msprime_name,
    build_samples_argument,
)
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import observed_reads
from bridge.statobs_parser import parse_statobs
from bridge.summary_statistics import compute_all_statistics_poolseq

REFERENCE_DIR = "reference/toy_example4"
SNP = f"{REFERENCE_DIR}/pseudo_observed_DATASET_SNP_POOLSEQ_4pops_Scenario3_MER.snp"

with open(f"{REFERENCE_DIR}/headerRF.txt") as f:
    header_text = f.read()
num_loci = parse_loci_description(header_text).total_loci["A"]

pool_sizes = build_samples_argument(SNP)
population_names = list(pool_sizes.keys())
reads = _reindex_reads_by_msprime_name(observed_reads(SNP), SNP)[:num_loci]

our_stats = compute_all_statistics_poolseq(reads, population_names, pool_sizes)
with open(f"{REFERENCE_DIR}/statobs.txt") as f:
    diyabc_stats = parse_statobs(f.read())

print(f"\n{'Statistique':<15} {'Nous':>15} {'DIYABC':>15} {'écart relatif':>15}")
print("-" * 65)
n_bad = 0
for key in sorted(our_stats):
    if key not in diyabc_stats:
        print(f"? {key:<15} {our_stats[key]:>15.6f} {'(absent)':>15}")
        continue
    ours, theirs = our_stats[key], diyabc_stats[key]
    rel_err = abs(ours - theirs) / abs(theirs) if theirs != 0 else abs(ours - theirs)
    ok = "✓" if rel_err < 0.01 else "✗"
    if rel_err >= 0.01:
        n_bad += 1
    print(f"{ok} {key:<15} {ours:>15.6f} {theirs:>15.6f} {rel_err:>15.2e}")

print(f"\n{n_bad} / {len(our_stats)} statistiques divergent de plus de 1%.")
