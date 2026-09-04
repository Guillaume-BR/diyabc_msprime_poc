"""
Constantes partagées entre les modules de bridge/ : offsets de graine
(pour qu'aucun tirage aléatoire indépendant ne collisionne avec un
autre -- voir feedback_seed_reuse_pattern, ce bug est déjà arrivé 5+
fois dans ce projet quand ces constantes étaient éparpillées) et
tailles de lot (batching de simulate_independent_loci).
"""

# ── Offsets de graine (chacun doit rester unique dans ce fichier) ──────────

_GROUP_PRIOR_SEED_OFFSET = 10_000_000  # parameter_sampling.py
_MAF_REJECTION_SEED_OFFSET = 2_000_000  # ancestry_simulation.py
_MRC_REJECTION_SEED_OFFSET = 3_000_000  # ancestry_simulation.py
_BINOMIAL_SEED_OFFSET = 4_000_000  # ancestry_simulation.py
_SCENARIO_DRAW_SEED_OFFSET = 50_000_000  # reftable_loop.py
_KAPPA1_SEED_OFFSET = 60_000_000  # ancestry_simulation.py
_KAPPA2_SEED_OFFSET = 70_000_000  # ancestry_simulation.py
_MUS_RATE_SEED_OFFSET = 80_000_000  # ancestry_simulation.py
_SITE_RATE_SEED_OFFSET = 90_000_000  # ancestry_simulation.py
_MUTATION_SEED_OFFSET = 100_000_000  # ancestry_simulation.py
_ANCESTRY_SEED_OFFSET = 110_000_000  # ancestry_simulation.py
_SHARED_M_ANCESTRY_SEED_OFFSET = 120_000_000  # ancestry_simulation.py
_SHARED_Y_ANCESTRY_SEED_OFFSET = 130_000_000  # ancestry_simulation.py

# pipeline.py -- un offset par type de locus SNP (<A>/<H>/<X>/<Y>/<M>),
# pour que _simulate_genotypes_for_all_locus_types dérive une graine
# distincte par type sans jamais réutiliser la même seed brute pour deux
# types différents (voir pipeline.py pour la justification empirique).
_LOCUS_TYPE_SEED_OFFSET = {
    "A": 0,
    "H": 10_000_000,
    "X": 20_000_000,
    "Y": 30_000_000,
    "M": 40_000_000,
}

# ── Tailles de lot ──────────────────────────────────────────────────────────

_MAF_BATCH_SIZE = 20  # ancestry_simulation.py (with_maf_filter)
_MRC_BATCH_SIZE = 20  # ancestry_simulation.py (with_mrc_filter)
