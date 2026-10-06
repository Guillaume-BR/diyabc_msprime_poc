"""
Constantes partagées entre les modules de bridge/ : offsets de graine
(pour qu'aucun tirage aléatoire indépendant ne collisionne avec un
autre -- voir feedback_seed_reuse_pattern, ce bug est déjà arrivé 5+
fois dans ce projet quand ces constantes étaient éparpillées) et
tailles de lot (batching de simulate_independent_loci), et dérivation des
graines (`_locus_seed`).
"""

import numpy as np

# ── Offsets de graine (chacun doit rester unique dans ce fichier) ──────────

_GROUP_PRIOR_SEED_OFFSET = 10_000_000  # parameter_sampling.py
_MAF_REJECTION_SEED_OFFSET = 20_000_000  # ancestry_simulation.py
_MRC_REJECTION_SEED_OFFSET = 30_000_000  # ancestry_simulation.py
_BINOMIAL_SEED_OFFSET = 40_000_000  # ancestry_simulation.py
_SCENARIO_DRAW_SEED_OFFSET = 50_000_000  # reftable_loop.py
_KAPPA1_SEED_OFFSET = 60_000_000  # ancestry_simulation.py
_KAPPA2_SEED_OFFSET = 70_000_000  # ancestry_simulation.py
_MUS_RATE_SEED_OFFSET = 80_000_000  # ancestry_simulation.py
_SITE_RATE_SEED_OFFSET = 90_000_000  # ancestry_simulation.py
_MUTATION_SEED_OFFSET = 100_000_000  # ancestry_simulation.py
_ANCESTRY_SEED_OFFSET = 110_000_000  # ancestry_simulation.py
_SHARED_M_ANCESTRY_SEED_OFFSET = 120_000_000  # ancestry_simulation.py
_SHARED_Y_ANCESTRY_SEED_OFFSET = 130_000_000  # ancestry_simulation.py
_MICROSAT_MUT_RATE_SEED_OFFSET = 140_000_000  # ancestry_simulation.py
_MICROSAT_PGEOM_SEED_OFFSET = 150_000_000  # ancestry_simulation.py
_MICROSAT_SNI_SEED_OFFSET = 160_000_000  # ancestry_simulation.py
_LIKELIHOOD_SEED_OFFSET = 170_000_000  # summary_statistics.py
_MRC_BATCH_SEED_OFFSET = 180_000_000  # ancestry_simulation.py (with_mrc_filter)
_MAF_BATCH_SEED_OFFSET = 190_000_000  # ancestry_simulation.py (with_maf_filter)
_GROUP_STAT_SEED_OFFSET = (
    200_000_000  # summary_statistics.py (AML, une graine par groupe)
)

# pipeline.py -- un offset par type de locus SNP (<A>/<H>/<X>/<Y>/<M>),
# pour que _simulate_genotypes_for_all_locus_types dérive une graine
# distincte par type sans jamais réutiliser la même seed brute pour deux
# types différents (voir pipeline.py pour la justification empirique).
_LOCUS_TYPE_SEED_OFFSET = {
    "A": 0,
    "H": 1_000_000,
    "X": 2_000_000,
    "Y": 3_000_000,
    "M": 4_000_000,
}

# ── Tailles de lot ──────────────────────────────────────────────────────────

_MAF_BATCH_SIZE = 20  # ancestry_simulation.py (with_maf_filter)
_MRC_BATCH_SIZE = 20  # ancestry_simulation.py (with_mrc_filter)


# ── Dérivation des graines ─────────────────────────────────────────────────


def _locus_seed(seed: int, offset: int, locus_index: int) -> int:
    """Dérive une graine indépendante à partir de celle de la particule.

    Remplace toute addition `seed + offset + i` : les graines de particule
    sont consécutives, donc `seed + offset + i` donnait au locus `i` de la
    particule `s` la même graine qu'au locus `i - 1` de la particule `s + 1`,
    et chaque flux aléatoire était réutilisé par ~nloci particules voisines
    qui n'étaient donc pas indépendantes (voir "`DTA_2_2` residual (05/10)"
    dans CLAUDE.md). `SeedSequence` mélange les trois entrées, si bien que
    deux triplets distincts ne partagent jamais de flux.

    Args:
        seed: La graine de la particule.
        offset: L'offset du type de tirage (arbre, mutation, lots...).
        locus_index: L'indice dans ce tirage (locus, lot, tentative, groupe ;
            0 pour un tirage unique par particule).

    Returns:
        Un entier dans [1, 2**32 - 1] (msprime rejette `seed=0`).
    """
    state = np.random.SeedSequence([seed, offset, locus_index]).generate_state(1)[0]
    return int(state) % (2**32 - 1) + 1
