"""
Simulation de l'ancestralité (coalescence) pour des loci SNP indépendants
(pas de recombinaison interne, pas de liaison génétique entre loci) --
cas du dataset human, dont les 51250 loci sont déclarés <A> (autosomal,
diploïde classique).

Chaque locus est un réplicat indépendant : msprime.sim_ancestry avec
num_replicates dérive correctement une graine distincte par réplicat à
partir d'une seule random_seed, de façon reproductible (vérifié
empiriquement, voir notes/exploration.md).

ploidy=2 (valeur par défaut de msprime) est cohérent avec le code
d'héritage <A> de human : chaque "sample individual" = 2 lignées
génomiques, et l'échelle de temps de la coalescence est calée en
générations diploïdes -- cohérent avec les bornes des priors de temps
(en générations) de header.txt.
"""

from collections.abc import Iterator
import random

import msprime

from bridge.observed_data import population_index_to_name, count_samples_per_population


def build_samples_argument(
    snp_file_path: str,
) -> dict[str, int]:
    """Construit l'argument `samples` attendu par msprime.sim_ancestry :
    {nom_population_msprime: nombre_d_individus}, où le nom de population
    msprime ("pop1", "pop2"...) correspond à l'indice utilisé dans
    header.txt, mappé sur le nombre réel d'individus observés pour la
    population correspondante (voir observed_data.py pour la
    justification du mapping par ordre d'apparition).
    """
    index_to_name = population_index_to_name(snp_file_path)
    counts_by_name = count_samples_per_population(snp_file_path)

    return {
        f"pop{index}": counts_by_name[name]
        for index, name in index_to_name.items()
    }


def simulate_independent_loci(
    demography: msprime.Demography,
    samples: dict[str, int],
    num_loci: int,
    seed: int,
) -> Iterator[msprime.TreeSequence]:
    """Simule num_loci généalogies indépendantes (un locus SNP = un
    réplicat, pas de recombinaison interne ni de liaison entre loci),
    sous la démographie donnée.

    Retourne un itérateur (pas une liste) : pour 51250 loci, matérialiser
    toutes les TreeSequence en mémoire simultanément serait coûteux --
    l'appelant doit consommer cet itérateur au fil de l'eau (ex: pour
    calculer des statistiques résumées locus par locus).
    """
    return msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=1,
        num_replicates=num_loci,
        random_seed=seed,
        ploidy=2,
    )


def mutate_independent_loci(
    tree_sequences: Iterator[msprime.TreeSequence],
    mutation_rate: float,
    seed: int,
) -> Iterator[msprime.TreeSequence]:
    """Applique une mutation binaire (alleles "0"/"1", ancestral toujours
    "0") sur chaque TreeSequence de l'itérateur, avec une graine DISTINCTE
    par locus -- dérivée de `seed` et de l'indice du locus, pour éviter
    toute corrélation artificielle entre loci tout en restant reproductible
    avec une seule seed globale.

    LIMITE CONNUE (voir notes/exploration.md) : le modèle de mutation réel
    utilisé par DIYABC pour les SNP de human n'a pas été élucidé
    précisément (aucun MEANMU/MEANSNI déclaré dans header.txt, aucune
    branche [M]/[S]/[P] activée). Ce modèle binaire à taux fixe est un
    choix simplifié et raisonnable pour le POC, pas une reproduction
    fidèle de l'algorithme DIYABC -- à revoir avant toute comparaison
    statistique fine avec le reftableRF.bin de référence.
    """
    model = msprime.BinaryMutationModel()
    rng = random.Random(seed)
    for ts in tree_sequences:
        locus_seed = rng.randrange(1, 2**31)
        yield msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            model=model,
            random_seed=locus_seed,
        )