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
import numpy as np

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
        f"pop{index}": counts_by_name[name] for index, name in index_to_name.items()
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


def _draw_single_mutation_edge_child(ts, rng: random.Random) -> int:
    """Tire le noeud portant la mutation unique, avec probabilité
    proportionnelle à la longueur de sa branche -- algorithme de Hudson,
    entièrement vectorisé via les tables (pas d'appel branch_length() par
    noeud). Valable pour un arbre unique (sequence_length=1).

    Chaque edge = une branche (couple parent-enfant) ; edges.child liste
    donc tous les noeuds ayant une branche au-dessus d'eux (tous sauf la
    racine). Longueur = time[parent] - time[child], calculé en numpy.

    Validé empiriquement (proportions observées vs attendues <1% ; valeurs
    de statistiques identiques à la version par branch_length() -- voir
    notes/exploration.md).
    """
    edges = ts.tables.edges
    node_times = ts.tables.nodes.time
    lengths = node_times[edges.parent] - node_times[edges.child]

    total = lengths.sum()
    target = rng.uniform(0, total)
    idx = np.searchsorted(np.cumsum(lengths), target)
    if idx >= len(edges.child):
        idx = len(edges.child) - 1
    return int(edges.child[idx])


def _draw_single_mutation_node_fast(tree, ts, rng: random.Random) -> int:
    """Version vectorisée : longueur de branche = temps(parent) - temps(noeud),
    calculé en numpy sur tous les noeuds d'un coup."""
    node_times = ts.tables.nodes.time  # numpy array, tous les temps
    nodes = np.fromiter((u for u in tree.nodes() if u != tree.root), dtype=np.int64)
    parents = np.fromiter((tree.parent(u) for u in nodes), dtype=np.int64)
    lengths = node_times[parents] - node_times[nodes]

    total = lengths.sum()
    target = rng.uniform(0, total)
    cumulative = np.cumsum(lengths)
    idx = np.searchsorted(cumulative, target)
    if idx >= len(nodes):
        idx = len(nodes) - 1
    return int(nodes[idx])


def _draw_single_mutation_node_vectorized(ts, rng: random.Random):
    """Tire le noeud portant la mutation, entièrement en numpy depuis les
    tables (pas d'appel branch_length() par noeud). Valable pour un arbre
    unique (sequence_length=1, une seule TreeSequence)."""
    edges = ts.tables.edges
    node_times = ts.tables.nodes.time

    children = edges.child  # array des noeuds enfants
    parents = edges.parent  # array des parents
    lengths = node_times[parents] - node_times[children]  # longueurs, vectorisé

    total = lengths.sum()
    target = rng.uniform(0, total)
    idx = np.searchsorted(np.cumsum(lengths), target)
    if idx >= len(children):
        idx = len(children) - 1
    return int(children[idx])


def simulate_snp_genotypes(
    tree_sequences: Iterator[msprime.TreeSequence],
    seed: int,
) -> Iterator[dict[str, list[int]]]:
    """Pour chaque TreeSequence (un locus = un arbre indépendant), tire
    une mutation UNIQUE selon l'algorithme de Hudson (vectorisé), et
    retourne les génotypes (0=ancestral, 1=dérivé) REGROUPÉS PAR
    POPULATION.

    Voir _draw_single_mutation_edge_child pour l'algorithme de tirage, et
    la docstring d'origine pour la justification du modèle (doc DIYABC
    section 2.4.3 : exactement une mutation par locus, locus toujours
    polymorphe).
    """
    rng = random.Random(seed)
    for ts in tree_sequences:
        tree = ts.first()
        mutated_node = _draw_single_mutation_edge_child(ts, rng)
        derived_samples = set(tree.samples(mutated_node))

        genotypes_by_population = {}
        for pop_index, population in enumerate(ts.tables.populations):
            pop_name = population.metadata.get("name") if population.metadata else None
            sample_ids = ts.samples(population=pop_index)
            if len(sample_ids) == 0:
                continue
            genotypes_by_population[pop_name] = [
                1 if s in derived_samples else 0 for s in sample_ids
            ]
        yield genotypes_by_population
