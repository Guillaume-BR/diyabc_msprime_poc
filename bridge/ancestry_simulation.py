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


def _draw_single_mutation_node(tree, rng: random.Random) -> int:
    """Tire un noeud de l'arbre (portant la mutation unique), avec une
    probabilité proportionnelle à la longueur de sa branche -- algorithme
    de Hudson (2002), tel que décrit dans la doc DIYABC (section 2.4.3) :
    "il est supposé qu'il y a eu une et une seule mutation dans l'arbre
    de coalescence des gènes échantillonnés".

    Validé empiriquement (20000 tirages sur un arbre test, proportions
    observées vs attendues alignées à <1% -- voir notes/exploration.md).
    """
    nodes = [u for u in tree.nodes() if u != tree.root]
    lengths = [tree.branch_length(u) for u in nodes]
    total = sum(lengths)
    target = rng.uniform(0, total)
    cumulative = 0.0
    for node, length in zip(nodes, lengths):
        cumulative += length
        if cumulative >= target:
            return node
    return nodes[-1]  # garde-fou contre un arrondi flottant en bord de plage


def simulate_snp_genotypes(
    tree_sequences: Iterator[msprime.TreeSequence],
    seed: int,
) -> Iterator[dict[str, list[int]]]:
    """Pour chaque TreeSequence (un locus = un arbre indépendant), tire
    une mutation UNIQUE selon l'algorithme de Hudson, et retourne les
    génotypes (0=ancestral, 1=dérivé) REGROUPÉS PAR POPULATION -- forme
    directement utilisable pour les formules q1/q2/HW/HB/FST1/ML1
    (sumstat.cpp), qui travaillent toujours par population, jamais par
    lignée individuelle isolée.

    Remplace mutate_independent_loci (modèle à taux fixe), qui était
    structurellement incorrect pour des SNP au sens DIYABC : la doc
    utilisateur (section 2.4.3) précise que les loci SNP sont (par
    construction) toujours polymorphes, avec exactement une mutation
    dans tout l'arbre -- pas un processus de Poisson à taux variable,
    qui pourrait produire des loci monomorphes ou multi-mutés.

    Retourne un itérateur de dict {nom_population: [génotypes...]}, où
    nom_population est le nom donné à add_population dans
    demography_builder.py (ex: "pop1", "pop2"...). L'indice interne
    msprime (0-based, dans l'ordre d'ajout des populations -- vérifié
    empiriquement) est utilisé pour interroger
    TreeSequence.samples(population=i), puis traduit vers ce nom.
    """
    rng = random.Random(seed)
    for ts in tree_sequences:
        tree = ts.first()
        mutated_node = _draw_single_mutation_node(tree, rng)
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