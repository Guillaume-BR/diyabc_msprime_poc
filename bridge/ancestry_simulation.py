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

import itertools
import random
from collections.abc import Iterator

import msprime
import numpy as np

from bridge.demography_builder import rescale_demography
from bridge.observed_data import (
    coalescence_coefficient,
    count_samples_per_population,
    individual_sexes_per_population,
    parse_sex_ratio,
    population_index_to_name,
)


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


def build_sex_stratified_samples_argument(
    snp_file_path: str,
) -> list[msprime.SampleSet]:
    """Construit l'argument `samples` pour un locus <X> : contrairement à
    build_samples_argument (un compte par population, ploidy uniforme),
    <X> a besoin d'une ploidy DIFFÉRENTE par individu selon son sexe
    (femelles=2 copies, mâles=1 -- voir ParticleC::calploidy,
    particuleC.cpp:220-233), donc une liste de msprime.SampleSet plutôt
    qu'un simple dict.

    Retourne 2 SampleSet par population : un pour les femelles
    (ploidy=2), un pour les mâles (ploidy=1) -- population= doit être le
    nom msprime ("pop1", "pop2"...), PAS le nom réel du fichier .snp
    (même traduction que build_samples_argument, via
    population_index_to_name).

    IMPORTANT -- le ploidy PAR SampleSet ne contrôle QUE le nombre de
    lignées regroupées par individu dans le résultat, PAS le taux de
    coalescence : cette liste doit être utilisée avec
    simulate_independent_loci(..., ploidy=1) et une `demography` déjà
    rescalée via rescale_demography(demography,
    coalescence_coefficient("X", sex_ratio) / 2) -- vérifié
    empiriquement avec le mentor que c'est le ploidy GLOBAL de
    sim_ancestry qui interprète initial_size, pas celui des SampleSet.

    Lève une erreur si un individu a le sexe "9" (inconnu) -- ex:
    human_snp_all22chr_maf5.snp, où AUCUN individu n'est sexé (dataset
    <A>-only) : on ne peut pas construire un échantillonnage <X> dessus,
    mieux vaut le signaler explicitement que de produire un résultat
    silencieusement faux (individual_sexes_per_population laisse ce
    choix à l'appelant, c'est ici qu'il se prend).
    """
    samples_by_population = []

    sexes_by_population = individual_sexes_per_population(snp_file_path)
    index_to_name = population_index_to_name(snp_file_path)
    name_to_index = {name: index for index, name in index_to_name.items()}
    for name in sexes_by_population:
        if "9" in sexes_by_population[name]:
            raise ValueError(
                f"Individu avec sexe inconnu trouvé dans la population {name}"
            )
        nb_femelles = sexes_by_population[name].count("F")
        nb_males = sexes_by_population[name].count("M")
        pop_index = name_to_index[name]

        # On ajoute les SampleSet pour les femelles et les mâles avec le ploidy approprié
        samples_by_population.append(
            msprime.SampleSet(
                num_samples=nb_femelles, population=f"pop{pop_index}", ploidy=2
            )
        )
        samples_by_population.append(
            msprime.SampleSet(
                num_samples=nb_males, population=f"pop{pop_index}", ploidy=1
            )
        )

    return samples_by_population


def build_male_only_samples_argument(snp_file_path: str) -> dict[str, int]:
    """Construit l'argument `samples` attendu par msprime.sim_ancestry pour
    un locus <Y> : {nom_population_msprime: nombre_d_individus_mâles}, où
    le nom de population msprime ("pop1", "pop2"...) correspond à
    l'indice utilisé dans header.txt, mappé sur le nombre réel
    d'individus MÂLES observés pour la population correspondante (voir
    observed_data.py pour la justification du mapping par ordre
    d'apparition).

    PAS pour <M> : le mitochondrial est transmis uniquement par les
    mères, mais présent et échantillonné chez TOUS les individus
    (mâles et femelles), contrairement à <Y> qui n'existe que chez les
    mâles -- <M> doit réutiliser build_samples_argument (tout le monde)
    avec ploidy=1, pas cette fonction.
    """
    samples_by_population = {}

    sexes_by_population = individual_sexes_per_population(snp_file_path)
    index_to_name = population_index_to_name(snp_file_path)
    name_to_index = {name: index for index, name in index_to_name.items()}
    for name in sexes_by_population:
        if "9" in sexes_by_population[name]:
            raise ValueError(
                f"Individu avec sexe inconnu trouvé dans la population {name}"
            )
        nb_males = sexes_by_population[name].count("M")
        pop_index = name_to_index[name]
        samples_by_population[f"pop{pop_index}"] = nb_males

    return samples_by_population


def simulate_independent_loci(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    seed: int,
    ploidy: int = 2,
) -> Iterator[msprime.TreeSequence]:
    """Simule num_loci généalogies indépendantes (un locus SNP = un
    réplicat, pas de recombinaison interne ni de liaison entre loci),
    sous la démographie donnée.

    Retourne un itérateur (pas une liste) : pour 51250 loci, matérialiser
    toutes les TreeSequence en mémoire simultanément serait coûteux --
    l'appelant doit consommer cet itérateur au fil de l'eau (ex: pour
    calculer des statistiques résumées locus par locus).

    samples : dict[str, int] (un compte par population, ploidy uniforme
    -- <A>/<M>, voir build_samples_argument) ou list[msprime.SampleSet]
    (ploidy hétérogène par sous-groupe au sein d'une population -- <X>,
    voir build_sex_stratified_samples_argument). Les deux formes sont
    transmises telles quelles à msprime.sim_ancestry, qui les accepte
    indifféremment.

    ploidy : 2 (défaut) pour <A>, cohérent avec une transmission diploïde
    classique -- chaque "sample individual" de `samples` compte pour 2
    lignées génomiques. Pour <Y>/<M>, passer ploidy=1 avec une
    `demography` déjà rescalée par rescale_demography (voir
    demography_builder.py) : ces loci sont haploïdes (une seule copie de
    gène transmise), et le facteur de rescaling de Ne
    (coalescence_coefficient, observed_data.py) suppose cette
    combinaison ploidy=1 + Ne rescalé, pas ploidy=2 + Ne d'origine. Pour
    <X>, passer aussi ploidy=1 (voir build_sex_stratified_samples_argument
    : c'est le ploidy PAR SampleSet, pas ce paramètre global, qui donne
    2 copies aux femelles et 1 aux mâles -- ce paramètre-ci ne fixe que
    le taux de coalescence, via la Demography déjà rescalée).
    """
    return msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=1,
        num_replicates=num_loci,
        random_seed=seed,
        ploidy=ploidy,
    )


def simulate_shared_ancestry_loci(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    seed: int,
    ploidy: int = 1,
) -> Iterator[msprime.TreeSequence]:
    """Simule UNE SEULE généalogie puis la retourne répétée num_loci fois
    -- pour <Y>/<M>, dont tous les loci d'un même type partagent la même
    généalogie réelle (non-recombinants, transmission uniparentale),
    contrairement à <A>/<X> qui tirent un arbre indépendant par locus
    (simulate_independent_loci). Reproduit le comportement de
    particuleC.cpp:2422-2435 (GeneTreeY/GeneTreeM : premier locus <Y> ou
    <M> tire un arbre normalement, tous les suivants COPIENT ce même
    arbre -- seule la mutation change d'un locus à l'autre).

    samples/ploidy : même contrat que simulate_independent_loci (voir sa
    docstring) -- cette fonction ne fait que réutiliser
    simulate_independent_loci avec num_loci=1, elle ne réinterprète pas
    ces paramètres.

    IMPORTANT -- ne PAS réimplémenter le tirage de mutation ici :
    simulate_snp_genotypes(tree_sequences, seed) lit déjà tree_sequences
    au fil de l'eau sans jamais modifier les TreeSequence qu'elle reçoit,
    et son rng avance à chaque itération -- lui donner le MÊME objet
    TreeSequence répété num_loci fois (au lieu de num_loci objets
    différents) suffit à obtenir num_loci mutations indépendantes sur
    UNE SEULE généalogie, sans aucune modification de cette fonction
    (vérifié empiriquement : 5 répétitions du même arbre -> 5 génotypes
    différents).
    """

    shared_genealogy = next(
        simulate_independent_loci(
            demography, samples, num_loci=1, seed=seed, ploidy=ploidy
        )
    )
    return itertools.repeat(shared_genealogy, num_loci)


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


# def _draw_single_mutation_node_fast(tree, ts, rng: random.Random) -> int:
#    """Version vectorisée : longueur de branche = temps(parent) - temps(noeud),
#    calculé en numpy sur tous les noeuds d'un coup."""
#    node_times = ts.tables.nodes.time  # numpy array, tous les temps
#    nodes = np.fromiter((u for u in tree.nodes() if u != tree.root), dtype=np.int64)
#    parents = np.fromiter((tree.parent(u) for u in nodes), dtype=np.int64)
#    lengths = node_times[parents] - node_times[nodes]
#
#    total = lengths.sum()
#    target = rng.uniform(0, total)
#    cumulative = np.cumsum(lengths)
#    idx = np.searchsorted(cumulative, target)
#    if idx >= len(nodes):
#        idx = len(nodes) - 1
#    return int(nodes[idx])
#
#
# def _draw_single_mutation_node_vectorized(ts, rng: random.Random):
#    """Tire le noeud portant la mutation, entièrement en numpy depuis les
#    tables (pas d'appel branch_length() par noeud). Valable pour un arbre
#    unique (sequence_length=1, une seule TreeSequence)."""
#    edges = ts.tables.edges
#    node_times = ts.tables.nodes.time
#
#    children = edges.child  # array des noeuds enfants
#    parents = edges.parent  # array des parents
#    lengths = node_times[parents] - node_times[children]  # longueurs, vectorisé
#
#    total = lengths.sum()
#    target = rng.uniform(0, total)
#    idx = np.searchsorted(np.cumsum(lengths), target)
#    if idx >= len(children):
#        idx = len(children) - 1
#    return int(children[idx])


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


def simulate_genotypes_for_locus_type(
    demography: msprime.Demography,
    snp_file_path: str,
    locus_type: str,
    num_loci: int,
    seed: int,
) -> Iterator[dict[str, list[int]]]:
    """Point d'entrée unique par type de locus : choisit la bonne
    combinaison samples/demography-rescalée-ou-non/ploidy/fonction de
    simulation-indépendante-ou-partagée selon locus_type, puis retourne
    les génotypes simulés (même contrat de sortie que
    simulate_snp_genotypes, qu'on appelle en dernière étape dans tous
    les cas -- elle ne dépend jamais de locus_type elle-même).

    `demography` : la démographie <A> "de base" (construite par
    build_demography, PAS encore rescalée) -- c'est CETTE fonction qui
    décide si/comment la rescaler selon locus_type, jamais l'appelant.

    sex_ratio n'est PAS un paramètre : il est dérivé automatiquement de
    snp_file_path via parse_sex_ratio, comme tout le reste (samples,
    sexes par individu) -- l'appelant n'a jamais besoin de le connaître.

    Dispatch (voir le tableau établi avec le mentor) :
      - "A" : build_samples_argument, demography TELLE QUELLE (pas de
        rescale_demography), ploidy=2, simulate_independent_loci.
      - "H" : build_samples_argument, demography rescalée par
        coalescence_coefficient("H", sex_ratio) / 2, ploidy=1,
        simulate_independent_loci.
      - "X" : build_sex_stratified_samples_argument, demography
        rescalée par coalescence_coefficient("X", sex_ratio) / 2,
        ploidy=1, simulate_independent_loci.
      - "Y" : build_male_only_samples_argument, demography rescalée par
        coalescence_coefficient("Y", sex_ratio) / 2, ploidy=1,
        simulate_shared_ancestry_loci (arbre unique partagé).
      - "M" : build_samples_argument (TOUT le monde, pas mâles seuls --
        voir la docstring de build_male_only_samples_argument sur ce
        point précis), demography rescalée par
        coalescence_coefficient("M", sex_ratio) / 2, ploidy=1,
        simulate_shared_ancestry_loci.
      - tout autre locus_type : lever NotImplementedError (même style
        que coalescence_coefficient pour un type inconnu).
    """

    sex_ratio = parse_sex_ratio(snp_file_path)

    if locus_type == "Y":
        samples = build_male_only_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        trees = simulate_shared_ancestry_loci(
            rescaled_demography, samples, num_loci, seed, ploidy=1
        )
        return simulate_snp_genotypes(trees, seed)
    elif locus_type == "M":
        samples = build_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        trees = simulate_shared_ancestry_loci(
            rescaled_demography, samples, num_loci, seed, ploidy=1
        )
        return simulate_snp_genotypes(trees, seed)
    elif locus_type == "A":
        samples = build_samples_argument(snp_file_path)
        trees = simulate_independent_loci(demography, samples, num_loci, seed, ploidy=2)
        return simulate_snp_genotypes(trees, seed)
    elif locus_type == "H":
        samples = build_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        trees = simulate_independent_loci(
            rescaled_demography, samples, num_loci, seed, ploidy=1
        )
        return simulate_snp_genotypes(trees, seed)
    elif locus_type == "X":
        samples = build_sex_stratified_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        trees = simulate_independent_loci(
            rescaled_demography, samples, num_loci, seed, ploidy=1
        )
        return simulate_snp_genotypes(trees, seed)
    else:
        raise NotImplementedError(f"Type de locus non supporté: {locus_type!r}")
