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
from pathlib import Path

import msprime
import numpy as np
import tskit

from bridge.demography_builder import rescale_demography
from bridge.header_dataclasses import LociDescriptionDetailed
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    base_frequency_by_locus,
    coalescence_coefficient,
    count_samples_per_population,
    individual_sexes_per_population,
    observed_count_population,
    observed_mrc,
    observed_reads,
    observed_sequences,
    parse_maf_ratio,
    parse_mrc_ratio,
    parse_sex_ratio,
    population_index_to_name,
)
from bridge.parameter_sampling import (
    draw_group_parameter_values,
    sample_site_rates,
    sampling_group_local_param,
)
from bridge.prior_parser import get_parameter_used_by_model, parse_group_priors

# Offset de graine dédié à la boucle de rejet MAF et du rejet MRC, distinct du +1_000_000
# déjà utilisé partout ailleurs dans le projet pour séparer la graine de
# mutation de la graine de généalogie (voir notebooks/scripts) -- ne
# jamais réutiliser 1_000_000 ici, ça collisionnerait avec cette
# convention existante plutôt qu'avec autre chose.
_MAF_REJECTION_SEED_OFFSET = 2_000_000
_MRC_REJECTION_SEED_OFFSET = 3_000_000

# Taille de lot PLANCHER pour with_maf_filter : un seul appel
# simulate_independent_loci (num_replicates=batch_size) au lieu d'un appel
# par tentative de rejet. Mesuré empiriquement (script jetable,
# toy_example3 scenario3, maf=0.05, 100 loci) : le setup Python/msprime
# (construction de la Demography, initialisation du simulateur) domine
# largement le coût réel de coalescence pour ce genre de petit échantillon
# -- répété à CHAQUE tentative sans batching, il fait payer ce setup une
# fois par rejet au lieu de le mutualiser. batch_size=20 mesuré à ~5.6x
# plus rapide que l'appel un-par-un (batch_size=1) sur ce cas.
#
# PLANCHER, pas la taille réelle : with_maf_filter calcule
# batch_size = max(_MAF_BATCH_SIZE, num_loci // 4) -- un lot fixe à 20 ne
# scale pas avec num_loci, donc le NOMBRE de lots (et le setup payé par
# lot) croît proportionnellement à num_loci, ce qui dégrade le facteur de
# gain sur les gros num_loci (mesuré empiriquement le 24/07/2026 :
# batch=20 fixe donne encore ~1.35-1.44x d'écart vs batch=num_loci//4 sur
# 500 loci, alors que le gain est négligeable sur 100 loci -- déjà sur le
# plateau de rendement décroissant à ce volume). num_loci//4 capture déjà
# ~94% du gain d'un lot égal à num_loci entier ; pas la peine d'aller
# jusque là.
_MAF_BATCH_SIZE = 20

# Taille de lot pour with_mrc_filter (PoolSeq) : contrairement à
# _MAF_BATCH_SIZE, ce lot est PARTAGÉ ENTRE TOUS LES LOCI (pas de scaling
# avec num_loci ici -- mesuré empiriquement le 24/07/2026, toy_example4,
# mrc=5 : le partage du pool entre loci apporte ~1.5x, mais faire varier
# la taille du lot une fois le pool partagé n'apporte quasiment rien en
# plus). Voir with_mrc_filter pour le détail du design.
_MRC_BATCH_SIZE = 20

# pour séparer le tirage binomial du tirage de mutation
_BINOMIAL_SEED_OFFSET = 4_000_000

# pour séparer les tirages des kappas (un tirage par locus, donc un tirage par réplicat) du tirage de mutation
_KAPPA1_SEED_OFFSET = 60_000_000
_KAPPA2_SEED_OFFSET = 70_000_000
_MUS_RATE_SEED_OFFSET = 80_000_000
_SITE_RATE_SEED_OFFSET = 90_000_000
_MUTATION_SEED_OFFSET = 100_000_000
_ANCESTRY_SEED_OFFSET = 110_000_000
_SHARED_M_ANCESTRY_SEED_OFFSET = 120_000_000


# ── Construction de l'argument samples (un builder par type de locus) ──────


def build_samples_argument(
    snp_file_path: str,
) -> dict[str, int]:
    """Construit l'argument `samples` de msprime.sim_ancestry pour un locus <A>.

    Le nom de population msprime ("pop1", "pop2"...) correspond à
    l'indice utilisé dans header.txt, mappé sur le nombre réel
    d'individus observés pour la population correspondante (voir
    observed_data.py pour la justification du mapping par ordre
    d'apparition).

    Un seul appel à count_samples_per_population (pas
    population_index_to_name EN PLUS, qui relit et rescanne tout le
    fichier .snp pour ne faire que redériver les mêmes clés dans le même
    ordre) -- l'indice 1-based se déduit directement de la position dans
    ce même dict, garanti dans l'ordre de première apparition (voir sa
    docstring).

    Args:
        snp_file_path: Chemin du fichier .snp.

    Returns:
        Un dict {nom_population_msprime: nombre_d_individus}.
    """
    counts_by_name = count_samples_per_population(snp_file_path)

    return {
        f"pop{index}": count
        for index, count in enumerate(counts_by_name.values(), start=1)
    }


def build_sex_stratified_samples_argument(
    snp_file_path: str,
) -> list[msprime.SampleSet]:
    """Construit l'argument `samples` de msprime.sim_ancestry pour un locus <X>.

    Contrairement à build_samples_argument (un compte par population,
    ploidy uniforme), <X> a besoin d'une ploidy DIFFÉRENTE par individu
    selon son sexe (femelles=2 copies, mâles=1 -- voir
    ParticleC::calploidy, particuleC.cpp:220-233), donc une liste de
    msprime.SampleSet plutôt qu'un simple dict : 2 SampleSet par
    population, un pour les femelles (ploidy=2), un pour les mâles
    (ploidy=1) -- population= est le nom msprime ("pop1", "pop2"...),
    PAS le nom réel du fichier .snp (même traduction que
    build_samples_argument, via population_index_to_name).

    IMPORTANT -- le ploidy PAR SampleSet ne contrôle QUE le nombre de
    lignées regroupées par individu dans le résultat, PAS le taux de
    coalescence : cette liste doit être utilisée avec
    simulate_independent_loci(..., ploidy=1) et une `demography` déjà
    rescalée via rescale_demography(demography,
    coalescence_coefficient("X", sex_ratio) / 2) -- vérifié
    empiriquement avec le mentor que c'est le ploidy GLOBAL de
    sim_ancestry qui interprète initial_size, pas celui des SampleSet.

    Args:
        snp_file_path: Chemin du fichier .snp.

    Returns:
        La liste des msprime.SampleSet (2 par population).

    Raises:
        ValueError: Si un individu a le sexe "9" (inconnu) -- ex:
            human_snp_all22chr_maf5.snp, où AUCUN individu n'est sexé
            (dataset <A>-only) : on ne peut pas construire un
            échantillonnage <X> dessus, mieux vaut le signaler
            explicitement que de produire un résultat silencieusement
            faux (individual_sexes_per_population laisse ce choix à
            l'appelant, c'est ici qu'il se prend).
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
    """Construit l'argument `samples` de msprime.sim_ancestry pour un locus <Y>.

    Le nom de population msprime ("pop1", "pop2"...) correspond à
    l'indice utilisé dans header.txt, mappé sur le nombre réel
    d'individus MÂLES observés pour la population correspondante (voir
    observed_data.py pour la justification du mapping par ordre
    d'apparition).

    PAS pour <M> : le mitochondrial est transmis uniquement par les
    mères, mais présent et échantillonné chez TOUS les individus
    (mâles et femelles), contrairement à <Y> qui n'existe que chez les
    mâles -- <M> doit réutiliser build_samples_argument (tout le monde)
    avec ploidy=1, pas cette fonction.

    Args:
        snp_file_path: Chemin du fichier .snp.

    Returns:
        Un dict {nom_population_msprime: nombre_d_individus_mâles}.

    Raises:
        ValueError: Si un individu a le sexe "9" (inconnu).
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


# ── Simulation des généalogies (arbres indépendants ou partagés) ───────────


def simulate_independent_loci(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    seed: int,
    ploidy: int = 2,
) -> Iterator[tskit.TreeSequence]:
    """Simule num_loci généalogies indépendantes sous la démographie donnée.

    Un locus SNP = un réplicat, pas de recombinaison interne ni de
    liaison entre loci. Retourne un itérateur (pas une liste) : pour
    51250 loci, matérialiser toutes les TreeSequence en mémoire
    simultanément serait coûteux -- l'appelant doit consommer cet
    itérateur au fil de l'eau (ex: pour calculer des statistiques
    résumées locus par locus).

    Args:
        demography: La démographie msprime.
        samples: dict[str, int] (un compte par population, ploidy
            uniforme -- <A>/<M>, voir build_samples_argument) ou
            list[msprime.SampleSet] (ploidy hétérogène par sous-groupe
            au sein d'une population -- <X>, voir
            build_sex_stratified_samples_argument). Les deux formes
            sont transmises telles quelles à msprime.sim_ancestry, qui
            les accepte indifféremment.
        num_loci: Le nombre de généalogies indépendantes à simuler.
        seed: La graine de la simulation.
        ploidy: 2 (défaut) pour <A>, cohérent avec une transmission
            diploïde classique -- chaque "sample individual" de
            `samples` compte pour 2 lignées génomiques. Pour <Y>/<M>,
            passer ploidy=1 avec une `demography` déjà rescalée par
            rescale_demography (voir demography_builder.py) : ces
            loci sont haploïdes (une seule copie de gène transmise),
            et le facteur de rescaling de Ne (coalescence_coefficient,
            observed_data.py) suppose cette combinaison ploidy=1 + Ne
            rescalé, pas ploidy=2 + Ne d'origine. Pour <X>, passer
            aussi ploidy=1 (voir
            build_sex_stratified_samples_argument : c'est le ploidy
            PAR SampleSet, pas ce paramètre global, qui donne 2
            copies aux femelles et 1 aux mâles -- ce paramètre-ci ne
            fixe que le taux de coalescence, via la Demography déjà
            rescalée).

    Returns:
        Un itérateur de num_loci TreeSequence indépendantes.
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
) -> Iterator[tskit.TreeSequence]:
    """Simule UNE SEULE généalogie puis la retourne répétée num_loci fois.

    Pour <Y>/<M>, dont tous les loci d'un même type partagent la même
    généalogie réelle (non-recombinants, transmission uniparentale),
    contrairement à <A>/<X> qui tirent un arbre indépendant par locus
    (simulate_independent_loci). Reproduit le comportement de
    particuleC.cpp:2422-2435 (GeneTreeY/GeneTreeM : premier locus <Y> ou
    <M> tire un arbre normalement, tous les suivants COPIENT ce même
    arbre -- seule la mutation change d'un locus à l'autre).

    IMPORTANT -- ne PAS réimplémenter le tirage de mutation ici :
    simulate_snp_genotypes(tree_sequences, seed) lit déjà tree_sequences
    au fil de l'eau sans jamais modifier les TreeSequence qu'elle reçoit,
    et son rng avance à chaque itération -- lui donner le MÊME objet
    TreeSequence répété num_loci fois (au lieu de num_loci objets
    différents) suffit à obtenir num_loci mutations indépendantes sur
    UNE SEULE généalogie, sans aucune modification de cette fonction
    (vérifié empiriquement : 5 répétitions du même arbre -> 5 génotypes
    différents).

    Args:
        demography: La démographie msprime.
        samples: Même contrat que simulate_independent_loci (voir sa
            docstring) -- cette fonction ne fait que réutiliser
            simulate_independent_loci avec num_loci=1, elle ne
            réinterprète pas ce paramètre.
        num_loci: Le nombre de fois où répéter la généalogie unique.
        seed: La graine de la simulation.
        ploidy: Pour <Y>/<M>, passer ploidy=1 (défaut) avec une `demography` déjà rescalée par
            rescale_demography (voir demography_builder.py) : ces
            loci sont haploïdes (une seule copie de gène transmise),
            et le facteur de rescaling de Ne (coalescence_coefficient,
            observed_data.py) suppose cette combinaison ploidy=1 + Ne
            rescalé, pas ploidy=2 + Ne d'origine.

    Returns:
        Un itérateur de num_loci TreeSequence, toutes identiques (le
        même objet Python répété).
    """

    shared_genealogy = next(
        simulate_independent_loci(
            demography, samples, num_loci=1, seed=seed, ploidy=ploidy
        )
    )
    return itertools.repeat(shared_genealogy, num_loci)


# ── Mutation (algorithme de Hudson) ────────────────────────────────────────


def _draw_single_mutation_edge_child(ts: tskit.TreeSequence, rng: random.Random) -> int:
    """Tire le noeud portant la mutation unique de l'algorithme de Hudson.

    Probabilité proportionnelle à la longueur de sa branche --
    entièrement vectorisé via les tables (pas d'appel branch_length()
    par noeud). Valable pour un arbre unique (sequence_length=1).

    Chaque edge = une branche (couple parent-enfant) ; edges.child liste
    donc tous les noeuds ayant une branche au-dessus d'eux (tous sauf la
    racine). Longueur = time[parent] - time[child], calculé en numpy.

    Validé empiriquement (proportions observées vs attendues <1% ; valeurs
    de statistiques identiques à la version par branch_length() -- voir
    notes/exploration.md).

    Args:
        ts: La TreeSequence (un seul arbre, sequence_length=1).
        rng: Le générateur aléatoire à utiliser.

    Returns:
        L'ID du noeud enfant portant la mutation.
    """
    edges = ts.tables.edges
    node_times = ts.tables.nodes.time

    children = edges.child  # array des noeuds enfants
    parents = edges.parent  # array des parents
    lengths = node_times[parents] - node_times[children]  # longueurs, vectorisé

    total = lengths.sum()
    target = rng.uniform(0, total)
    idx = np.searchsorted(
        np.cumsum(lengths), target
    )  # premier index où la somme cumulée dépasse target
    if idx >= len(edges.child):
        idx = len(edges.child) - 1
    return int(edges.child[idx])


def compute_population_layout(
    ts: tskit.TreeSequence,
) -> list[tuple[str | None, np.ndarray]]:
    """Calcule le layout (nom de population, IDs des noeuds échantillons) d'une TreeSequence.

    Factorisé pour pouvoir être calculé UNE SEULE FOIS et réutilisé sur
    plusieurs loci/tentatives qui partagent la même `demography`/
    `samples` d'origine -- seule la topologie coalescente varie d'un
    réplicat à l'autre, jamais l'assignation des noeuds échantillons aux
    populations (vérifié empiriquement, y compris entre réplicats tirés
    avec des graines différentes). Voir `simulate_snp_genotypes` (cache
    par défaut sur un flux de plusieurs loci) et les boucles de rejet MAF
    de `with_maf_filter`/`with_maf_filter_shared_ancestry` (cache
    explicite à travers les tentatives, voir notes/exploration.md,
    entrée du 20/07/2026).

    Args:
        ts: La TreeSequence à inspecter.

    Returns:
        La liste des (nom_population, IDs des noeuds échantillons de
        cette population), une entrée par population non vide.
    """
    layout = []
    for pop_index, population in enumerate(ts.tables.populations):
        sample_ids = ts.samples(population=pop_index)
        if len(sample_ids) == 0:
            continue
        pop_name = population.metadata.get("name") if population.metadata else None
        layout.append((pop_name, sample_ids))
    return layout


def simulate_snp_genotypes(
    tree_sequences: Iterator[tskit.TreeSequence],
    seed: int,
    population_layout: list[tuple[str | None, np.ndarray]] | None = None,
) -> Iterator[dict[str, list[int]]]:
    """Tire une mutation par locus (Hudson) et retourne les génotypes par population.

    Pour chaque TreeSequence (un locus = un arbre indépendant), tire
    une mutation UNIQUE selon l'algorithme de Hudson (vectorisé), et
    retourne les génotypes (0=ancestral, 1=dérivé) REGROUPÉS PAR
    POPULATION. Voir _draw_single_mutation_edge_child pour l'algorithme
    de tirage, et la docstring d'origine pour la justification du
    modèle (doc DIYABC section 2.4.3 : exactement une mutation par
    locus, locus toujours polymorphe).

    Args:
        tree_sequences: Un itérateur de TreeSequence, un arbre
            indépendant par locus.
        seed: La graine du tirage de mutation.
        population_layout: Voir `compute_population_layout`. Si
            `None` (cas d'un appel unique sur tout un flux de loci,
            ex: chemin `maf=0.0`), calculé UNE SEULE FOIS ici même, au
            premier locus, et réutilisé pour tous les suivants --
            valable car tous les `tree_sequences` d'un même appel
            partagent la même `demography`/`samples` d'origine (mêmes
            réplicats d'un seul appel à simulate_independent_loci/
            simulate_shared_ancestry_loci) : seule la topologie
            coalescente varie d'un locus à l'autre, jamais
            l'assignation des noeuds échantillons aux populations
            (vérifié empiriquement). Si fourni par l'appelant (ex:
            boucles de rejet MAF de `with_maf_filter`/`with_maf_
            filter_shared_ancestry`, qui appellent cette fonction une
            fois PAR TENTATIVE et calculent donc leur propre cache à
            travers les tentatives), utilisé tel quel sans jamais être
            recalculé. Sans ce cache, le redécodage du metadata des
            populations et le refiltrage de ts.samples(population=...)
            à CHAQUE locus représentaient à eux seuls ~20% du temps
            d'une particule sur 5000 loci (voir notes/exploration.md,
            entrée du 20/07/2026) -- le plus gros poste évitable du
            surcoût tskit par locus identifié dans cette investigation.

    Returns:
        Un itérateur de dicts {nom_population: [génotype, ...]} (un
        dict par locus).
    """
    rng = random.Random(seed)

    for ts in tree_sequences:
        tree = ts.first()
        mutated_node = _draw_single_mutation_edge_child(ts, rng)
        derived_samples = set(tree.samples(mutated_node))

        if population_layout is None:
            population_layout = compute_population_layout(ts)

        genotypes_by_population = {
            pop_name: [1 if s in derived_samples else 0 for s in sample_ids]
            for pop_name, sample_ids in population_layout
        }
        yield genotypes_by_population


# Passage des différents filtres pour les SNP


def observed_maf(locus_genotypes: dict[str, list[int]]) -> float:
    """Calcule la MAF poolée sur toutes les populations pour un locus.

    Comme ParticleC::mafreached (min(dérivé, ancestral) / total) -- pas
    juste la fréquence dérivée.

    Args:
        locus_genotypes: Dict {nom_population: [génotype, ...]} (0 ou
            1) pour un seul locus.

    Returns:
        La MAF (fréquence de l'allèle minoritaire).
    """
    all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
    n1 = sum(all_genotypes)
    n0 = len(all_genotypes) - n1
    return min(n0, n1) / len(all_genotypes)


def with_maf_filter(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    maf: float,
    seed: int,
    ploidy: int = 2,
) -> Iterator[dict[str, list[int]]]:
    """Simule des loci SNP indépendants avec filtre MAF.

    MAF = minor allele frequency, cf. doc DIYABC section 2.4.3 : si la
    fréquence de l'allèle MINORITAIRE (le moins fréquent des deux,
    dérivé ou ancestral -- pas forcément le dérivé) est strictement
    inférieure à `maf`, on rejette ce locus et on en resimule un
    nouveau (nouvelle généalogie + nouvelle mutation, jamais de
    recyclage de l'arbre rejeté) jusqu'à obtenir `num_loci` loci
    acceptés. Reproduit `ParticleC::mafreached`
    (particuleC.cpp:2194-2210).

    `maf=0.0` (équivalent DIYABC de `<MAF=hudson>` ou d'un tag absent)
    délègue directement à `simulate_independent_loci` +
    `simulate_snp_genotypes` avec la même graine pour les deux (comme
    le fait déjà chaque branche de `simulate_genotypes_for_locus_type`)
    -- comportement et résultats identiques à un appel direct de ces
    deux fonctions, pour ne rien changer aux datasets déjà validés qui
    n'ont pas de filtre MAF actif (human, toy_example5, ...).

    `maf>0.0` : les tentatives sont tirées PAR LOT de
    `max(_MAF_BATCH_SIZE, num_loci // 4)` (un seul appel
    `simulate_independent_loci(num_replicates=batch_size)` au lieu d'un
    appel par tentative individuelle) -- mesuré empiriquement ~5.6x plus
    rapide qu'un appel un-par-un sur toy_example3/scenario3/maf=0.05 (voir
    _MAF_BATCH_SIZE), le coalescent restant identique (nouvelle
    généalogie à chaque rejet, jamais de recyclage d'un arbre rejeté,
    comme avant). Le lot scale avec `num_loci` plutôt que d'être fixe :
    voir _MAF_BATCH_SIZE pour le détail. La structure population/
    échantillons (`population_layout`) ne dépend que de `demography`/
    `samples`, jamais de la graine tirée -- elle est donc calculée une
    seule fois, à la toute première tentative, et réutilisée pour toutes
    les suivantes (voir notes/exploration.md, entrée du 20/07/2026).

    Args:
        demography: La démographie msprime.
        samples: Même contrat que simulate_independent_loci.
        num_loci: Le nombre de loci acceptés à produire.
        maf: Le seuil MAF, déjà extrait (ex: via `parse_maf_ratio` sur
            le fichier .snp observé) -- cette fonction ne lit aucun
            fichier, à l'appelant de décider d'où vient le seuil.
        seed: La graine de la simulation.
        ploidy: Transmis tel quel à `simulate_independent_loci` (même
            contrat -- 2 pour <A>, 1 pour <H>/<X> avec une
            `demography` déjà rescalée, voir sa docstring).

    Returns:
        Un itérateur de `num_loci` dicts {nom_population:
        [génotype, ...]}, tous au-dessus du seuil MAF.
    """
    if maf == 0.0:
        tree_sequences = simulate_independent_loci(
            demography, samples, num_loci=num_loci, seed=seed, ploidy=ploidy
        )
        yield from simulate_snp_genotypes(tree_sequences, seed=seed)
        return

    # batch_size dérive num_loci -- voir _MAF_BATCH_SIZE pour la
    # justification empirique. La graine de lot (batch_seed = seed +
    # batch_index * batch_size) utilise directement batch_size comme
    # marge entre lots : contrairement à un stride fixe, elle ne peut
    # jamais collisionner avec les graines de mutation du lot précédent
    # (attempt_in_batch va de 0 à batch_size-1), quelle que soit la
    # taille du lot.
    batch_size = max(_MAF_BATCH_SIZE, num_loci // 4)
    accepted_loci = 0
    population_layout = None
    batch_index = 0
    while accepted_loci < num_loci:
        batch_seed = seed + batch_index * batch_size
        tree_sequences = simulate_independent_loci(
            demography,
            samples,
            num_loci=batch_size,
            seed=batch_seed,
            ploidy=ploidy,
        )
        for attempt_in_batch, ts in enumerate(tree_sequences):
            if population_layout is None:
                population_layout = compute_population_layout(ts)
            genotypes_by_population = next(
                simulate_snp_genotypes(
                    [ts],
                    seed=batch_seed + attempt_in_batch + _MAF_REJECTION_SEED_OFFSET,
                    population_layout=population_layout,
                )
            )

            maf_observed = observed_maf(genotypes_by_population)

            if maf_observed >= maf:
                yield genotypes_by_population
                accepted_loci += 1
                if accepted_loci >= num_loci:
                    return

        batch_index += 1


def with_maf_filter_shared_ancestry(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    maf: float,
    seed: int,
    ploidy: int = 1,
) -> Iterator[dict[str, list[int]]]:
    """Variante de with_maf_filter pour <Y>/<M> (généalogie partagée).

    Contrairement aux loci <A>/<H>/<X> (chaque locus = sa propre
    généalogie indépendante), tous les loci <Y> (resp. <M>) d'une même
    particule PARTAGENT UNE SEULE généalogie (voir
    simulate_shared_ancestry_loci) -- seule la mutation diffère d'un
    locus à l'autre.

    Reproduit exactement `particuleC.cpp:2424-2495` : le cache
    GeneTreeY/GeneTreeM est rempli AVANT le test MAF, donc indépendamment
    de son résultat -- la généalogie est tirée UNE SEULE FOIS (au tout
    premier appel), et un rejet MAF ne fait jamais redessiner l'arbre,
    seulement retirer une nouvelle mutation SUR CE MÊME ARBRE, jusqu'à
    obtenir `num_loci` loci acceptés. Voir aussi with_maf_filter (loci
    <A>/<H>/<X>), qui redessine au contraire une toute nouvelle
    généalogie à chaque rejet -- les deux mécanismes sont réellement
    différents côté DIYABC, pas juste une simplification.

    Args:
        demography: La démographie msprime (déjà rescalée par
            l'appelant si nécessaire).
        samples: Même contrat que simulate_independent_loci.
        num_loci: Le nombre de loci acceptés à produire.
        maf: Le seuil MAF, déjà extrait (voir with_maf_filter).
            `maf=0.0` délègue directement à
            simulate_shared_ancestry_loci + simulate_snp_genotypes
            avec la même graine pour les deux, comportement identique
            à un appel direct de ces deux fonctions.
        seed: La graine de la simulation.
        ploidy: Transmis tel quel à simulate_independent_loci/
            simulate_shared_ancestry_loci.

    Returns:
        Un itérateur de `num_loci` dicts {nom_population:
        [génotype, ...]}, tous au-dessus du seuil MAF.
    """
    if maf == 0.0:
        tree_sequences = simulate_shared_ancestry_loci(
            demography, samples, num_loci, seed, ploidy=ploidy
        )
        yield from simulate_snp_genotypes(tree_sequences, seed=seed)
        return

    shared_tree = next(
        simulate_independent_loci(
            demography, samples, num_loci=1, seed=seed, ploidy=ploidy
        )
    )
    # Calculée une seule fois : même généalogie PARTAGÉE à chaque
    # tentative, donc même structure population/échantillons -- voir
    # population_layout.
    population_layout = compute_population_layout(shared_tree)

    attempt = 0
    accepted_loci = 0
    while accepted_loci < num_loci:
        genotypes_by_population = next(
            simulate_snp_genotypes(
                [shared_tree],
                seed=seed + attempt + _MAF_REJECTION_SEED_OFFSET,
                population_layout=population_layout,
            )
        )

        maf_observed = observed_maf(genotypes_by_population)

        if maf_observed >= maf:
            yield genotypes_by_population
            accepted_loci += 1

        attempt += 1


# ── Dispatch par type de locus (compose tout ce qui précède) ──────────────


def simulate_genotypes_for_locus_type(
    demography: msprime.Demography,
    snp_file_path: str,
    locus_type: str,
    num_loci: int,
    seed: int,
) -> Iterator[dict[str, list[int]]]:
    """Point d'entrée unique de simulation de génotypes SNP, par type de locus.

    Choisit la bonne combinaison samples/demography-rescalée-ou-non/
    ploidy/fonction de simulation-indépendante-ou-partagée selon
    locus_type, puis retourne les génotypes simulés (même contrat de
    sortie que simulate_snp_genotypes, qu'on appelle en dernière étape
    dans tous les cas -- elle ne dépend jamais de locus_type
    elle-même).

    sex_ratio n'est PAS un paramètre : il est dérivé automatiquement de
    snp_file_path via parse_sex_ratio, comme tout le reste (samples,
    sexes par individu) -- l'appelant n'a jamais besoin de le connaître.

    Dispatch :
      - "A" : build_samples_argument, demography TELLE QUELLE (pas de
        rescale_demography), ploidy=2, with_maf_filter.
      - "H" : build_samples_argument, demography rescalée par
        coalescence_coefficient("H", sex_ratio) / 2, ploidy=1,
        with_maf_filter.
      - "X" : build_sex_stratified_samples_argument, demography
        rescalée par coalescence_coefficient("X", sex_ratio) / 2,
        ploidy=1, with_maf_filter.
      - "Y" : build_male_only_samples_argument, demography rescalée par
        coalescence_coefficient("Y", sex_ratio) / 2, ploidy=1,
        with_maf_filter_shared_ancestry (arbre unique partagé).
      - "M" : build_samples_argument (TOUT le monde, pas mâles seuls --
        voir la docstring de build_male_only_samples_argument sur ce
        point précis), demography rescalée par
        coalescence_coefficient("M", sex_ratio) / 2, ploidy=1,
        with_maf_filter_shared_ancestry.
      - tout autre locus_type : lever NotImplementedError (même style
        que coalescence_coefficient pour un type inconnu).

    MAF : le seuil est lu une fois via parse_maf_ratio(snp_file_path) et
    délégué à with_maf_filter ("A"/"H"/"X", généalogie indépendante par
    locus) ou with_maf_filter_shared_ancestry ("Y"/"M", généalogie
    partagée) -- les deux gèrent elles-mêmes le cas maf=0.0 (pas de
    filtre, comportement identique à un appel direct des fonctions
    sous-jacentes) et le cas maf>0.0 (boucle de rejet).

    "Y"/"M" (MAF quelconque) utilisent with_maf_filter_shared_ancestry,
    pas with_maf_filter : ces deux types partagent UNE SEULE généalogie
    entre tous leurs loci (simulate_shared_ancestry_loci) -- un rejet MAF
    ne redessine jamais l'arbre, seulement la mutation (voir la docstring
    de with_maf_filter_shared_ancestry, qui reproduit exactement
    particuleC.cpp:2424-2495).

    Args:
        demography: La démographie <A> "de base" (construite par
            build_demography, PAS encore rescalée) -- c'est CETTE
            fonction qui décide si/comment la rescaler selon
            locus_type, jamais l'appelant.
        snp_file_path: Chemin du fichier .snp observé.
        locus_type: "A", "H", "X", "Y" ou "M".
        num_loci: Le nombre de loci à simuler.
        seed: La graine de la simulation.

    Returns:
        Un itérateur de `num_loci` dicts {nom_population:
        [génotype, ...]}.

    Raises:
        NotImplementedError: Si locus_type est inconnu.
    """

    sex_ratio = parse_sex_ratio(snp_file_path)
    maf_ratio = parse_maf_ratio(snp_file_path)

    if locus_type == "Y":
        samples = build_male_only_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        return with_maf_filter_shared_ancestry(
            rescaled_demography, samples, num_loci, maf_ratio, seed, ploidy=1
        )
    elif locus_type == "M":
        samples = build_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        return with_maf_filter_shared_ancestry(
            rescaled_demography, samples, num_loci, maf_ratio, seed, ploidy=1
        )
    elif locus_type == "A":
        samples = build_samples_argument(snp_file_path)
        return with_maf_filter(demography, samples, num_loci, maf_ratio, seed, ploidy=2)
    elif locus_type == "H":
        samples = build_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        return with_maf_filter(
            rescaled_demography, samples, num_loci, maf_ratio, seed, ploidy=1
        )
    elif locus_type == "X":
        samples = build_sex_stratified_samples_argument(snp_file_path)
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(locus_type, sex_ratio) / 2
        )
        return with_maf_filter(
            rescaled_demography, samples, num_loci, maf_ratio, seed, ploidy=1
        )
    else:
        raise NotImplementedError(f"Type de locus non supporté: {locus_type!r}")


# -----------------------------------------------
# Simulation des lectures PoolSeq
# -----------------------------------------------


def simulate_poolseq_reads(
    tree_sequences: Iterator[tskit.TreeSequence],
    observed_reads_per_locus: list[dict[str, tuple[int, int]]],
    seed: int,
    population_layout: list[tuple[str | None, np.ndarray]] | None = None,
) -> Iterator[dict[str, tuple[int, int]]]:
    """Simule les lectures PoolSeq de chaque locus.

    Tire une mutation par locus (même algorithme de Hudson que
    simulate_snp_genotypes), puis convertit la proportion de lignées
    dérivées de chaque population en un tirage binomial de lectures,
    calé sur la profondeur totale RÉELLEMENT observée à ce locus/cette
    population (`observed_reads_per_locus`) -- seule la répartition
    allèle1/allèle2 est simulée, jamais la couverture elle-même.

    Args:
        tree_sequences: Un itérateur de TreeSequence simulées, une par
            locus.
        observed_reads_per_locus: Une liste de dicts {nom_population:
            (nreads_dérivé, nreads_total)} observés, un par locus.
        seed: La graine du tirage de mutation.
        population_layout: Voir `compute_population_layout`. Si
            `None` (cas d'un appel unique sur tout un flux de loci,
            ex: chemin `mrc<=0`), calculé UNE SEULE FOIS ici même, au
            premier locus, et réutilisé pour tous les suivants --
            même principe que `simulate_snp_genotypes`. Si fourni par
            l'appelant (ex: boucle de rejet MRC de `with_mrc_filter`,
            qui appelle cette fonction une fois PAR TENTATIVE et
            calcule donc son propre cache à travers les tentatives),
            utilisé tel quel sans jamais être recalculé.

    Returns:
        Un itérateur de dicts {nom_population: (nreads_dérivé,
        nreads_total)} simulés, un par locus.
    """

    rng = random.Random(seed)
    binom_rng = np.random.default_rng(
        seed + _BINOMIAL_SEED_OFFSET
    )  # graine aléatoire séparée pour le tirage binomial, pour ne pas interférer avec le tirage de mutation

    for ts, reads_observed in zip(
        tree_sequences, observed_reads_per_locus, strict=False
    ):
        tree = ts.first()
        mutated_node = _draw_single_mutation_edge_child(ts, rng)
        # set(...) impératif ICI (pas dans la boucle plus bas) : tree.samples()
        # renvoie un générateur, épuisé après la 1ère population -- sans ce
        # set() immédiat, pop_derived_count tombait silencieusement à 0 pour
        # TOUTES les populations sauf la première de population_layout
        # (confirmé empiriquement le 22/07/2026 : seule pop1 montrait jamais
        # de variation dans tout un reftable simulé). simulate_snp_genotypes
        # fait déjà ce set() immédiat, c'est le bon modèle à suivre.
        derived_samples = set(tree.samples(mutated_node))
        if population_layout is None:
            population_layout = compute_population_layout(ts)

        reads_by_population = {}
        for pop_name, sample_ids in population_layout:
            total_reads = reads_observed[pop_name][1]
            pop_derived_count = len(derived_samples.intersection(sample_ids))
            p = pop_derived_count / len(sample_ids) if len(sample_ids) > 0 else 0.0
            if total_reads > 0:
                derived_reads = binom_rng.binomial(total_reads, p)
                reads_by_population[pop_name] = (derived_reads, total_reads)
            else:
                reads_by_population[pop_name] = (0, 0)
        yield reads_by_population


def _reindex_reads_by_msprime_name(
    observed_reads_per_locus: list[dict[str, tuple[int, int]]],
    snp_file_path: str,
) -> list[dict[str, tuple[int, int]]]:
    """Reindexe les lectures observées pour utiliser les noms de population msprime.

    Args:
        observed_reads_per_locus: Une liste de dicts {nom_population
            réel: (nreads_dérivé, nreads_total)}, un par locus.
        snp_file_path: Chemin du fichier .snp, pour obtenir la
            correspondance des noms de population.

    Returns:
        La même liste, avec les clés remplacées par les noms de
        population msprime ("pop1", "pop2"...).
    """
    index_to_name = population_index_to_name(
        snp_file_path
    )  # {1: "POP1", 2: "POP2", 3: "POP3", 4: "POP4"} pour toy_example4.
    real_name_to_msprime_name = {
        name: f"pop{index}" for index, name in index_to_name.items()
    }
    return [
        {
            real_name_to_msprime_name[pop_name]: reads
            for pop_name, reads in locus_reads.items()
        }
        for locus_reads in observed_reads_per_locus
    ]


def with_mrc_filter(
    demography: msprime.Demography,
    samples: dict[str, int] | list[msprime.SampleSet],
    num_loci: int,
    mrc: float,
    observed_reads_per_locus: list[dict[str, tuple[int, int]]],
    seed: int,
    ploidy: int = 2,
) -> Iterator[dict[str, tuple[int, int]]]:
    """Simule des loci SNP indépendants avec filtre MRC.

    MRC = minimum read count. Si le nombre de lectures dérivées est
    strictement inférieur à `mrc`, on rejette ce locus et on en
    resimule un nouveau. Reproduit le comportement de
    `ParticleC::mrc_reached`.

    `mrc>0` : les tentatives sont tirées depuis un POOL PARTAGÉ ENTRE TOUS
    LES LOCI, pas un pool privé par locus -- contrairement à un batching
    naïf "par locus" (un nouveau lot de `_MRC_BATCH_SIZE` généalogies à
    chaque `locus_index`, même si ce locus n'a besoin que d'UNE seule
    tentative), qui paie un plancher d'un appel `simulate_independent_loci`
    PAR LOCUS quel que soit le taux d'acceptation. Ici, un seul flux
    continu de généalogies (régénéré par lot de `_MRC_BATCH_SIZE`
    uniquement quand épuisé) est consommé par n'importe quel locus qui a
    besoin d'une nouvelle tentative -- si la plupart des loci passent dès
    le premier tirage (cas courant), un seul lot peut servir des dizaines
    de loci au lieu d'un lot par locus. Gain mesuré empiriquement (script
    jetable, toy_example4, mrc=5) : ~1.5x supplémentaire par rapport au
    batching par-locus, quelle que soit la taille du lot (le partage
    compte, pas la taille).

    Le compteur `attempt` est GLOBAL et n'est jamais remis à zéro par
    locus -- ça élimine par construction le risque de corrélation qui
    existait avec l'ancien design par-locus (deux loci ayant besoin du
    même nombre de tentatives tiraient alors le même arbre/mutation,
    confirmé empiriquement le 22/07/2026) : chaque tentative, tous loci
    confondus, consomme une position distincte dans un flux continu,
    jamais réutilisée.

    Args:
        demography: La démographie msprime.
        samples: Même contrat que simulate_independent_loci.
        num_loci: Le nombre de loci acceptés à produire.
        mrc: Le seuil MRC, déjà extrait via `parse_mrc_ratio` sur le
            fichier .snp observé.
        observed_reads_per_locus: Une liste de dicts {nom_population:
            (nreads_dérivé, nreads_total)} observés, un par locus.
        seed: La graine de la simulation.
        ploidy: Transmis tel quel à `simulate_independent_loci`.

    Returns:
        Un itérateur de `num_loci` dicts {nom_population:
        (nreads_dérivé, nreads_total)}, tous au-dessus du seuil MRC.
    """

    if mrc <= 0:
        tree_sequences = simulate_independent_loci(
            demography, samples, num_loci=num_loci, seed=seed, ploidy=ploidy
        )
        yield from simulate_poolseq_reads(
            tree_sequences, observed_reads_per_locus, seed=seed
        )  # liste de dictionnaires contenant le nombre de lectures dérivées et ancestrales par population pour chaque locus
        return
    # Calculée une seule fois, à la première tentative, et réutilisée pour
    # toutes les suivantes (tous les loci/tentatives partagent la même
    # demography/samples) -- voir population_layout et with_maf_filter
    # (même principe côté IndSeq).
    population_layout = None
    # tree_sequences_iter : itérateur du lot COURANT de _MRC_BATCH_SIZE
    # généalogies, PARTAGÉ entre tous les locus_index -- régénéré
    # (nouveau lot, nouvelle graine) uniquement quand épuisé, jamais
    # réinitialisé au passage à un nouveau locus (voir docstring).
    tree_sequences_iter = None
    batch_index = 0
    attempt = 0
    for locus_index in range(num_loci):
        while True:
            if tree_sequences_iter is None:
                batch_seed = seed + batch_index * _MRC_BATCH_SIZE
                tree_sequences_iter = simulate_independent_loci(
                    demography,
                    samples,
                    num_loci=_MRC_BATCH_SIZE,
                    seed=batch_seed,
                    ploidy=ploidy,
                )
                batch_index += 1
            ts = next(tree_sequences_iter, None)
            if ts is None:
                tree_sequences_iter = None
                continue
            if population_layout is None:
                population_layout = compute_population_layout(ts)
            reads_by_population = next(
                simulate_poolseq_reads(
                    [ts],
                    observed_reads_per_locus[locus_index : locus_index + 1],
                    seed=seed + attempt + _MRC_REJECTION_SEED_OFFSET,
                    population_layout=population_layout,
                )
            )
            attempt += 1

            # calcul du mrc observé
            mrc_observed = observed_mrc(reads_by_population)

            if mrc_observed >= mrc:
                yield reads_by_population
                break


def prepare_poolseq_observed_reads(
    snp_file_path: str, num_loci: int
) -> list[dict[str, tuple[int, int]]]:
    """Prépare les lectures observées pour la simulation PoolSeq.

    Lit le fichier .snp et tronque aux `num_loci` premières entrées.

    Args:
        snp_file_path: Chemin du fichier .snp (doit être POOLSEQ).
        num_loci: Le nombre de loci à conserver.

    Returns:
        Une liste de dicts {nom_population_msprime: (nreads_dérivé,
        nreads_total)}, un par locus.
    """
    raw_reads = observed_reads(snp_file_path, num_loci=num_loci)
    reindexed_reads = _reindex_reads_by_msprime_name(raw_reads, snp_file_path)
    return reindexed_reads


def simulate_poolseq_reads_with_mrc_filter(
    demography: msprime.Demography,
    snp_file_path: str,
    seed: int,
    num_loci: int,
    observed_reads_per_locus: list[dict[str, tuple[int, int]]] = None,
) -> Iterator[dict[str, tuple[int, int]]]:
    """Point d'entrée unique de simulation de lectures pour un fichier PoolSeq.

    Pendant de simulate_genotypes_for_locus_type (IndSeq), mais sans
    dispatch multi-type : un fichier PoolSeq n'a jamais qu'un seul type
    de locus déclaré (`<A>`, cf. `data.cpp:529` -- seule la classe de
    locus autosomale diploïde est supportée pour PoolSeq côté DIYABC).
    Aucun rescale de Ne n'est nécessaire (même `coeffcoal` que l'IndSeq
    `<A>` standard, cf. `data.cpp:1589-1603` -- PoolSeq a `type=15`,
    `15 % 5 == 0`, donc tombe dans exactement la même branche que le
    cas autosomal diploïde standard).

    ploidy=2, PAS ploidy=1 (corrigé le 22/07/2026 -- voir
    notes/exploration.md) : DIYABC construit l'arbre de généalogie
    PoolSeq en réutilisant tel quel le chemin `<A>` standard --
    `HAPLOID_SAMPLE_SIZE` (déclaré dans le fichier .snp, `POOL
    pop:N`) est le nombre de COPIES DE GÈNES à échantillonner, mais le
    Ne (N1/N2/...) reste interprété en individus DIPLOÏDES, exactement
    comme `<A>` IndSeq (`particuleC.cpp:1185`, `data.cpp:970-974`,
    confirmé par exploration du code source réel -- aucune branche de
    simulation généalogique spécifique à PoolSeq n'existe, seule la
    lecture du fichier et le tirage des reads le sont). Passer
    `samples=build_samples_argument(...)` (comptes haploïdes) avec
    `ploidy=1` -- ce qui a été fait par erreur avant cette correction --
    revient à traiter le Ne comme s'il était HAPLOÏDE : la coalescence
    devient deux fois plus rapide en unités de générations réelles pour
    le MÊME Ne déclaré, ce qui gonfle artificiellement toute
    différenciation entre populations (FST/F3/F4 ~60-140% trop élevés,
    confirmé empiriquement sur toy_example4 en comparaison appariée
    contre un vrai reftable DIYABC) sans affecter les statistiques
    mono-population (HW/ML1, qui ne dépendent pas de la vitesse relative
    de coalescence entre populations). D'où le //2 ci-dessous : on donne
    à msprime un compte d'INDIVIDUS (`ploidy=2` double automatiquement
    en lignées), pas un compte de lignées déjà doublé -- le nombre total
    de lignées échantillonnées (donc la taille de l'arbre) reste
    identique (`HAPLOID_SAMPLE_SIZE`), seule la vitesse de coalescence
    relative au Ne change, pour retomber sur la même convention que
    l'IndSeq `<A>`.

    Compose, dans l'ordre :
      - `parse_mrc_ratio(snp_file_path)` -- seuil MRC (défaut 1 si
        `<MRC=...>` absent, PAS 0 comme pour MAF -- voir
        `parse_mrc_ratio`).
      - `build_samples_argument(snp_file_path)` -- retourne la taille
        HAPLOÏDE du pool par population (cf.
        `count_samples_per_population`/`_parse_pool_header_line`) --
        divisée par 2 ici pour obtenir un compte d'INDIVIDUS diploïdes
        (voir ci-dessus) ; utilisée TELLE QUELLE (non divisée) partout
        ailleurs, notamment comme `pool_sizes` dans
        `summary_statistics.py` (la correction de biais de lecture Q1
        a besoin du vrai `HAPLOID_SAMPLE_SIZE`, pas de sa moitié).
      - `observed_reads(snp_file_path)` -- les lectures RÉELLEMENT
        observées par locus/population, ensuite retraduites vers les
        noms de population msprime (`"pop1"`, `"pop2"`...) via
        `_reindex_reads_by_msprime_name` (les noms réels du fichier .snp
        n'ont aucune raison de coïncider avec cette convention -- voir
        son docstring). Tronquées aux `num_loci` premières entrées :
        c'est CETTE couverture réelle, fixe par emplacement de locus,
        qui sert de paramètre `n` au tirage binomial dans
        `simulate_poolseq_reads`, jamais retirée au hasard (voir
        `with_mrc_filter`/`simulate_poolseq_reads`).
      - `with_mrc_filter(..., ploidy=2)` -- simulation + rejet-et-
        resimule si le critère MRC (min des reads dérivés/ancestraux,
        toutes populations combinées) n'est pas atteint.

    Args:
        demography: La démographie de base (PAS encore rescalée --
            comme pour `<A>` en IndSeq, aucun rescale n'est
            nécessaire ici).
        snp_file_path: Chemin du fichier .snp (doit être POOLSEQ).
        seed: La graine de la simulation.
        num_loci: Le nombre de loci à simuler.
        observed_reads_per_locus: Si `None` (défaut), calculé via
            `prepare_poolseq_observed_reads(snp_file_path, num_loci)`.
            Sinon, utilisé tel quel.

    Returns:
        Un itérateur de `num_loci` dicts {nom_population:
        (nreads_dérivé, nreads_total)}, tous au-dessus du seuil MRC.
    """
    mrc = parse_mrc_ratio(snp_file_path)
    haploid_pool_sizes = build_samples_argument(snp_file_path)
    samples = {pop: count // 2 for pop, count in haploid_pool_sizes.items()}
    if observed_reads_per_locus is None:
        observed_reads_per_locus = prepare_poolseq_observed_reads(
            snp_file_path, num_loci
        )

    return with_mrc_filter(
        demography, samples, num_loci, mrc, observed_reads_per_locus, seed, ploidy=2
    )


# -------------------------------------------------------
# Séquence ADN
# -------------------------------------------------------


# Simulation des mutations pour les séquences ADN
def build_transition_matrix(
    name_model: str, kappas: tuple[float, float], frequences_by_locus: dict[str, float]
) -> np.ndarray:
    """Calcule la matrice de transition (matQ) pour un modèle donné et un locus.

    Args:
        name_model: "JK", "K2P", "HKY" ou "TN".
        kappas: (k1, k2) -- k2 ignoré pour "JK"/"K2P"/"HKY".
        frequences_by_locus: Dict {"pi_A": ..., "pi_C": ..., "pi_G":
            ..., "pi_T": ...}, les fréquences de bases observées de ce
            locus (voir observed_data.base_frequency_by_locus).

    Returns:
        La matrice 4x4 (ordre A/C/G/T) des probabilités de
        transition, chaque ligne normalisée à somme 1, diagonale nulle.

    Raises:
        NotImplementedError: Si name_model est inconnu.
    """
    pi_A, pi_C, pi_G, pi_T = (
        frequences_by_locus["pi_A"],
        frequences_by_locus["pi_C"],
        frequences_by_locus["pi_G"],
        frequences_by_locus["pi_T"],
    )
    k1, k2 = kappas[0], kappas[1]
    # Initialisation de la matrice de transition
    transition_matrix = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if i == j:
                transition_matrix[i, j] = 0
            else:
                transition_matrix[i, j] = 1
    if name_model == "JK":
        pass
    elif name_model == "K2P":
        (
            transition_matrix[0, 2],
            transition_matrix[1, 3],
            transition_matrix[2, 0],
            transition_matrix[3, 1],
        ) = k1, k1, k1, k1
    elif name_model == "HKY":
        # ligne1
        transition_matrix[0, 1] = pi_C
        transition_matrix[0, 2] = k1 * pi_G
        transition_matrix[0, 3] = pi_T
        # ligne2
        transition_matrix[1, 0] = pi_A
        transition_matrix[1, 2] = pi_G
        transition_matrix[1, 3] = k1 * pi_T
        # ligne3
        transition_matrix[2, 0] = k1 * pi_A
        transition_matrix[2, 1] = pi_C
        transition_matrix[2, 3] = pi_T
        # ligne4
        transition_matrix[3, 0] = pi_A
        transition_matrix[3, 1] = k1 * pi_C
        transition_matrix[3, 2] = pi_G
    elif name_model == "TN":
        # ligne1
        transition_matrix[0, 1] = pi_C
        transition_matrix[0, 2] = k1 * pi_G
        transition_matrix[0, 3] = pi_T
        # ligne2
        transition_matrix[1, 0] = pi_A
        transition_matrix[1, 2] = pi_G
        transition_matrix[1, 3] = k2 * pi_T
        # ligne3
        transition_matrix[2, 0] = k1 * pi_A
        transition_matrix[2, 1] = pi_C
        transition_matrix[2, 3] = pi_T
        # ligne4
        transition_matrix[3, 0] = pi_A
        transition_matrix[3, 1] = k2 * pi_C
        transition_matrix[3, 2] = pi_G
    else:
        raise NotImplementedError(
            f"Modèle de substitution non supporté: {name_model!r}"
        )

    # Normalisation de la matrice de transition
    transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
    return transition_matrix


def build_rate_map(
    mutsit: list[float], mus_rate: float, dnalength: int
) -> msprime.RateMap:
    """Construit le profil de taux de mutation par site d'un locus (msprime.RateMap).

    `mus_rate` est un taux moyen PAR SITE (donc le taux total attendu
    sur tout le locus est `mus_rate * dnalength`) ; `mutsit` répartit
    ce budget total entre les sites (poids relatifs qui somment à 1,
    voir parameter_sampling.sample_site_rates -- un site invariant a un
    poids de 0). Le taux absolu d'un site donné est donc `mus_rate *
    dnalength * mutsit[site]` : on distribue le taux total du locus
    site par site, proportionnellement à son poids relatif.

    Args:
        mutsit: Le poids de mutation relatif par site, longueur
            dnalength, normalisé à somme 1.
        mus_rate: Le taux de mutation moyen PAR SITE du locus.
        dnalength: La longueur du locus.

    Returns:
        La msprime.RateMap correspondante (un taux absolu par site).

    Raises:
        ValueError: Si len(mutsit) != dnalength.
    """
    if len(mutsit) != dnalength:
        raise ValueError(
            f"Le nombre de sites de mutation ({len(mutsit)}) ne correspond pas à la longueur de la séquence ({dnalength})."
        )
    # Crée le profil de taux (msprime.RateMap), un taux absolu par site
    rate_map = msprime.RateMap(
        position=[i for i in range(dnalength + 1)],
        rate=[mus_rate * dnalength * mutsit[i] for i in range(dnalength)],
    )
    return rate_map


def count_loci_per_group(list_loci: list[LociDescriptionDetailed]) -> dict[str, int]:
    """Compte le nombre de loci par groupe.

    Args:
        list_loci: La liste des loci détaillés.

    Returns:
        Un dict {nom_groupe: nombre_de_loci}.

    Raises:
        ValueError: Si un même groupe mélange plusieurs types de loci
            ("M" et "S"), non supporté.
    """
    loci_count = {}
    loci_type_per_group = {}
    for locus in list_loci:
        group_name = locus.group
        loci_type = locus.ms_or_seq
        if group_name not in loci_count:
            loci_count[group_name] = 0
            loci_type_per_group[group_name] = set()
        loci_count[group_name] += 1
        loci_type_per_group[group_name].add(loci_type)

    if any(len(types) > 1 for types in loci_type_per_group.values()):
        raise ValueError(
            "Différents types de loci dans le même groupe, ce qui n'est pas supporté."
        )

    return loci_count


def build_group_local_param_per_locus(
    header_text: str, seed: int
) -> dict[str, tuple[float, float, float]]:
    """Tire k1/k2/mus_rate par locus (hiérarchie à deux niveaux, groupe puis locus).

    Pour chaque groupe `[S]` de header.txt : `draw_group_parameter_values`
    donne la valeur moyenne par groupe (premier niveau), puis
    `sampling_group_local_param` en dérive une valeur par locus (second
    niveau, dispersion optionnelle autour de la moyenne du groupe).

    Args:
        header_text: Texte complet de header.txt.
        seed: La graine du tirage.

    Returns:
        Un dict {nom_locus: (k1, k2, mus_rate)} -- triplet à arité
        fixe pour chaque locus `[S]`, quel que soit le modèle de
        substitution utilisé par son groupe (0.0 pour les kappas non
        utilisés).
    """
    params_per_locus = {}
    group_priors = parse_group_priors(header_text)
    list_loci = parse_loci_description(header_text)

    list_loci_seq = [locus for locus in list_loci if locus.ms_or_seq == "S"]
    nloc_per_group = count_loci_per_group(list_loci_seq)

    # Un seul appel pour tous les groupes -- draw_group_parameter_values gère
    # déjà en interne son propre décalage de graine (_GROUP_PRIOR_SEED_OFFSET),
    # distinct de _KAPPA1_SEED_OFFSET/_KAPPA2_SEED_OFFSET utilisés plus bas
    # pour le tirage par locus -- pas besoin (et pas correct) de la rappeler
    # une fois par groupe/par kappa avec une graine décalée différente.
    values = draw_group_parameter_values(group_priors, seed)

    # Une seule instance de RNG par paramètre, créée AVANT la boucle sur les
    # groupes et réutilisée/avancée à chaque groupe -- même motif que le rng
    # de build_rate_map_per_locus (déjà correct). Créer un random.Random(seed
    # + OFFSET) FRAIS À CHAQUE ITÉRATION du groupe (comme avant ce correctif)
    # fait que tous les groupes utilisant le même paramètre (ex: G2 et G3,
    # tous deux K2P donc tous deux kappa1) rejouent EXACTEMENT la même
    # séquence de tirages -- une corrélation inter-groupes invisible tant
    # qu'un seul groupe de chaque modèle existe, révélée seulement en
    # explorant la piste <M> avec toy_example2_ms_dna_50loci (G2/G3 tous
    # deux K2P). Voir feedback_seed_reuse_pattern (même classe de bug déjà
    # vue 4 fois dans ce projet).
    mus_rate_rng = random.Random(seed + _MUS_RATE_SEED_OFFSET)
    kappa1_rng = random.Random(seed + _KAPPA1_SEED_OFFSET)
    kappa2_rng = random.Random(seed + _KAPPA2_SEED_OFFSET)

    for group in nloc_per_group:
        list_locus_in_group = [locus for locus in list_loci_seq if locus.group == group]
        # calcul du mus_rate par locus
        mus_rate = sampling_group_local_param(
            next(gp for gp in group_priors[group] if gp.name == "GAMMU"),
            k_moy=values[group]["MEANMU"],
            n_loci=nloc_per_group[group],
            check_nloc=True,
            list_loci=list_locus_in_group,
            rng=mus_rate_rng,
        )
        # calcul des kappas par locus selon le modèle de substitution du groupe
        gp_model = next(gp for gp in group_priors[group] if gp.model)
        verif_kappas = get_parameter_used_by_model(gp_model)
        if not verif_kappas[0] and not verif_kappas[1]:
            for locus in list_locus_in_group:
                # Handle the case where neither k1 nor k2 is used by the model
                params_per_locus[locus.name] = (0.0, 0.0, mus_rate[locus.name])
        elif verif_kappas[0] and not verif_kappas[1]:  # Model K2P et HKY
            gp_gamk1 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK1"), None
            )
            if gp_gamk1 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa1, mais aucun GAMK1 n'a été trouvé pour le groupe {group}."
                )
            k1_moy = values[group]["MEANK1"]
            kappa1_values = sampling_group_local_param(
                gp_gamk1,
                k_moy=k1_moy,
                n_loci=nloc_per_group[group],
                check_nloc=True,
                list_loci=list_locus_in_group,
                rng=kappa1_rng,
            )
            for locus in list_locus_in_group:
                params_per_locus[locus.name] = (
                    kappa1_values[locus.name],
                    0.0,
                    mus_rate[locus.name],
                )
        elif verif_kappas[0] and verif_kappas[1]:  # Modèle TN
            gp_gamk1 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK1"), None
            )
            if gp_gamk1 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa1, mais aucun GAMK1 n'a été trouvé pour le groupe {group}."
                )
            k1_moy = values[group]["MEANK1"]
            gp_gamk2 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK2"), None
            )
            if gp_gamk2 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa2, mais aucun GAMK2 n'a été trouvé pour le groupe {group}."
                )
            k2_moy = values[group]["MEANK2"]
            kappa1_values = sampling_group_local_param(
                gp_gamk1,
                k_moy=k1_moy,
                n_loci=nloc_per_group[group],
                check_nloc=True,
                list_loci=list_locus_in_group,
                rng=kappa1_rng,
            )
            kappa2_values = sampling_group_local_param(
                gp_gamk2,
                k_moy=k2_moy,
                n_loci=nloc_per_group[group],
                check_nloc=False,
                list_loci=list_locus_in_group,
                rng=kappa2_rng,
            )
            for locus in list_locus_in_group:
                params_per_locus[locus.name] = (
                    kappa1_values[locus.name],
                    kappa2_values[locus.name],
                    mus_rate[locus.name],
                )
    return params_per_locus


def build_matrix_per_locus(
    header_text: str, mss_file_path: str | Path, seed: int
) -> dict[str, np.ndarray]:
    """Construit la matrice de transition (matQ) de chaque locus [S].

    Pipeline complet header.txt + .mss + seed -> {nom_locus: matQ},
    en composant build_group_local_param_per_locus (k1/k2 par locus),
    base_frequency_by_locus (pi par locus) et build_transition_matrix.

    Args:
        header_text: Texte complet de header.txt.
        mss_file_path: Chemin du fichier .mss.
        seed: La graine du tirage.

    Returns:
        Un dict {nom_locus: matQ} pour chaque locus [S].
    """
    list_loci = parse_loci_description(header_text)
    params_per_locus = build_group_local_param_per_locus(header_text, seed)
    sequences_by_indiv = observed_sequences(mss_file_path, list_loci)
    frequencies_by_locus = base_frequency_by_locus(sequences_by_indiv)
    group_priors = parse_group_priors(header_text)

    matrix_per_locus = {}
    for locus in list_loci:
        if locus.ms_or_seq == "S":
            group_locus = locus.group
            gp_model = next(gp for gp in group_priors[group_locus] if gp.model)
            name_model = gp_model.name_model
            kappas = params_per_locus[locus.name][0], params_per_locus[locus.name][1]
            frequencies = frequencies_by_locus[locus.name]
            matrix_per_locus[locus.name] = build_transition_matrix(
                name_model, kappas, frequencies
            )
    return matrix_per_locus


def build_rate_map_per_locus(header_text: str, seed: int) -> dict[str, msprime.RateMap]:
    """Construit le profil de taux de mutation (msprime.RateMap) de chaque locus [S].

    Args:
        header_text: Texte complet de header.txt.
        seed: La graine du tirage.

    Returns:
        Un dict {nom_locus: msprime.RateMap} pour chaque locus [S].
    """
    list_loci = parse_loci_description(header_text)
    params_per_locus = build_group_local_param_per_locus(header_text, seed)
    mus_rate_per_locus = {
        locus.name: params_per_locus[locus.name][2]
        for locus in list_loci
        if locus.ms_or_seq == "S"
    }
    group_priors = parse_group_priors(header_text)

    rate_map_per_locus = {}
    rng = random.Random(seed + _SITE_RATE_SEED_OFFSET)
    for locus in list_loci:
        if locus.ms_or_seq == "S":
            gp_model = next(gp for gp in group_priors[locus.group] if gp.model)
            mutsit = sample_site_rates(
                gp_model.p_fixe, gp_model.gams, locus.dnalength, rng=rng
            )
            mus_rate = mus_rate_per_locus[locus.name]
            rate_map_per_locus[locus.name] = build_rate_map(
                mus_rate=mus_rate, mutsit=mutsit, dnalength=locus.dnalength
            )
    return rate_map_per_locus


def simulate_dna_mutations(
    tree_sequence: tskit.TreeSequence,
    transition_matrix: np.ndarray,
    frequencies: dict[str, float],
    rate_map: msprime.RateMap,
    seed: int,
) -> tskit.TreeSequence:
    """Place les mutations sur une généalogie ADN via msprime.sim_mutations.

    Args:
        tree_sequence: La généalogie non mutée (un locus).
        transition_matrix: La matrice de transition (matQ) du locus.
        frequencies: Dict {"pi_A": ..., "pi_C": ..., "pi_G": ...,
            "pi_T": ...}, distribution ancestrale à la racine.
        rate_map: Le profil de taux de mutation par site du locus.
        seed: La graine du tirage de mutation.

    Returns:
        La TreeSequence mutée.
    """
    alleles = ["A", "C", "G", "T"]
    list_frequencies = [
        frequencies["pi_A"],
        frequencies["pi_C"],
        frequencies["pi_G"],
        frequencies["pi_T"],
    ]
    model = msprime.MatrixMutationModel(
        alleles=alleles,
        root_distribution=list_frequencies,
        transition_matrix=transition_matrix,
    )

    mutated_ts = msprime.sim_mutations(
        tree_sequence,
        rate=rate_map,
        model=model,
        random_seed=seed,
    )

    return mutated_ts


def dna_ancestry_parameters_for_heritage(
    heritage: str, demography: msprime.Demography, sex_ratio: float
) -> tuple[msprime.Demography, int]:
    """Détermine la démographie (rescalée ou non) et la ploïdie pour un locus ADN.

    Reproduit le même dispatch que `simulate_genotypes_for_locus_type`
    côté SNP : "A" utilise la démographie <A> telle quelle en
    ploidy=2 ; "H"/"M" la rescalent par
    `coalescence_coefficient(heritage, sex_ratio) / 2` en ploidy=1
    (mêmes coefficients, mêmes formules, aucune raison structurelle
    qu'ils diffèrent entre SNP et séquences ADN -- comp_matQ/
    put_mutations opèrent sur le même arbre msprime que le Hudson SNP,
    seule la mutation diffère).

    Args:
        heritage: "A", "H", "M", "X" ou "Y".
        demography: La démographie <A> de base (PAS encore rescalée).
        sex_ratio: Fraction de mâles (voir observed_data.parse_sex_ratio).

    Returns:
        Le tuple (demography, ploidy) à utiliser pour ce locus.

    Raises:
        NotImplementedError: Pour "X"/"Y" -- contrairement au .snp
            (colonnes IND/SEX/POP), le format .mss (genepop) ne porte
            aucun sexe par individu -- `build_sex_stratified_samples_
            argument`/`build_male_only_samples_argument` n'ont pas
            d'équivalent exploitable sur ce format, et on ne devine
            pas un sexe par individu qui n'existe pas dans le fichier
            observé. Et pour tout autre type d'héritage inconnu.
    """
    if heritage == "A":
        return demography, 2
    elif heritage in ("H", "M"):
        rescaled_demography = rescale_demography(
            demography, coalescence_coefficient(heritage, sex_ratio) / 2
        )
        return rescaled_demography, 1
    elif heritage in ("X", "Y"):
        raise NotImplementedError(
            f"Locus ADN de type <{heritage}> non supporté : le format .mss "
            "ne porte pas de sexe par individu, nécessaire pour la "
            "stratification par sexe qu'exige ce type d'héritage."
        )
    else:
        raise NotImplementedError(
            f"Type d'héritage de locus non supporté: {heritage!r}"
        )


def dna_mutation_simulation_per_locus(
    header_text: str,
    mss_file_path: str | Path,
    demography: msprime.Demography,
    seed: int,
) -> dict[str, tskit.TreeSequence]:
    """Assemble le pipeline complet de simulation ADN, par locus.

    Pour chaque locus [S] : généalogie (msprime.sim_ancestry direct,
    pas simulate_independent_loci -- sequence_length variable par
    locus) + mutation (matQ/RateMap déjà construits).

    La démographie et la ploïdie utilisées pour l'arbre de coalescence de
    chaque locus dépendent de son type d'héritage (<A>/<H>/<M>, voir
    `dna_ancestry_parameters_for_heritage`) -- un groupe peut mélanger des
    loci de types différents (ex. toy_example2_ms_dna : G2 <A>, G3 <M>),
    donc ce dispatch se fait par locus, jamais une fois pour tout le
    dataset.

    Pour les loci [S] de type <A>, on tire une graine différente pour chaque
    locus (seed + _ANCESTRY_SEED_OFFSET + i),     pour que chaque locus <A>
    ait sa propre généalogie indépendante tandis que pour les loci mitochondriaux,
    on utilise la graine (seed + _SHARED_M_ANCESTRY_SEED_OFFSET) pour que tous les
    loci mitochondriaux partagent la même généalogie.
    Il n'y pas pour l'instant de support pour les loci X/Y, car le format .mss
    ne porte pas de sexe par individu et donc ne permet pas de stratification par sexe.

    Args:
        header_text: Texte complet de header.txt.
        mss_file_path: Chemin du fichier .mss.
        demography: La démographie <A> de base (PAS encore rescalée).
        seed: La graine de la simulation.

    Returns:
        Un dict {nom_locus: TreeSequence mutée} pour chaque locus [S].
    """
    rate_map_per_locus = build_rate_map_per_locus(header_text, seed)
    matrix_per_locus = build_matrix_per_locus(header_text, mss_file_path, seed)
    list_loci = parse_loci_description(header_text)
    frequencies_by_locus = base_frequency_by_locus(
        observed_sequences(mss_file_path, list_loci)
    )
    samples = observed_count_population(mss_file_path=mss_file_path)
    sex_ratio = parse_sex_ratio(mss_file_path)
    mutated_tree_sequences = {}

    for i, locus in enumerate(list_loci):
        if locus.ms_or_seq != "S":
            continue
        else:
            locus_demography, ploidy = dna_ancestry_parameters_for_heritage(
                locus.heritage, demography, sex_ratio
            )
            seed_offset = (
                seed + _SHARED_M_ANCESTRY_SEED_OFFSET
                if locus.heritage == "M"
                else seed + _ANCESTRY_SEED_OFFSET + i
            )
            tree_sequences = msprime.sim_ancestry(
                samples=samples,
                demography=locus_demography,
                sequence_length=locus.dnalength,
                random_seed=seed_offset,
                ploidy=ploidy,
            )
            transition_matrix = matrix_per_locus[locus.name]
            rate_map = rate_map_per_locus[locus.name]
            frequencies = frequencies_by_locus[locus.name]
            mutated_ts = simulate_dna_mutations(
                tree_sequences,
                transition_matrix,
                frequencies,
                rate_map,
                seed + _MUTATION_SEED_OFFSET + i,
            )
            mutated_tree_sequences[locus.name] = mutated_ts
    return mutated_tree_sequences


# equivalent à partir des valeurs tirées via diyabc


def _group_prior_values_from_columns(
    group_priors_values: dict[str, float], group_priors: dict
) -> dict[str, dict[str, float]]:
    """Reconstruit le dict nested {groupe: {prior: valeur}} depuis les colonnes du reftable réel.

    Reshape les colonnes plates du vrai reftable DIYABC (ex:
    "µseq_2", "k1seq_2") dans la forme nested que
    draw_group_parameter_values produit normalement, pour que
    build_group_local_param_per_locus_from_values puisse réutiliser
    tel quel le corps de build_group_local_param_per_locus.

    Args:
        group_priors_values: Dict {nom_colonne: valeur} tel que lu
            dans le vrai reftable (voir
            reftable_loop.parse_real_reftable_params_with_group_priors).
        group_priors: Dict {nom_groupe: [GroupPrior, ...]} (voir
            prior_parser.parse_group_priors).

    Returns:
        Un dict {nom_groupe: {nom_prior: valeur}}, pour les groupes de
        loci ADN ([S]) uniquement -- les groupes MicroSat sont ignorés.
    """
    group_values = {}
    for group, priors in group_priors.items():
        if next(prior for prior in priors).ms_or_seq != "S":
            continue  # Ignore les groupes de loci non ADN
        group_number = group[1:]  # Extrait le numéro du groupe (ex: "G1" -> "1")
        group_values[group] = {}
        gp_model = next(gp for gp in group_priors[group] if gp.model)
        k1_used, k2_used = get_parameter_used_by_model(gp_model)
        for prior in priors:
            if prior.name == "MEANMU":
                group_values[group][prior.name] = group_priors_values[
                    f"µseq_{group_number}"
                ]
            elif prior.name == "MEANK1" and k1_used:
                group_values[group][prior.name] = group_priors_values[
                    f"k1seq_{group_number}"
                ]
            elif prior.name == "MEANK2" and k2_used:
                group_values[group][prior.name] = group_priors_values[
                    f"k2seq_{group_number}"
                ]
            else:
                continue
        else:
            continue
    return group_values


def build_group_local_param_per_locus_from_values(
    header_text: str, group_priors_values: dict[str, float], seed: int
) -> dict[str, tuple[float, float, float]]:
    """Variante replay de build_group_local_param_per_locus (tier 1 = valeurs réelles).

    Le tirage par-groupe (premier niveau) est remplacé par les valeurs
    réellement tirées par DIYABC (`group_priors_values`, via
    `_group_prior_values_from_columns`) ; le tirage par-locus (second
    niveau, dispersion autour de la moyenne) N'EST PAS remplacé -- il
    continue de dépendre de `seed`, car DIYABC n'enregistre pas cette
    dispersion dans le reftable, il n'y a donc rien à rejouer pour elle.

    Args:
        header_text: Texte complet de header.txt.
        group_priors_values: Dict {nom_colonne: valeur} tel que lu
            dans le vrai reftable.
        seed: La graine du tirage par-locus (second niveau).

    Returns:
        Un dict {nom_locus: (k1, k2, mus_rate)}, même contrat que
        build_group_local_param_per_locus.
    """
    params_per_locus = {}
    group_priors = parse_group_priors(header_text)
    list_loci = parse_loci_description(header_text)

    list_loci_seq = [locus for locus in list_loci if locus.ms_or_seq == "S"]
    nloc_per_group = count_loci_per_group(list_loci_seq)

    # Un seul appel pour tous les groupes -- draw_group_parameter_values gère
    # déjà en interne son propre décalage de graine (_GROUP_PRIOR_SEED_OFFSET),
    # distinct de _KAPPA1_SEED_OFFSET/_KAPPA2_SEED_OFFSET utilisés plus bas
    # pour le tirage par locus -- pas besoin (et pas correct) de la rappeler
    # une fois par groupe/par kappa avec une graine décalée différente.
    values = _group_prior_values_from_columns(
        group_priors_values=group_priors_values, group_priors=group_priors
    )

    # Voir le commentaire équivalent dans build_group_local_param_per_locus :
    # un random.Random(seed + OFFSET) créé À CHAQUE ITÉRATION du groupe fait
    # rejouer la même séquence de tirages à tous les groupes partageant le
    # même modèle (ex: G2/G3 tous deux K2P sur toy_example2_ms_dna) --
    # corrigé en créant chaque rng UNE SEULE FOIS avant la boucle.
    mus_rate_rng = random.Random(seed + _MUS_RATE_SEED_OFFSET)
    kappa1_rng = random.Random(seed + _KAPPA1_SEED_OFFSET)
    kappa2_rng = random.Random(seed + _KAPPA2_SEED_OFFSET)

    for group in nloc_per_group:
        list_locus_in_group = [locus for locus in list_loci_seq if locus.group == group]
        # calcul du mus_rate par locus
        mus_rate = sampling_group_local_param(
            next(gp for gp in group_priors[group] if gp.name == "GAMMU"),
            k_moy=values[group]["MEANMU"],
            n_loci=nloc_per_group[group],
            check_nloc=True,
            list_loci=list_locus_in_group,
            rng=mus_rate_rng,
        )
        # calcul des kappas par locus selon le modèle de substitution du groupe
        gp_model = next(gp for gp in group_priors[group] if gp.model)
        verif_kappas = get_parameter_used_by_model(gp_model)
        if not verif_kappas[0] and not verif_kappas[1]:
            for locus in list_locus_in_group:
                # Handle the case where neither k1 nor k2 is used by the model
                params_per_locus[locus.name] = (0.0, 0.0, mus_rate[locus.name])
        elif verif_kappas[0] and not verif_kappas[1]:  # Model K2P et HKY
            gp_gamk1 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK1"), None
            )
            if gp_gamk1 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa1, mais aucun GAMK1 n'a été trouvé pour le groupe {group}."
                )
            k1_moy = values[group]["MEANK1"]
            kappa1_values = sampling_group_local_param(
                gp_gamk1,
                k_moy=k1_moy,
                n_loci=nloc_per_group[group],
                check_nloc=True,
                list_loci=list_locus_in_group,
                rng=kappa1_rng,
            )
            for locus in list_locus_in_group:
                params_per_locus[locus.name] = (
                    kappa1_values[locus.name],
                    0.0,
                    mus_rate[locus.name],
                )
        elif verif_kappas[0] and verif_kappas[1]:  # Modèle TN
            gp_gamk1 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK1"), None
            )
            if gp_gamk1 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa1, mais aucun GAMK1 n'a été trouvé pour le groupe {group}."
                )
            k1_moy = values[group]["MEANK1"]
            gp_gamk2 = next(
                (gp for gp in group_priors[group] if gp.name == "GAMK2"), None
            )
            if gp_gamk2 is None:
                raise ValueError(
                    f"Le modèle {gp_model} nécessite kappa2, mais aucun GAMK2 n'a été trouvé pour le groupe {group}."
                )
            k2_moy = values[group]["MEANK2"]
            kappa1_values = sampling_group_local_param(
                gp_gamk1,
                k_moy=k1_moy,
                n_loci=nloc_per_group[group],
                check_nloc=True,
                list_loci=list_locus_in_group,
                rng=kappa1_rng,
            )
            kappa2_values = sampling_group_local_param(
                gp_gamk2,
                k_moy=k2_moy,
                n_loci=nloc_per_group[group],
                check_nloc=False,
                list_loci=list_locus_in_group,
                rng=kappa2_rng,
            )
            for locus in list_locus_in_group:
                params_per_locus[locus.name] = (
                    kappa1_values[locus.name],
                    kappa2_values[locus.name],
                    mus_rate[locus.name],
                )
    return params_per_locus


def build_matrix_per_locus_from_values(
    header_text: str,
    mss_file_path: str | Path,
    group_priors_values: dict[str, float],
    seed: int,
) -> dict[str, np.ndarray]:
    """Variante replay de build_matrix_per_locus (premier niveau = valeurs réelles).

    Ne tire PAS les k1/k2 moyens par groupe (premier niveau, voir
    build_group_local_param_per_locus_from_values) : elle réutilise des
    valeurs déjà connues, typiquement les tirages RÉELS d'un reftable
    DIYABC existant (voir reftable_loop.
    parse_real_reftable_params_with_group_priors) -- pour comparer DIYABC
    et msprime sur EXACTEMENT les mêmes valeurs de k1/k2 par groupe, sans
    le biais possible de deux tirages indépendants.

    Le tirage par-locus (second niveau, la dispersion de k1/k2 autour de
    la moyenne du groupe, via sampling_group_local_param à l'intérieur de
    build_group_local_param_per_locus_from_values) N'EST PAS remplacé --
    il continue d'être tiré depuis `seed`, car cette valeur par-locus
    n'est jamais enregistrée dans le vrai reftable DIYABC (seule la
    moyenne de groupe l'est) : il n'y a donc rien à rejouer pour lui.

    Sinon identique à build_matrix_per_locus (même construction de
    matrice de transition via build_transition_matrix, mêmes fréquences
    de bases observées).

    Args:
        header_text: Texte complet de header.txt.
        mss_file_path: Chemin du fichier .mss.
        group_priors_values: Dict {nom_colonne: valeur} tel que lu
            dans le vrai reftable.
        seed: La graine du tirage par-locus (second niveau).

    Returns:
        Un dict {nom_locus: matQ} pour chaque locus [S], même contrat
        que build_matrix_per_locus.
    """
    list_loci = parse_loci_description(header_text)
    params_per_locus = build_group_local_param_per_locus_from_values(
        header_text, group_priors_values, seed
    )
    sequences_by_indiv = observed_sequences(mss_file_path, list_loci)
    frequencies_by_locus = base_frequency_by_locus(sequences_by_indiv)
    group_priors = parse_group_priors(header_text)

    matrix_per_locus = {}
    for locus in list_loci:
        if locus.ms_or_seq == "S":
            group_locus = locus.group
            gp_model = next(gp for gp in group_priors[group_locus] if gp.model)
            name_model = gp_model.name_model
            kappas = params_per_locus[locus.name][0], params_per_locus[locus.name][1]
            frequencies = frequencies_by_locus[locus.name]
            matrix_per_locus[locus.name] = build_transition_matrix(
                name_model, kappas, frequencies
            )
    return matrix_per_locus


def build_rate_map_per_locus_from_values(
    header_text: str, group_priors_values: dict[str, float], seed: int
) -> dict[str, msprime.RateMap]:
    """Variante replay de build_rate_map_per_locus (premier niveau = valeurs réelles).

    Ne tire PAS le mus_rate moyen par groupe (premier niveau) :
    réutilise group_priors_values, comme
    build_matrix_per_locus_from_values -- même principe, voir sa
    docstring pour le détail complet.

    Le tirage de `mutsit` (sample_site_rates, la dispersion du taux de
    mutation par SITE au sein d'un locus) N'EST PAS remplacé -- il
    continue de dépendre de `seed` (rng = random.Random(seed +
    _SITE_RATE_SEED_OFFSET)), pour la même raison que le tirage
    par-locus de k1/k2 : `mutsit` n'est jamais enregistré dans le vrai
    reftable DIYABC, il n'y a donc rien à rejouer pour lui non plus.

    Args:
        header_text: Texte complet de header.txt.
        group_priors_values: Dict {nom_colonne: valeur} tel que lu
            dans le vrai reftable.
        seed: La graine du tirage par-locus (second niveau) et de
            `mutsit`.

    Returns:
        Un dict {nom_locus: msprime.RateMap} pour chaque locus [S],
        même contrat que build_rate_map_per_locus.
    """
    list_loci = parse_loci_description(header_text)
    params_per_locus = build_group_local_param_per_locus_from_values(
        header_text, group_priors_values, seed
    )
    mus_rate_per_locus = {
        locus.name: params_per_locus[locus.name][2]
        for locus in list_loci
        if locus.ms_or_seq == "S"
    }
    group_priors = parse_group_priors(header_text)

    rate_map_per_locus = {}
    rng = random.Random(seed + _SITE_RATE_SEED_OFFSET)
    for locus in list_loci:
        if locus.ms_or_seq == "S":
            gp_model = next(gp for gp in group_priors[locus.group] if gp.model)
            mutsit = sample_site_rates(
                gp_model.p_fixe, gp_model.gams, locus.dnalength, rng=rng
            )
            mus_rate = mus_rate_per_locus[locus.name]
            rate_map_per_locus[locus.name] = build_rate_map(
                mus_rate=mus_rate, mutsit=mutsit, dnalength=locus.dnalength
            )
    return rate_map_per_locus


def dna_mutation_simulation_per_locus_from_values(
    header_text: str,
    mss_file_path: str | Path,
    demography: msprime.Demography,
    group_priors_values: dict[str, float],
    seed: int,
) -> dict[str, tskit.TreeSequence]:
    """Variante replay de dna_mutation_simulation_per_locus (rejeu apparié DIYABC/msprime).

    Voir build_matrix_per_locus_from_values pour le principe général :
    appelle build_rate_map_per_locus_from_values/
    build_matrix_per_locus_from_values au lieu des originales, pour que
    les k1/k2/mus_rate moyens par groupe soient ceux RÉELLEMENT tirés
    par DIYABC (group_priors_values) plutôt que tirés à nouveau depuis
    `seed`.

    `demography` reste un paramètre déjà construit par l'appelant, comme
    dans la version d'origine -- rien à changer ici pour les paramètres
    historiques (N1, ta, ts...), leur propre rejeu "from values" est géré
    en amont par pipeline.build_demography_for_scenario_index, réutilisée
    telle quelle (générique, ne sait rien de SNP vs ADN).

    Pour les loci [S] de type <A>, on tire une graine différente pour chaque
    locus (seed + _ANCESTRY_SEED_OFFSET + i), pour que chaque locus <A>
    ait sa propre généalogie indépendante tandis que pour les loci mitochondriaux,
    on utilise la graine (seed + _SHARED_M_ANCESTRY_SEED_OFFSET) pour que tous
    les loci mitochondriaux partagent la même généalogie.
    Il n'y pas pour l'instant de support pour les loci X/Y, car le format
    .mss ne porte pas de sexe par individu et donc ne permet pas de stratification
    par sexe.

    Args:
        header_text: Texte complet de header.txt.
        mss_file_path: Chemin du fichier .mss.
        demography: La démographie <A> de base (déjà construite par
            l'appelant à partir des valeurs historiques réelles).
        group_priors_values: Dict {nom_colonne: valeur} tel que lu
            dans le vrai reftable.
        seed: La graine du tirage par-locus (second niveau), de la
            généalogie et de la mutation.

    Returns:
        Un dict {nom_locus: TreeSequence mutée} pour chaque locus [S],
        même contrat que dna_mutation_simulation_per_locus.
    """
    rate_map_per_locus = build_rate_map_per_locus_from_values(
        header_text, group_priors_values, seed
    )
    matrix_per_locus = build_matrix_per_locus_from_values(
        header_text, mss_file_path, group_priors_values, seed
    )
    list_loci = parse_loci_description(header_text)
    frequencies_by_locus = base_frequency_by_locus(
        observed_sequences(mss_file_path, list_loci)
    )
    samples = observed_count_population(mss_file_path=mss_file_path)
    sex_ratio = parse_sex_ratio(mss_file_path)
    mutated_tree_sequences = {}

    for i, locus in enumerate(list_loci):
        if locus.ms_or_seq != "S":
            continue
        else:
            locus_demography, ploidy = dna_ancestry_parameters_for_heritage(
                locus.heritage, demography, sex_ratio
            )
            seed_offset = (
                seed + _SHARED_M_ANCESTRY_SEED_OFFSET
                if locus.heritage == "M"
                else seed + _ANCESTRY_SEED_OFFSET + i
            )
            tree_sequences = msprime.sim_ancestry(
                samples=samples,
                demography=locus_demography,
                sequence_length=locus.dnalength,
                random_seed=seed_offset,
                ploidy=ploidy,
            )
            transition_matrix = matrix_per_locus[locus.name]
            rate_map = rate_map_per_locus[locus.name]
            frequencies = frequencies_by_locus[locus.name]
            mutated_ts = simulate_dna_mutations(
                tree_sequences,
                transition_matrix,
                frequencies,
                rate_map,
                seed + _MUTATION_SEED_OFFSET + i,
            )
            mutated_tree_sequences[locus.name] = mutated_ts
    return mutated_tree_sequences
