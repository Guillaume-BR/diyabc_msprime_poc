"""
Point d'entrée de haut niveau du pont DIYABC -> msprime.

Compose les briques indépendantes déjà testées (scenario_parser,
prior_parser, parameter_sampling, demography_builder) pour aller du texte
brut de header.txt jusqu'à une msprime.Demography prête à simuler.

Ce module ne contient aucune nouvelle logique de parsing ou de
construction : il orchestre uniquement.

Deux familles de points d'entrée, chacune de bout en bout dans sa propre
section ci-dessous :
  - tirage ALÉATOIRE des paramètres (build_random_demography,
    run_poc_for_directory, compute_summary_statistics) ;
  - valeurs de paramètres DÉJÀ CONNUES / rejeu (build_demography_for_
    scenario_index, run_poc_for_directory_with_values,
    compute_summary_statistics_from_values) -- voir
    reftable_loop.replay_reftable_simulation.
Les deux partagent les mêmes helpers (section "Fondations" ci-dessous).
"""

import itertools
from collections.abc import Iterator
from pathlib import Path

import msprime

from bridge.ancestry_simulation import (
    build_samples_argument,
    dna_mutation_simulation_per_locus,
    dna_mutation_simulation_per_locus_from_values,
    simulate_genotypes_for_locus_type,
    simulate_poolseq_reads_with_mrc_filter,
)
from bridge.demography_builder import build_demography
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import detect_snp_file_type, observed_count_population
from bridge.parameter_sampling import draw_parameter_values
from bridge.prior_parser import parse_priors
from bridge.scenario_parser import parse_header_scenarios
from bridge.scenario_types import Scenario
from bridge.stats_group_parser import parse_requested_statistic_names
from bridge.summary_statistics import (
    compute_all_statistics,
    compute_all_statistics_dna,
    compute_all_statistics_poolseq,
)

# Décalage fixe et bien séparé par type de locus, à ajouter à la seed de
# base avant d'appeler simulate_genotypes_for_locus_type : réutiliser la
# MÊME seed pour tous les types corrèle artificiellement la position
# relative de leur première mutation (random.Random(seed) fraîchement
# recréé à chaque appel de simulate_snp_genotypes tire la même fraction
# relative target/total au premier tirage, quel que soit le total réel
# de branches -- vérifié empiriquement). Des offsets petits/positionnels
# (0,1,2,3...) ne suffisent PAS : ils entreraient en collision avec les
# seeds d'autres particules (seed = particle_index + 1 dans
# reftable_loop.py, donc de petits entiers eux aussi). Des offsets
# grands et bien séparés rendent une collision (particule, type) quasi
# impossible même sur des dizaines de milliers de particules -- msprime
# accepte des seeds jusqu'à 2**32-1.
LOCUS_TYPE_SEED_OFFSET = {
    "A": 0,
    "H": 10_000_000,
    "X": 20_000_000,
    "Y": 30_000_000,
    "M": 40_000_000,
}


# ── Fondations et helpers partagés par les deux familles ──────────────────


def read_header_text(directory: Path) -> str:
    """Lit header.txt si présent, sinon headerRF.txt en repli -- les deux
    noms coexistent selon les jeux de données (header.txt = config
    initiale fournie par l'utilisateur, headerRF.txt = variante produite
    par un run DIYABC réel ; nos sous-dossiers de test n'auront au
    départ que l'un des deux)."""
    header_path = directory / "header.txt"
    if not header_path.exists():
        header_path = directory / "headerRF.txt"
    return header_path.read_text()


def _simulate_genotypes_for_all_locus_types(
    demography: msprime.Demography,
    header_text: str,
    snp_path: Path,
    *,
    num_loci: int | None = None,
    seed: int,
) -> Iterator[dict[str, list[int]]]:
    """Boucle sur TOUS les types de locus déclarés dans la section 'loci
    description' de header_text (parse_loci_description(header_text).
    total_loci -- dict[str, int], ex: {"A": 5000} pour human, {"A": 70,
    "X": 10, "M": 10, "Y": 10} pour toy_example5), et concatène les
    génotypes simulés pour chacun via simulate_genotypes_for_locus_type.

    IMPORTANT -- num_loci est un compte PAR TYPE, pas un total : pour un
    dataset <A>-only comme human (un seul type déclaré), c'est
    rigoureusement identique au comportement actuel (num_loci loci de
    type <A>, point). Pour un dataset multi-type comme toy_example5,
    si l'on rentre une valeur précise de num_loci,
    num_loci loci sont simulés pour CHAQUE type déclaré -- pas les vrais
    comptes du header.txt (70/10/10/10) pour lesquels il faut passer
    num_loci=None. Le comportement par défaut (num_loci=None) est donc de
    simuler le nombre exact de loci déclaré dans header.txt pour chaque
    type, ce qui est le plus souvent ce que l'on veut pour un POC ou
    un test de validation.

    Un seed DISTINCT est dérivé par type de locus via
    LOCUS_TYPE_SEED_OFFSET (voir sa docstring en tête de fichier pour la
    justification empirique) -- ne JAMAIS appeler
    simulate_genotypes_for_locus_type avec la même seed brute pour
    plusieurs types dans cette boucle.

    Retourne un itérateur (pas une liste) : voir simulate_independent_loci
    pour la justification (ne pas matérialiser 51250 TreeSequence en
    mémoire simultanément).
    """

    total_loci = parse_loci_description(header_text).total_loci

    liste_iterateurs_par_type = []

    for locus_type, declared_count in total_loci.items():
        seed_for_type = seed + LOCUS_TYPE_SEED_OFFSET[locus_type]
        loci_count = num_loci if num_loci is not None else declared_count
        liste_iterateurs_par_type.append(
            simulate_genotypes_for_locus_type(
                demography, snp_path, locus_type, loci_count, seed_for_type
            )
        )
    return itertools.chain(*liste_iterateurs_par_type)


def _population_names(
    genotypes_list: list[dict[str, list[int]]], snp_path: Path
) -> list[str]:
    """Noms de population ("pop1", "pop2"...) dans le même ordre que
    build_samples_argument -- dérivés GRATUITEMENT des clés du premier
    locus déjà simulé (simulate_snp_genotypes construit ce dict avec
    exactement les mêmes noms, voir ancestry_simulation._population_layout)
    plutôt que de rescanner le fichier .snp une deuxième fois par
    particule (mesuré : ~4% du temps d'une particule sur human, voir
    notes/exploration.md, entrée du 20/07/2026). Repli sur
    build_samples_argument si genotypes_list est vide (num_loci=0, cas
    dégénéré qui n'arrive pas en pratique)."""
    if genotypes_list:
        return list(genotypes_list[0].keys())
    return list(build_samples_argument(snp_path).keys())


def _filter_statistics(
    summary_stats: dict[str, float],
    header_text: str,
    stats_filter: str,
) -> dict[str, float]:
    """Applique stats_filter ('ALL' ou 'HEADER') à un dict de statistiques
    déjà calculé -- factorisé entre compute_summary_statistics et
    compute_summary_statistics_from_values (même logique de filtrage,
    seule la source des valeurs de paramètres diffère entre les deux)."""
    if stats_filter == "ALL":
        return summary_stats
    elif stats_filter == "HEADER":
        requested_names = parse_requested_statistic_names(header_text)
        missing = [name for name in requested_names if name not in summary_stats]
        if missing:
            raise ValueError(
                f"header.txt déclare des statistiques non calculées par "
                f"compute_all_statistics (vocabulaire obsolète ou non "
                f"implémenté) : {missing}"
            )
        return {name: summary_stats[name] for name in requested_names}
    else:
        raise NotImplementedError(
            f"stats_filter={stats_filter!r} non géré (valeurs connues : "
            f"'ALL', 'HEADER')"
        )


# ── Tirage ALÉATOIRE des paramètres ────────────────────────────────────────


def build_random_demography(
    scenario: Scenario,
    header_text: str,
    seed: int,
) -> tuple[msprime.Demography, dict[str, float]]:
    """Tire des valeurs de paramètres à partir des priors déclarés dans
    header_text, puis construit la Demography msprime correspondant à
    scenario avec ces valeurs.

    Toutes les valeurs de priors du fichier sont tirées (pas seulement
    celles utilisées par ce scenario précis) : plus simple, et évite de
    casser des contraintes d'ordre qui pourraient porter sur des
    paramètres d'autres scénarios.

    Retourne (demography, values) -- les valeurs tirées sont renvoyées en
    plus de la Demography, car elles seront nécessaires plus tard pour
    écrire le reftable.bin (colonnes de paramètres).
    """
    priors, constraints = parse_priors(header_text)
    values = draw_parameter_values(priors, constraints, seed)
    demography = build_demography(scenario, values)
    return demography, values


def build_random_demography_for_scenario_index(
    header_text: str,
    scenario_index: int,
    seed: int,
) -> tuple[msprime.Demography, dict[str, float]]:
    """Variante pratique : sélectionne le scénario par son index (1-indexed,
    comme dans header.txt) plutôt que de demander un objet Scenario déjà
    parsé. Utile pour les tests et l'utilisation interactive.
    """
    scenarios = parse_header_scenarios(header_text)
    scenario = next((s for s in scenarios if s.index == scenario_index), None)
    if scenario is None:
        raise ValueError(
            f"Scénario {scenario_index} non trouvé ou non géré par le parser "
            f"(scénarios disponibles : {sorted(s.index for s in scenarios)})"
        )
    return build_random_demography(scenario, header_text, seed)


def run_poc_for_directory(
    directory: str | Path,
    scenario_index: int,
    *,
    num_loci: int | None = None,
    seed: int,
):
    """Point d'entrée de haut niveau : équivalent du `-p ./` de DIYABC.

    Prend un dossier contenant header.txt et le fichier de données
    observées (.snp), et produit num_loci TreeSequence mutées, simulées
    sous le scénario demandé.

    Le nom du fichier de données est lu sur la PREMIÈRE LIGNE de
    header.txt (ex: "human_snp_all22chr_maf5.snp"), pas deviné par
    extension -- c'est le contrat du format DIYABC.

    Retourne (mutated_tree_sequences, values) : l'itérateur des
    TreeSequence mutées, et le dict des valeurs de paramètres tirées
    (nécessaires plus tard pour écrire le reftable.bin).
    """
    directory = Path(directory)
    header_text = read_header_text(directory)

    snp_filename = header_text.splitlines()[0].strip()
    snp_path = directory / snp_filename

    demography, values = build_random_demography_for_scenario_index(
        header_text, scenario_index, seed
    )

    mutated = _simulate_genotypes_for_all_locus_types(
        demography, header_text, snp_path, num_loci=num_loci, seed=seed
    )

    return mutated, values


def compute_summary_statistics(
    reference_directory: str | Path,
    scenario_index: int,
    *,
    num_loci: int | None = None,
    seed: int,
    work_directory: str | Path = None,  # gardé pour compatibilité, ignoré
    general_binary_path: str | Path = None,  # gardé pour compatibilité, ignoré
    stats_filter: str = "ALL",
    observed_reads_per_locus: list[dict[str, tuple[int, int]]] = None,
) -> tuple[dict[str, float], dict[str, float]]:
    """Calcule les statistiques résumées sur des données SIMULÉES par
    notre pipeline, en utilisant nos formules Python validées
    (summary_statistics.py) -- remplace la délégation au binaire C++
    (subprocess + fichier .snp intermédiaire).

    stats_filter :
      - "ALL" (défaut) : retourne toutes les statistiques implémentées
        (compute_all_statistics), sans filtrage.
      - "HEADER" : ne garde, dans l'ordre de déclaration, que les
        statistiques listées dans la section 'group summary statistics'
        de header.txt (voir stats_group_parser.
        parse_requested_statistic_names) -- nécessaire pour que
        reftable_msprime.txt/.bin aient EXACTEMENT les mêmes colonnes
        que le vrai reftable DIYABC (sinon toute comparaison
        colonne-par-nom entre les deux pipelines est faussée, comme
        découvert sur toy_example5_modif : 'ML3p_1.2.3' calculé par
        nous mais absent du vrai DIYABC). Lève ValueError si header.txt
        déclare une statistique qu'on ne sait pas calculer (vocabulaire
        obsolète, ex: human/header.txt -- voir notes/exploration.md).

    Retourne (summary_statistics, parameter_values).
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    snp_filename = header_text.splitlines()[0].strip()
    snp_path = reference_directory / snp_filename

    if detect_snp_file_type(snp_path) == "IND":
        genotypes_per_locus, values = run_poc_for_directory(
            reference_directory,
            scenario_index=scenario_index,
            num_loci=num_loci,
            seed=seed,
        )
        genotypes_list = list(genotypes_per_locus)

        population_names = _population_names(genotypes_list, snp_path)
        summary_stats = compute_all_statistics(genotypes_list, population_names)
        summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)
    else:
        total_loci_poolseq = parse_loci_description(header_text).total_loci["A"]
        demography, values = build_random_demography_for_scenario_index(
            header_text, scenario_index, seed
        )

        reads_list = list(
            simulate_poolseq_reads_with_mrc_filter(
                demography,
                snp_path,
                seed,
                num_loci=total_loci_poolseq,
                observed_reads_per_locus=observed_reads_per_locus,
            )
        )
        pool_sizes = build_samples_argument(snp_path)
        population_names = list(pool_sizes.keys())
        summary_stats = compute_all_statistics_poolseq(
            reads_list, population_names, pool_sizes
        )
        summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)
    return summary_stats, values


# ── Valeurs de paramètres DÉJÀ CONNUES (rejeu, voir reftable_loop) ─────────


def build_demography_for_scenario_index(
    header_text: str,
    scenario_index: int,
    values: dict[str, float],
) -> msprime.Demography:
    """Variante de build_random_demography_for_scenario_index qui NE TIRE
    AUCUNE valeur : construit la Demography directement à partir de
    valeurs de paramètres déjà connues (ex: reprises telles quelles d'un
    reftable DIYABC réel pour servir d'oracle-- voir reftable_loop.replay_reftable_simulation).
    """
    scenarios = parse_header_scenarios(header_text)
    scenario = next((s for s in scenarios if s.index == scenario_index), None)
    if scenario is None:
        raise ValueError(
            f"Scénario {scenario_index} non trouvé ou non géré par le parser "
            f"(scénarios disponibles : {sorted(s.index for s in scenarios)})"
        )
    return build_demography(scenario, values)


def run_poc_for_directory_with_values(
    directory: str | Path,
    scenario_index: int,
    values: dict[str, float],
    *,
    num_loci: int | None = None,
    seed: int,
):
    """Variante de run_poc_for_directory qui prend des valeurs de
    paramètres déjà connues au lieu d'en tirer de nouvelles -- même
    contrat par ailleurs (lecture du nom de fichier .snp sur la première
    ligne de header.txt, échantillonnage, simulation)."""
    directory = Path(directory)
    header_text = read_header_text(directory)

    snp_filename = header_text.splitlines()[0].strip()
    snp_path = directory / snp_filename

    demography = build_demography_for_scenario_index(
        header_text, scenario_index, values
    )

    return _simulate_genotypes_for_all_locus_types(
        demography, header_text, snp_path, num_loci=num_loci, seed=seed
    )


def compute_summary_statistics_from_values(
    reference_directory: str | Path,
    scenario_index: int,
    values: dict[str, float],
    *,
    num_loci: int | None = None,
    seed: int,
    stats_filter: str = "ALL",
    observed_reads_per_locus: list[dict[str, tuple[int, int]]] = None,
) -> dict[str, float]:
    """Variante de compute_summary_statistics qui NE TIRE AUCUNE valeur de
    prior : reprend telles quelles des valeurs de paramètres déjà
    connues, typiquement les tirages RÉELS d'un reftable DIYABC existant
    (voir reftable_loop.replay_reftable_simulation) -- permet de comparer
    DIYABC et msprime sur EXACTEMENT les mêmes tirages de priors, sans le
    biais possible de deux tirages indépendants.
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    snp_filename = header_text.splitlines()[0].strip()
    snp_path = reference_directory / snp_filename

    if detect_snp_file_type(snp_path) == "IND":
        genotypes_per_locus = run_poc_for_directory_with_values(
            reference_directory, scenario_index, values, num_loci=num_loci, seed=seed
        )
        genotypes_list = list(genotypes_per_locus)

        population_names = _population_names(genotypes_list, snp_path)
        summary_stats = compute_all_statistics(genotypes_list, population_names)
        summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)
    else:
        # l'argument num_loci est ignoré ici, côté PoolSeq on simule tous les loci déclarés dans header.txt
        # alors qu'en Indseq num_loci sert à limiter le nombre de loci simulés

        total_loci_poolseq = parse_loci_description(header_text).total_loci["A"]
        demography = build_demography_for_scenario_index(
            header_text, scenario_index, values
        )

        reads_list = list(
            simulate_poolseq_reads_with_mrc_filter(
                demography,
                snp_path,
                seed,
                num_loci=total_loci_poolseq,
                observed_reads_per_locus=observed_reads_per_locus,
            )
        )
        pool_sizes = build_samples_argument(snp_path)
        population_names = list(pool_sizes.keys())
        summary_stats = compute_all_statistics_poolseq(
            reads_list, population_names, pool_sizes
        )
        summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)
    return summary_stats


def compute_summary_statistics_dna(
    reference_directory: str | Path,
    scenario_index: int,
    *,
    seed: int,
    stats_filter: str = "ALL",
) -> tuple[dict[str, float], dict[str, float]]:
    """Calcule les 13 statistiques résumées ADN (compute_all_statistics_dna)
    sur des données SIMULÉES par msprime -- équivalent ADN de
    compute_summary_statistics (chemin IND/PoolSeq), pour les datasets
    qui déclarent des loci séquence (`[S]`, groupes `G2`/`G3`... de
    header.txt) plutôt que des SNP.

    Tire les paramètres historiques (N1, ta, ts...) ET les priors de
    groupe (k1/k2/mus_rate par groupe ADN, en interne à
    dna_mutation_simulation_per_locus) depuis `seed` -- voir
    compute_summary_statistics_dna_from_values pour la variante qui
    rejoue des valeurs déjà connues plutôt que d'en tirer de nouvelles
    (paired comparison avec un vrai reftable DIYABC).

    `values` (le second élément du tuple retourné) ne contient QUE les
    paramètres historiques, pas les priors de groupe -- dna_mutation_
    simulation_per_locus ne renvoie nulle part les valeurs de k1/k2/
    mus_rate qu'elle a tirées en interne, donc ce `values` seul ne
    suffirait pas à rejouer exactement cette même particule (contrairement
    au chemin SNP, où `values` capture tout ce qui a été tiré).

    Args:
        reference_directory: dossier contenant header.txt/headerRF.txt
            et le fichier .mss observé (son nom lu sur la première ligne
            du header).
        scenario_index: le scénario à utiliser pour construire la
            démographie (pas de tirage pondéré multi-scénario ici,
            contrairement à reftable_loop.run_reftable_simulation).
        seed: graine de la particule -- dérive toutes les graines
            internes (tirage des paramètres historiques, des priors de
            groupe, des généalogies et mutations par locus).
        stats_filter: "ALL" (toutes les stats implémentées) ou "HEADER"
            (seulement celles déclarées dans header.txt, voir
            compute_summary_statistics pour le détail).

    Returns:
        (summary_stats, values) -- summary_stats est le dict {nom_
        colonne_diyabc: valeur} de compute_all_statistics_dna (ex.
        "NSS_2_1"), values est {nom_paramètre_historique: valeur}.
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename

    demography, values = build_random_demography_for_scenario_index(
        header_text, scenario_index, seed
    )

    mutated = dna_mutation_simulation_per_locus(
        header_text,
        mss_path,
        demography,
        seed,
    )

    population_names = list(observed_count_population(mss_path).keys())
    summary_stats = compute_all_statistics_dna(header_text, mutated, population_names)
    summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)

    return summary_stats, values


def compute_summary_statistics_dna_from_values(
    reference_directory: str | Path,
    scenario_index: int,
    values: dict[str, float],
    group_priors_values: dict[str, float],
    *,
    seed: int,
    stats_filter: str = "ALL",
) -> dict[str, float]:
    """Variante de compute_summary_statistics_dna qui ne tire AUCUNE valeur
    de prior : reprend telles quelles des valeurs de paramètres déjà
    connues, typiquement les tirages RÉELS d'un reftable DIYABC existant
    (voir reftable_loop.replay_reftable_simulation) -- permet de comparer
    DIYABC et msprime sur EXACTEMENT les mêmes tirages de priors, sans le
    biais possible de deux tirages indépendants.

    group_priors_values : dict[str, float] -- {nom_param: valeur}
    pour tous les groupes ADN déclaré dans header.txt. Ces valeurs sont celles que dna_mutation
    simulation_per_locus aurait tirées en interne si on avait appelé la
    variante "random" (compute_summary_statistics_dna) -- elles ne sont
    pas capturées par le dict `values` retourné par cette fonction.
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename

    demography = build_demography_for_scenario_index(
        header_text, scenario_index, values
    )

    mutated = dna_mutation_simulation_per_locus_from_values(
        header_text,
        mss_path,
        demography,
        group_priors_values=group_priors_values,
        seed=seed,
    )

    population_names = list(observed_count_population(mss_path).keys())
    summary_stats = compute_all_statistics_dna(header_text, mutated, population_names)
    summary_stats = _filter_statistics(summary_stats, header_text, stats_filter)

    return summary_stats
