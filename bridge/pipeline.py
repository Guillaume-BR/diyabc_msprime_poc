"""
Point d'entrée de haut niveau du pont DIYABC -> msprime.

Compose les briques indépendantes déjà testées (scenario_parser,
prior_parser, parameter_sampling, demography_builder) pour aller du texte
brut de header.txt jusqu'à une msprime.Demography prête à simuler.

Ce module ne contient aucune nouvelle logique de parsing ou de
construction : il orchestre uniquement.
"""

from pathlib import Path

import msprime

from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_independent_loci,
    simulate_snp_genotypes,
)
from bridge.demography_builder import build_demography
from bridge.parameter_sampling import draw_parameter_values
from bridge.prior_parser import parse_priors
from bridge.scenario_parser import parse_header_scenarios
from bridge.scenario_types import Scenario
from bridge.stats_group_parser import parse_requested_statistic_names
from bridge.summary_statistics import compute_all_statistics


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
    num_loci: int,
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
    samples = build_samples_argument(snp_path)

    tree_sequences = simulate_independent_loci(demography, samples, num_loci, seed)
    mutated = simulate_snp_genotypes(tree_sequences, seed)

    return mutated, values


def compute_summary_statistics(
    reference_directory: str | Path,
    scenario_index: int,
    num_loci: int,
    seed: int,
    work_directory: str | Path = None,  # gardé pour compatibilité, ignoré
    general_binary_path: str | Path = None,  # gardé pour compatibilité, ignoré
    stats_filter: str = "ALL",
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

    genotypes_per_locus, values = run_poc_for_directory(
        reference_directory, scenario_index, num_loci, seed
    )
    genotypes_list = list(genotypes_per_locus)

    samples = build_samples_argument(snp_path)
    population_names = list(samples.keys())

    summary_stats = compute_all_statistics(genotypes_list, population_names)

    if stats_filter == "ALL":
        pass
    elif stats_filter == "HEADER":
        requested_names = parse_requested_statistic_names(header_text)
        missing = [name for name in requested_names if name not in summary_stats]
        if missing:
            raise ValueError(
                f"header.txt déclare des statistiques non calculées par "
                f"compute_all_statistics (vocabulaire obsolète ou non "
                f"implémenté) : {missing}"
            )
        summary_stats = {name: summary_stats[name] for name in requested_names}
    else:
        raise NotImplementedError(
            f"stats_filter={stats_filter!r} non géré (valeurs connues : "
            f"'ALL', 'HEADER')"
        )

    return summary_stats, values
