"""
Point d'entrée de haut niveau du pont DIYABC -> msprime.

Compose les briques indépendantes déjà testées (scenario_parser,
prior_parser, parameter_sampling, demography_builder) pour aller du texte
brut de header.txt jusqu'à une msprime.Demography prête à simuler.

Ce module ne contient aucune nouvelle logique de parsing ou de
construction : il orchestre uniquement.
"""

import subprocess
import shutil
from pathlib import Path

import msprime

from bridge.scenario_parser import parse_header_scenarios
from bridge.prior_parser import parse_priors
from bridge.parameter_sampling import draw_parameter_values
from bridge.demography_builder import build_demography
from bridge.scenario_types import Scenario
from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_independent_loci,
    simulate_snp_genotypes,
)
from bridge.snp_writer import write_snp_file
from bridge.loci_parser import rewrite_loci_count
from bridge.statobs_parser import parse_statobs


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
    header_text = (directory / "header.txt").read_text()

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
    general_binary_path: str | Path,
    work_directory: str | Path,
    stats_filter: str = "ALL",
) -> tuple[dict[str, float], dict[str, float]]:
    """Calcule les statistiques résumées sur des données SIMULÉES par
    notre pipeline, en délégant le calcul au binaire C++ `general`
    existant (HeaderC::calstatobs) -- plutôt que de réimplémenter
    FST1/ML1/F3/F4/NEI/AML en Python.

    Étapes :
    1. Simule num_loci génotypes SNP sous un tirage de paramètres pour
       scenario_index (via run_poc_for_directory)
    2. Écrit ces génotypes dans un fichier .snp (write_snp_file)
    3. Copie header.txt (en adaptant son nombre de loci déclaré à
       num_loci -- rewrite_loci_count) et RNG_state_0000.bin dans
       work_directory
    4. Appelle `general -p <work_directory> -R <stats_filter> -r 1 ...`
       en sous-processus -- -r 1 suffit, calstatobs() s'exécute avant
       même la boucle de simulation (voir notes/exploration.md)
    5. Parse statobsRF.txt en dict {nom_stat: valeur}

    work_directory doit déjà exister (n'est PAS créé automatiquement).

    Retourne (summary_statistics, parameter_values) : les statistiques
    calculées par le C++, et les valeurs de paramètres tirées par notre
    pipeline (pour constituer plus tard une ligne complète du
    reftable.bin : param[] + stat[]).
    """
    reference_directory = Path(reference_directory)
    work_directory = Path(work_directory)
    general_binary_path = Path(general_binary_path)

    original_header_text = (reference_directory / "header.txt").read_text()
    snp_filename = original_header_text.splitlines()[0].strip()

    genotypes_per_locus, values = run_poc_for_directory(
        reference_directory, scenario_index, num_loci, seed
    )

    write_snp_file(list(genotypes_per_locus), work_directory / snp_filename)

    adapted_header_text = rewrite_loci_count(original_header_text, num_loci)
    (work_directory / "header.txt").write_text(adapted_header_text)

    rng_state_source = reference_directory / "RNG_state_0000.bin"
    if rng_state_source.exists():
        shutil.copy(rng_state_source, work_directory / "RNG_state_0000.bin")
        print(f"Copié {rng_state_source} -> {work_directory / 'RNG_state_0000.bin'}")

    subprocess.run(
        [
            str(general_binary_path),
            "-p", "./",
            "-R", stats_filter,
            "-r", "1",
            "-g", "50",
            "-m",
            "-t", "1",
        ],
        cwd=work_directory,
        check=True,
        capture_output=True,
    )

    statobs_text = (work_directory / "statobsRF.txt").read_text()
    summary_statistics = parse_statobs(statobs_text)

    return summary_statistics, values