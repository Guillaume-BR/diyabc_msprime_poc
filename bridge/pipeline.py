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

from bridge.scenario_parser import parse_header_scenarios
from bridge.prior_parser import parse_priors
from bridge.parameter_sampling import draw_parameter_values
from bridge.demography_builder import build_demography
from bridge.scenario_types import Scenario
from bridge.observed_data import count_samples_per_population, population_index_to_name
from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_independent_loci,
    simulate_snp_genotypes,
)



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
    observées (.snp), tire un jeu de paramètres démographiques, et
    simule les génotypes SNP de num_loci loci indépendants sous ce
    tirage -- un seul tirage de paramètres pour toute la particule,
    cohérent avec le fonctionnement de DIYABC (drawscenario/
    setHistParamValue appelés une fois, avant la boucle sur les loci).

    Le nom du fichier de données est lu sur la PREMIÈRE LIGNE de
    header.txt (ex: "human_snp_all22chr_maf5.snp"), pas deviné par
    extension -- c'est le contrat du format DIYABC.

    Chaque locus est garanti polymorphe par construction (algorithme de
    Hudson, une mutation unique par locus -- voir
    ancestry_simulation.simulate_snp_genotypes).

    Retourne (genotypes_per_locus, values) :
    - genotypes_per_locus : itérateur de num_loci dicts
      {id_lignée: 0 ou 1}, un par locus (PAS des TreeSequence)
    - values : dict des valeurs de paramètres démographiques tirées,
      nécessaires plus tard pour écrire le reftable.bin (colonnes de
      paramètres)
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


