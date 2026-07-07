"""
Boucle d'itération produisant les nrec "particules" (lignes) d'un futur
reftable.bin : pour chaque particule, un tirage de scénario et de
paramètres distinct, et une simulation msprime complète (calcul des
statistiques résumées via compute_summary_statistics, 100% Python --
plus de subprocess ni de fichier intermédiaire sur disque).

Parallélisé via ProcessPoolExecutor : chaque particule est indépendante
des autres (son propre tirage, sa propre simulation), donc embarrassingly
parallel.
"""

import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.parameter_sampling import draw_scenario
from bridge.pipeline import compute_summary_statistics, read_header_text
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.scenario_parser import parse_header_scenarios
from bridge.scenario_types import Scenario


@dataclass
class ParticleResult:
    """Le résultat d'une particule : une future ligne du reftable.bin."""

    particle_index: int
    scenario_index: int
    parameter_values: dict[str, float]
    summary_statistics: dict[str, float]


def _run_single_particle(
    particle_index: int,
    reference_directory: Path,
    scenarios: list[Scenario],
    num_loci: int,
    stats_filter: str,
) -> ParticleResult:
    """Calcule une seule particule -- fonction top-level (picklable),
    appelée par chaque worker du ProcessPoolExecutor.

    La seed utilisée est dérivée de particle_index, garantissant un
    tirage distinct et reproductible par particule (même particle_index
    -> même résultat, peu importe l'ordre d'exécution des workers).

    IMPORTANT : seed = particle_index + 1, jamais particle_index seul.
    msprime.sim_ancestry rejette explicitement seed=0 (ValueError "seeds
    must be greater than 0 and less than 2^32") -- vérifié empiriquement.
    Donc particle_index=0 (le cas le plus probable, première particule)
    utilise seed=1, pas seed=0.
    """
    seed = particle_index + 1
    drawn_scenario = draw_scenario(scenarios, seed)

    summary_statistics, parameter_values = compute_summary_statistics(
        reference_directory=reference_directory,
        scenario_index=drawn_scenario.index,
        num_loci=num_loci,
        seed=seed,
        stats_filter=stats_filter,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=drawn_scenario.index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
    )


def run_reftable_simulation(
    reference_directory: str | Path,
    scenarios: list[Scenario],
    num_loci: int,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Produit nrec particules (lignes de reftable.bin) en parallèle.

    `scenarios` est la liste des scénarios candidats (typiquement TOUS
    les scénarios déclarés dans header.txt) : chaque particule tire le
    SIEN au hasard, pondéré par son `weight` (voir
    parameter_sampling.draw_scenario, sémantique vérifiée contre
    particuleC.cpp::ParticleC::drawscenario) -- une même particule peut
    donc finir sur n'importe lequel des scénarios de la liste, pas
    forcément le même pour toutes.

    N'écrit rien sur disque par particule (compute_summary_statistics
    est 100% Python, en mémoire) -- pas de dossier de travail à créer
    ni à nettoyer, contrairement à l'ancienne architecture par
    subprocess+.snp intermédiaire.

    Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
    nrec-1), pas dans l'ordre de complétion des workers -- important
    pour la reproductibilité de l'ordre des lignes du reftable final.

    max_workers : nombre de process en parallèle (défaut : laissé à
    ProcessPoolExecutor, généralement le nombre de cœurs disponibles).
    """
    reference_directory = Path(reference_directory)

    results_by_index: dict[int, ParticleResult] = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle,
                particle_index,
                reference_directory,
                scenarios,
                num_loci,
                stats_filter,
            ): particle_index
            for particle_index in range(nrec)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()

    return [results_by_index[i] for i in range(nrec)]


def _kept_param_names_by_scenario(
    priors: list, scenarios: list[Scenario]
) -> dict[int, list[str]]:
    """Pour chaque scénario, la liste (dans l'ordre de déclaration des
    priors) des noms de paramètres à garder : non constants
    (is_constant_prior) ET référencés par CE scénario précis
    (get_parameter_names_used_by_scenario) -- même filtre à deux
    critères que l'ancienne version single-scenario, appliqué
    séparément par scénario."""
    result = {}
    for scenario in scenarios:
        used_param_names = get_parameter_names_used_by_scenario(scenario)
        result[scenario.index] = [
            p.name
            for p in priors
            if not is_constant_prior(p) and p.name in used_param_names
        ]
    return result


def write_reftable_bin(
    results: list[ParticleResult],
    priors: list,
    scenarios: list[Scenario],
    output_path: str | Path,
) -> None:
    """Écrit un reftable.bin au format binaire DIYABC (vérifié contre
    reftable.cpp et un vrai reftableRF.bin multi-scénario -- voir
    docs/synthese_diyabc_msprime.docx section 5).

    `scenarios` est la liste de TOUS les scénarios candidats déclarés
    dans header.txt (pas seulement ceux effectivement tirés dans
    `results`) : nscen = len(scenarios), et le numéro de scénario écrit
    par ligne est le numéro 1-indexed du header.txt (scenario.index),
    jamais renuméroté localement -- vérifié dans particuleC.cpp::
    drawscenario et reftable.cpp.

    IMPORTANT -- format à LONGUEUR VARIABLE par ligne, PAS d'union de
    colonnes ni de valeur NA écrite sur disque : chaque ligne écrit
    SEULEMENT nparam[scenario-1] floats de paramètres, ceux de son
    propre scénario -- vérifié empiriquement contre un vrai
    reftableRF.bin multi-scénario (dataset MER modelchoice/PoolSeq) et
    contre reftable.cpp (boucle d'écriture indexée par
    nparam[numscen-1]). La reconstruction en matrice rectangulaire avec
    NA pour les colonnes non concernées est une responsabilité du
    LECTEUR (readReftable.R), jamais de l'écrivain -- ne PAS essayer de
    remplir les colonnes manquantes ici.

    nparam[i] : nombre de paramètres (non constants, référencés) du
    i-ème scénario de `scenarios` -- pilote directement la taille de
    chaque enregistrement, comme dans reftable.cpp.

    Ne gère PAS les paramètres de mutation (absents de human) -- à
    ajouter (toujours en dernière position, après les paramètres
    démographiques -- voir readReftable.R) si un dataset avec
    microsatellites/séquences est traité plus tard.
    """
    if not results:
        raise ValueError("results est vide : au moins une particule est requise")

    known_indices = {s.index for s in scenarios}
    unknown = {r.scenario_index for r in results} - known_indices
    if unknown:
        raise ValueError(
            f"results contient des scenario_index absents de `scenarios` : {unknown}"
        )

    kept_param_names_by_scenario = _kept_param_names_by_scenario(priors, scenarios)
    stat_names = list(results[0].summary_statistics.keys())

    nrec = len(results)
    nscen = len(scenarios)
    nrecscen = [
        sum(1 for r in results if r.scenario_index == s.index) for s in scenarios
    ]
    nparam = [len(kept_param_names_by_scenario[s.index]) for s in scenarios]
    nstat = len(stat_names)

    with open(output_path, "wb") as f:
        f.write(struct.pack("<i", nrec))
        f.write(struct.pack("<i", nscen))
        for n in nrecscen:
            f.write(struct.pack("<i", n))
        for n in nparam:
            f.write(struct.pack("<i", n))
        f.write(struct.pack("<i", nstat))

        for result in results:
            f.write(struct.pack("<i", result.scenario_index))
            for name in kept_param_names_by_scenario[result.scenario_index]:
                f.write(struct.pack("<f", result.parameter_values[name]))
            for name in stat_names:
                f.write(struct.pack("<f", result.summary_statistics[name]))


def write_reftable_txt(
    results: list[ParticleResult],
    priors: list,
    scenarios: list[Scenario],
    output_path: str | Path,
) -> None:
    """Écrit les résultats au format texte de DIYABC :
    first_records_of_the_reference_table_0.txt.

    Format reproduit depuis particleset.cpp / header.cpp :
      - Ligne 1 : noms de colonnes, chacun centré sur 14 caractères
        (fonction C++ centre(s1, 14)).
      - Lignes suivantes : "%3d  " pour le numéro de scénario, puis
        "  %12.6f" pour chaque paramètre et statistique.

    Note sur le format des paramètres : le C++ distingue categ<2 (%12.0f,
    entiers) vs categ>=2 (%12.3f, flottants), mais cette catégorie n'est
    pas exposée dans les priors Python. On utilise %12.6f uniformément --
    suffisant pour la comparaison statistique des distributions.

    Le texte utilise un jeu de colonnes de paramètres FIXE : l'UNION
    (dans l'ordre de déclaration des priors) des paramètres utilisés par
    au moins un des `scenarios`. Pour une ligne dont le scénario tiré
    n'utilise pas tel paramètre, on écrit quand même sa valeur RÉELLEMENT
    TIRÉE (r.parameter_values[name] -- toujours présente, puisque
    draw_parameter_values tire TOUS les priors déclarés, indépendamment
    du scénario, voir build_random_demography), jamais une case vide.

    IMPORTANT : ne JAMAIS laisser de case vide ici, même si elle
    représente une valeur non pertinente pour le scénario de la ligne --
    un parseur par espaces (ex: pandas read_csv(sep=r'\\s+'), ou même un
    simple line.split()) ne produit AUCUN token pour une case vide,
    ce qui décale d'une colonne TOUTES les valeurs suivantes sur la
    ligne. Bug découvert empiriquement (comparaison DIYABC/msprime sur
    toy_example5_modif, colonne 'r' non utilisée par le scénario actif
    laissée en blanc -> décalage systématique des 51 statistiques sur
    les 1000 lignes, provoquant des "écarts" massifs et incohérents qui
    n'avaient rien à voir avec un vrai écart de simulation).
    """
    if not results:
        raise ValueError("results est vide : au moins une particule est requise")

    kept_param_names_by_scenario = _kept_param_names_by_scenario(priors, scenarios)
    used_by_any = {
        name for names in kept_param_names_by_scenario.values() for name in names
    }
    all_param_names = [p.name for p in priors if p.name in used_by_any]
    stat_names = list(results[0].summary_statistics.keys())

    def _centre(s: str, width: int = 14) -> str:
        """Reproduit centre() de mesutils.cpp : centrage sur width chars."""
        return s.center(width)

    with open(output_path, "w", encoding="utf-8") as f:
        # En-tête : "scenario" + noms des paramètres + noms des stats
        header = (
            _centre("scenario")
            + "".join(_centre(n) for n in all_param_names)
            + "".join(_centre(n) for n in stat_names)
        )
        f.write(header + "\n")

        # Une ligne par particule
        for r in results:
            line = f"{r.scenario_index:3d}  "
            for name in all_param_names:
                line += f"  {r.parameter_values[name]:12.6f}"
            for name in stat_names:
                line += f"  {r.summary_statistics[name]:12.6f}"
            f.write(line + "\n")


def simulate_reference_directory(
    test_directory: str | Path,
    num_loci: int,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Point d'entrée pour un sous-dossier de test sous reference/ (ex:
    reference/mon_test/), qui ne contient au départ qu'un header.txt (ou
    headerRF.txt, repli si absent -- voir pipeline.read_header_text) et
    le fichier .snp observé (nommé sur la première ligne du header, pas
    un nom fixe -- voir pipeline.run_poc_for_directory).

    Tire les scénarios candidats pondérés par leur `weight` parmi TOUS
    ceux déclarés dans le header (voir parameter_sampling.draw_scenario),
    simule nrec particules, et écrit le résultat dans
    test_directory/reftable_msprime.txt -- jamais "reftable.txt", pour
    ne pas être confondu avec le first_records_of_the_reference_table_0.
    txt qu'un vrai run DIYABC produirait dans le même dossier.

    Retourne les ParticleResult (utile pour appeler write_reftable_bin
    en plus, si besoin -- pas fait ici, ce point d'entrée ne produit que
    le texte).
    """
    test_directory = Path(test_directory)
    header_text = read_header_text(test_directory)

    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)

    results = run_reftable_simulation(
        reference_directory=test_directory,
        scenarios=scenarios,
        num_loci=num_loci,
        nrec=nrec,
        stats_filter=stats_filter,
        max_workers=max_workers,
    )

    write_reftable_txt(
        results, priors, scenarios, test_directory / "reftable_msprime.txt"
    )

    write_reftable_bin(
        results, priors, scenarios, test_directory / "reftable_msprime.bin"
    )

    return results
