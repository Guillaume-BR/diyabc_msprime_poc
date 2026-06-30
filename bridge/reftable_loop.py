"""
Boucle d'itération produisant les nrec "particules" (lignes) d'un futur
reftable.bin : pour chaque particule, un tirage de paramètres distinct,
une simulation msprime complète, et un calcul de statistiques résumées
délégué au binaire C++ (compute_summary_statistics).

Parallélisé via ProcessPoolExecutor : chaque particule est indépendante
des autres (son propre tirage, sa propre simulation), donc embarrassingly
parallel. Chaque worker utilise un work_directory DISTINCT (basé sur
l'index de la particule), pour éviter toute collision d'écriture entre
processus concurrents sur les mêmes fichiers (.snp, statobsRF.txt...).
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import struct

from bridge.pipeline import compute_summary_statistics
from bridge.prior_parser import is_constant_prior
from bridge.demography_builder import get_parameter_names_used_by_scenario


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
    scenario_index: int,
    num_loci: int,
    general_binary_path: Path,
    base_work_directory: Path,
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
    work_directory = base_work_directory / f"particle_{particle_index}"
    work_directory.mkdir(parents=True, exist_ok=True)

    summary_statistics, parameter_values = compute_summary_statistics(
        reference_directory=reference_directory,
        scenario_index=scenario_index,
        num_loci=num_loci,
        seed=particle_index + 1,
        general_binary_path=general_binary_path,
        work_directory=work_directory,
        stats_filter=stats_filter,
    )

    return ParticleResult(
        particle_index=particle_index,
        scenario_index=scenario_index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
    )


def run_reftable_simulation(
    reference_directory: str | Path,
    scenario_index: int,
    num_loci: int,
    nrec: int,
    general_binary_path: str | Path,
    base_work_directory: str | Path,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Produit nrec particules (lignes de reftable.bin) en parallèle.

    base_work_directory doit déjà exister ; un sous-dossier
    "particle_<i>" y est créé pour chacune des nrec particules (donc
    nrec sous-dossiers au total -- à nettoyer par l'appelant si besoin,
    pas fait automatiquement ici).

    Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
    nrec-1), pas dans l'ordre de complétion des workers -- important
    pour la reproductibilité de l'ordre des lignes du reftable final.

    max_workers : nombre de process en parallèle (défaut : laissé à
    ProcessPoolExecutor, généralement le nombre de cœurs disponibles).
    """
    reference_directory = Path(reference_directory)
    general_binary_path = Path(general_binary_path)
    base_work_directory = Path(base_work_directory)

    results_by_index: dict[int, ParticleResult] = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle,
                particle_index,
                reference_directory,
                scenario_index,
                num_loci,
                general_binary_path,
                base_work_directory,
                stats_filter,
            ): particle_index
            for particle_index in range(nrec)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()

    return [results_by_index[i] for i in range(nrec)]


def write_reftable_bin(
    results: list[ParticleResult],
    priors: list,
    scenario,
    output_path: str | Path,
) -> None:
    """Écrit un reftable.bin au format binaire DIYABC (vérifié contre
    reftable.cpp, readReftable.R, et abcranger/readreftable.cpp -- voir
    docs/synthese_diyabc_msprime.docx section 5).

    Limité à un SEUL scénario actif (cohérent avec ce POC, qui ne traite
    que le scénario 1 de human) : nscen=1, donc tous les résultats
    doivent partager le même scenario_index -- vérifié, lève ValueError
    sinon.

    Filtre les colonnes de paramètres sur DEUX critères, dans cet ordre :
    1. is_constant_prior : exclut les priors quasi-dégénérés (comme
       readReftable.R / abcranger)
    2. get_parameter_names_used_by_scenario(scenario) : exclut les
       priors non référencés par CE scénario précis -- correction d'un
       bug découvert empiriquement (readReftable.R levait "indice hors
       limites" : notre code gardait les 21 priors du header.txt entier,
       alors que le scénario 1 n'en référence que 16 -- voir notes/
       exploration.md).

    Ne gère PAS les paramètres de mutation (absents de human) -- à
    ajouter (toujours en dernière position, après les paramètres
    démographiques -- voir readReftable.R) si un dataset avec
    microsatellites/séquences est traité plus tard.
    """
    if not results:
        raise ValueError("results est vide : au moins une particule est requise")

    scenario_indices = {r.scenario_index for r in results}
    if len(scenario_indices) != 1:
        raise NotImplementedError(
            f"write_reftable_bin ne gère qu'un seul scénario actif à la "
            f"fois -- scénarios trouvés dans results : {scenario_indices}"
        )
    scenario_index = scenario_indices.pop()

    used_param_names = get_parameter_names_used_by_scenario(scenario)
    kept_param_names = [
        p.name
        for p in priors
        if not is_constant_prior(p) and p.name in used_param_names
    ]
    stat_names = sorted(results[0].summary_statistics.keys())

    nrec = len(results)
    nscen = 1
    nrecscen = [nrec]
    nparam = [len(kept_param_names)]
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
            f.write(struct.pack("<i", scenario_index))
            for name in kept_param_names:
                f.write(struct.pack("<f", result.parameter_values[name]))
            for name in stat_names:
                f.write(struct.pack("<f", result.summary_statistics[name]))
