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

from bridge.pipeline import compute_summary_statistics


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