"""
Écrit un fichier .snp DIYABC (format IND SEX POP <génotypes>) à partir de
génotypes simulés par msprime, pour pouvoir déléguer le calcul des
statistiques résumées (FST1, ML1, F3, F4, NEI, AML...) au binaire C++
existant (HeaderC::calstatobs), plutôt que de réimplémenter ces formules
en Python -- décision prise après une tentative ratée de réimplémentation
de FST1 (Fst négatif obtenu, voir notes/exploration.md), et au vu de
l'ampleur du corpus de statistiques à reproduire (doc DIYABC section 2.6.3).

Transformation effectuée : genotypes_per_locus contient, pour chaque
locus, des génotypes HAPLOÏDES (une valeur 0/1 par lignée génomique),
regroupés par population (voir ancestry_simulation.simulate_snp_genotypes).
Ce module les agrège par PAIRES DE LIGNÉES CONSÉCUTIVES en génotypes
DIPLOÏDES (0/1/2), conformément au format réel du fichier .snp (vérifié
empiriquement : msprime associe les lignées [2i, 2i+1] au même individu
diploïde i -- voir notes/exploration.md).
"""

from pathlib import Path


def _genotypes_to_diploid(haploid_genotypes: list[int]) -> list[int]:
    """Agrège une liste de génotypes haploïdes (une valeur par lignée) en
    génotypes diploïdes (0/1/2), en sommant les paires de lignées
    consécutives [2i, 2i+1] -- correspondant aux deux copies
    chromosomiques d'un même individu (vérifié empiriquement avec
    ts.individuals()[i].nodes).

    Lève ValueError si le nombre de lignées est impair (incohérent avec
    une simulation en ploidy=2).
    """
    if len(haploid_genotypes) % 2 != 0:
        raise ValueError(
            f"Nombre de lignées impair ({len(haploid_genotypes)}) -- "
            f"incohérent avec une simulation diploïde (ploidy=2)."
        )
    return [
        haploid_genotypes[2 * i] + haploid_genotypes[2 * i + 1]
        for i in range(len(haploid_genotypes) // 2)
    ]


def write_snp_file(
    genotypes_per_locus: list[dict[str, list[int]]],
    output_path: str | Path,
) -> None:
    """Écrit un fichier .snp DIYABC à partir de num_loci dicts
    {nom_population: [génotypes haploïdes...]}, un par locus (la forme
    produite par ancestry_simulation.simulate_snp_genotypes).

    Le nom de chaque individu simulé est généré comme "sim_<pop>_<n>"
    (ex: "sim_pop1_1", "sim_pop1_2"...). La colonne SEX est fixée à "9"
    pour tous les individus -- valeur arbitraire, non confirmée comme
    sans impact pour des loci autosomaux <A> (voir notes/exploration.md
    pour la justification de cette hypothèse).

    genotypes_per_locus doit contenir AU MOINS un locus, et toutes les
    populations doivent être présentes et avoir le même nombre de lignées
    à chaque locus (cohérence vérifiée par la simulation elle-même, pas
    revérifiée ici).
    """
    if not genotypes_per_locus:
        raise ValueError("genotypes_per_locus est vide : au moins un locus est requis")

    population_names = list(genotypes_per_locus[0].keys())

    # Construit, pour chaque population, la matrice diploïde
    # [individu][locus] -- nécessaire car le fichier organise les données
    # par individu (toutes ses lignes de loci sur une seule ligne).
    diploid_matrix_per_population: dict[str, list[list[int]]] = {}
    for pop_name in population_names:
        per_locus_diploid = [
            _genotypes_to_diploid(locus_genotypes[pop_name])
            for locus_genotypes in genotypes_per_locus
        ]
        # Transpose : de [locus][individu] vers [individu][locus]
        num_individuals = len(per_locus_diploid[0])
        diploid_matrix_per_population[pop_name] = [
            [per_locus_diploid[loc][ind] for loc in range(len(per_locus_diploid))]
            for ind in range(num_individuals)
        ]

    num_loci = len(genotypes_per_locus)
    header_cols = ["IND", "SEX", "POP"] + ["A"] * num_loci

    lines = [" ".join(header_cols)]
    for pop_name, individuals_matrix in diploid_matrix_per_population.items():
        for ind_index, diploid_genotypes in enumerate(individuals_matrix, start=1):
            ind_name = f"sim_{pop_name}_{ind_index}"
            row = [ind_name, "9", pop_name] + [str(g) for g in diploid_genotypes]
            lines.append(" ".join(row))

    Path(output_path).write_text("\n".join(lines) + "\n")