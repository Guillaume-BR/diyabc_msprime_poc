"""
Lecture des données observées DIYABC (fichiers .snp, format individu par
ligne) -- pour l'instant, uniquement le comptage du nombre d'individus par
population, nécessaire pour savoir combien d'échantillons demander à
msprime.sim_ancestry().

Référence : src-JMC-C++/data.cpp (détection du format "IND SEX POP").
Ce module ne lit PAS les génotypes eux-mêmes : on simule des données
artificielles avec msprime, on ne réutilise jamais les données observées
réelles dans le pipeline de simulation.
"""

from collections import Counter
from pathlib import Path


def count_samples_per_population(snp_file_path: str | Path) -> dict[str, int]:
    """Compte le nombre d'individus par population dans un fichier .snp
    DIYABC au format 'IND SEX POP <génotypes...>'.

    Ex: pour human_snp_all22chr_maf5.snp -> {"ASW": 30, "YRI": 30, ...}

    L'en-tête 'IND SEX POP' peut être précédé ou non d'un commentaire libre
    en première ligne (comportement observé dans data.cpp, qui teste les
    deux cas) : on recherche son index plutôt que de supposer sa position,
    pour ne perdre aucune ligne de données quel que soit le cas.

    Lève ValueError si l'en-tête n'est trouvé dans aucune des deux
    premières lignes.
    """
    path = Path(snp_file_path)
    lines = path.read_text().splitlines()

    header_index = next(
        (i for i in range(min(2, len(lines))) if lines[i].split()[:3] == ["IND", "SEX", "POP"]),
        None,
    )
    if header_index is None:
        raise ValueError(
            f"En-tête 'IND SEX POP' non trouvé dans les deux premières "
            f"lignes de {path}. Lignes lues : {lines[:2]!r}"
        )

    pop_index = lines[header_index].split().index("POP")

    counts = Counter()
    for line in lines[header_index + 1:]:
        fields = line.split()
        if not fields:
            continue
        counts[fields[pop_index]] += 1

    return dict(counts)