"""
Parseur de la section 'group summary statistics' de header.txt : la
liste des statistiques que DIYABC calcule réellement pour ce dataset.

Sert à filtrer summary_statistics.compute_all_statistics (qui calcule
TOUTES les statistiques implémentées) pour ne garder que celles
réellement déclarées par header.txt -- sinon reftable_msprime.txt/.bin
a des colonnes en trop (ou en moins, pour un header.txt au vocabulaire
obsolète, voir notes/exploration.md) par rapport au vrai reftable
DIYABC, ce qui fausse toute comparaison colonne-par-nom entre les deux
pipelines (découvert empiriquement sur toy_example5_modif : 'ML3p_1.2.3'
présent côté msprime, absent côté DIYABC).

Format géré, un seul groupe (suffisant pour human/toy_example5) :
    group summary statistics (N)
    group G1 (N)
    ML1p 1 2 3
    ML2p 1.2 1.3 2.3
    ...
-> noms de colonnes "STAT_index" (ex: "ML1p_1", "ML2p_1.2"), même
convention que summary_statistics.py.
"""

import re

_SECTION_HEADER_RE = re.compile(r"^group summary statistics\s*\((\d+)\)\s*$")
_GROUP_LINE_RE = re.compile(r"^group\s+(\S+)\s*\((\d+)\)\s*$")


def parse_requested_statistic_names(header_text: str) -> list[str]:
    """Extrait, dans l'ordre de déclaration, les noms de colonnes de
    statistiques attendues par header.txt (section 'group summary
    statistics').

    Limité à un seul groupe de statistiques (ex: "group G1 (N)") --
    lève NotImplementedError si plusieurs groupes sont déclarés (non
    rencontré sur human/toy_example5).
    """
    lines = header_text.splitlines()

    section_index = next(
        (i for i, line in enumerate(lines) if _SECTION_HEADER_RE.match(line.strip())),
        None,
    )
    if section_index is None:
        raise ValueError(
            "Section 'group summary statistics' non trouvée dans header.txt"
        )

    group_line = lines[section_index + 1].strip()
    group_match = _GROUP_LINE_RE.match(group_line)
    if not group_match:
        raise ValueError(f"Ligne de groupe inattendue : {group_line!r}")
    expected_count = int(group_match.group(2))

    content_lines = []
    for line in lines[section_index + 2 :]:
        stripped = line.strip()
        if not stripped:
            break
        if stripped.startswith("group "):
            raise NotImplementedError(
                "Plusieurs groupes de statistiques déclarés -- non géré par "
                "ce parser (un seul groupe attendu, ex: 'group G1 (N)')"
            )
        content_lines.append(stripped)

    names = []
    for stripped in content_lines:
        tokens = stripped.split()
        stat_name, indices = tokens[0], tokens[1:]
        names.extend(f"{stat_name}_{index}" for index in indices)

    if len(names) != expected_count:
        raise ValueError(
            f"'{group_line}' annonce {expected_count} statistiques mais "
            f"{len(names)} ont été trouvées en parsant les lignes suivantes"
        )

    return names
