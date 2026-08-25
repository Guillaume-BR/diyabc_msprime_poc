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

Formats gérés :
    - un seul groupe (suffisant pour human/toy_example5) :
    group summary statistics (N)
    group G1 (N)
    ML1p 1 2 3
    ML2p 1.2 1.3 2.3
    ...
-> noms de colonnes "STAT_index" (ex: "ML1p_1", "ML2p_1.2"), même
convention que summary_statistics.py.

    - ou plusieurs groupes (ex: toy_example2) :
    group summary statistics (N)
    group G1 (N1)
    ML1p 1 2 3
    ML2p 1.2 1.3 2.3
    ...
    group G2 (N2)
    ML1p 1 2 3
    ML2p 1.2 1.3 2.3
    ...
-> noms de colonnes "STAT_group_index" (ex: "ML1p_1_1", "ML2p_1_1.2"), même
convention que summary_statistics.py.

"""

import re

_SECTION_HEADER_RE = re.compile(r"^group summary statistics\s*\((\d+)\)\s*$")
_GROUP_LINE_RE = re.compile(r"^group\s+(\S+)\s*\((\d+)\)\s*$")


def _split_stats_blocks(header_text: str) -> list[str]:
    """Découpe le texte complet de header.txt en blocs bruts, un par
    groupe, chaque bloc commençant par sa ligne d'en-tête
    'group G1 (N)' et s'arrêtant juste avant le bloc suivant
    (ou la fin de la section, ex: 'scenario')."""

    lines = header_text.splitlines()
    # repère la ligne d'index où démarre la section "group summary statistics (N)".
    section_start_index = next(
        (i for i, line in enumerate(lines) if _SECTION_HEADER_RE.match(line.strip())),
        None,
    )
    if section_start_index is None:
        raise ValueError(
            "Section 'group summary statistics' non trouvée dans header.txt"
        )

    # Repère les lignes d'index où démarre chaque groupe de statistiques (ex: "group G1 (N)").
    start_indices = [
        i
        for i, line in enumerate(lines)
        if _GROUP_LINE_RE.match(line.strip()) and i > section_start_index
    ]
    if not start_indices:
        raise ValueError("Aucun bloc 'group Gx (N)' trouvé dans le texte fourni")

    # Borne de fin pour le tout dernier groupe : la dernière ligne du fichier "scenario N1 N2 ...",
    # si elle est présente, sinon la fin du fichier.
    last_line_section = next(
        (
            i
            for i, line in enumerate(lines)
            if line.strip().startswith("scenario") and i > start_indices[-1]
        ),
        len(lines),
    )

    blocks = []
    for k, start in enumerate(start_indices):
        end = start_indices[k + 1] if k + 1 < len(start_indices) else last_line_section
        block = "\n".join(lines[start:end]).strip()
        blocks.append(block)
    return blocks


def parse_requested_statistic_names(header_text: str) -> list[str]:
    """Extrait, dans l'ordre de déclaration, les noms de colonnes de
    statistiques attendues par header.txt (section 'group summary
    statistics').
    """
    blocks = _split_stats_blocks(header_text)
    nb_blocks = len(blocks)

    stats_blocks = {}
    for block in blocks:
        group_line = block.splitlines()[0].strip()
        group_match = _GROUP_LINE_RE.match(group_line)
        if not group_match:
            raise ValueError(f"Ligne de groupe inattendue : {group_line!r}")
        numero_group = group_match.group(1)[1:]  # ex: "G1" -> "1"

        expected_count = int(group_match.group(2))

        content_lines = []
        for line in block.splitlines()[1:]:
            stripped = line.strip()
            if not stripped:
                break
            content_lines.append(stripped)

        names = []
        if nb_blocks == 1:
            for stripped in content_lines:
                tokens = stripped.split()
                stat_name, indices = tokens[0], tokens[1:]
                names.extend(f"{stat_name}_{index}" for index in indices)
        elif nb_blocks > 1:
            for stripped in content_lines:
                tokens = stripped.split()
                stat_name, indices = tokens[0], tokens[1:]
                names.extend(f"{stat_name}_{numero_group}_{index}" for index in indices)

        if len(names) != expected_count:
            raise ValueError(
                f"'{group_line}' annonce {expected_count} statistiques mais "
                f"{len(names)} ont été trouvées en parsant les lignes suivantes"
            )

        stats_blocks[group_match.group(1)] = names

    return [name for group_names in stats_blocks.values() for name in group_names]
