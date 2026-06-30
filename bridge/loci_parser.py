"""
Parseur de la section 'loci description' de header.txt -- limité au
format CONDENSÉ observé dans human ("5000 <A> G1 from 1"), pas au format
détaillé un-locus-par-ligne observé dans sequences-mut ("Lep04 <A> [M] G1
2 40"), qui est un format différent (header.cpp distingue ces deux cas
selon dataobs.filetype).

Référence : header.cpp::readHeaderLoci, branche "else" (fichier SNP).

Format condensé géré : "<n1> <type1> <groupe> from <indice>"
Un seul type d'héritage par ligne -- le format multi-types
("70 <A> 10 <X> ...") n'est PAS géré (non nécessaire pour human).
"""

import re

from bridge.scenario_types import LociDescription

_SECTION_HEADER_RE = re.compile(r"^loci description\s*\((\d+)\)\s*$")

# "5000 <A> G1 from 1" -> total_loci=5000, heritage="<A>", group="G1", start=1
_CONDENSED_SINGLE_TYPE_RE = re.compile(r"^(\d+)\s+<[AHXYM]>\s+(\S+)\s+from\s+(\d+)\s*$")


def parse_loci_description(header_text: str) -> LociDescription:
    """Extrait la description des loci à partir de header.txt, pour le
    format condensé à un seul type d'héritage (cas de human).

    Lève NotImplementedError si le format détecté est le format détaillé
    (plusieurs lignes, un locus nommé par ligne) ou le format condensé
    multi-types -- non nécessaires pour human, à implémenter si on
    généralise à un autre dataset.
    """
    lines = header_text.splitlines()

    section_index = next(
        (i for i, line in enumerate(lines) if _SECTION_HEADER_RE.match(line.strip())),
        None,
    )
    if section_index is None:
        raise ValueError("Section 'loci description' non trouvée dans header.txt")

    num_lines = int(_SECTION_HEADER_RE.match(lines[section_index].strip()).group(1))
    if num_lines != 1:
        raise NotImplementedError(
            f"Format détaillé (loci description ({num_lines}), {num_lines} "
            f"lignes attendues) non géré par ce parser -- limité au format "
            f"condensé à une seule ligne (cas de human)."
        )

    content_line = lines[section_index + 1].strip()
    match = _CONDENSED_SINGLE_TYPE_RE.match(content_line)
    if not match:
        raise NotImplementedError(
            f"Format de ligne non géré (probablement multi-types, ex: "
            f"'70 <A> 10 <X> ...') : {content_line!r}"
        )

    total_loci, group, start_1based = match.groups()
    return LociDescription(
        total_loci=int(total_loci),
        group=group,
        start_index=int(start_1based)
        - 1,  # conversion 1-based -> 0-based, comme prem = N-1 en C++
    )


def rewrite_loci_count(header_text: str, new_total_loci: int) -> str:
    """Retourne une copie de header_text où le nombre de loci déclaré dans
    'loci description' est remplacé par new_total_loci -- nécessaire pour
    tester avec un nombre de loci réduit sans avoir à maintenir un
    header.txt séparé à la main.

    Limité au même format condensé à un seul type que parse_loci_description
    (lève NotImplementedError dans les mêmes cas).
    """
    lines = header_text.splitlines()

    section_index = next(
        (i for i, line in enumerate(lines) if _SECTION_HEADER_RE.match(line.strip())),
        None,
    )
    if section_index is None:
        raise ValueError("Section 'loci description' non trouvée dans header.txt")

    content_index = section_index + 1
    content_line = lines[content_index].strip()
    match = _CONDENSED_SINGLE_TYPE_RE.match(content_line)
    if not match:
        raise NotImplementedError(
            f"Format de ligne non géré par rewrite_loci_count : {content_line!r}"
        )

    _, group, start_1based = match.groups()
    heritage_match = re.search(r"<[AHXYM]>", content_line)
    heritage = heritage_match.group(0)

    lines[content_index] = f"{new_total_loci} {heritage} {group} from {start_1based}"
    return "\n".join(lines)
