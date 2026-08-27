"""
Parseur de la section 'loci description' de header.txt. Gère le format
condenséobservé dans human et toy_example5 pour les snp: ("5000 <A> G1 from
1", ou "70 <A> 10 <X> 10 <M> 10 <Y> G1 from 1"), ainsi que le format détaillé
un-locus-par-ligne observé dans sequences-mut ("Lep04 <A> [M] G1 2 40"),
qui est un format différent (header.cpp distingue ces deux cas selon
dataobs.filetype).
"""

import re

from bridge.header_dataclasses import LociDescription, LociDescriptionDetailed

_SECTION_HEADER_RE = re.compile(r"^loci description\s*\((\d+)\)\s*$")

# "5000 <A> G1 from 1" -> total_loci=5000, heritage="<A>", group="G1", start=1
_CONDENSED_SINGLE_TYPE_RE = re.compile(r"^(\d+)\s+<[AHXYM]>\s+(\S+)\s+from\s+(\d+)\s*$")

# Locus_M_A_1_ <A> [M] G1 2 40
_DETAILED_TYPE_RE = re.compile(
    r"^(\S+)\s+<([AHXYM])>\s+\[([MS])\]\s+(\S+)\s+(\d+)(?:\s+(\d+))?\s*$"
)

_LOCI_PAIR_RE = re.compile(r"(\d+)\s*<([^>]+)>")
_LOCI_TRAILER_RE = re.compile(r"(\S+)\s+from\s+(\d+)\s*$")


def _find_loci_description_section_index(lines: list[str]) -> tuple[int, int]:
    """Repère l'index de la ligne d'en-tête 'loci description (N)' --
    factorisé entre parse_loci_description et rewrite_loci_count, qui en
    ont toutes deux besoin.

    Retourne l'index de la locus description
    Lève ValueError si la section n'est trouvée nulle part.
    """
    section_index = next(
        (i for i, line in enumerate(lines) if _SECTION_HEADER_RE.match(line.strip())),
        None,
    )

    if section_index is None:
        raise ValueError("Section 'loci description' non trouvée dans header.txt")

    return section_index


def _extract_loci_info_condensed(text: str) -> tuple[dict[str, int], str, str]:
    """Découpe une ligne de contenu condensée en (total_loci, groupe,
    indice_1based) -- une ou plusieurs paires '<n> <type>' suivies de
    '<groupe> from <indice>' (single-type ou multi-type indifféremment).

    Lève NotImplementedError si le format ne correspond pas.
    """
    paires = _LOCI_PAIR_RE.findall(text)
    if not paires:
        raise NotImplementedError(
            f"Format de ligne non géré (aucune paire '<n> <type>' trouvée) : {text!r}"
        )
    total_loci = {heritage: int(count) for count, heritage in paires}

    reste = _LOCI_PAIR_RE.sub("", text).strip()
    trailer_match = _LOCI_TRAILER_RE.match(reste)
    if trailer_match is None:
        raise NotImplementedError(
            f"Format de ligne non géré (attendu '<groupe> from <indice>' "
            f"après les paires '<n> <type>', trouvé {reste!r}) : {text!r}"
        )
    group, start_1based = trailer_match.groups()
    return total_loci, group, start_1based


def parse_loci_description(
    header_text: str,
) -> LociDescription | list[LociDescriptionDetailed]:
    """Extrait la description des loci à partir de header.txt, pour le
    format condensé (single-type comme human, ou multi-type comme
    toy_example5) et le format description locus par locus
    détaillé observé dans sequences-mut.

    Lève NotImplementedError si le format détecté n'est pas d'un de ces type là.
    """
    lines = header_text.splitlines()
    section_index = _find_loci_description_section_index(lines)
    num_lines = int(_SECTION_HEADER_RE.match(lines[section_index].strip()).group(1))

    if num_lines != 1:
        list_loci = []
        for i in range(section_index + 1, section_index + 1 + num_lines):
            content_line = lines[i].strip()
            match = _DETAILED_TYPE_RE.match(content_line)
            if match is None:
                raise NotImplementedError(
                    f"Format de ligne non géré par parse_loci_description : {content_line!r}"
                )
            if match.groups()[2] == "M":
                name, heritage, ms_or_seq, group, motif_size, motif_range = (
                    match.groups()
                )
                list_loci.append(
                    LociDescriptionDetailed(
                        name=name,
                        heritage=heritage,
                        ms_or_seq=ms_or_seq,
                        group=group,
                        motif_size=int(motif_size) if motif_size else None,
                        motif_range=int(motif_range) if motif_range else None,
                        dnalength=None,
                    )
                )
            else:
                name, heritage, ms_or_seq, group, dnalength, _ = match.groups()
                list_loci.append(
                    LociDescriptionDetailed(
                        name=name,
                        heritage=heritage,
                        ms_or_seq=ms_or_seq,
                        group=group,
                        motif_size=None,
                        motif_range=None,
                        dnalength=int(dnalength) if dnalength else None,
                    )
                )
        return list_loci
    else:
        content_line = lines[section_index + 1].strip()
        total_loci, group, start_1based = _extract_loci_info_condensed(content_line)

        return LociDescription(
            total_loci=total_loci,
            group=group,
            start_index=int(start_1based)
            - 1,  # conversion 1-based -> 0-based, comme prem = N-1 en C++
        )


def rewrite_loci_count(header_text: str, new_total_loci: int) -> str:
    """Retourne une copie de header_text où le nombre de loci déclaré dans
    'loci description' est remplacé par new_total_loci -- nécessaire pour
    tester avec un nombre de loci réduit sans avoir à maintenir un
    header.txt séparé à la main.

    Contrairement à parse_loci_description, limité au format condensé à
    un SEUL type d'héritage (lève NotImplementedError sur une ligne
    multi-type comme celle de toy_example5) -- pas encore mis à jour
    pour ce cas, non nécessaire pour human.
    """
    lines = header_text.splitlines()
    section_index = _find_loci_description_section_index(lines)

    num_lines = int(_SECTION_HEADER_RE.match(lines[section_index].strip()).group(1))

    if num_lines != 1:
        raise NotImplementedError(
            f"Format détaillé (loci description ({num_lines}), {num_lines} "
            f"lignes attendues) non géré par ce parser -- limité au format "
            f"condensé à une seule ligne (cas de human/toy_example5/toy_example3)."
        )

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
