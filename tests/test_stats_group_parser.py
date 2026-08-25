"""Vérifie stats_group_parser : extraction des noms de statistiques
demandées dans la section 'group summary statistics' de header.txt,
vocabulaire ancien (obsolète) et moderne."""

import pytest

from bridge.stats_group_parser import (
    _split_stats_blocks,
    parse_requested_statistic_names,
)


def test_split_stats_blocks(header_text_te2):
    """Vérifie que split_stats_blocks() découpe correctement la section
    'group summary statistics' en blocs de lignes, un bloc par groupe
    de statistiques (ex: "group G1 (N)"), et que le dernier bloc est
    bien limité à la fin de la section (avant la ligne vide ou le
    début d'une autre section)."""
    blocks = _split_stats_blocks(header_text_te2)

    assert len(blocks) == 3

    # Vérifie que le premier bloc contient bien les lignes du groupe G1
    assert len(blocks[0].splitlines()) == 12
    assert blocks[0].startswith("group G1 (16)")


def test_parse_requested_statistic_names_human_old_vocabulary(header_text):
    """Le parseur est purement syntaxique -- il fonctionne aussi sur le
    vocabulaire ANCIEN (obsolète) de human/header.txt (HP0/HM1/HV1/HMO/
    FP0.../AP0...), même si ces noms ne correspondent à aucune
    statistique calculée par summary_statistics.py (incohérence connue
    de ce header.txt, voir notes/exploration.md)."""
    names = parse_requested_statistic_names(header_text)

    assert len(names) == 112
    assert names[:5] == ["HP0_1", "HP0_2", "HP0_3", "HP0_4", "HM1_1"]
    assert names[-1] == "AMO_4.1.2"


def test_parse_requested_statistic_names_modern_vocabulary(header_text_te2):
    """Format condensé moderne (ex: toy_example5_modif ou toy_example2) : noms de
    colonnes 'STAT_index' ou 'STAT_group_index' selon le nombre de groupes, même convention que summary_statistics.py."""
    header_text = (
        "group summary statistics (9)\n"
        "group G1 (9)\n"
        "ML1p 1 2 3\n"
        "ML2p 1.2 1.3 2.3\n"
        "HWm 1 2 3\n"
    )

    names = parse_requested_statistic_names(header_text)

    assert names == [
        "ML1p_1",
        "ML1p_2",
        "ML1p_3",
        "ML2p_1.2",
        "ML2p_1.3",
        "ML2p_2.3",
        "HWm_1",
        "HWm_2",
        "HWm_3",
    ]

    # On test si il y a plusieurs groupes, les noms de colonnes sont "STAT_group_index"
    names_te2 = parse_requested_statistic_names(header_text_te2)
    assert names_te2[:2] == ["NAL_1_1", "NAL_1_2"]


def test_parse_requested_statistic_names_missing_section():
    with pytest.raises(ValueError, match="group summary statistics"):
        parse_requested_statistic_names("pas de section ici\n")


def test_parse_requested_statistic_names_count_mismatch():
    header_text = "group summary statistics (5)\ngroup G1 (5)\nML1p 1 2 3\n"

    with pytest.raises(ValueError, match="annonce 5"):
        parse_requested_statistic_names(header_text)
