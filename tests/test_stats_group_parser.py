"""Vérifie stats_group_parser : extraction des noms de statistiques
demandées dans la section 'group summary statistics' de header.txt,
vocabulaire ancien (obsolète) et moderne."""

import pytest

from bridge.stats_group_parser import parse_requested_statistic_names


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


def test_parse_requested_statistic_names_modern_vocabulary():
    """Format condensé moderne (ex: toy_example5_modif) : noms de
    colonnes 'STAT_index', même convention que summary_statistics.py."""
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


def test_parse_requested_statistic_names_missing_section():
    with pytest.raises(ValueError, match="group summary statistics"):
        parse_requested_statistic_names("pas de section ici\n")


def test_parse_requested_statistic_names_count_mismatch():
    header_text = "group summary statistics (5)\ngroup G1 (5)\nML1p 1 2 3\n"

    with pytest.raises(ValueError, match="annonce 5"):
        parse_requested_statistic_names(header_text)
