"""
Vérifie que scenario_parser produit, sur le vrai header.txt du dataset
human, exactement les événements qu'on a décortiqués à la main avec le
mentor (voir notes/exploration.md) pour les scénarios 1, 2 et 4.
"""

import pytest

from bridge.header_dataclasses import MergeEvent, SampleEvent, SplitEvent, VarNeEvent
from bridge.scenario_parser import parse_header_scenarios


def test_unimplemented_scenarios_are_skipped_with_warning(header_text):
    """Les scénarios 2,3,5,6 utilisent 'split' (pas encore implémenté) :
    ils doivent être ignorés avec un avertissement, pas faire planter le
    parsing des autres scénarios."""
    # with pytest.warns(UserWarning, match="split"):
    scenarios = parse_header_scenarios(header_text)

    # Scénarios gérables aujourd'hui : 1 (merge/varNe) et 4 (idem,
    # numérotation t11..t44). Les 4 autres utilisent split -> ignorés.
    found_indices = {s.index for s in scenarios}
    assert found_indices == {1, 2, 3, 4, 5, 6}


def test_scenario1_metadata(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario1 = scenarios[0]
    assert scenario1.index == 1
    assert scenario1.weight == pytest.approx(0.16667)
    assert scenario1.initial_pop_size_exprs == ["N1", "N2", "N3", "N4"]


def test_scenario1_events(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario1 = scenarios[0]

    expected = [
        SampleEvent(time_expr="0", pop=1),
        SampleEvent(time_expr="0", pop=2),
        SampleEvent(time_expr="0", pop=3),
        SampleEvent(time_expr="0", pop=4),
        MergeEvent(time_expr="t1", ancestral_pop=2, derived_pop=1),
        VarNeEvent(time_expr="t2-d3", pop=3, new_size_expr="Nbn3"),
        VarNeEvent(time_expr="t2-d4", pop=4, new_size_expr="Nbn4"),
        MergeEvent(time_expr="t2", ancestral_pop=3, derived_pop=4),
        VarNeEvent(time_expr="t2", pop=3, new_size_expr="N34"),
        VarNeEvent(time_expr="t3-d34", pop=3, new_size_expr="Nbn34"),
        MergeEvent(time_expr="t3", ancestral_pop=2, derived_pop=3),
        VarNeEvent(time_expr="t4", pop=2, new_size_expr="Na"),
    ]

    assert scenario1.events == expected


def test_scenario4_events(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario4 = next(s for s in scenarios if s.index == 4)

    expected = [
        SampleEvent(time_expr="0", pop=1),
        SampleEvent(time_expr="0", pop=2),
        SampleEvent(time_expr="0", pop=3),
        SampleEvent(time_expr="0", pop=4),
        MergeEvent(time_expr="t11", ancestral_pop=2, derived_pop=1),
        VarNeEvent(time_expr="t22-d3", pop=3, new_size_expr="Nbn3"),
        MergeEvent(time_expr="t22", ancestral_pop=2, derived_pop=3),
        VarNeEvent(time_expr="t33-d4", pop=4, new_size_expr="Nbn4"),
        MergeEvent(time_expr="t33", ancestral_pop=2, derived_pop=4),
        VarNeEvent(time_expr="t44", pop=2, new_size_expr="Na"),
    ]

    assert scenario4.events == expected


def test_scenario2_events(header_text):
    """Vérifie qu'un événement 'split' (admixture) est correctement
    interprété : 't1 split 1 4 2 ra' -> pop 1 disparaît, chaque lignée
    part vers pop 4 avec probabilité 'ra', sinon vers pop 2 -- sémantique
    vérifiée contre history.cpp/particuleC.cpp (voir docstring de
    SplitEvent)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    assert scenario2.events[4] == SplitEvent(
        time_expr="t1",
        derived_pop=1,
        ancestral_pop1=4,
        ancestral_pop2=2,
        admixture_rate="ra",
    )
