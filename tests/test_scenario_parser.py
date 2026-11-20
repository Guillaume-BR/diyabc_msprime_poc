"""Vérifie que scenario_parser produit, sur le vrai header.txt du dataset
human, exactement les événements qu'on a décortiqués à la main avec le mentor
(voir notes/exploration.md) pour les scénarios 1, 2 et 4."""

import pytest

from bridge.header_dataclasses import (
    BottleneckEvent,
    MergeEvent,
    MigrationEvent,
    SampleEvent,
    SplitEvent,
    VarNeEvent,
)
from bridge.scenario_parser import (
    _parse_event_line,
    check_scenario_logic,
    parse_header_scenarios,
    parse_scenario_block,
    split_scenario_blocks,
)
from tests.conftest import REFERENCE_DIR


def test_unimplemented_scenarios_are_skipped_with_warning(header_text):
    """Les scénarios 2,3,5,6 utilisent 'split' (pas encore implémenté) :
    ils doivent être ignorés avec un avertissement, pas faire planter le
    parsing des autres scénarios.
    """
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
    """Vérifie qu'un événement 'split' (admixture) est correctement interprété
    : 't1 split 1 4 2 ra' -> pop 1 disparaît, chaque lignée part vers pop 4
    avec probabilité 'ra', sinon vers pop 2 -- sémantique vérifiée contre
    history.cpp/particuleC.cpp (voir docstring de SplitEvent)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    assert scenario2.events[4] == SplitEvent(
        time_expr="t1",
        derived_pop=1,
        ancestral_pops=[4, 2],
        admixture_rates=["ra"],
    )


def test_parse_event_line_three_population_split():
    line = "t1 split 5 1 2 3 r1 r2"
    event = _parse_event_line(line)
    assert isinstance(event, SplitEvent)
    assert event.time_expr == "t1"
    assert event.derived_pop == 5
    assert event.ancestral_pops == [1, 2, 3]
    assert event.admixture_rates == ["r1", "r2"]


def test_parse_event_line_wrong_format():
    line = "t1 split 5 1 2 3 r1"
    with pytest.raises(ValueError, match="SplitEvent mal formé"):
        _parse_event_line(line)

    line = "t1 split 5 1"
    with pytest.raises(
        ValueError,
        match="Au moins deux populations ancestrales sont attendues pour un événement split",
    ):
        _parse_event_line(line)

    line = line = "t1 split 5 a 2 3 r1 r2"
    with pytest.raises(ValueError, match="SplitEvent mal formé"):
        _parse_event_line(line)


def test_parse_event_line_with_mig():
    line1 = "t1 mig 0 1 m01"
    event = _parse_event_line(line1)
    assert isinstance(event, MigrationEvent)
    assert event.dest == 1
    assert event.source == 0
    assert not event.symmetric
    assert event.migration_rate == "m01"

    line2 = "t1 smig 0 1 m01"
    event = _parse_event_line(line2)
    assert isinstance(event, MigrationEvent)
    assert event.symmetric


def test_parse_event_line_with_wrong_mig():
    line1 = "t1 mig 0 1 2 m01"
    with pytest.raises(
        ValueError, match="MigrationEvent mal formé : 't1 mig 0 1 2 m01'"
    ):
        _parse_event_line(line1)

    line2 = "t1 mig pop0 1 m01"
    with pytest.raises(ValueError, match="Les populations concernées"):
        _parse_event_line(line2)


def test_parse_event_line_with_btn():
    line1 = "t1 bn 0 pbn"
    event = _parse_event_line(line1)
    assert isinstance(event, BottleneckEvent)
    assert event.pop == 0
    assert not event.instantaneous
    assert event.strength is None
    assert event.proportion == "pbn"

    line2 = "t1 ibn 0 sbn"
    event = _parse_event_line(line2)
    assert isinstance(event, BottleneckEvent)
    assert event.pop == 0
    assert event.instantaneous
    assert event.strength == "sbn"
    assert event.proportion is None


def test_parse_event_line_with_wrong_btn():
    line1 = "t1 bn 0 1 pbn"
    with pytest.raises(ValueError, match="BottleneckEvent mal formé : 't1 bn 0 1 pbn'"):
        _parse_event_line(line1)

    line2 = "t1 bn pop0 pbn"
    with pytest.raises(
        ValueError, match="La population concernée doit être un entier."
    ):
        _parse_event_line(line2)


def test_parse_header_scenarios(header_text_te1):
    results = parse_header_scenarios(header_text_te1)
    assert len(results) == 2


def test_parse_header_scenario_with_wrong_key_word():
    header_text = """
    human_snp_all22chr_maf5.snp
    21 parameters and 130 summary statistics

    3 scenarios: 2 2 2 
    scenario 1 [0.16667] (1)
    N1
    0 sample 1
    scenario 2 [0.16667] (1)
    N1
    0 toto 1
    scenario 3 [0.16667] (1)
    N1
    0 sample 1
    """

    with pytest.warns(UserWarning, match="Bloc"):
        parse_header_scenarios(header_text)


def test_parse_header_scenarios_with_wrong_split_format(header_text):
    # Ajout d'un scénario avec un split mal formé
    header_text_with_wrong_split = header_text.replace(
        "t11 split 1 3 2 ra", "t11 split 1 3 2"
    )
    with pytest.raises(ValueError, match="SplitEvent mal formé"):
        parse_header_scenarios(header_text_with_wrong_split)


def test_split_scenario_blocks(header_text_te1):
    """Vérifie que le parser de scénarios découpe correctement le bloc header
    en sous-blocs par scénario, même si les scénarios sont séparés par des
    lignes vides ou des commentaires."""
    results = split_scenario_blocks(header_text_te1)
    assert len(results) == 2


def test_check_scenario_logic():
    for dir in REFERENCE_DIR.iterdir():
        if dir.is_dir():
            path_header = REFERENCE_DIR / dir / "headerRF.txt"
            header_text = path_header.read_text()
            scenarios = parse_header_scenarios(header_text)
            for s in scenarios:
                assert check_scenario_logic(s) is None


def test_check_scenario_logic_with_wrong_merge():
    # Pop dérivée non présentes dans le sampling
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 merge 2 15
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="La population dérivée 15"):
        check_scenario_logic(scenario)

    # pop dérivée = pop_ancestrale
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 merge 2 2
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(
        ValueError, match="Les populations ancestrale et dérivée de l'évènement"
    ):
        check_scenario_logic(scenario)


def test_check_scenario_logic_with_wrong_varne():
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 varne 14 Na
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="La population 14 de l'évènement"):
        check_scenario_logic(scenario)


def test_check_scenario_logic_with_wrong_bottleneck():
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 bn 12 Na
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="La population 12 de l'évènement"):
        check_scenario_logic(scenario)


def test_check_scenario_logic_with_wrong_split():
    # Pop dérivée non échantillonnée
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 split 7 4 2 ra
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="La population dérivée 7"):
        check_scenario_logic(scenario)

    # Pop ancestrales égales
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 split 1 4 4 ra
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="Les populations ancestrales"):
        check_scenario_logic(scenario)


def test_check_scenario_logic_with_wrong_mig():
    # Pop dérivée non échantillonnée
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 mig 3 2 ra
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="Les populations source et dérivé"):
        check_scenario_logic(scenario)

    # Pop ancestrales égales
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        t1 mig 1 1 ra
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="Populations source et dérivé"):
        check_scenario_logic(scenario)


def test_check_scenario_logic_with_wrong_final_ancestral_pop():
    # Pop dérivée non échantillonnée
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 N3 N4
        0 sample 1
        0 sample 2
        0 sample 3
        0 sample 4
        t1 merge 1 2
        t2 merge 3 4
        """
    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="Nombre de populations ancestrales :"):
        check_scenario_logic(scenario)
