"""Vérifie demography_builder : traduction d'un Scenario + valeurs de
paramètres tirées en msprime.Demography (merges, admixture, extraction des noms
de paramètres référencés)."""

import msprime
import numpy as np
import pytest
from conftest import _check_replace_modification
from msprime.demography import InstantaneousBottleneck, SimpleBottleneck
from scipy import stats

from bridge.demography_builder import (
    _founding_time,
    _size_in_force,
    build_demography,
    evaluate_expression,
    extract_referenced_names,
    get_parameter_names_used_by_scenario,
    rescale_demography,
)
from bridge.header_dataclasses import VarNeLogEvent
from bridge.observed_data import coalescence_coefficient
from bridge.scenario_parser import (
    check_scenario_logic,
    parse_header_scenarios,
    parse_scenario_block,
)
from bridge.varnelog import add_varnelog


@pytest.mark.parametrize(
    "scenario_index,founding_event, values, expected",
    [
        (1, "t2 merge 1 2", {"t2": 50}, 50),  # nom
        (2, "50 merge 1 2", {}, 50),  # littéral
        (3, "t2-d3 merge 1 2", {"t2": 50, "d3": 8}, 42),  # expression
        (
            4,
            "t2 split 2 1 3 r\n t3 merge 1 3",
            {"t2": 50, "r": 0.5, "t3": 80},
            50,
        ),  # split
    ],
)
def test_founding_time(scenario_index, founding_event, values, expected):
    block = f"""
        scenario {scenario_index} [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t1 VARNELOG 2 NF r 
        {founding_event}
        """
    scenario = parse_scenario_block(block)

    # Garde-fou : un scenario invalide fausserait le test
    check_scenario_logic(scenario)

    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    assert _founding_time(scenario, idx, values) == expected


def test_founding_time_with_merge_before_foundation():
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        0 sample 3
        t1 VARNELOG 2 NF r 
        t2 merge 1 3
        t3 merge 1 2
        """

    values = {
        "N1": 1000,
        "N2": 1000,
        "t1": 10,
        "NF": 500,
        "r": 0.2,
        "t2": 50,
        "t3": 100,
    }
    scenario = parse_scenario_block(block)

    # Garde-fou : un scenario invalide fausserait le test
    check_scenario_logic(scenario)

    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    assert _founding_time(scenario, idx, values) == values["t3"]


def test_founding_time_without_foundation():
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t1 VARNELOG 2 NF r 
        t2 merge 2 1
        """

    values = {
        "N1": 1000,
        "N2": 1000,
        "t1": 10,
        "NF": 500,
        "r": 0.2,
        "t2": 50,
    }
    scenario = parse_scenario_block(block)
    # Pas de garde fou check_logic_scenario pour voir si l'erreur de _founding_time est levée
    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    with pytest.raises(ValueError, match="La population 2 n'est pas fondée"):
        _founding_time(scenario, idx, values)


def test_size_in_force():
    block = """
    scenario 1 [0.16667] (6)
    N1 N2 N3 
    0   sample 1
    0   sample 2
    0   sample 3
    t1  merge 1 2        
    t2  VARNELOG 1 NF r  
    t3  merge 3 1
    """
    values = {
        "N1": 1000,
        "N2": 2000,
        "N3": 3000,
        "t1": 10,
        "NF": 500,
        "r": 0.2,
        "t2": 50,
        "t3": 100,
    }
    scenario = parse_scenario_block(block)
    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    assert _size_in_force(scenario, idx, values) == values["N1"]


def test_size_in_force_with_varne_before_varnelog():
    block = """
    scenario 1 [0.16667] (6)
    N1 N2
    0   sample 1
    0   sample 2
    t1  varne 1 Nbn1  
    t1  varne 1 Nbn2     
    t2  VARNELOG 1 NF r  
    t3  merge 2 1
    """
    values = {
        "N1": 1000,
        "N2": 2000,
        "t1": 10,
        "Nbn1": 500,
        "Nbn2": 800,
        "NF": 500,
        "r": 0.2,
        "t2": 50,
        "t3": 100,
    }
    scenario = parse_scenario_block(block)
    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    assert _size_in_force(scenario, idx, values) == values["Nbn2"]


def test_size_in_force_with_missing_value():
    block = """
    scenario 1 [0.16667] (6)
    N1 N2
    0   sample 1
    0   sample 2     
    t2  VARNELOG 1 NF r  
    t3  merge 2 1
    """
    values = {"N2": 2000, "NF": 500, "r": 0.2, "t2": 50, "t3": 100}
    scenario = parse_scenario_block(block)
    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    with pytest.raises(ValueError, match="Expression 'N1' non évaluable"):
        _size_in_force(scenario, idx, values)


def test_size_in_force_with_wrong_population():
    block = """
    scenario 1 [0.16667] (6)
    N1 N2
    0   sample 1
    0   sample 2     
    t2  VARNELOG 3 NF r  
    t3  merge 2 1
    """
    values = {"N2": 2000, "NF": 500, "r": 0.2, "t2": 50, "t3": 100}
    scenario = parse_scenario_block(block)
    idx = next(i for i, e in enumerate(scenario.events) if isinstance(e, VarNeLogEvent))
    with pytest.raises(ValueError, match="La population 3 n'a pas de taille"):
        _size_in_force(scenario, idx, values)


def test_build_demography_scenario2_admixture(header_text):
    """Vérifie que build_demography traduit un SplitEvent en événement msprime
    Admixture avec les bonnes populations et proportions (rate, 1-rate)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "ra": 0.3,
    }

    demography = build_demography(scenario2, values)

    admixtures = [
        event
        for event in demography.events
        if isinstance(event, msprime.demography.Admixture)
    ]
    assert len(admixtures) == 1
    admixture = admixtures[0]
    assert admixture.time == 10
    assert admixture.derived == "pop1"
    assert admixture.ancestral == ["pop4", "pop2"]
    assert admixture.proportions == pytest.approx([0.3, 0.7])


def test_build_demography_scenario2_modified(header_text):
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "t1 split 1 4 2 3 r1 r2"
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "r1": 0.3,
        "r2": 0.2,
    }

    demography = build_demography(scenario2, values)

    admixtures = [
        event
        for event in demography.events
        if isinstance(event, msprime.demography.Admixture)
    ]
    assert len(admixtures) == 1
    admixture = admixtures[0]
    assert admixture.time == 10
    assert admixture.derived == "pop1"
    assert admixture.ancestral == ["pop4", "pop2", "pop3"]
    assert admixture.proportions == pytest.approx([0.3, 0.2, 0.5])


def test_build_demography_with_sum_proportions_greater_than_one(header_text):
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "t1 split 1 4 2 3 r1 r2"
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    # Garde avec

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "r1": 0.3,
        "r2": 0.8,  # somme > 1
    }
    with pytest.raises(ValueError, match="La somme des taux d'admixture"):
        build_demography(scenario2, values)


def test_build_demography_with_negative_proportions(header_text):
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "t1 split 1 4 2 3 r1 r2"
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "r1": -0.3,  # négatif
        "r2": 0.5,
    }
    with pytest.raises(ValueError, match="Un des taux d'admixture est négatif"):
        build_demography(scenario2, values)


def test_build_demography_with_three_ancestral_populations(header_text):
    """Vérifie que build_demography traduit un SplitEvent en événement msprime
    Admixture avec les bonnes populations et proportions (rate1, rate2, 1 - rate1 - rate2)."""
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "t1 split 1 4 2 3 r1 r2"
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "r1": 0.3,
        "r2": 0.2,
    }

    demography = build_demography(scenario2, values)

    samples = {"pop1": 150, "pop2": 150, "pop3": 150, "pop4": 150}
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=1,
        random_seed=42,
        record_migrations=True,
    )

    mig = ts.tables.migrations
    mask = (mig.time == values["t1"]) & (mig.source == demography["pop1"].id)
    dest = mig.dest[mask]
    _, counts = np.unique(
        dest, return_counts=True
    )  # Attention, counts dans l'ordre des populations déclarées dans l'ordre de msprime pop2,pop3,pop4

    # Bon ordre
    _, p_value = stats.chisquare(
        f_obs=counts,
        f_exp=np.array([values["r2"], 1 - values["r1"] - values["r2"], values["r1"]])
        * np.sum(counts),
    )
    assert p_value > 0.001

    # Ordre permuté
    _, p_value = stats.chisquare(
        f_obs=counts,
        f_exp=np.array([values["r1"], 1 - values["r1"] - values["r2"], values["r2"]])
        * np.sum(counts),
    )
    assert p_value < 0.001


def test_demography_builder_with_migration(header_text):
    header_text_modified = _check_replace_modification(
        header_text,
        "t1 split 1 4 2 ra",
        "tm1 mig 1 4 r1\n tm2 mig 1 4 r2\n t1 split 1 4 2 ra",
    )

    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 2000,
        "t2": 5000,
        "tm1": 50,
        "tm2": 1000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "ra": 0.3,
        "r1": 0.2,
        "r2": 0.0,
    }

    demography = build_demography(scenario2, values)
    assert demography.debug().epochs[1].migration_matrix[0, 3] == 0.2
    assert demography.debug().epochs[1].migration_matrix[3, 0] == 0.0
    assert demography.debug().epochs[2].migration_matrix[0, 3] == 0.0


def test_demography_builder_with_bottleneck():
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    values = {
        "N1": 1000,
        "N2": 1000,
        "tbn": 10,
        "p": 0.6,
        "tibn": 20,
        "s": 100,
        "t1": 500,
    }

    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values)
    assert (
        SimpleBottleneck(time=10, population="pop1", proportion=0.6)
        in demography.events
    )
    assert (
        InstantaneousBottleneck(time=20, population="pop2", strength=100)
        in demography.events
    )


def test_demography_builder_with_wrong_simple_bottleneck():
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    values = {
        "N1": 1000,
        "N2": 1000,
        "tbn": 10,
        "p": -0.6,
        "tibn": 20,
        "s": 100,
        "t1": 500,
    }

    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="comprise entre 0 et 1 : -0.6."):
        build_demography(scenario, values)

    # Temps de départ du bottleneck à 0
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    values = {
        "N1": 1000,
        "N2": 1000,
        "tbn": 0,
        "p": 0.6,
        "tibn": 20,
        "s": 100,
        "t1": 500,
    }

    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="Bottleneck silencieux et sans effet"):
        build_demography(scenario, values)


def test_demography_builder_with_wrong_instantaneous_bottleneck():
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    values = {
        "N1": 1000,
        "N2": 1000,
        "tbn": 10,
        "p": 0.6,
        "tibn": 20,
        "s": -100,
        "t1": 500,
    }

    scenario = parse_scenario_block(block)
    with pytest.raises(ValueError, match="instantané doit être positive : -100"):
        build_demography(scenario, values)


@pytest.mark.parametrize("mode", ["harmonic", "exact"])
def test_demography_builder_with_varnelog(mode):
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t1 VARNELOG 1 NF r 
        t2 merge 2 1
        """
    values = {"N1": 1000, "N2": 5000, "t1": 10, "NF": 50, "r": 0.2, "t2": 50}

    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values, varnelog_mode=mode)

    # Démographie construite à la main
    time_start = values["t1"]
    time_found = values["t2"]
    NS = values["N1"]
    NF = values["NF"]
    r = values["r"]
    oracle_demography = msprime.Demography()
    oracle_demography.add_population(
        initial_size=values["N1"], name="pop1", initially_active=True
    )
    oracle_demography.add_population(
        initial_size=values["N2"], name="pop2", initially_active=True
    )
    add_varnelog(
        oracle_demography, "pop1", time_start, time_found, NF, NS, r, mode=mode
    )
    oracle_demography.add_population_split(
        time=time_found, derived=["pop1"], ancestral="pop2"
    )
    oracle_demography.sort_events()

    assert len(demography.debug().epochs) == len(oracle_demography.debug().epochs)
    for e, oracle_e in zip(
        demography.debug().epochs, oracle_demography.debug().epochs, strict=True
    ):
        assert e.start_time == oracle_e.start_time
        assert e.end_time == oracle_e.end_time

        for name in ("pop1", "pop2"):
            pid = demography[name].id
            assert e.populations[pid].start_size == oracle_e.populations[pid].start_size


@pytest.mark.parametrize("mode", ["harmonic", "exact"])
def test_demography_builder_with_varne_before_varnelog(mode):
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t0 varne 1 Nbn0
        t1 varne 1 Nbn1
        t2 VARNELOG 1 NF r 
        t3 merge 2 1
        """
    values = {
        "N1": 1000,
        "N2": 5000,
        "t0": 5,
        "Nbn0": 1000,
        "t1": 10,
        "Nbn1": 500,
        "t2": 30,
        "NF": 50,
        "r": 0.2,
        "t3": 50,
    }

    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values, varnelog_mode=mode)

    # Démographie construite à la main
    time_start = values["t2"]
    time_found = values["t3"]
    NS = values["Nbn1"]
    NF = values["NF"]
    r = values["r"]
    oracle_demography = msprime.Demography()
    oracle_demography.add_population(
        initial_size=values["N1"], name="pop1", initially_active=True
    )
    oracle_demography.add_population(
        initial_size=values["N2"], name="pop2", initially_active=True
    )
    oracle_demography.add_population_parameters_change(
        time=values["t0"], population="pop1", initial_size=values["Nbn0"]
    )
    oracle_demography.add_population_parameters_change(
        time=values["t1"], population="pop1", initial_size=values["Nbn1"]
    )
    add_varnelog(
        oracle_demography, "pop1", time_start, time_found, NF, NS, r, mode=mode
    )
    oracle_demography.add_population_split(
        time=time_found, derived=["pop1"], ancestral="pop2"
    )
    oracle_demography.sort_events()

    assert len(demography.debug().epochs) == len(oracle_demography.debug().epochs)
    assert len(demography.debug().epochs) > 3
    for e, oracle_e in zip(
        demography.debug().epochs, oracle_demography.debug().epochs, strict=True
    ):
        assert e.start_time == oracle_e.start_time
        assert e.end_time == oracle_e.end_time

        for name in ("pop1", "pop2"):
            pid = demography[name].id
            assert e.populations[pid].start_size == oracle_e.populations[pid].start_size


@pytest.mark.parametrize("mode", ["harmonic", "exact"])
def test_demography_builder_with_varne_on_other_pop(mode):
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t0 varne 1 Nbn0
        t1 varne 2 Nbn1
        t2 VARNELOG 1 NF r 
        t3 merge 2 1
        """
    values = {
        "N1": 1000,
        "N2": 5000,
        "t0": 5,
        "Nbn0": 1000,
        "t1": 10,
        "Nbn1": 500,
        "t2": 30,
        "NF": 50,
        "r": 0.2,
        "t3": 50,
    }

    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values, varnelog_mode=mode)

    # Démographie construite à la main
    time_start = values["t2"]
    time_found = values["t3"]
    NS = values["Nbn0"]
    NF = values["NF"]
    r = values["r"]
    oracle_demography = msprime.Demography()
    oracle_demography.add_population(
        initial_size=values["N1"], name="pop1", initially_active=True
    )
    oracle_demography.add_population(
        initial_size=values["N2"], name="pop2", initially_active=True
    )
    oracle_demography.add_population_parameters_change(
        time=values["t0"], population="pop1", initial_size=values["Nbn0"]
    )
    oracle_demography.add_population_parameters_change(
        time=values["t1"], population="pop2", initial_size=values["Nbn1"]
    )
    add_varnelog(
        oracle_demography, "pop1", time_start, time_found, NF, NS, r, mode=mode
    )
    oracle_demography.add_population_split(
        time=time_found, derived=["pop1"], ancestral="pop2"
    )
    oracle_demography.sort_events()

    assert len(demography.debug().epochs) == len(oracle_demography.debug().epochs)
    assert len(demography.debug().epochs) > 3
    for e, oracle_e in zip(
        demography.debug().epochs, oracle_demography.debug().epochs, strict=True
    ):
        assert e.start_time == oracle_e.start_time
        assert e.end_time == oracle_e.end_time

        for name in ("pop1", "pop2"):
            pid = demography[name].id
            assert e.populations[pid].start_size == oracle_e.populations[pid].start_size


@pytest.mark.parametrize("mode", ["harmonic", "exact"])
def test_demography_builder_with_rescale(mode):
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t1 VARNELOG 1 NF r 
        t2 merge 2 1
        """
    values = {"N1": 1000, "N2": 5000, "t1": 10, "NF": 50, "r": 0.2, "t2": 50}

    factor = 0.5
    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values, varnelog_mode=mode)
    rescaled = rescale_demography(demography, factor=factor)
    # Démographie construite à la main
    time_start = values["t1"]
    time_found = values["t2"]
    NS = values["N1"]
    NF = values["NF"]
    r = values["r"]
    oracle_demography = msprime.Demography()
    oracle_demography.add_population(
        initial_size=values["N1"], name="pop1", initially_active=True
    )
    oracle_demography.add_population(
        initial_size=values["N2"], name="pop2", initially_active=True
    )
    add_varnelog(
        oracle_demography, "pop1", time_start, time_found, NF, NS, r, mode=mode
    )
    oracle_demography.add_population_split(
        time=time_found, derived=["pop1"], ancestral="pop2"
    )
    oracle_demography.sort_events()

    assert len(rescaled.debug().epochs) == len(oracle_demography.debug().epochs)
    assert len(rescaled.debug().epochs)
    for e, oracle_e in zip(
        rescaled.debug().epochs, oracle_demography.debug().epochs, strict=True
    ):
        assert e.start_time == oracle_e.start_time
        assert e.end_time == oracle_e.end_time

        for name in ("pop1", "pop2"):
            pid = rescaled[name].id
            assert e.populations[pid].start_size == pytest.approx(
                factor * oracle_e.populations[pid].start_size
            )


def test_demography_builder_with_symmetric_migration(header_text):
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "tm1 smig 1 4 r1\n t1 split 1 4 2 ra"
    )

    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "tm1": 5,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "ra": 0.1,
        "r1": 0.3,
    }

    demography = build_demography(scenario2, values)
    assert demography.debug().epochs[1].migration_matrix[0, 3] == 0.3
    assert demography.debug().epochs[1].migration_matrix[3, 0] == 0.3


def test_demography_builder_with_wrong_migration_rate(header_text):
    header_text_modified = _check_replace_modification(
        header_text, "t1 split 1 4 2 ra", "tm1 smig 1 4 r1\n t1 split 1 4 2 ra"
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "tm1": 5,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "ra": 0.2,
        "r1": -0.3,
    }

    with pytest.raises(
        ValueError, match="Le taux de migration doit être supérieur ou égal à 0"
    ):
        build_demography(scenario2, values)


def test_get_parameter_names_used_by_scenario2(header_text):
    """Vérifie que le taux d'admixture 'ra' est bien inclus dans les paramètres
    référencés par un scénario qui contient un SplitEvent -- sinon il serait
    exclu à tort des colonnes du reftable.bin pour ce scénario (même bug que
    celui déjà corrigé pour t11..t44)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    used_names = get_parameter_names_used_by_scenario(scenario2)

    assert "ra" in used_names
    expected = {
        "N1",
        "N2",
        "N3",
        "N4",
        "t1",
        "ra",
        "t2",
        "d3",
        "Nbn3",
        "d4",
        "Nbn4",
        "N34",
        "t3",
        "d34",
        "Nbn34",
        "t4",
        "Na",
    }
    assert used_names == expected


def test_get_parameter_names_used_by_scenario_with_migration(header_text):
    """Vérifie que le taux demigration 'r1' est bien inclus dans les paramètres
    référencés par un scénario qui contient un MigrationEvent -- sinon il serait
    exclu à tort des colonnes du reftable.bin pour ce scénario (même bug que
    celui déjà corrigé pour t11..t44). Par contre le second taux littéral 0.0,
    lui n'est pas présent."""
    header_text_modified = _check_replace_modification(
        header_text,
        "t1 split 1 4 2 ra",
        "tm1 mig 1 4 r1\n tm2 mig 1 4 0.0\n t1 split 1 4 2 ra",
    )
    scenarios = parse_header_scenarios(header_text_modified)
    scenario2 = next(s for s in scenarios if s.index == 2)

    used_names = get_parameter_names_used_by_scenario(scenario2)

    assert "r1" in used_names
    expected = {
        "N1",
        "N2",
        "N3",
        "N4",
        "t1",
        "ra",
        "t2",
        "d3",
        "Nbn3",
        "d4",
        "Nbn4",
        "N34",
        "t3",
        "d34",
        "Nbn34",
        "t4",
        "Na",
        "tm1",
        "tm2",
        "r1",
    }
    assert used_names == expected


def test_get_parameter_used_by_scenario_with_bottleneck():
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    scenario = parse_scenario_block(block)
    used_names = get_parameter_names_used_by_scenario(scenario)
    assert {"tbn", "p", "tibn", "s"} <= used_names


def test_get_parameter_used_by_scenario_with_varnelog():
    block = """
        scenario 1 [0.16667] (4)
        N1 N2 
        0 sample 1
        0 sample 2
        t1  VARNELOG 2 NF r 
        t2  merge 1 2
        """

    scenario = parse_scenario_block(block)
    used_names = get_parameter_names_used_by_scenario(scenario)
    assert used_names == {"N1", "N2", "t1", "t2", "NF", "r"}


def test_evaluate_expression():
    values = {"t1": 12.3, "t2": 4881.0, "d3": 35.0}

    assert evaluate_expression("t1", values) == 12.3
    assert evaluate_expression("0", values) == 0.0
    assert evaluate_expression("t2-d3", values) == 4881.0 - 35.0
    assert evaluate_expression("t2+d3", values) == 4881.0 + 35.0

    with pytest.raises(ValueError):
        evaluate_expression("inconnu", values)


def test_build_demography_scenario1(header_text):
    """Vérifie que build_demography produit la bonne structure d'événements
    pour le scénario 1 de human, avec des valeurs de paramètres fixées à la
    main (pas de tirage aléatoire ici, pour isoler le test de la logique de
    construction de la démographie)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    # Valeurs choisies à la main, cohérentes avec les contraintes
    # (t4 > t3 > t2 > t2-d3, t2-d4 ; t3 > t3-d34)
    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
    }

    demography = build_demography(scenario1, values)

    # 4 populations créées
    assert len(demography.populations) == 4
    assert {p.name for p in demography.populations} == {"pop1", "pop2", "pop3", "pop4"}

    # Les événements de fusion sont bien présents, avec les bons temps
    splits = [
        e
        for e in demography.events
        if isinstance(e, msprime.demography.PopulationSplit)
    ]
    assert len(splits) == 3

    split_by_time = {s.time: s for s in splits}
    assert (
        split_by_time[10].derived == ["pop1"] and split_by_time[10].ancestral == "pop2"
    )
    assert (
        split_by_time[5000].derived == ["pop4"]
        and split_by_time[5000].ancestral == "pop3"
    )
    assert (
        split_by_time[8000].derived == ["pop3"]
        and split_by_time[8000].ancestral == "pop2"
    )


def test_extract_referenced_names():
    """Vérifie l'extraction de noms sur des cas simples."""
    assert extract_referenced_names("t1") == {"t1"}
    assert extract_referenced_names("0") == set()
    assert extract_referenced_names("t2-d3") == {"t2", "d3"}
    assert extract_referenced_names("t2+d3") == {"t2", "d3"}


def test_get_parameter_names_used_by_scenario1(header_text):
    """Vérifie que le scénario 1 référence bien exactement les 16 paramètres
    attendus (21 priors déclarés au total dans header.txt, moins
    ra/t11/t22/t33/t44 qui appartiennent aux scénarios 2-6)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    used_names = get_parameter_names_used_by_scenario(scenario1)

    expected = {
        "N1",
        "N2",
        "N3",
        "N4",
        "t1",
        "t2",
        "d3",
        "Nbn3",
        "d4",
        "Nbn4",
        "N34",
        "t3",
        "d34",
        "Nbn34",
        "t4",
        "Na",
    }
    assert used_names == expected
    assert len(used_names) == 16


def test_rescale_demography(header_text):
    """Vérifie que la mise à l'échelle de la démographie fonctionne
    correctement."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
    }

    # Vérifie que les tailles de population sont bien mises à l'échelle
    demography = build_demography(scenario1, values)
    scaled_demography = rescale_demography(demography, 2.0)
    for pop_scaled, pop in zip(
        scaled_demography.populations, demography.populations, strict=True
    ):
        assert pop_scaled.initial_size == pop.initial_size * 2.0

    # Vérifie que les tailles de population sont bien mises à l'échelle pour les variations de population hors bottleneck
    change_events = [
        e
        for e in demography.events
        if hasattr(e, "initial_size") and e.initial_size is not None
    ]
    scaled_change_events = [
        e
        for e in scaled_demography.events
        if hasattr(e, "initial_size") and e.initial_size is not None
    ]
    for event_scaled, event in zip(scaled_change_events, change_events, strict=True):
        assert event_scaled.initial_size == event.initial_size * 2.0

    # Vérification que l'original n'est pas modifié
    assert demography.populations[0].initial_size == 50000


def test_rescale_demography_with_migration(header_text):
    for text in ["mig", "smig"]:
        header_text_modified = _check_replace_modification(
            header_text, "t1 split 1 4 2 ra", f"tm1 {text} 1 4 r1\n t1 split 1 4 2 ra"
        )

        scenarios = parse_header_scenarios(header_text_modified)
        scenario2 = next(s for s in scenarios if s.index == 2)

        values = {
            "N1": 1000,
            "N2": 1000,
            "N3": 1000,
            "N4": 1000,
            "tm1": 5,
            "t1": 10,
            "t2": 5000,
            "d3": 30,
            "Nbn3": 200,
            "d4": 20,
            "Nbn4": 300,
            "N34": 60000,
            "t3": 8000,
            "d34": 25,
            "Nbn34": 250,
            "t4": 9000,
            "Na": 40000,
            "ra": 0.1,
            "r1": 0.3,
        }
        factor = 0.4
        demography = build_demography(scenario2, values)
        demography_rescaled = rescale_demography(demography, factor=factor)
        for epoch, epoch_rescaled in zip(
            demography.debug().epochs, demography_rescaled.debug().epochs, strict=True
        ):
            assert np.array_equal(
                epoch.migration_matrix, epoch_rescaled.migration_matrix
            )
        for pop, pop_rescaled in zip(
            demography.populations, demography_rescaled.populations, strict=True
        ):
            assert pop.initial_size * factor == pop_rescaled.initial_size


def test_rescale_demography_with_bottleneck():
    block = """
    scenario 1 [0.16667] (4)
    N1 N2
    0 sample 1
    0 sample 2
    tbn bn 1 p
    tibn ibn 2 s
    t1 merge 1 2
    """

    values = {
        "N1": 1000,
        "N2": 1000,
        "tbn": 10,
        "p": 0.6,
        "tibn": 20,
        "s": 100,
        "t1": 500,
    }
    scenario = parse_scenario_block(block)
    demography = build_demography(scenario, values)
    demography_rescaled = rescale_demography(demography, 0.5)

    for pop, pop_rescaled in zip(
        demography.populations, demography_rescaled.populations, strict=True
    ):
        assert pop.initial_size * 0.5 == pop_rescaled.initial_size

    for event, event_rescaled in zip(
        demography.events, demography_rescaled.events, strict=True
    ):
        assert event == event_rescaled


# (type de locus, ploïdie msprime, facteur de rescale appliqué par le pipeline) à sex ratio 0.5 :
# <A> reste en ploïdie 2 sans rescale (facteur 1, cf. CLAUDE.md), les autres types passent
# en ploïdie 1 avec rescale_demography(coalescence_coefficient(type, 0.5) / 2).
_BOTTLENECK_LOCUS_CASES = [
    ("A", 2, 1.0),
    ("X", 1, coalescence_coefficient("X", 0.5) / 2),
    ("Y", 1, coalescence_coefficient("Y", 0.5) / 2),
    ("M", 1, coalescence_coefficient("M", 0.5) / 2),
]


def _fraction_coalesced_at_bottleneck(block, values, ploidy, factor, time, seed):
    """Fraction de paires de lignées dont la coalescence a lieu exactement à `time`.

    Args:
        block: Bloc de scénario (une seule population, un seul bottleneck).
        values: Valeurs des paramètres du scénario.
        ploidy: Ploïdie passée à msprime (2 pour <A>, 1 sinon).
        factor: Facteur de rescale des tailles (1.0 = pas de rescale).
        time: Date du bottleneck (générations).
        seed: Graine msprime.

    Returns:
        (nombre de paires coalescées à `time`, nombre de réplicats).
    """
    n_replicates = 4000
    demography = build_demography(parse_scenario_block(block), values)
    if factor != 1.0:
        demography = rescale_demography(demography, factor)
    n_individuals = 1 if ploidy == 2 else 2  # toujours 2 lignées échantillonnées
    replicates = msprime.sim_ancestry(
        samples={"pop1": n_individuals},
        demography=demography,
        ploidy=ploidy,
        sequence_length=1,
        num_replicates=n_replicates,
        random_seed=seed,
    )
    tmrca = np.array([ts.first().tmrca(0, 1) for ts in replicates])
    # La date n'est pas rescalée : les coalescences du bottleneck tombent pile à `time`.
    at_bottleneck = np.isclose(tmrca, time)
    return int(at_bottleneck.sum()), n_replicates


@pytest.mark.parametrize("locus_type, ploidy, factor", _BOTTLENECK_LOCUS_CASES)
def test_simple_bottleneck_effect_is_not_rescaled(locus_type, ploidy, factor):
    """La date et la proportion d'un `bn` ne sont jamais rescalées, quel que soit le
    type de locus : P(deux lignées coalescent au bottleneck) = proportion^2,
    à la date exacte du bottleneck, en ploïdie 1 comme en ploïdie 2."""
    block = """
    scenario 1 [1.0] (3)
    N1
    0 sample 1
    tbn bn 1 p
    """
    time, proportion = 10.0, 0.6
    # N très grand : la coalescence "de fond" avant `time` est négligeable (~1e-4).
    values = {"N1": 1e6, "tbn": time, "p": proportion}

    k, n = _fraction_coalesced_at_bottleneck(block, values, ploidy, factor, time, 101)

    result = stats.binomtest(k, n, proportion**2)
    assert result.pvalue > 1e-3, (
        f"<{locus_type}> : {k}/{n} coalescences à t={time}, attendu {proportion**2:.3f}"
    )


@pytest.mark.parametrize("locus_type, ploidy, factor", _BOTTLENECK_LOCUS_CASES)
def test_instantaneous_bottleneck_effect_scales_with_gene_copies(
    locus_type, ploidy, factor
):
    """La force d'un `ibn` n'est jamais rescalée, mais son effet suit le type de locus
    par la taille rescalée : P(coalescence au bottleneck) = 1 - exp(-s / (ploïdie * N'))
    avec N' = N * facteur (vérifié pour msprime 1.4.2). C'est ce qui rend un même
    `strength` plus fort pour <M>/<Y> que pour <A>."""
    block = """
    scenario 1 [1.0] (3)
    N1
    0 sample 1
    tibn ibn 1 s
    """
    time, strength, size = 10.0, 50_000.0, 100_000.0
    values = {"N1": size, "tibn": time, "s": strength}
    expected = 1 - np.exp(-strength / (ploidy * size * factor))

    k, n = _fraction_coalesced_at_bottleneck(block, values, ploidy, factor, time, 202)

    result = stats.binomtest(k, n, expected)
    assert result.pvalue > 1e-3, (
        f"<{locus_type}> : {k}/{n} = {k / n:.3f} coalescences à t={time}, attendu {expected:.3f}"
    )
