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
    build_demography,
    evaluate_expression,
    extract_referenced_names,
    get_parameter_names_used_by_scenario,
    rescale_demography,
)
from bridge.observed_data import coalescence_coefficient
from bridge.scenario_parser import parse_header_scenarios, parse_scenario_block


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
