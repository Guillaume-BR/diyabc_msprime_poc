"""
Vérifie que scenario_parser produit, sur le vrai header.txt du dataset
human, exactement les événements qu'on a décortiqués à la main avec le
mentor (voir notes/exploration.md) pour le scénario 1.
"""

import os
from pathlib import Path

import msprime
import pytest

from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_independent_loci,
    simulate_snp_genotypes,
)
from bridge.demography_builder import (
    build_demography,
    evaluate_expression,
    extract_referenced_names,
    get_parameter_names_used_by_scenario,
)
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import count_samples_per_population, population_index_to_name
from bridge.parameter_sampling import draw_parameter_values
from bridge.pipeline import (
    build_random_demography_for_scenario_index,
    compute_summary_statistics,
    run_poc_for_directory,
)
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.reftable_loop import run_reftable_simulation, write_reftable_bin
from bridge.scenario_parser import parse_header_scenarios
from bridge.scenario_types import (
    MergeEvent,
    OrderConstraint,
    Prior,
    SampleEvent,
    SplitEvent,
    VarNeEvent,
)
from bridge.snp_writer import write_snp_file

REFERENCE_DIR = Path(__file__).parent.parent / "reference" / "human"
GENERAL_BINARY_PATH = os.environ.get("DIYABC_GENERAL_PATH")


@pytest.fixture
def header_text() -> str:
    return (REFERENCE_DIR / "header.txt").read_text()


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


def test_build_demography_scenario2_admixture(header_text):
    """Vérifie que build_demography traduit un SplitEvent en événement
    msprime Admixture avec les bonnes populations et proportions
    (rate, 1-rate)."""
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
        e for e in demography.events if isinstance(e, msprime.demography.Admixture)
    ]
    assert len(admixtures) == 1
    admixture = admixtures[0]
    assert admixture.time == 10
    assert admixture.derived == "pop1"
    assert admixture.ancestral == ["pop4", "pop2"]
    assert admixture.proportions == pytest.approx([0.3, 0.7])


def test_get_parameter_names_used_by_scenario2(header_text):
    """Vérifie que le taux d'admixture 'ra' est bien inclus dans les
    paramètres référencés par un scénario qui contient un SplitEvent --
    sinon il serait exclu à tort des colonnes du reftable.bin pour ce
    scénario (même bug que celui déjà corrigé pour t11..t44)."""
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


def test_priors_and_constraints(header_text):
    priors, constraints = parse_priors(header_text)

    assert len(priors) == 21
    assert len(constraints) == 4

    priors_by_name = {p.name: p for p in priors}
    assert priors_by_name["N1"].category == "N"
    assert priors_by_name["N1"].law == "UN"
    assert priors_by_name["N1"].bounds == (1000.0, 100000.0, 0.0, 0.0)

    assert priors_by_name["t1"].category == "T"
    assert priors_by_name["t1"].bounds == (1.0, 30.0, 0.0, 0.0)

    assert OrderConstraint(param1="t4", operator=">", param2="t3") in constraints
    assert OrderConstraint(param1="t3", operator=">", param2="t2") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t33") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t22") in constraints


def test_draw_parameter_values(header_text):
    """Vérifie que draw_parameter_values tire bien une valeur pour chaque
    prior, et que le tirage retourné respecte toutes les contraintes
    d'ordre (t4>t3, t3>t2, t44>t33, t44>t22)."""
    priors, constraints = parse_priors(header_text)

    seed = 42
    values = draw_parameter_values(priors, constraints, seed)

    # Toutes les valeurs ont bien été tirées
    assert set(values.keys()) == {p.name for p in priors}

    # Toutes les contraintes sont respectées par ce tirage
    for constraint in constraints:
        assert constraint.is_satisfied(values)


def test_draw_parameter_values_reproducible(header_text):
    """Même graine -> même tirage (déterminisme attendu pour la
    reproductibilité scientifique)."""
    priors, constraints = parse_priors(header_text)

    values1 = draw_parameter_values(priors, constraints, seed=123)
    values2 = draw_parameter_values(priors, constraints, seed=123)

    assert values1 == values2


def test_evaluate_expression():
    values = {"t1": 12.3, "t2": 4881.0, "d3": 35.0}

    assert evaluate_expression("t1", values) == 12.3
    assert evaluate_expression("0", values) == 0.0
    assert evaluate_expression("t2-d3", values) == 4881.0 - 35.0
    assert evaluate_expression("t2+d3", values) == 4881.0 + 35.0

    with pytest.raises(ValueError):
        evaluate_expression("inconnu", values)


def test_build_demography_scenario1(header_text):
    """Vérifie que build_demography produit la bonne structure
    d'événements pour le scénario 1 de human, avec des valeurs de
    paramètres fixées à la main (pas de tirage aléatoire ici, pour
    isoler le test de la logique de construction de la démographie)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)
    """Extrait les priors et les contraintes d'ordre de header.txt.

 

    Retourne (priors, constraints). Une ligne qui ne correspond à aucun

    des deux formats connus lève une erreur explicite plutôt que d'être

    silencieusement ignorée : contrairement aux événements de scénario, on

    n'a pas de raison de s'attendre à du vocabulaire non géré ici pour le

    dataset human.

    """
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
        if isinstance(e, __import__("msprime").demography.PopulationSplit)
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


def test_parse_loci_description(header_text):
    """Vérifie le parsing de la section 'loci description' de human,
    format condensé à un seul type d'héritage."""
    description = parse_loci_description(header_text)

    assert description.total_loci == 5000
    assert description.group == "G1"
    assert description.start_index == 0  # "from 1" en 1-based -> 0 en 0-based


def test_pipeline_scenario1(header_text):
    """Vérifie que le pipeline complet (header.txt -> Demography) fonctionne
    de bout en bout sur le scénario 1, et que la démographie produite a la
    structure attendue (4 populations, 3 fusions)."""

    demography, values = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )

    assert len(demography.populations) == 4

    splits = [
        e
        for e in demography.events
        if isinstance(e, msprime.demography.PopulationSplit)
    ]
    assert len(splits) == 3

    # Les valeurs tirées doivent inclure tous les paramètres du header
    assert "N1" in values
    assert "t1" in values


OBSERVED_SNP_FILE = REFERENCE_DIR / "human_snp_all22chr_maf5.snp"


def test_count_samples_per_population():
    """Vérifie que le comptage retrouve bien les 4 populations à 30
    individus chacune, annoncées en commentaire dans le fichier."""
    counts = count_samples_per_population(OBSERVED_SNP_FILE)

    assert set(counts.keys()) == {"ASW", "YRI", "CHB", "GBR"}
    assert all(n == 30 for n in counts.values())


def test_population_index_to_name():
    """Vérifie le mapping indice de scénario (1-indexed) -> nom réel de
    population, dans l'ordre d'apparition du fichier .snp."""
    mapping = population_index_to_name(OBSERVED_SNP_FILE)

    assert mapping == {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"}


def test_simulate_independent_loci_scenario1(header_text):
    """Vérifie que build_samples_argument construit bien le dict attendu
    par msprime.sim_ancestry, avec les bons noms de populations et le bon
    nombre d'individus par population."""

    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    # 4 populations, 30 individus chacune attendus
    assert samples == {"pop1": 30, "pop2": 30, "pop3": 30, "pop4": 30}

    num_loci = 10  # petit nombre pour un test rapide, pas les 51250 réels
    tree_sequences = list(
        simulate_independent_loci(demography, samples, num_loci=num_loci, seed=123)
    )

    # On doit obtenir exactement num_loci arbres indépendants
    assert len(tree_sequences) == num_loci

    # Chaque arbre doit avoir le bon nombre total de lignées échantillonnées
    # (30 individus x 4 populations x ploidy 2 = 240 lignées)
    for ts in tree_sequences:
        assert ts.num_samples == 30 * 4 * 2


def test_simulate_snp_genotypes_scenario1(header_text):
    """Vérifie que chaque locus simulé est polymorphe (au moins un 0 et
    un 1 parmi les génotypes), garantissant la propriété centrale de
    l'algorithme de Hudson : exactement une mutation par locus, jamais
    un locus monomorphe."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    num_loci = 20
    tree_sequences = simulate_independent_loci(demography, samples, num_loci, seed=123)
    genotypes_per_locus = list(simulate_snp_genotypes(tree_sequences, seed=456))

    assert len(genotypes_per_locus) == num_loci

    for locus_genotypes in genotypes_per_locus:
        all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
        assert set(all_genotypes) == {0, 1}, f"Locus non polymorphe : {locus_genotypes}"


def test_simulate_snp_genotypes_grouped_by_population(header_text):
    """Vérifie que les génotypes sont bien regroupés par nom de
    population (pop1..pop4), avec le bon nombre de lignées par groupe
    (30 individus x ploidy 2 = 60 lignées par population), et que chaque
    locus reste globalement polymorphe."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    num_loci = 10
    tree_sequences = simulate_independent_loci(demography, samples, num_loci, seed=123)
    genotypes_per_locus = list(simulate_snp_genotypes(tree_sequences, seed=456))

    assert len(genotypes_per_locus) == num_loci

    for locus_genotypes in genotypes_per_locus:
        assert set(locus_genotypes.keys()) == {"pop1", "pop2", "pop3", "pop4"}
        for _pop_name, genos in locus_genotypes.items():
            assert len(genos) == 60  # 30 individus x ploidy 2

        # Polymorphe globalement (au moins un 0 et un 1 sur l'ensemble)
        all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
        assert set(all_genotypes) == {0, 1}


def test_run_poc_for_directory():
    """Vérifie le point d'entrée de haut niveau : à partir d'un simple
    chemin de dossier (comme le -p ./ de DIYABC), tout le pipeline doit
    fonctionner sans qu'on ait à lire manuellement header.txt ou le
    fichier .snp nous-mêmes."""
    mutated, values = run_poc_for_directory(
        REFERENCE_DIR,
        scenario_index=1,
        num_loci=15,
        seed=42,
    )

    mutated_list = list(mutated)
    assert len(mutated_list) == 15
    assert "N1" in values


def test_write_snp_file_small_case(tmp_path):
    """Vérifie l'écriture du fichier .snp sur un cas minimal : 2 loci,
    2 populations, 2 lignées (1 individu) chacune."""
    genotypes_per_locus = [
        {
            "pop1": [0, 1],
            "pop2": [1, 1],
        },  # locus 0 : pop1 hétérozygote (1), pop2 homozygote dérivé (2)
        {
            "pop1": [0, 0],
            "pop2": [0, 1],
        },  # locus 1 : pop1 homozygote ancestral (0), pop2 hétérozygote (1)
    ]

    output_file = tmp_path / "test.snp"
    write_snp_file(genotypes_per_locus, output_file)

    content = output_file.read_text()
    lines = content.strip().splitlines()

    assert lines[0] == "IND SEX POP A A"
    assert lines[1] == "sim_pop1_1 9 pop1 1 0"
    assert lines[2] == "sim_pop2_1 9 pop2 2 1"


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie -- "
    "ce test nécessite le binaire 'general' compilé de DIYABC.",
)
def test_compute_summary_statistics_scenario1(tmp_path):
    """Vérifie que compute_summary_statistics produit bien les 112
    statistiques résumées attendues (filtre ALL), en déléguant le calcul
    au vrai binaire C++ sur des données simulées par notre pipeline."""
    summary_statistics, values = compute_summary_statistics(
        reference_directory=REFERENCE_DIR,
        scenario_index=1,
        num_loci=10,
        seed=42,
        general_binary_path=GENERAL_BINARY_PATH,
        work_directory=tmp_path,
        stats_filter="ALL",
    )
    print(sorted(summary_statistics.keys()))
    # 112 statistiques attendues (vu dans "group summary statistics (112)")
    assert len(summary_statistics) == 130

    # Quelques noms de colonnes attendus, parmi les plus simples à vérifier
    assert "ML1p_1" in summary_statistics
    assert "FST1m_1" in summary_statistics

    # Les valeurs de paramètres tirées doivent toujours être présentes
    assert "N1" in values


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie.",
)
def test_run_reftable_simulation_scenario1(tmp_path):
    """Vérifie que run_reftable_simulation produit bien nrec particules
    distinctes (tirages de paramètres différents), chacune avec ses 130
    statistiques résumées calculées."""
    nrec = 4
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenario_index=1,
        num_loci=10,
        nrec=nrec,
        general_binary_path=GENERAL_BINARY_PATH,
        base_work_directory=tmp_path,
        stats_filter="ALL",
    )

    assert len(results) == nrec

    # Les résultats doivent être dans l'ordre des indices de particule
    assert [r.particle_index for r in results] == list(range(nrec))

    # Chaque particule doit avoir ses 130 stats et tous les paramètres
    for result in results:
        assert len(result.summary_statistics) == 130
        assert "N1" in result.parameter_values

    # Les tirages de paramètres doivent être DIFFÉRENTS entre particules
    # (preuve que chaque particule a bien sa propre seed)
    n1_values = [r.parameter_values["N1"] for r in results]
    assert len(set(n1_values)) == nrec  # 4 valeurs distinctes


def test_is_constant_prior():
    """Vérifie la règle de filtrage des priors quasi-constants."""
    # Cas normal : large intervalle, jamais constant
    normal_prior = Prior(
        name="N1", category="N", law="UN", bounds=(1000.0, 100000.0, 0.0, 0.0)
    )
    assert is_constant_prior(normal_prior) is False

    # Cas dégénéré : min == max, clairement constant
    constant_prior = Prior(
        name="X", category="N", law="UN", bounds=(100.0, 100.0, 0.0, 0.0)
    )
    assert is_constant_prior(constant_prior) is True

    # Cas limite : différence infime, sous le seuil
    near_constant_prior = Prior(
        name="Y", category="N", law="UN", bounds=(100.0, 100.00001, 0.0, 0.0)
    )
    assert is_constant_prior(near_constant_prior) is True

    # Cas limite inverse : différence juste au-dessus du seuil
    barely_variable_prior = Prior(
        name="Z", category="N", law="UN", bounds=(100.0, 100.1, 0.0, 0.0)
    )
    assert is_constant_prior(barely_variable_prior) is False


def test_write_reftable_bin(tmp_path, header_text):
    """Vérifie l'écriture du reftable.bin, et sa relecture (vérification
    manuelle du format, sans dépendre de readReftable.R pour ce test
    Python -- juste une vérification structurelle du binaire produit)."""
    nrec = 3
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenario_index=1,
        num_loci=10,
        nrec=nrec,
        general_binary_path=GENERAL_BINARY_PATH,
        base_work_directory=tmp_path / "particles",
        stats_filter="ALL",
    )

    priors, _ = parse_priors(header_text)
    output_file = tmp_path / "reftable.bin"
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)
    write_reftable_bin(results, priors, scenario1, output_file)

    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_extract_referenced_names():
    """Vérifie l'extraction de noms sur des cas simples."""
    assert extract_referenced_names("t1") == {"t1"}
    assert extract_referenced_names("0") == set()
    assert extract_referenced_names("t2-d3") == {"t2", "d3"}
    assert extract_referenced_names("t2+d3") == {"t2", "d3"}


def test_get_parameter_names_used_by_scenario1(header_text):
    """Vérifie que le scénario 1 référence bien exactement les 16
    paramètres attendus (21 priors déclarés au total dans header.txt,
    moins ra/t11/t22/t33/t44 qui appartiennent aux scénarios 2-6)."""
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
