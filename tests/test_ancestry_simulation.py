"""Vérifie ancestry_simulation : construction de l'argument samples pour
msprime.sim_ancestry, simulation de généalogies indépendantes par locus,
et mutation Hudson (exactement une mutation par locus, toujours
polymorphe)."""

import msprime
import numpy as np
import pytest
from conftest import (
    OBSERVED_MSS_FILE_TE2_XY,
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE4,
    OBSERVED_SNP_FILE_TE5,
)

from bridge.ancestry_simulation import (
    _distribution_from_position,
    _group_prior_values_from_columns,
    _place_gsm_row_on_dense_grid,
    _reindex_reads_by_msprime_name,
    _sni_row_on_dense_grid,
    build_group_local_param_per_locus,
    build_group_local_param_per_locus_from_values,
    build_male_only_samples_argument,
    build_male_only_samples_argument_ms_dna,
    build_matrix_microsat_per_locus,
    build_matrix_per_locus,
    build_microsat_local_param_per_locus,
    build_microsat_transition_matrix,
    build_microsat_transition_matrix_with_sni,
    build_rate_map,
    build_rate_map_per_locus,
    build_sample_sets_from_scenario,
    build_samples_argument,
    build_sex_stratified_samples_argument,
    build_sex_stratified_samples_argument_ms_dna,
    build_transition_matrix,
    compute_population_layout,
    compute_sample_layout,
    count_loci_per_group,
    dna_mutation_simulation_per_locus,
    microsat_mutation_simulation_per_locus,
    ms_dna_ancestry_parameters_for_heritage,
    observed_maf,
    simulate_genotypes_for_locus_type,
    simulate_independent_loci,
    simulate_poolseq_reads,
    simulate_poolseq_reads_with_mrc_filter,
    simulate_shared_ancestry_loci,
    simulate_snp_genotypes,
    with_maf_filter,
    with_maf_filter_shared_ancestry,
    with_mrc_filter,
)
from bridge.demography_builder import rescale_demography
from bridge.header_dataclasses import LociDescriptionDetailed, SampleEvent, Scenario
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    coalescence_coefficient,
    observed_reads,
    parse_sex_ratio,
)
from bridge.pipeline import build_random_demography_for_scenario_index
from bridge.prior_parser import parse_group_priors
from bridge.scenario_parser import parse_header_scenarios


# Petit helper pour créer un TreeSequence minimal avec une population et un individu, pour tester compute_sample_layout
def _ts_from_sample_sets(sample_sets):
    """TreeSequence minimale pour tester un layout : la topologie n'importe pas."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    return msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=3,
    )


def test_build_sample_sets_from_scenario(microsat_context_te1):
    """Sur le vrai header sériel : 4 événements `sample`, tous sur pop1.

    Vérifie surtout que les effectifs sont pris DANS L'ORDRE des
    événements et non par nom de population -- d'où des effectifs
    volontairement tous différents : une indexation par
    `f"pop{event.pop}"` renverrait 20 quatre fois, les quatre `sample`
    de ce scénario pointant tous la population 1.
    """
    scenarios = parse_header_scenarios(microsat_context_te1.header_text)
    counts = {"pop1": 20, "pop2": 15, "pop3": 30, "pop4": 10}

    sample_sets = build_sample_sets_from_scenario(
        scenario=scenarios[0], values={}, counts_by_samples=counts
    )

    assert len(sample_sets) == 4
    assert [s.population for s in sample_sets] == ["pop1"] * 4
    assert [s.num_samples for s in sample_sets] == [20, 15, 30, 10]
    assert [s.time for s in sample_sets] == [0, 50, 200, 500]


def test_build_sample_sets_from_scenario_evaluates_time_expressions():
    """Un temps littéral et un temps paramétré dans le même scénario.

    Le chemin "temps exprimé par un nom de paramètre"
    (particuleC.cpp:599-605) n'est exercé par AUCUN jeu de données réel
    du dépôt -- toy_example1_ms n'a que des littéraux -- d'où ce
    scénario synthétique.
    """
    fake_scenario = Scenario(
        index=1,
        weight=1.0,
        initial_pop_size_exprs=["100"],
        events=[
            SampleEvent(time_expr="0", pop=1),
            SampleEvent(time_expr="tbn", pop=1),
        ],
    )
    # 317.0 : valeur non confondable avec un littéral du scénario, pour
    # qu'une assertion ne puisse pas réussir par coïncidence.
    values = {"tbn": 317.0}
    counts = {"pop1": 20, "pop2": 15}

    sample_sets = build_sample_sets_from_scenario(
        scenario=fake_scenario, values=values, counts_by_samples=counts
    )

    assert len(sample_sets) == 2
    assert [s.population for s in sample_sets] == ["pop1", "pop1"]
    assert [s.num_samples for s in sample_sets] == [20, 15]
    assert [s.time for s in sample_sets] == [0.0, 317.0]


def test_build_sample_sets_from_scenario_raises_on_count_mismatch():
    """Vérifie que build_sample_sets_from_scenario lève une ValueError si le
    nombre de comptes d'échantillons ne correspond pas au nombre d'événements sample.
    """
    # test pour lever l'erreur
    fake_scenario = Scenario(
        index=1,
        weight=1.0,
        initial_pop_size_exprs=["100"],
        events=[SampleEvent(time_expr="0", pop=1), SampleEvent(time_expr="50", pop=2)],
    )
    values = {}
    counts = {"pop1": 20, "pop2": 15, "pop3": 30}
    with pytest.raises(
        ValueError, match="2 événements sample mais 3 échantillons observés"
    ):
        build_sample_sets_from_scenario(
            scenario=fake_scenario, values=values, counts_by_samples=counts
        )


def test_simulate_independent_loci_scenario1(header_text):
    """Vérifie que build_samples_argument construit bien le dict attendu
    par msprime.sim_ancestry, avec les bons noms de populations et le bon
    nombre d'individus par population."""

    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE_HUMAN)

    # 4 populations, 30 individus chacune attendus
    assert samples == {"pop1": 30, "pop2": 30, "pop3": 30, "pop4": 30}

    num_loci = 10  # petit nombre pour un test rapide, pas les 51250 réels
    tree_sequences = list(
        simulate_independent_loci(
            demography, samples, num_loci=num_loci, seed=123, ploidy=2
        )
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
    samples = build_samples_argument(OBSERVED_SNP_FILE_HUMAN)

    num_loci = 20
    tree_sequences = simulate_independent_loci(
        demography, samples, num_loci, seed=123, ploidy=2
    )
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
    samples = build_samples_argument(OBSERVED_SNP_FILE_HUMAN)

    num_loci = 10
    tree_sequences = simulate_independent_loci(
        demography, samples, num_loci, seed=123, ploidy=2
    )
    genotypes_per_locus = list(simulate_snp_genotypes(tree_sequences, seed=456))

    assert len(genotypes_per_locus) == num_loci

    for locus_genotypes in genotypes_per_locus:
        assert set(locus_genotypes.keys()) == {"pop1", "pop2", "pop3", "pop4"}
        for _pop_name, genos in locus_genotypes.items():
            assert len(genos) == 60  # 30 individus x ploidy 2

        # Polymorphe globalement (au moins un 0 et un 1 sur l'ensemble)
        all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
        assert set(all_genotypes) == {0, 1}


def test_build_sex_stratified_samples_argument():
    """Vérifie que build_samples_argument construit bien le dict attendu
    par msprime.sim_ancestry, avec les bons noms de populations et le bon
    nombre d'individus par population, en tenant compte du sexe des individus (pour les loci <X>/<Y>/<M>).
    """

    with pytest.raises(ValueError, match="sexe inconnu"):
        build_sex_stratified_samples_argument(
            OBSERVED_SNP_FILE_HUMAN
        )  # sexe non renseigné

    liste_samples = build_sex_stratified_samples_argument(OBSERVED_SNP_FILE_TE5)
    assert len(liste_samples) == 6  # 3 populations x 2 sexes (M/F)
    for sample_set in liste_samples:
        assert sample_set.population in {"pop1", "pop2", "pop3"}
        assert sample_set.num_samples in {10}  # 10 M ou 10 F par population
        assert sample_set.ploidy in {1, 2}  # M=1, F=2


def test_build_male_only_samples_argument():
    """Vérifie que build_male_only_samples_argument construit bien un
    dict {population: nombre_de_mâles} (PAS une liste de SampleSet,
    contrairement à build_sex_stratified_samples_argument) -- <Y> n'a
    besoin que d'un ploidy uniforme=1 parmi les mâles, pas d'hétérogénéité
    au sein d'une population."""
    with pytest.raises(ValueError, match="sexe inconnu"):
        build_male_only_samples_argument(OBSERVED_SNP_FILE_HUMAN)  # sexe non renseigné

    samples = build_male_only_samples_argument(OBSERVED_SNP_FILE_TE5)
    assert samples == {"pop1": 10, "pop2": 10, "pop3": 10}


def test_simulate_shared_ancestry_loci(header_text):
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE_HUMAN)
    shared_ancestry = simulate_shared_ancestry_loci(
        demography=demography, samples=samples, num_loci=5, seed=42, ploidy=1
    )

    trees = list(shared_ancestry)
    assert len(trees) == 5
    assert all(t is trees[0] for t in trees)


def test_simulate_genotypes_for_locus_type(snp_context_te5):
    demography, _ = build_random_demography_for_scenario_index(
        snp_context_te5.header_text, scenario_index=1, seed=42
    )
    for locus_type in ["A", "X", "Y", "M", "H"]:
        genotypes = simulate_genotypes_for_locus_type(
            demography=demography,
            locus_type=locus_type,
            context=snp_context_te5,
            num_loci=5,
            seed=42,
        )
        assert len(list(genotypes)) == 5
        for locus_genotypes in genotypes:
            all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
            assert set(all_genotypes) == {
                0,
                1,
            }, f"Locus non polymorphe : {locus_genotypes}"


def test_compute_sample_layout_different_count():
    """Vérifie que compute_sample_layout renvoie bien une liste de tuples
    (nom_population, sample_ids) avec les bons noms de populations et le
    bon nombre d'individus par population.
    """

    # test effectifs inégaux
    ts = _ts_from_sample_sets(
        [msprime.SampleSet(7, "pop1"), msprime.SampleSet(3, "pop2")]
    )
    samples = {"pop1": 7, "pop2": 3}
    computed_layout = compute_sample_layout(ts, samples)

    assert [pop_name for pop_name, _ in computed_layout] == ["pop1", "pop2"]
    assert np.array_equal(computed_layout[0][1], np.arange(14))
    assert np.array_equal(computed_layout[1][1], np.arange(14, 20))


def test_compute_sample_layout_with_mixed_ploidy():
    # test sur la ploidie
    ts = _ts_from_sample_sets(
        [
            msprime.SampleSet(3, "pop1", ploidy=2),
            msprime.SampleSet(2, "pop1", ploidy=1),
            msprime.SampleSet(4, "pop2", ploidy=2),
            msprime.SampleSet(1, "pop2", ploidy=1),
        ]
    )
    samples = {"pop1": 5, "pop2": 5}
    computed_layout = compute_sample_layout(ts, samples)

    assert [pop_name for pop_name, _ in computed_layout] == ["pop1", "pop2"]
    assert np.array_equal(computed_layout[0][1], np.arange(8))
    assert np.array_equal(computed_layout[1][1], np.arange(8, 17))


def test_compute_sample_layout_with_null_effectif():
    # test avec effectif nul pour une population
    ts = _ts_from_sample_sets(
        [
            msprime.SampleSet(5, "pop1"),
            msprime.SampleSet(0, "pop2"),
            msprime.SampleSet(3, "anc"),
        ]
    )
    samples = {"pop1": 5, "pop2": 0, "anc": 3}
    computed_layout = compute_sample_layout(ts, samples)

    assert [pop_name for pop_name, _ in computed_layout] == ["pop1", "anc"]
    assert np.array_equal(computed_layout[0][1], np.arange(10))
    assert np.array_equal(computed_layout[1][1], np.arange(10, 16))


def test_compute_sample_layout_with_mismatched_count():
    # test sur le garde
    ts = _ts_from_sample_sets(
        [
            msprime.SampleSet(5, "pop1"),
            msprime.SampleSet(0, "pop2"),
            msprime.SampleSet(3, "anc"),
        ]
    )
    samples = {"pop1": 5, "pop2": 0, "anc": 2}
    with pytest.raises(ValueError, match="Le nombre total d'individus"):
        compute_sample_layout(ts, samples)


def test_compute_sample_layout_matches_population_layout_on_non_serial_data(
    dna_context_te2,
):
    """Sur un jeu NON sériel, le découpage par échantillon doit être
    STRICTEMENT identique au découpage par population. C'est cette
    propriété qui rend sûre la substitution dans summary_statistics.py :
    sans elle, remplacer compute_population_layout changerait les
    résultats sur tous les jeux déjà validés.
    """
    scenarios = parse_header_scenarios(dna_context_te2.header_text)
    demography, values = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=7
    )
    counts = dna_context_te2.samples_default

    ts = msprime.sim_ancestry(
        samples=build_sample_sets_from_scenario(scenarios[0], values, counts),
        demography=demography,
        sequence_length=1,
        random_seed=7,
        ploidy=2,
    )

    by_population = compute_population_layout(ts)
    by_sample = compute_sample_layout(ts, counts)

    assert [n for n, _ in by_population] == [n for n, _ in by_sample]
    assert all(
        np.array_equal(a, b)
        for (_, a), (_, b) in zip(by_population, by_sample, strict=True)
    )


# tests relatifs aux filtres MAF
def test_with_maf_filter_no_filter_matches_direct_call(header_text):
    """maf=0.0 doit produire EXACTEMENT le même résultat qu'un appel
    direct à simulate_independent_loci + simulate_snp_genotypes (même
    graine pour les deux, comme fait déjà chaque branche de
    simulate_genotypes_for_locus_type) -- garantit que with_maf_filter
    ne change rien aux datasets déjà validés sans filtre MAF actif
    (human, toy_example5)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE_HUMAN)
    num_loci = 10

    direct = list(
        simulate_snp_genotypes(
            simulate_independent_loci(
                demography, samples, num_loci, seed=123, ploidy=1
            ),
            seed=123,
        )
    )
    via_filter = list(
        with_maf_filter(demography, samples, num_loci, maf=0.0, seed=123, ploidy=1)
    )

    assert via_filter == direct


def test_with_maf_filter_with_same_layout_matches_null_maf():
    """Vérifie que l'argument si l'argument layout correspond à compute_population_layout(ts) ne change pas le résultat de with_maf_filter."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    result_with_layout = with_maf_filter(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.0,
        seed=123,
        ploidy=2,
        counts_by_samples={"pop1": 7, "pop2": 3},
    )
    result_without_layout = with_maf_filter(
        demography, sample_sets, num_loci=10, maf=0.0, seed=123, ploidy=2
    )

    assert list(result_with_layout) == list(result_without_layout)


def test_with_maf_filter_with_same_layout_matches_nonnull_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter si le MAF est nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    result_with_layout = with_maf_filter(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.05,
        seed=123,
        ploidy=2,
        counts_by_samples={"pop1": 7, "pop2": 3},
    )
    result_without_layout = with_maf_filter(
        demography, sample_sets, num_loci=10, maf=0.05, seed=123, ploidy=2
    )

    assert list(result_with_layout) == list(result_without_layout)


def test_with_maf_filter_with_different_layout_null_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter si le MAF est nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=3,
    )

    layout = compute_population_layout(ts)
    split_layout = []
    for pop_name, node_ids in layout:
        half = len(node_ids) // 2
        split_layout.append((f"{pop_name}_a", node_ids[:half]))
        split_layout.append((f"{pop_name}_b", node_ids[half:]))

    with_split = with_maf_filter(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.0,
        seed=123,
        ploidy=2,
        counts_by_samples={"pop1_a": 4, "pop1_b": 3, "pop2_a": 2, "pop2_b": 1},
    )

    without_layout = with_maf_filter(
        demography, sample_sets, num_loci=10, maf=0.0, seed=123, ploidy=2
    )

    for split, plain in zip(with_split, without_layout, strict=True):
        assert list(split) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]  # discrimine
        assert split["pop1_a"] + split["pop1_b"] == plain["pop1"]  # recollement
        assert split["pop2_a"] + split["pop2_b"] == plain["pop2"]


def test_with_maf_filter_with_different_layout_nonnull_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter si le MAF est non nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=3,
    )

    layout = compute_population_layout(ts)
    split_layout = []
    for pop_name, node_ids in layout:
        half = len(node_ids) // 2
        split_layout.append((f"{pop_name}_a", node_ids[:half]))
        split_layout.append((f"{pop_name}_b", node_ids[half:]))

    with_split = with_maf_filter(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.05,
        seed=123,
        ploidy=2,
        counts_by_samples={"pop1_a": 4, "pop1_b": 3, "pop2_a": 2, "pop2_b": 1},
    )

    without_layout = with_maf_filter(
        demography, sample_sets, num_loci=10, maf=0.05, seed=123, ploidy=2
    )

    for split, plain in zip(with_split, without_layout, strict=True):
        assert list(split) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]  # discrimine
        assert split["pop1_a"] + split["pop1_b"] == plain["pop1"]  # recollement
        assert split["pop2_a"] + split["pop2_b"] == plain["pop2"]


# test relatifs à with_maf_filter_shared_ancestry


def test_with_maf_filter_shared_ancestry_no_filter_matches_direct_call(
    header_text_te5,
):
    """maf=0.0 doit produire EXACTEMENT le même résultat qu'un appel
    direct à simulate_shared_ancestry_loci + simulate_snp_genotypes
    (même graine pour les deux) -- garantit que
    with_maf_filter_shared_ancestry ne change rien à toy_example5, qui
    n'a pas de filtre MAF actif."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te5, scenario_index=1, seed=42
    )
    sex_ratio = parse_sex_ratio(OBSERVED_SNP_FILE_TE5)
    samples = build_male_only_samples_argument(OBSERVED_SNP_FILE_TE5)
    rescaled_demography = rescale_demography(
        demography, coalescence_coefficient("Y", sex_ratio) / 2
    )
    num_loci = 8

    direct = list(
        simulate_snp_genotypes(
            simulate_shared_ancestry_loci(
                rescaled_demography, samples, num_loci, seed=99, ploidy=2
            ),
            seed=99,
        )
    )
    via_filter = list(
        with_maf_filter_shared_ancestry(
            rescaled_demography, samples, num_loci, maf=0.0, seed=99, ploidy=2
        )
    )

    assert via_filter == direct


def test_with_maf_filter_shared_ancestry_rejects_low_maf_loci(header_text_te5):
    """Avec maf>0, chaque locus retourné doit respecter le seuil -- même
    contrat que with_maf_filter, mais ici sur une généalogie UNIQUE
    partagée entre tous les loci (reproduit particuleC.cpp:2424-2495 :
    le cache GeneTreeY est rempli avant le test MAF, donc un rejet ne
    redessine jamais l'arbre, seulement la mutation)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te5, scenario_index=1, seed=42
    )
    sex_ratio = parse_sex_ratio(OBSERVED_SNP_FILE_TE5)
    samples = build_male_only_samples_argument(OBSERVED_SNP_FILE_TE5)
    rescaled_demography = rescale_demography(
        demography, coalescence_coefficient("Y", sex_ratio) / 2
    )

    num_loci = 15
    maf = 0.2
    loci = list(
        with_maf_filter_shared_ancestry(
            rescaled_demography, samples, num_loci, maf=maf, seed=7, ploidy=1
        )
    )

    assert len(loci) == num_loci
    for locus_genotypes in loci:
        assert observed_maf(locus_genotypes) >= maf


def test_with_maf_filter_shared_ancestry_with_same_layout_matches_null_maf():
    """Vérifie que l'argument si l'argument layout correspond à compute_population_layout(ts) ne change pas le résultat de with_maf_filter."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    result_with_layout = with_maf_filter_shared_ancestry(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.0,
        seed=123,
        ploidy=1,
        counts_by_samples={"pop1": 7, "pop2": 3},
    )
    result_without_layout = with_maf_filter_shared_ancestry(
        demography, sample_sets, num_loci=10, maf=0.0, seed=123, ploidy=1
    )

    assert list(result_with_layout) == list(result_without_layout)


def test_with_maf_filter_shared_ancestry_with_same_layout_matches_nonnull_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter si le MAF est nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    result_with_layout = with_maf_filter_shared_ancestry(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.05,
        seed=123,
        ploidy=1,
        counts_by_samples={"pop1": 7, "pop2": 3},
    )
    result_without_layout = with_maf_filter_shared_ancestry(
        demography, sample_sets, num_loci=10, maf=0.05, seed=123, ploidy=1
    )

    assert list(result_with_layout) == list(result_without_layout)


def test_with_maf_filter_shared_ancestry_with_different_layout_null_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter_shared_ancestry si le MAF est nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=3,
        ploidy=1,
    )

    layout = compute_population_layout(ts)
    split_layout = []
    for pop_name, node_ids in layout:
        half = len(node_ids) // 2
        split_layout.append((f"{pop_name}_a", node_ids[:half]))
        split_layout.append((f"{pop_name}_b", node_ids[half:]))

    with_split = with_maf_filter_shared_ancestry(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.0,
        seed=123,
        ploidy=1,
        counts_by_samples={"pop1_a": 4, "pop1_b": 3, "pop2_a": 2, "pop2_b": 1},
    )

    without_layout = with_maf_filter_shared_ancestry(
        demography, sample_sets, num_loci=10, maf=0.0, seed=123, ploidy=1
    )

    for split, plain in zip(with_split, without_layout, strict=True):
        assert list(split) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]  # discrimine
        assert split["pop1_a"] + split["pop1_b"] == plain["pop1"]  # recollement
        assert split["pop2_a"] + split["pop2_b"] == plain["pop2"]


def test_with_maf_filter_shared_ancestry_with_different_layout_nonnull_maf():
    """Vérifie que l'argument layout ne change pas le résultat de with_maf_filter_shared_ancestry si le MAF est non nul, même si le layout est différent de compute_population_layout(ts)."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    sample_sets = [
        msprime.SampleSet(7, "pop1"),
        msprime.SampleSet(3, "pop2"),
    ]

    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=3,
        ploidy=1,
    )

    layout = compute_population_layout(ts)
    split_layout = []
    for pop_name, node_ids in layout:
        half = len(node_ids) // 2
        split_layout.append((f"{pop_name}_a", node_ids[:half]))
        split_layout.append((f"{pop_name}_b", node_ids[half:]))

    with_split = with_maf_filter_shared_ancestry(
        demography,
        sample_sets,
        num_loci=10,
        maf=0.05,
        seed=123,
        ploidy=1,
        counts_by_samples={"pop1_a": 4, "pop1_b": 3, "pop2_a": 2, "pop2_b": 1},
    )

    without_layout = with_maf_filter_shared_ancestry(
        demography, sample_sets, num_loci=10, maf=0.05, seed=123, ploidy=1
    )

    for split, plain in zip(with_split, without_layout, strict=True):
        assert list(split) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]  # discrimine
        assert split["pop1_a"] + split["pop1_b"] == plain["pop1"]  # recollement
        assert split["pop2_a"] + split["pop2_b"] == plain["pop2"]


# ------------------------------------------------------
# Tests sur les fonctions relatives aux reads Pool-seq
# ------------------------------------------------------


def test_reindex_reads_by_msprime_name():
    """Vérifie que reindex_reads_by_msprime_name renvoie bien un dict
    {nom_population: (derived_reads, total_reads)} avec les bons noms de
    populations, et que le nombre total de reads est correct."""
    observed_reads_te4 = observed_reads(OBSERVED_SNP_FILE_TE4)
    reindexed = _reindex_reads_by_msprime_name(
        observed_reads_te4, OBSERVED_SNP_FILE_TE4
    )[0]
    assert set(reindexed.keys()) == {"pop1", "pop2", "pop3", "pop4"}


def test_with_mrc_filter(header_text_te4):
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te4, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE_TE4)
    num_loci = 10
    mrc = 5
    observed_reads_per_locus = observed_reads(OBSERVED_SNP_FILE_TE4)
    observed_reads_per_locus = _reindex_reads_by_msprime_name(
        observed_reads_per_locus, OBSERVED_SNP_FILE_TE4
    )
    loci = list(
        with_mrc_filter(
            demography,
            samples,
            num_loci,
            mrc,
            observed_reads_per_locus,
            seed=7,
            ploidy=2,
        )
    )
    assert len(loci) == num_loci
    for reads_by_population in loci:
        sum_derived = sum(
            derived_reads for derived_reads, _ in reads_by_population.values()
        )
        sum_total = sum(total_reads for _, total_reads in reads_by_population.values())
        mrc_observed = (
            min(sum_derived, sum_total - sum_derived) if sum_total > 0 else 0.0
        )
        assert mrc_observed >= mrc


def test_simulate_poolseq_reads(header_text_te4):
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te4, scenario_index=1, seed=42
    )
    observed_reads_per_locus = observed_reads(OBSERVED_SNP_FILE_TE4)
    observed_reads_per_locus = _reindex_reads_by_msprime_name(
        observed_reads_per_locus, OBSERVED_SNP_FILE_TE4
    )

    n = 30

    def run():
        tree_sequences = simulate_independent_loci(
            demography,
            build_samples_argument(OBSERVED_SNP_FILE_TE4),
            num_loci=n,
            seed=123,
            ploidy=2,
        )
        return list(
            simulate_poolseq_reads(
                tree_sequences, observed_reads_per_locus[:n], seed=12
            )
        )

    results1 = run()
    results2 = run()
    assert results1 == results2, (
        "simulate_poolseq_reads should be deterministic with the same seed"
    )

    valeurs_pop1 = {r["pop1"] for r in results1}
    assert len(valeurs_pop1) > 1, (
        "simulate_poolseq_reads should produce different read counts for different loci"
    )


def test_simulate_poolseq_reads_with_mrc_filter(snp_context_te4):
    demography, _ = build_random_demography_for_scenario_index(
        snp_context_te4.header_text, scenario_index=1, seed=42
    )

    num_loci = snp_context_te4.loci_description.loci_counts_by_heritage["A"]
    mrc = snp_context_te4.mrc_ratio

    results = list(
        simulate_poolseq_reads_with_mrc_filter(
            demography,
            snp_context_te4,
            num_loci=num_loci,
            seed=12,
        )
    )

    assert len(results) == 100
    assert results[0].keys() == {"pop1", "pop2", "pop3", "pop4"}
    for reads_by_population in results:
        sum_derived = sum(
            derived_reads for derived_reads, _ in reads_by_population.values()
        )
        sum_total = sum(total_reads for _, total_reads in reads_by_population.values())
        mrc_observed = (
            min(sum_derived, sum_total - sum_derived) if sum_total > 0 else 0.0
        )
        assert mrc_observed >= mrc


# ------------------------------------------------------
# Tests sur les fonctions relatives aux séquences d'ADN
# ------------------------------------------------------


def test_transition_matrix_jk():
    """Vérifie que la matrice de transition est bien contruite pour les différents
    modèles de mutation (JK, K2P, HKY, TN) et que les paramètres sont corrects."""
    kappas = (2, 3)
    frequences_by_locus = {"pi_A": 0.1, "pi_C": 0.2, "pi_G": 0.3, "pi_T": 0.4}
    # test pour le modèle JK
    name_model = "JK"
    transition_matrix = build_transition_matrix(name_model, kappas, frequences_by_locus)
    expected_matrix = np.array(
        [
            [0, 1 / 3, 1 / 3, 1 / 3],
            [1 / 3, 0, 1 / 3, 1 / 3],
            [1 / 3, 1 / 3, 0, 1 / 3],
            [1 / 3, 1 / 3, 1 / 3, 0],
        ]
    )
    assert np.allclose(transition_matrix, expected_matrix)


def test_transition_matrix_k2p():
    """Vérifie que la matrice de transition est bien construite pour le modèle K2P."""
    kappas = (2, 3)
    frequences_by_locus = {"pi_A": 0.1, "pi_C": 0.2, "pi_G": 0.3, "pi_T": 0.4}
    # test pour le modèle K2P
    name_model = "K2P"
    transition_matrix = build_transition_matrix(name_model, kappas, frequences_by_locus)
    before_normalisation = np.array(
        [[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]]
    )
    expected_matrix = before_normalisation / before_normalisation.sum(
        axis=1, keepdims=True
    )
    assert np.allclose(transition_matrix, expected_matrix)


def test_transition_matrix_hky():
    """Vérifie que la matrice de transition est bien construite pour le modèle HKY."""
    kappas = (2, 3)
    frequences_by_locus = {"pi_A": 0.1, "pi_C": 0.2, "pi_G": 0.3, "pi_T": 0.4}
    # test pour le modèle HKY
    name_model = "HKY"
    before_normalisation = np.array(
        [[0, 0.2, 0.6, 0.4], [0.1, 0, 0.3, 0.8], [0.2, 0.2, 0, 0.4], [0.1, 0.4, 0.3, 0]]
    )
    expected_matrix = before_normalisation / before_normalisation.sum(
        axis=1, keepdims=True
    )
    transition_matrix = build_transition_matrix(name_model, kappas, frequences_by_locus)
    assert np.allclose(transition_matrix, expected_matrix)


def test_transition_matrix_tn():
    """Vérifie que la matrice de transition est bien construite pour le modèle TN."""
    kappas = (2, 3)
    frequences_by_locus = {"pi_A": 0.1, "pi_C": 0.2, "pi_G": 0.3, "pi_T": 0.4}
    # test pour le modèle TN
    name_model = "TN"
    before_normalisation = np.array(
        [[0, 0.2, 0.6, 0.4], [0.1, 0, 0.3, 1.2], [0.2, 0.2, 0, 0.4], [0.1, 0.6, 0.3, 0]]
    )
    expected_matrix = before_normalisation / before_normalisation.sum(
        axis=1, keepdims=True
    )
    transition_matrix = build_transition_matrix(name_model, kappas, frequences_by_locus)
    assert np.allclose(transition_matrix, expected_matrix)


def test_transition_matrix_invalid_model():
    """Vérifie que la fonction build_transition_matrix lève une exception pour un modèle invalide."""
    kappas = (2, 3)
    frequences_by_locus = {"pi_A": 0.1, "pi_C": 0.2, "pi_G": 0.3, "pi_T": 0.4}
    with pytest.raises(NotImplementedError, match="Modèle de"):
        build_transition_matrix("INVALID_MODEL", kappas, frequences_by_locus)


def test_count_loci_per_group(header_text_te2):
    """Vérifie que la fonction count_loci_per_group retourne le bon nombre de loci par groupe."""
    list_loci = parse_loci_description(header_text_te2)
    counts = count_loci_per_group(list_loci)
    assert counts == {"G1": 10, "G2": 5, "G3": 5}

    # Test avec un autre jeu de loci pour vérifier que l'erreur est bien levée.
    list_loci_invalid = [
        LociDescriptionDetailed(
            name="locus1",
            heritage="A",
            ms_or_seq="S",
            group="G1",
            motif_size=None,
            motif_range=None,
            dnalength=4,
        ),
        LociDescriptionDetailed(
            name="locus2",
            heritage="A",
            ms_or_seq="M",
            group="G1",
            motif_size=None,
            motif_range=None,
            dnalength=4,
        ),
    ]

    with pytest.raises(ValueError, match="Différents types de loci"):
        count_loci_per_group(list_loci_invalid)


def test_build_group_local_param_per_locus(header_text_te2):
    """Vérifie que la fonction build_group_local_param_per_locus retourne le bon dictionnaire
    de kappa1 et kappa2 par locus pour le fichier toy_example2 (dataset <A>+<M>
    avec 3 populations).
    Test de reproductibilité avec la même graine.
    Il manque un test pour vérifier lorsuqe le model est JK ou TN
    """
    params_per_locus = build_group_local_param_per_locus(header_text_te2, seed=42)

    assert len(params_per_locus) == 10
    assert len(params_per_locus["Locus_S_A_11_"]) == 3
    assert params_per_locus["Locus_S_A_11_"][1] == 0.0
    assert params_per_locus["Locus_S_A_11_"][0] > 0.0
    all_k1_values = [k[0] for k in params_per_locus.values()]
    assert len(set(all_k1_values)) == 10  # Tous les kappa1 sont différents

    all_mus_rate_values = [k[2] for k in params_per_locus.values()]
    assert len(set(all_mus_rate_values)) == 10  # Tous les mus_rate sont différents

    # test de reproductibilité avec la même graine
    params_per_locus_2 = build_group_local_param_per_locus(header_text_te2, seed=42)
    assert params_per_locus == params_per_locus_2


def test_build_matrix_per_locus(dna_context_te2):
    """Vérifie que la fonction build_matrix_per_locus retourne le bon dictionnaire
    de matrices de transition par locus pour le fichier toy_example2 (dataset <A>+<M>
    avec 3 populations).
    Test de reproductibilité avec la même graine.
    """
    matrix_per_locus = build_matrix_per_locus(dna_context_te2, seed=42)

    assert len(matrix_per_locus) == 10
    for matrix in matrix_per_locus.values():
        assert matrix.shape == (4, 4)
        assert np.allclose(matrix.sum(axis=1), 1.0)  # Chaque ligne doit sommer à 1

    # test de reproductibilité avec la même graine
    matrix_per_locus_2 = build_matrix_per_locus(dna_context_te2, seed=42)
    for locus in matrix_per_locus:
        assert np.allclose(matrix_per_locus[locus], matrix_per_locus_2[locus])


def test_build_rate_map():
    """Vérifie que la fonction build_rate_map retourne le bon dictionnaire
    de cartes de taux par locus pour le fichier toy_example2 (dataset <A>+<M>
    avec 3 populations).
    Test de reproductibilité avec la même graine.
    """
    # Test avec un exemple qui ne passerait pas
    mutsit = [0.1, 0.2, 0.3, 0.4]
    dnalength = 3
    with pytest.raises(ValueError, match="Le nombre de sites de mutation"):
        build_rate_map(mutsit, mus_rate=0.01, dnalength=dnalength)

    # Test avec un exemple correct
    dnalength = 4
    rate_map = build_rate_map(mutsit, mus_rate=0.01, dnalength=dnalength)
    assert all(
        rate_map.rate[i] == 0.01 * dnalength * mutsit[i] for i in range(dnalength)
    )


def test_build_rate_map_per_locus(header_text_te2):
    """Vérifie que la fonction build_rate_map_per_locus retourne le bon dictionnaire
    de cartes de taux par locus pour le fichier toy_example2 (dataset <A>+<M>
    avec 3 populations).
    Test de reproductibilité avec la même graine.
    """
    rate_map_per_locus = build_rate_map_per_locus(header_text_te2, seed=42)

    assert len(rate_map_per_locus) == 10
    for rate_map in rate_map_per_locus.values():
        assert isinstance(rate_map, msprime.RateMap)

    # le nombre de sites a taux nul doit correspondre a p_fixe du groupe
    # (Locus_S_A_11_ est dans G2, p_fixe=10, dnalength=100 -> nsv=90,
    # donc 10 sites fixes)
    group_priors = parse_group_priors(header_text_te2)
    gp_model_g2 = next(gp for gp in group_priors["G2"] if gp.model)
    dnalength = 100
    nsv = int(dnalength * (1 - 0.01 * gp_model_g2.p_fixe) + 0.5)
    nb_sites_fixes = dnalength - nsv
    rate_g2_locus1 = rate_map_per_locus["Locus_S_A_11_"].rate
    assert sum(1 for r in rate_g2_locus1 if r == 0.0) == nb_sites_fixes

    # deux loci differents ne doivent pas avoir le meme motif de taux
    rate_g2_locus2 = rate_map_per_locus["Locus_S_A_12_"].rate
    assert list(rate_g2_locus1) != list(rate_g2_locus2)

    # test de reproductibilité avec la même graine
    rate_map_per_locus_2 = build_rate_map_per_locus(header_text_te2, seed=42)
    for locus in rate_map_per_locus:
        assert rate_map_per_locus[locus] == rate_map_per_locus_2[locus]


def test_dna_mutation_simulation_per_locus(dna_context_te2, dna_context_te2_xy):
    """Vérifie que dna_mutation_simulation_per_locus produit bien une
    TreeSequence mutée par locus séquence (pas les loci microsat), avec
    une généalogie ET des mutations indépendantes d'un locus à l'autre
    (pas la même graine réutilisée partout), et reproductible avec la
    même graine de particule."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )

    mutated_tree_sequences = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    # 10 loci séquences (5 <A> + 5 <M>), pas les 10 loci microsat du même header
    assert len(mutated_tree_sequences) == 10
    assert set(mutated_tree_sequences.keys()) == {
        f"Locus_S_A_{11 + i}_" for i in range(5)
    } | {f"Locus_S_M_{16 + i}_" for i in range(5)}

    # Deux loci différents ne doivent pas partager la même généalogie ni
    # les mêmes positions de mutation -- sinon la graine par locus serait
    # réutilisée telle quelle (bug qu'on a corrigé plus tôt).
    ts1 = mutated_tree_sequences["Locus_S_A_11_"]
    ts2 = mutated_tree_sequences["Locus_S_A_12_"]
    assert ts1.tables.edges != ts2.tables.edges
    assert list(ts1.tables.sites.position) != list(ts2.tables.sites.position)

    # Reproductibilité : même graine de particule -> même résultat pour
    # chaque locus.
    mutated_tree_sequences_2 = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )
    for locus_name in mutated_tree_sequences:
        assert (
            mutated_tree_sequences[locus_name].tables.edges
            == mutated_tree_sequences_2[locus_name].tables.edges
        )
        assert list(mutated_tree_sequences[locus_name].tables.sites.position) == list(
            mutated_tree_sequences_2[locus_name].tables.sites.position
        )

    # test mais à partir de header_text_te2_xy pour vérifier que les loci <X> et <Y> sont bien supportés
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated_tree_sequences = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2_xy,
        seed=42,
    )

    # Vérifier que les loci <X> et <Y> sont bien simulés
    assert "Locus_S_A_11_" in mutated_tree_sequences
    assert "Locus_S_M_16_" in mutated_tree_sequences

    ts1 = mutated_tree_sequences["Locus_S_A_11_"]
    ts2 = mutated_tree_sequences["Locus_S_M_16_"]
    ts3 = mutated_tree_sequences["Locus_S_M_20_"]

    assert ts1.num_samples == 60
    assert ts2.num_samples == 20
    assert ts2.tables.edges == ts3.tables.edges


def test_ms_dna_ancestry_parameters_for_heritage(header_text_te2):
    """Vérifie le dispatch démographie/ploïdie par type d'héritage pour les
    loci ADN, exactement le même que celui de simulate_genotypes_for_locus_type
    côté SNP : <A> -- démographie inchangée, ploidy=2 ; <H>/<M> -- démographie
    rescalée par coalescence_coefficient/2, ploidy=1 ; <X>/<Y> -- pas
    supportés sur .mss (pas de sexe par individu dans ce format)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    sex_ratio = 0.4

    resolved_demography, ploidy = ms_dna_ancestry_parameters_for_heritage(
        "A", demography, sex_ratio
    )
    assert ploidy == 2
    assert resolved_demography is demography

    for heritage in ("H", "M", "X", "Y"):
        rescaled_demography, ploidy = ms_dna_ancestry_parameters_for_heritage(
            heritage, demography, sex_ratio
        )
        assert ploidy == 1
        factor = coalescence_coefficient(heritage, sex_ratio) / 2
        for pop, rescaled_pop in zip(
            demography.populations, rescaled_demography.populations, strict=True
        ):
            assert rescaled_pop.initial_size == pytest.approx(pop.initial_size * factor)


def test_dna_mutation_simulation_per_locus_ploidy_matches_heritage(dna_context_te2):
    """Vérifie que le nombre de lignées échantillonnées reflète bien la
    ploïdie attendue par héritage : un locus <A> (G2) doit avoir 2x plus de
    "samples" msprime qu'un locus <M> (G3) pour la même population -- avant
    la correction, les deux étaient simulés en ploidy=2 sans distinction."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )

    mutated_tree_sequences = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    ts_a = mutated_tree_sequences["Locus_S_A_11_"]
    ts_m = mutated_tree_sequences["Locus_S_M_16_"]
    assert ts_a.num_samples == 2 * ts_m.num_samples


# tests sur les valeurs de prior par groupe tirées par diyabc


def test_group_prior_values_from_columns(header_text_te2):
    group_priors_values = {
        "µmic_1": 0.0007375,
        "pmic_1": 0.2029,
        "snimic_1": 2.32e-07,
        "µseq_2": 4.068e-07,
        "k1seq_2": 2.684,
        "µseq_3": 8.112e-06,
        "k1seq_3": 13.42,
    }

    group_priors = parse_group_priors(header_text_te2)
    result = _group_prior_values_from_columns(group_priors_values, group_priors)
    assert result == {
        "G2": {"MEANMU": 4.068e-07, "MEANK1": 2.684},
        "G3": {"MEANMU": 8.112e-06, "MEANK1": 13.42},
    }


def test_build_group_local_param_per_locus_from_values(header_text_te2):
    group_priors_values = {
        "µmic_1": 0.0007375,
        "pmic_1": 0.2029,
        "snimic_1": 2.32e-07,
        "µseq_2": 4.068e-07,
        "k1seq_2": 2.684,
        "µseq_3": 8.112e-06,
        "k1seq_3": 13.42,
    }

    result = build_group_local_param_per_locus_from_values(
        header_text_te2, group_priors_values, seed=42
    )
    assert len(result) == 10
    assert len(result["Locus_S_A_11_"]) == 3
    assert result["Locus_S_A_11_"][1] == 0.0
    assert result["Locus_S_A_11_"][0] > 0.0


def test_build_sex_stratified_samples_argument_ms_dna(header_text_te2_XY):
    """Vérifie que build_sex_stratified_samples construit bien la liste de SampleSet attendue par msprime.sim_ancestry,"""
    liste_loci = parse_loci_description(header_text_te2_XY)
    locus_name = "Locus_M_A_1_"
    samples = build_sex_stratified_samples_argument_ms_dna(
        OBSERVED_MSS_FILE_TE2_XY, liste_loci, locus_name
    )

    assert len(samples) == 4
    for sample_set in samples:
        assert sample_set.population in {"pop1", "pop2"}
        assert sample_set.num_samples in {1, 19, 20, 0}  # 1 M ou 19 F par population
        assert sample_set.ploidy in {1, 2}  # M=1, F=2


def test_build_male_only_samples_argument_ms_dna(header_text_te2_XY):
    """Vérifie que build_male_only_samples_argument construit bien un dict {population: nombre_de_mâles} (PAS une liste de SampleSet, contrairement à build_sex_stratified_samples"""

    liste_loci = parse_loci_description(header_text_te2_XY)
    locus_name = "Locus_M_A_1_"
    samples = build_male_only_samples_argument_ms_dna(
        OBSERVED_MSS_FILE_TE2_XY, liste_loci, locus_name
    )
    assert samples == {"pop1": 1, "pop2": 0}


# ------------------------------------------------------
# Tests sur les fonctions relatives aux microsatellites
# ------------------------------------------------------
def test_distribution_from_position():
    """Vérifie que la fonction distribution_from_position retourne la bonne distribution de mutation pour un locus microsatellite donné."""
    # Test avec un exemple simple
    position = 4
    kmin = 1
    motif_size = 2

    matrix = np.array(
        np.arange(0, 12).reshape(3, 4)
    )  # matrice de transition simple pour le test

    distribution = _distribution_from_position(position, matrix, kmin, motif_size)

    assert isinstance(distribution, np.ndarray)
    assert distribution.shape == (4,)
    assert np.array_equal(distribution, np.array([4, 5, 6, 7]))


def test_place_gsm_row_on_dense_grid():
    """Vérifie que la fonction place_gsm_row_on_dense_grid retourne la bonne ligne de matrice de transition pour un locus microsatellite donné."""
    local_distribution = np.array([0.05, 0.05, 0.2, 0.3, 0.4])
    position = 3
    kmin = 0
    kmax = 10
    motif_size = 2
    result = _place_gsm_row_on_dense_grid(
        local_distribution, position, kmin, kmax, motif_size
    )

    assert np.allclose(
        result, np.array([0.0, 0.05, 0.0, 0.05, 0.0, 0.2, 0, 0.3, 0, 0.4, 0])
    )


def test_sni_row_on_dense_grid():
    """Vérifie que la fonction sni_row_on_dense_grid retourne la bonne ligne de matrice de transition pour un locus microsatellite donné."""

    kmin = 0
    kmax = 5

    # test sur une position intermédiaire
    position1 = 3
    result1 = _sni_row_on_dense_grid(position1, kmin, kmax)
    assert np.allclose(result1, np.array([0, 0, 0.5, 0, 0.5, 0]))

    # test sur une position à l'extrémité
    position2 = kmin
    result2 = _sni_row_on_dense_grid(position2, kmin, kmax)
    assert np.allclose(result2, np.array([0.5, 0.5, 0, 0, 0, 0]))

    position3 = kmax
    result3 = _sni_row_on_dense_grid(position3, kmin, kmax)
    assert np.allclose(result3, np.array([0, 0, 0, 0, 0.5, 0.5]))


def test_build_microsat_transition_matrix_with_sni():
    """Vérifie que la fonction build_microsat_transition_matrix_with_sni retourne la bonne matrice de transition pour les loci microsatellites."""
    kmin = 0
    kmax = 6
    motif_size = 2
    Pgeom = 0.2
    sni_rate = 0.1
    mut_rate = 0.5

    matrix = build_microsat_transition_matrix_with_sni(
        kmin=kmin,
        kmax=kmax,
        motif_size=motif_size,
        Pgeom=Pgeom,
        sni_rate=sni_rate,
        mut_rate=mut_rate,
    ).transition_matrix

    assert matrix.shape == (kmax - kmin + 1, kmax - kmin + 1)
    assert np.allclose(matrix.sum(axis=1), 1.0)  # Chaque ligne doit sommer à 1

    # test pour savoir si on obtient la mêm chose sans sni avec sni_rate =0
    model_without_sni = build_microsat_transition_matrix(
        kmin=kmin,
        kmax=kmax,
        motif_size=motif_size,
        Pgeom=Pgeom,
        epsilon=1e-16,
    )
    matrix_without_sni = model_without_sni.transition_matrix

    sni_rate = 0
    matrix_with_sni = build_microsat_transition_matrix_with_sni(
        kmin=kmin,
        kmax=kmax,
        motif_size=motif_size,
        Pgeom=Pgeom,
        sni_rate=sni_rate,
        mut_rate=mut_rate,
    ).transition_matrix

    alleles = model_without_sni.alleles
    numeric_allele = [int(a) for a in alleles]
    expected_index = [numeric_allele[i] - kmin for i in range(len(numeric_allele))]
    ixgrid = np.ix_(expected_index, expected_index)

    extracted_matrix_with_sni = matrix_with_sni[ixgrid]

    assert np.allclose(matrix_without_sni, extracted_matrix_with_sni)


def test_build_microsat_transition_matrix():
    """Vérifie que la fonction build_microsat_transition_matrix retourne la bonne matrice de transition pour les loci microsatellites."""

    model = build_microsat_transition_matrix(
        kmin=162,
        kmax=240,
        motif_size=2,
        Pgeom=0.2,
        epsilon=1e-16,
    )

    assert len(model.alleles) == 39
    assert np.allclose(
        model.transition_matrix.sum(axis=1), 1.0
    )  # Chaque ligne doit sommer à 1

    # root=201, n_minus=19 pour kmin=162, kmax=240, motif_size=2 -- calculé à la main,
    assert model.alleles[19] == "201"

    # On va tester avec un Pgeom limite
    model2 = build_microsat_transition_matrix(
        kmin=162,
        kmax=240,
        motif_size=2,
        Pgeom=0,
        epsilon=1e-16,
    )

    assert np.isclose(model2.transition_matrix[19, 18], 0.5)
    assert np.isclose(model2.transition_matrix[19, 20], 0.5)
    assert np.isclose(model2.transition_matrix.sum(axis=1)[19], 1.0)

    model3 = build_microsat_transition_matrix(
        kmin=162,
        kmax=240,
        motif_size=2,
        Pgeom=1,
        epsilon=1e-16,
    )

    assert np.isclose(model3.transition_matrix[19, 19], 0.0)
    row_without_root = np.delete(model3.transition_matrix[19, :], 19)
    assert np.allclose(row_without_root, 1.0 / 38)


def test_build_microsat_local_param_per_locus(header_text_te2_XY):
    """Vérifie que la fonction buiomparaison directe de forme), le raleld_microsat_local_param_per_locus retourne le bon dictionnaire
    Test de reproductibilité avec la même graine.
    """
    params_per_locus = build_microsat_local_param_per_locus(header_text_te2_XY, seed=42)

    list_loci = parse_loci_description(header_text_te2_XY)

    # Vérification qu'il n'y a que des loci microsatellites dans le dictionnaire
    for locus in params_per_locus:
        locus_type = next(
            (
                locus_desc.ms_or_seq
                for locus_desc in list_loci
                if locus_desc.name == locus
            ),
            None,
        )
        assert locus_type == "M", f"Locus {locus} n'est pas un microsatellite"

    assert len(params_per_locus) == 10
    all_mus_rate_values = [k[0] for k in params_per_locus.values()]
    assert len(set(all_mus_rate_values)) == 10  # Tous les mus_rate sont différents
    all_Pgeom_values = [k[1] for k in params_per_locus.values()]
    assert len(set(all_Pgeom_values)) == 10  # Tous les Pgeom sont différents
    all_sni_rate_values = [k[2] for k in params_per_locus.values()]
    assert len(set(all_sni_rate_values)) == 10  # Tous les sni_rate sont différents

    for locus in params_per_locus:
        assert len(params_per_locus[locus]) == 3
        assert params_per_locus[locus][0] > 0.0
        assert params_per_locus[locus][1] >= 0.0 and params_per_locus[locus][1] <= 1.0
        assert params_per_locus[locus][2] >= 0.0

    assert params_per_locus["Locus_M_A_1_"][0] == 0.00010313460804143706
    assert params_per_locus["Locus_M_A_1_"][1] == 0.4820860729990509
    assert params_per_locus["Locus_M_A_1_"][2] == 3.1999081859003997e-07

    # test de reproductibilité avec la même graine
    params_per_locus_2 = build_microsat_local_param_per_locus(
        header_text_te2_XY, seed=42
    )
    assert params_per_locus == params_per_locus_2


def test_build_matrix_microsat_per_locus(microsat_context_te2_xy):
    """Vérifie que la fonction build_matrix_microsat_per_locus retourne le bon dictionnaire"""
    list_loci = microsat_context_te2_xy.list_loci
    matrix_per_locus = build_matrix_microsat_per_locus(microsat_context_te2_xy, seed=42)

    for locus in matrix_per_locus:
        locus_type = next(
            (
                locus_desc.ms_or_seq
                for locus_desc in list_loci
                if locus_desc.name == locus
            ),
            None,
        )
        assert locus_type == "M", f"Locus {locus} n'est pas un microsatellite"
    assert len(matrix_per_locus) == 10
    assert matrix_per_locus["Locus_M_A_1_"].transition_matrix.shape == (79, 79)

    matrix_per_locus2 = build_matrix_microsat_per_locus(
        microsat_context_te2_xy, seed=42
    )
    for locus in matrix_per_locus:
        assert np.allclose(
            matrix_per_locus[locus].transition_matrix,
            matrix_per_locus2[locus].transition_matrix,
        )
        assert np.array_equal(
            matrix_per_locus[locus].alleles, matrix_per_locus2[locus].alleles
        )
        assert np.allclose(
            matrix_per_locus[locus].root_distribution,
            matrix_per_locus2[locus].root_distribution,
        )


def test_microsat_mutation_simulation_per_locus(microsat_context_te2_xy):
    """Vérifie que microsat_mutation_simulation_per_locus produit bien une
    TreeSequence mutée par locus microsat, avec une généalogie ET des
    mutations indépendantes d'un locus à l'autre (pas la même graine
    réutilisée partout), et reproductible avec la même graine de particule."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )

    mutated_tree_sequences = microsat_mutation_simulation_per_locus(
        context=microsat_context_te2_xy,
        demography=demography,
        seed=42,
    )

    # 10 loci microsatellites (5 <A> + 5 <M>), pas les 10 loci séquences du même header
    assert len(mutated_tree_sequences) == 10
    assert set(mutated_tree_sequences.keys()) == {
        f"Locus_M_A_{1 + i}_" for i in range(10)
    }

    # Deux loci différents ne doivent pas partager la même généalogie ni
    # les mêmes positions de mutation -- sinon la graine par locus serait
    # réutilisée telle quelle (bug qu'on a corrigé plus tôt).
    ts3 = mutated_tree_sequences["Locus_M_A_3_"]
    ts4 = mutated_tree_sequences["Locus_M_A_4_"]
    assert ts3.tables.edges != ts4.tables.edges
    assert np.array_equal(ts3.genotype_matrix(), ts4.genotype_matrix()) is False

    # Reproductibilité : même graine de particule -> même résultat pour
    # chaque locus.
    mutated_tree_sequences_2 = microsat_mutation_simulation_per_locus(
        context=microsat_context_te2_xy,
        demography=demography,
        seed=42,
    )
    for locus_name in mutated_tree_sequences:
        assert (
            mutated_tree_sequences[locus_name].tables.edges
            == mutated_tree_sequences_2[locus_name].tables.edges
        )
        assert np.array_equal(
            mutated_tree_sequences[locus_name].genotype_matrix(),
            mutated_tree_sequences_2[locus_name].genotype_matrix(),
        )


def test_microsat_mutation_simulation_per_locus_ploidy_matches_heritage(
    microsat_context_te2_xy,
):
    """Vérifie que le nombre de lignées échantillonnées reflète bien la
    ploïdie attendue par héritage : un locus <A> (G2) doit avoir 2x plus de
    "samples" msprime qu'un locus <M> (G3) pour la même population -- avant
    la correction, les deux étaient simulés en ploidy=2 sans distinction."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )

    mutated_tree_sequences = microsat_mutation_simulation_per_locus(
        context=microsat_context_te2_xy,
        demography=demography,
        seed=42,
    )

    list_loci = microsat_context_te2_xy.list_loci

    ts_A = mutated_tree_sequences["Locus_M_A_3_"]
    ts_M = mutated_tree_sequences["Locus_M_A_2_"]
    ts_Y = mutated_tree_sequences["Locus_M_A_10_"]

    num_samples_male_only = sum(
        build_male_only_samples_argument_ms_dna(
            microsat_context_te2_xy.mss_path, list_loci, "Locus_M_A_10_"
        ).values()
    )
    assert ts_A.num_samples == 2 * ts_M.num_samples
    assert ts_Y.num_samples == num_samples_male_only
