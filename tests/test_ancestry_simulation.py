"""Vérifie ancestry_simulation : construction de l'argument samples pour
msprime.sim_ancestry, simulation de généalogies indépendantes par locus,
et mutation Hudson (exactement une mutation par locus, toujours
polymorphe)."""

import msprime
import numpy as np
import pytest
from conftest import (
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE2,
    OBSERVED_SNP_FILE_TE4,
    OBSERVED_SNP_FILE_TE5,
)

from bridge.ancestry_simulation import (
    _reindex_reads_by_msprime_name,
    build_group_local_param_per_locus,
    build_male_only_samples_argument,
    build_matrix_per_locus,
    build_rate_map,
    build_rate_map_per_locus,
    build_samples_argument,
    build_sex_stratified_samples_argument,
    build_transition_matrix,
    count_loci_per_group,
    dna_mutation_simulation_per_locus,
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
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    coalescence_coefficient,
    observed_reads,
    parse_mrc_ratio,
    parse_sex_ratio,
)
from bridge.pipeline import build_random_demography_for_scenario_index
from bridge.prior_parser import parse_group_priors
from bridge.scenario_types import LociDescriptionDetailed


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


def test_simulate_genotypes_for_locus_type(header_text_te5):
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te5, scenario_index=1, seed=42
    )
    for locus_type in ["A", "X", "Y", "M", "H"]:
        genotypes = simulate_genotypes_for_locus_type(
            demography=demography,
            locus_type=locus_type,
            snp_file_path=OBSERVED_SNP_FILE_TE5,
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


def test_simulate_poolseq_reads_with_mrc_filter(header_text_te4):
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te4, scenario_index=1, seed=42
    )

    num_loci = parse_loci_description(header_text_te4).total_loci["A"]
    mrc = parse_mrc_ratio(OBSERVED_SNP_FILE_TE4)

    results = list(
        simulate_poolseq_reads_with_mrc_filter(
            demography,
            OBSERVED_SNP_FILE_TE4,
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


def test_build_matrix_per_locus(header_text_te2):
    """Vérifie que la fonction build_matrix_per_locus retourne le bon dictionnaire
    de matrices de transition par locus pour le fichier toy_example2 (dataset <A>+<M>
    avec 3 populations).
    Test de reproductibilité avec la même graine.
    """
    matrix_per_locus = build_matrix_per_locus(
        header_text_te2, OBSERVED_SNP_FILE_TE2, seed=42
    )

    assert len(matrix_per_locus) == 10
    for matrix in matrix_per_locus.values():
        assert matrix.shape == (4, 4)
        assert np.allclose(matrix.sum(axis=1), 1.0)  # Chaque ligne doit sommer à 1

    # test de reproductibilité avec la même graine
    matrix_per_locus_2 = build_matrix_per_locus(
        header_text_te2, OBSERVED_SNP_FILE_TE2, seed=42
    )
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


def test_dna_mutation_simulation_per_locus(header_text_te2):
    """Vérifie que dna_mutation_simulation_per_locus produit bien une
    TreeSequence mutée par locus séquence (pas les loci microsat), avec
    une généalogie ET des mutations indépendantes d'un locus à l'autre
    (pas la même graine réutilisée partout), et reproductible avec la
    même graine de particule."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )

    mutated_tree_sequences = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
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
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
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
