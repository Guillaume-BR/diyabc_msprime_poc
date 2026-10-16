import dataclasses
import itertools

import msprime
import numpy as np
import pytest
import tskit

from bridge.ancestry_simulation import (
    build_transition_matrix_microsat,
    compute_population_layout,
    compute_sample_layout,
    dna_mutation_simulation_per_locus,
    microsat_mutation_simulation_per_locus,
)
from bridge.pipeline import build_random_demography_for_scenario_index
from bridge.summary_statistics import (
    FocalGenotypeFrequencies,
    _allele_counts_by_sample,
    _compute_DM2_for_one_locus,
    _compute_FST_constants_for_two_samples_combined,
    _compute_FST_constants_on_all_alleles_for_two_samples,
    _compute_H2P_for_one_pair,
    _compute_identical_pair_for_one_pair,
    _compute_LIK_for_one_locus,
    _compute_MGW_by_locus,
    _compute_N2P_for_one_locus,
    _compute_N2P_for_one_pair,
    _compute_num_den_lik_for_one_individual,
    _compute_V2P_constants,
    _compute_VAR_constants,
    _compute_VAR_for_one_sample,
    _count_distinct_haplotypes,
    _focal_genotype_frequencies,
    _genotype_matrix_by_sample,
    _genotypes_by_sample_and_individuals,
    _individuals_by_sample,
    _length_by_sample,
    _length_by_sample_and_individuals,
    _log_likelihood_admixture,
    _pairwise_distances_per_locus,
    _pairwise_hamming_distances,
    _parental_frequencies_with_cache,
    _pool_allele_counts_for_two_samples,
    _prepare_matrices_poolseq,
    _prepared_parental_frequencies,
    _vectorize_parent_for_focal,
    compute_all_statistics_dna,
    compute_all_statistics_microsat,
    compute_all_statistics_poolseq,
    compute_DAS,
    compute_DM2,
    compute_DTA,
    compute_FST,
    compute_H2P,
    compute_HET,
    compute_HST,
    compute_LIK,
    compute_MGW,
    compute_MNS,
    compute_MP2,
    compute_MPB,
    compute_MPD,
    compute_N2P,
    compute_NAL,
    compute_NH2,
    compute_NHA,
    compute_NS2,
    compute_NSS,
    compute_PSS,
    compute_V2P,
    compute_VAR,
    compute_VNS,
    compute_VPD,
    count_alleles_per_sample,
    total_genes_copies_per_sample,
)


def test_prepare_matrices_poolseq():
    reads_per_locus = [
        {"POP1": (0, 93), "POP2": (0, 100), "POP3": (1, 116), "POP4": (0, 139)},
        {"POP1": (1, 80), "POP2": (0, 90), "POP3": (0, 110), "POP4": (1, 120)},
    ]
    sample_names = ["POP1", "POP2", "POP3", "POP4"]
    counts, ns, freq0, freq1 = _prepare_matrices_poolseq(reads_per_locus, sample_names)

    assert counts.shape == (4, 2)
    assert np.array_equal(counts, np.array([[0, 1], [0, 0], [1, 0], [0, 1]]))


def test_compute_all_statistics_poolseq():
    """Vérifie que compute_all_statistics_poolseq renvoie bien un dictionnaire
    de statistiques pour les fichiers POOLSEQ toy_example4 et toy_example5."""
    # Préparer des données fictives pour le test
    reads_per_locus = [
        {"POP1": (0, 93), "POP2": (0, 100), "POP3": (1, 116), "POP4": (0, 139)},
        {"POP1": (1, 80), "POP2": (0, 90), "POP3": (0, 110), "POP4": (1, 120)},
    ]
    sample_names = ["POP1", "POP2", "POP3", "POP4"]
    pool_sizes = {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}

    results = compute_all_statistics_poolseq(reads_per_locus, sample_names, pool_sizes)

    # Vérifier que les résultats contiennent les clés attendues
    expected_keys = {
        "HWm_1",
        "FST1m_1",
        "AMLm_1.2.3",
        # Ajouter d'autres statistiques attendues ici si nécessaire
    }
    assert expected_keys.issubset(results.keys())

    assert len(results) == 130
    assert abs(results["FST1m_1"]) <= 1, "FST1m_1 should be between 0 and 1"
    assert abs(results["HWm_1"] - 0.0126) < 1e-4, "HWm_1 should be approximately 0.0126"
    assert abs(results["HWv_1"] - 0.0003) < 1e-4, "HWv_1 should be approximately 0.0003"


def test_genotype_matrix_by_sample(dna_context_te2):
    """Vérifie _genotype_matrix_by_sample sur un locus <A> (diploïde) et un
    locus <M> (haploïde) du même dataset : la forme retournée doit respecter le
    nombre de sites/samples réels de la TreeSequence, sans perte ni doublon de
    sample entre populations, et le nombre de samples par population doit
    refléter la ploïdie du locus (rapport 2:1 entre <A> et <M> -- couvre le bug
    de ploïdie corrigé le 2026-08-24 dans
    dna_ancestry_parameters_for_heritage)."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    def samples_nodes_per_population(ts):
        return {
            population.metadata["name"]: ts.samples(population=population.id)
            for population in ts.populations()
        }

    for locus_name in ("Locus_S_A_11_", "Locus_S_M_16_"):
        ts = mutated[locus_name]
        result = _genotype_matrix_by_sample(ts)
        expected_samples = samples_nodes_per_population(ts)

        assert result.keys() == {"pop1", "pop2"}
        total_samples = 0
        for sample_name, matrix in result.items():
            assert matrix.shape[0] == ts.num_sites
            assert matrix.shape[1] == len(expected_samples[sample_name])
            total_samples += matrix.shape[1]
        assert total_samples == ts.num_samples

    samples_a = _genotype_matrix_by_sample(mutated["Locus_S_A_11_"])
    samples_m = _genotype_matrix_by_sample(mutated["Locus_S_M_16_"])
    for sample_name in ("pop1", "pop2"):
        assert samples_a[sample_name].shape[1] == 2 * samples_m[sample_name].shape[1]


def test_genotype_matrix_by_sample_with_same_layout(dna_context_te2):
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    first_key = next(iter(mutated))
    ts = mutated[first_key]
    layout = compute_population_layout(ts)

    result_with_layout = _genotype_matrix_by_sample(ts, layout=layout)
    result_without_layout = _genotype_matrix_by_sample(ts)

    assert [sample_name for sample_name in result_with_layout] == [
        sample_name for sample_name in result_without_layout
    ]
    for sample_name in result_with_layout:
        assert np.array_equal(
            result_with_layout[sample_name], result_without_layout[sample_name]
        )


def _two_population_demography():
    """Démographie minimale : pop1 et pop2 fusionnent à t = 500 dans `anc`."""
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    return demography


def _mutated_two_population_ts(seed=3):
    """TreeSequence mutée minimale : deux populations de 10 individus
    diploïdes.

    Bien moins coûteuse qu'un `dna_mutation_simulation_per_locus`
    complet (qui simule les 10 loci du header) quand le test n'a besoin
    que d'UNE TreeSequence. La population ancestrale est indispensable :
    sans elle, les lignées des deux populations ne peuvent jamais
    coalescer et msprime lève `Infinite waiting time until next
    simulation event`.
    """
    return msprime.sim_mutations(
        msprime.sim_ancestry(
            samples=[msprime.SampleSet(10, "pop1"), msprime.SampleSet(10, "pop2")],
            demography=_two_population_demography(),
            sequence_length=100,
            random_seed=seed,
        ),
        rate=1e-3,
        random_seed=seed,
    )


def test_genotype_matrix_by_sample_with_different_layout():
    """Un layout qui change le REGROUPEMENT, pas seulement l'ordre des clés.

    Chaque population est coupée en deux : la sortie doit contenir
    quatre groupes de taille moitié, ce que le découpage par population
    ne peut structurellement pas produire. Une helper qui ignorerait son
    paramètre `layout` échouerait dès le nombre de clés -- contrairement
    à un layout simplement réordonné, qui laisse l'association nom ->
    noeuds intacte et ne discrimine donc presque rien.

    C'est la forme exacte du cas sériel, où plusieurs échantillons
    partagent une seule population msprime.
    """
    ts = _mutated_two_population_ts()
    default_layout = compute_population_layout(ts)

    # Les bornes de coupe doivent tomber sur des frontières d'individus
    # (ici multiples de 2, ploïdie diploïde) : une coupe au milieu d'un
    # individu l'attribuerait en entier au groupe de son premier noeud,
    # silencieusement -- voir _length_by_sample_and_individuals, qui retrouve
    # la population d'un individu par `nodes[0] in inds`.
    split_layout = []
    for sample_name, node_ids in default_layout:
        half = len(node_ids) // 2
        split_layout.append((f"{sample_name}_a", node_ids[:half]))
        split_layout.append((f"{sample_name}_b", node_ids[half:]))

    without_layout = _genotype_matrix_by_sample(ts)
    with_split_layout = _genotype_matrix_by_sample(ts, layout=split_layout)

    assert list(without_layout) == ["pop1", "pop2"]
    assert list(with_split_layout) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]

    for sample_name, node_ids in split_layout:
        assert with_split_layout[sample_name].shape == (ts.num_sites, len(node_ids))

    # Recoller les deux moitiés d'une population redonne exactement la
    # matrice du découpage par population : le layout change le groupement,
    # jamais le contenu.
    for sample_name in ("pop1", "pop2"):
        rejoined = np.hstack(
            [
                with_split_layout[f"{sample_name}_a"],
                with_split_layout[f"{sample_name}_b"],
            ]
        )
        assert np.array_equal(rejoined, without_layout[sample_name])


def test_mean_segregating_sites_per_group(dna_context_te2):
    """Vérifie compute_NSS séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre NSS_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NSS(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"pop1": pytest.approx(9.6), "pop2": pytest.approx(9.0)}

    mean_g3 = compute_NSS(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"pop1": pytest.approx(2.6), "pop2": pytest.approx(3.2)}


def test_mean_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    sample_names = ["pop1", "pop2"]
    assert compute_NSS([], sample_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group(dna_context_te2):
    """Vérifie compute_NHA séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre NH_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NHA(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"pop1": pytest.approx(7.4), "pop2": pytest.approx(7.2)}

    mean_g3 = compute_NHA(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"pop1": pytest.approx(3.2), "pop2": pytest.approx(3.8)}


def test_mean_distinct_haplotypes_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    sample_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = compute_NHA(tree_sequences_empty, sample_names)
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def _reference_count_distinct_haplotypes(matrix):
    """Version d'origine de `_count_distinct_haplotypes` (déduplication des
    colonnes par `np.unique(axis=1)`), gardée comme oracle : la version
    actuelle code chaque colonne par un entier en base K, plus rapide, et doit
    donner exactement les mêmes décomptes."""
    return np.unique(matrix, axis=1).shape[1]


def test_count_distinct_haplotypes_matches_reference_on_random_matrices():
    """La version rapide doit égaler l'oracle sur des matrices aléatoires (graine
    fixe) qui couvrent les tailles du pipeline (20, 40 et 80 colonnes), des
    alphabets avec trous (`{0, 3}`), des états négatifs (donnée manquante), les
    trois types d'entiers, des colonnes dupliquées, et des nombres de sites
    de part et d'autre du seuil de repli (30 à 33 sites avec 4 états)."""
    rng = np.random.default_rng(12345)
    n_sites_choices = [0, 1, 2, 3, 5, 11, 20, 30, 31, 32, 33, 40, 60, 120]
    n_samples_choices = [1, 2, 5, 20, 40, 80]
    dtypes = [np.int8, np.int32, np.int64]
    for _ in range(1500):
        n_sites = int(rng.choice(n_sites_choices))
        n_samples = int(rng.choice(n_samples_choices))
        n_states = int(rng.integers(1, 5))
        states = rng.choice(4, size=n_states, replace=False)
        if rng.random() < 0.15:
            states = states - 1  # inclut -1 (donnée manquante)
        dtype = dtypes[int(rng.integers(len(dtypes)))]
        matrix = rng.choice(states, size=(n_sites, n_samples)).astype(dtype)
        if n_samples > 2 and rng.random() < 0.3:
            matrix[:, 1] = matrix[:, 0]  # deux haplotypes identiques

        result = _count_distinct_haplotypes(matrix)
        expected = _reference_count_distinct_haplotypes(matrix)
        assert result == expected, (
            f"forme {matrix.shape}, type {matrix.dtype}, états {sorted(states)} : "
            f"{result} au lieu de {expected}"
        )
        assert type(result) is int


def test_count_distinct_haplotypes_known_values():
    """Valeurs dérivées à la main."""
    # Toutes les colonnes identiques : 1 haplotype.
    assert _count_distinct_haplotypes(np.ones((5, 4), dtype=np.int32)) == 1
    # Toutes les colonnes différentes : autant d'haplotypes que de colonnes.
    assert _count_distinct_haplotypes(np.eye(4, dtype=np.int32)) == 4
    # Cas qui a révélé une erreur de codage : les colonnes (2, 0) et (0, 1)
    # sont distinctes, mais coderaient le même nombre si la base était le
    # nombre de sites (2) au lieu du nombre d'états (3).
    assert _count_distinct_haplotypes(np.array([[2, 0], [0, 1]])) == 2
    # Une seule colonne : un seul haplotype.
    assert _count_distinct_haplotypes(np.array([[0], [3], [1]])) == 1


def test_count_distinct_haplotypes_without_site_or_sample():
    """Un locus sans mutation (aucun site) a un seul haplotype, vide ; une
    matrice sans échantillon est une erreur."""
    assert _count_distinct_haplotypes(np.zeros((0, 7), dtype=np.int32)) == 1
    with pytest.raises(ValueError, match="vide"):
        _count_distinct_haplotypes(np.zeros((5, 0), dtype=np.int32))


def test_count_distinct_haplotypes_around_the_overflow_threshold():
    """Avec 4 états, 4**30 tient dans un entier 64 bits mais 4**32 non : la
    fonction bascule alors sur la déduplication des colonnes. Les deux côtés
    du seuil doivent donner le même résultat que l'oracle (un dépassement non
    détecté confondrait des haplotypes sans lever d'erreur)."""
    assert 4**30 <= 2**62 < 4**32  # le test encadre bien le seuil
    rng = np.random.default_rng(7)
    for n_sites in (30, 31, 32, 33):
        matrix = rng.integers(0, 4, size=(n_sites, 40)).astype(np.int32)
        matrix[:, 5] = matrix[:, 4]  # au moins un doublon
        assert _count_distinct_haplotypes(matrix) == (
            _reference_count_distinct_haplotypes(matrix)
        ), f"{n_sites} sites"


def test_count_distinct_haplotypes_int8_extremes():
    """Soustraire le minimum dans le type d'origine déborderait pour `int8` :
    les bornes -128 et 127 doivent donner le bon décompte."""
    matrix = np.array([[127, -128, 127], [-128, 127, -128]], dtype=np.int8)
    assert _count_distinct_haplotypes(matrix) == 2
    assert _count_distinct_haplotypes(np.full((3, 4), -128, dtype=np.int8)) == 1


def _reference_pairwise_hamming_distances(matrix: np.ndarray) -> np.ndarray:
    """Calcule les distances de Hamming par paire d'échantillons.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Un vecteur 1D de longueur C(n_samples, 2) -- une valeur par
        paire (i, j) avec i < j, pas la matrice carrée
        (n_samples, n_samples) complète.

    Raises:
        ValueError: Si matrix.shape[1] == 0 (échantillon sans échantillon).
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")

    # Oracle volontairement naïf : comparaison de toutes les paires de colonnes
    # par diffusion, et indices du triangle supérieur recalculés à chaque appel
    # (sans le cache ni le produit matriciel de la version rapide).
    matrix_hamming = (matrix[:, :, None] != matrix[:, None, :]).sum(axis=0)
    triangle_sup = np.triu_indices(matrix_hamming.shape[0], k=1)
    return matrix_hamming[triangle_sup]


def compare_pairwise_hamming_distances(matrix: np.ndarray) -> None:
    """Compare les distances de Hamming par paire calculées par la fonction
    de référence et la fonction optimisée.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Raises:
        AssertionError: Si les résultats des deux fonctions diffèrent.
    """
    reference_distances = _reference_pairwise_hamming_distances(matrix)
    optimized_distances = _pairwise_hamming_distances(matrix)
    assert np.array_equal(reference_distances, optimized_distances), (
        f"Les distances de Hamming par paire diffèrent pour la matrice :\n{matrix}\n"
        f"Référence : {reference_distances}\nOptimisé : {optimized_distances}"
    )
    # Même forme et même type : les statistiques en aval (somme, moyenne,
    # variance) dépendent du type entier de ce vecteur.
    assert optimized_distances.shape == reference_distances.shape
    assert optimized_distances.dtype == reference_distances.dtype, (
        f"type {optimized_distances.dtype} au lieu de {reference_distances.dtype}"
    )


def test_pairwise_hamming_distances_known_values():
    """Vérifie que la fonction optimisée _pairwise_hamming_distances renvoie les
    mêmes résultats que la fonction de référence pour des matrices de génotypes
    connues."""
    # Cas simple : 3 échantillons avec 2 sites
    matrix1 = np.array([[0, 1, 0], [1, 0, 1]])
    compare_pairwise_hamming_distances(matrix1)

    # Cas avec des colonnes identiques : 4 échantillons avec 3 sites
    matrix2 = np.array([[0, 0, 1, 1], [1, 1, 0, 0], [0, 0, 1, 1]])
    compare_pairwise_hamming_distances(matrix2)

    # Cas sans aucun échantillon (zéro colonne) : doit lever une exception
    with pytest.raises(ValueError, match="La matrice de génotypes est vide"):
        _pairwise_hamming_distances(np.array([[]]))


def test_pairwise_hamming_distances_hand_derived_values():
    """Valeurs calculées à la main, sans passer par l'oracle."""
    # 3 sites, 3 échantillons : les haplotypes sont (0,2,0), (1,1,0) et (1,2,0).
    # Paires (0,1), (0,2), (1,2) : ils diffèrent aux sites {0,1}, {0}, {1}.
    matrix = np.array([[0, 1, 1], [2, 1, 2], [0, 0, 0]])
    assert _pairwise_hamming_distances(matrix).tolist() == [2, 1, 1]
    # Aucun site (locus sans mutation) : les 4 gènes sont identiques, donc les
    # 4*3/2 = 6 paires sont à distance 0 (et non un vecteur vide).
    no_site = _pairwise_hamming_distances(np.zeros((0, 4), dtype=np.int32))
    assert no_site.tolist() == [0, 0, 0, 0, 0, 0]
    assert no_site.dtype == np.int64
    # Un seul échantillon : aucune paire.
    single = _pairwise_hamming_distances(np.ones((5, 1), dtype=np.int32))
    assert single.shape == (0,)
    assert single.dtype == np.int64


def test_pairwise_hamming_distances_states_not_starting_at_zero():
    """Les états n'ont pas à commencer à 0 (ni à être consécutifs) : les
    distances ne dépendent que de l'égalité entre états. Un calcul qui
    comparerait les états décalés aux indices 0, 1, ... de la matrice
    d'origine ne verrait aucun accord."""
    # Colonnes (5,7), (7,7), (5,5) : paires (0,1), (0,2), (1,2) -> 1, 1 et 2
    # (la dernière paire diffère aux deux sites).
    matrix = np.array([[5, 7, 5], [7, 7, 5]])
    assert _pairwise_hamming_distances(matrix).tolist() == [1, 1, 2]
    # Mêmes distances avec des états négatifs (donnée manquante -1).
    matrix = np.array([[-1, 1, -1], [1, 1, -1]])
    assert _pairwise_hamming_distances(matrix).tolist() == [1, 1, 2]
    compare_pairwise_hamming_distances(matrix)
    compare_pairwise_hamming_distances(matrix + 3)


def test_pairwise_hamming_distances_matches_reference_on_random_matrices():
    """La version rapide doit égaler l'oracle (valeurs, forme et type) sur des
    matrices aléatoires à graine fixe : tailles du pipeline (20, 40 et 80
    colonnes), 1 à 12 états (de part et d'autre du seuil de 8 états qui change
    de méthode), alphabets à trous, états négatifs, et les trois types
    d'entiers."""
    rng = np.random.default_rng(2024)
    n_sites_choices = [0, 1, 2, 5, 11, 30, 60, 120]
    n_samples_choices = [1, 2, 5, 20, 40, 80]
    n_states_choices = [1, 2, 3, 4, 6, 8, 9, 12]
    dtypes = [np.int8, np.int32, np.int64]
    for _ in range(1500):
        n_sites = int(rng.choice(n_sites_choices))
        n_samples = int(rng.choice(n_samples_choices))
        n_states = int(rng.choice(n_states_choices))
        states = rng.choice(12, size=n_states, replace=False)
        if rng.random() < 0.15:
            states = states - 1  # inclut -1 (donnée manquante)
        dtype = dtypes[int(rng.integers(len(dtypes)))]
        matrix = rng.choice(states, size=(n_sites, n_samples)).astype(dtype)
        compare_pairwise_hamming_distances(matrix)


def test_pairwise_hamming_distances_around_the_state_threshold():
    """Jusqu'à 8 états la fonction passe par un produit matriciel par état ;
    au-delà, par la comparaison directe de toutes les paires. Les deux côtés
    du seuil doivent donner le même résultat que l'oracle."""
    rng = np.random.default_rng(11)
    for max_state in (6, 7, 8, 9, 10):  # 8 états = états 0..7, 9 états = 0..8
        matrix = rng.integers(0, max_state + 1, size=(20, 40)).astype(np.int32)
        matrix[0, :] = np.arange(40) % (max_state + 1)  # force tous les états présents
        compare_pairwise_hamming_distances(matrix)


def test_pairwise_distances_per_locus(dna_context_te2):
    """Le helper renvoie une liste alignée sur les loci du groupe (et non un
    dictionnaire par nom de locus), avec pour chaque locus un vecteur par
    échantillon égal à `_pairwise_hamming_distances` de sa matrice."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography, context=dna_context_te2, seed=42
    )
    sample_names = ["pop1", "pop2"]
    names = [locus.name for locus in dna_context_te2.list_loci if locus.group == "G2"]
    tree_sequences = [mutated[name] for name in names]
    genotype_matrices_per_locus = [
        _genotype_matrix_by_sample(ts) for ts in tree_sequences
    ]

    distances_per_locus = _pairwise_distances_per_locus(
        genotype_matrices_per_locus, sample_names
    )

    assert isinstance(distances_per_locus, list)
    assert len(distances_per_locus) == len(tree_sequences) > 0
    for matrices, distances in zip(
        genotype_matrices_per_locus, distances_per_locus, strict=True
    ):
        assert list(distances) == sample_names
        for samp_name in sample_names:
            expected = _pairwise_hamming_distances(matrices[samp_name])
            assert np.array_equal(distances[samp_name], expected)
            assert distances[samp_name].dtype == np.int64


def test_pairwise_distances_per_locus_empty_group():
    """Un groupe sans locus donne une liste vide, sans exception."""
    assert _pairwise_distances_per_locus([], ["pop1", "pop2"]) == []


def test_pairwise_hamming_distances_int8_extremes():
    """Soustraire le minimum dans le type d'origine déborderait pour `int8` :
    les bornes -128 et 127 doivent donner les bonnes distances."""
    matrix = np.array([[127, -128, 127], [-128, 127, -128]], dtype=np.int8)
    assert _pairwise_hamming_distances(matrix).tolist() == [2, 0, 2]
    compare_pairwise_hamming_distances(matrix)


def test_mean_pairwise_differences_per_group(dna_context_te2):
    """Vérifie compute_MPD séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre VPD_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    # Valeurs du 02/10 (reproduction du bug DIYABC sur le dénominateur, voir
    # compute_MPD/compute_VPD et "G2 `MPD`/`VPD` residual" dans CLAUDE.md),
    # régénérées le 05/10 avec la dérivation de graines par _locus_seed.
    # ATTENTION : avec ces graines G2 et G3 ont chacun leurs 5 loci
    # polymorphes, donc ce fixture n'exerce PLUS le garde `num_sites == 0` --
    # ces valeurs sont celles du calcul sans loci écartés. Le garde est
    # vérifié par un test dédié, indépendant des graines :
    # test_mean_pairwise_differences_excludes_loci_without_mutation.
    mean_g2 = compute_MPD(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {
        "pop1": pytest.approx(2.4733333333),
        "pop2": pytest.approx(2.03),
    }

    mean_g3 = compute_MPD(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {
        "pop1": pytest.approx(0.5631578947368421),
        "pop2": pytest.approx(0.8400000000000001),
    }


def test_mean_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    sample_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = compute_MPD(tree_sequences_empty, sample_names)
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_pairwise_differences_excludes_loci_without_mutation():
    """Vérifie le garde `num_sites == 0` de compute_MPD/compute_VPD, qui
    reproduit un bug de DIYABC (voir "G2 `MPD`/`VPD` residual" dans
    CLAUDE.md) : un locus où la simulation n'a tiré aucune mutation sort
    du dénominateur de MPD et de VPD (sa séquence vide est le marqueur de
    donnée manquante `SEQMISSING`), alors que NSS, dont le `OK` vient de
    l'échantillon OBSERVÉ, le compte toujours avec 0 site.

    Indépendant des graines : le même arbre sert deux fois, muté et non
    muté, donc la différence ne tient qu'à l'absence de mutation. Les
    fixtures `te2` n'exercent ce garde que si une réalisation donnée laisse
    un locus monomorphe -- ce n'est plus le cas depuis le 05/10.
    """
    sample_names = ["pop1", "pop2"]
    ts_unmutated = msprime.sim_ancestry(
        samples=[msprime.SampleSet(10, "pop1"), msprime.SampleSet(10, "pop2")],
        demography=_two_population_demography(),
        sequence_length=100,
        random_seed=3,
    )
    ts_mutated = msprime.sim_mutations(ts_unmutated, rate=1e-3, random_seed=3)
    assert ts_unmutated.num_sites == 0
    assert ts_mutated.num_sites > 0

    # Un locus sans mutation s'ajoute : MPD et VPD ne bougent pas.
    mpd_alone = compute_MPD([ts_mutated], sample_names)
    assert compute_MPD([ts_mutated, ts_unmutated], sample_names) == mpd_alone
    vpd_alone = compute_VPD([ts_mutated], sample_names)
    assert compute_VPD([ts_mutated, ts_unmutated], sample_names) == vpd_alone

    # NSS, lui, divise par les 2 loci : la moitié du cas à un seul locus.
    nss_alone = compute_NSS([ts_mutated], sample_names)
    nss_both = compute_NSS([ts_mutated, ts_unmutated], sample_names)
    for name in sample_names:
        assert nss_both[name] == pytest.approx(nss_alone[name] / 2)

    # Et un jeu où AUCUN locus n'est muté donne 0.0, comme `res = 0.0` du C++.
    assert compute_MPD([ts_unmutated], sample_names) == {"pop1": 0.0, "pop2": 0.0}


def test_variance_pairwise_differences_per_group(dna_context_te2):
    """Vérifie compute_VPD séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre VPD_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    # Valeurs du 02/10 (reproduction du bug DIYABC sur le dénominateur, voir
    # compute_MPD/compute_VPD et "G2 `MPD`/`VPD` residual" dans CLAUDE.md),
    # régénérées le 05/10 avec la dérivation de graines par _locus_seed.
    # ATTENTION : avec ces graines G2 et G3 ont chacun leurs 5 loci
    # polymorphes, donc ce fixture n'exerce PLUS le garde `num_sites == 0` --
    # ces valeurs sont celles du calcul sans loci écartés. Le garde est
    # vérifié par un test dédié, indépendant des graines :
    # test_mean_pairwise_differences_excludes_loci_without_mutation.
    variance_g2 = compute_VPD(tree_sequences_for_group("G2"), sample_names)
    assert variance_g2 == {
        "pop1": pytest.approx(3.9747631744),
        "pop2": pytest.approx(3.179673151),
    }

    variance_g3 = compute_VPD(tree_sequences_for_group("G3"), sample_names)
    assert variance_g3 == {
        "pop1": pytest.approx(0.4454079643553328),
        "pop2": pytest.approx(0.9032804232804231),
    }


def test_variance_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    sample_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    variance_empty = compute_VPD(tree_sequences_empty, sample_names)
    assert variance_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_tajima_d_per_group(dna_context_te2):
    """Vérifie compute_DTA séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre DTA_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    dta_g2 = compute_DTA(tree_sequences_for_group("G2"), sample_names)
    assert dta_g2 == {
        "pop1": pytest.approx(0.1808813166),
        "pop2": pytest.approx(-0.162658095),
    }

    dta_g3 = compute_DTA(tree_sequences_for_group("G3"), sample_names)
    assert dta_g3 == {
        "pop1": pytest.approx(-0.4111990563677869),
        "pop2": pytest.approx(-0.414526208396265),
    }


def test_mean_tajima_d_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    sample_names = ["pop1", "pop2"]
    assert compute_DTA([], sample_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_private_segregating_sites_per_group(dna_context_te2):
    """Vérifie compute_PSS séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre PSS_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    pss_g2 = compute_PSS(tree_sequences_for_group("G2"), sample_names)
    assert pss_g2 == {"pop1": pytest.approx(2.2), "pop2": pytest.approx(1.6)}

    pss_g3 = compute_PSS(tree_sequences_for_group("G3"), sample_names)
    assert pss_g3 == {"pop1": pytest.approx(1.4), "pop2": pytest.approx(2.0)}


def test_mean_private_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_PSS([], sample_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_minor_allele_count_per_group(dna_context_te2):
    """Vérifie compute_MNS séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre MNS_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mns_g2 = compute_MNS(tree_sequences_for_group("G2"), sample_names)
    assert mns_g2 == {
        "pop1": pytest.approx(6.4578571429),
        "pop2": pytest.approx(5.6341658342),
    }

    mns_g3 = compute_MNS(tree_sequences_for_group("G3"), sample_names)
    assert mns_g3 == {
        "pop1": pytest.approx(1.8),
        "pop2": pytest.approx(2.7833333333333337),
    }


def test_mean_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_MNS([], sample_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_variance_minor_allele_count_per_group(dna_context_te2):
    """Vérifie compute_VNS séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre VNS_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    vns_g2 = compute_VNS(tree_sequences_for_group("G2"), sample_names)
    assert vns_g2 == {
        "pop1": pytest.approx(30.5631437075),
        "pop2": pytest.approx(23.57072631),
    }

    vns_g3 = compute_VNS(tree_sequences_for_group("G3"), sample_names)
    assert vns_g3 == {
        "pop1": pytest.approx(5.9399999999999995),
        "pop2": pytest.approx(3.381944444444444),
    }


def test_variance_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_VNS([], sample_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group_pairwize(dna_context_te2):
    """Vérifie compute_NH2 séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre NH_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NH2(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"1.2": pytest.approx(9.6)}

    mean_g3 = compute_NH2(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"1.2": pytest.approx(5.0)}


def test_mean_distinct_haplotypes_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_NH2([], sample_names) == {
        "1.2": 0.0,
    }


def test_mean_segregating_sites_per_group_pairwize(dna_context_te2):
    """Vérifie compute_NS2 séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre NSS_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NS2(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"1.2": pytest.approx(11.2)}

    mean_g3 = compute_NS2(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"1.2": pytest.approx(4.6)}


def test_mean_segregating_sites_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_NS2([], sample_names) == {
        "1.2": 0.0,
    }


def test_mean_pairwise_differences_per_group_pairwize(dna_context_te2):
    """Vérifie compute_MP2 séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre VPD_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_MP2(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"1.2": pytest.approx(2.2516666667)}

    mean_g3 = compute_MP2(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"1.2": pytest.approx(0.7015789473684211)}


def test_mean_pairwise_differences_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_MP2([], sample_names) == {
        "1.2": 0.0,
    }


def test_mean_pairwise_differences_between_per_group_pairwize(dna_context_te2):
    """Vérifie compute_MPB séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre VPD_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_MPB(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"1.2": pytest.approx(2.30225)}

    mean_g3 = compute_MPB(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"1.2": pytest.approx(0.742)}


def test_mean_pairwise_differences_between_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population du
    résultat ni lever d'exception -- chaque population attendue garde une
    valeur (0.0)."""
    sample_names = ["pop1", "pop2"]
    assert compute_MPB([], sample_names) == {
        "1.2": 0.0,
    }


def test_mean_hst_per_group_pairwize(dna_context_te2):
    """Vérifie compute_HST séparément sur G2 (<A>) et G3 (<M>) de
    toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque chaque
    `group Gx` du header calcule son propre HST_i à partir de ses seuls
    loci."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    loci_description = dna_context_te2.list_loci
    sample_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_HST(tree_sequences_for_group("G2"), sample_names)
    assert mean_g2 == {"1.2": pytest.approx(0.02197126)}

    mean_g3 = compute_HST(tree_sequences_for_group("G3"), sample_names)
    assert mean_g3 == {"1.2": pytest.approx(0.05447581217193929)}


def test_compute_all_statistics_dna(dna_context_te2):
    """Vérifie compute_all_statistics_dna sur toy_example2_ms_dna."""
    demography, _ = build_random_demography_for_scenario_index(
        dna_context_te2.header_text, scenario_index=1, seed=42
    )
    mutated, _ = dna_mutation_simulation_per_locus(
        demography=demography,
        context=dna_context_te2,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_all_statistics_dna(
        dna_context_te2.header_text, mutated, sample_names
    )

    # Vérifier que les résultats contiennent les clés attendues
    expected_keys = {
        "NHA_2_1",
        "NHA_2_2",
        "NSS_2_1",
        "NSS_2_2",
        "MPD_2_1",
        "MPD_2_2",
        "VPD_2_1",
    }
    assert expected_keys.issubset(results.keys())

    assert len(results) == 42
    assert pytest.approx(results["NSS_2_1"]) == 9.6
    assert pytest.approx(results["HST_2_1.2"]) == 0.02197126
    assert pytest.approx(results["NH2_3_1.2"]) == 5.0


# -----------------------------------------------------------------------
# Pour les microsat
# -----------------------------------------------------------------------


def _mutated_microsat_ts(seed=3):
    """TreeSequence mutée minimale pour les helpers MICROSAT.

    `sequence_length=1` : un locus microsat est UN site à grand
    alphabet, jamais plusieurs sites indépendants.
    """
    demography = msprime.Demography()
    for name in ("pop1", "pop2", "anc"):
        demography.add_population(name=name, initial_size=1000)
    demography.add_population_split(time=500, derived=["pop1", "pop2"], ancestral="anc")
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(10, "pop1"), msprime.SampleSet(10, "pop2")],
        demography=demography,
        sequence_length=1,
        random_seed=seed,
    )
    return msprime.sim_mutations(
        ts,
        rate=1e-3,
        random_seed=seed,
        model=build_transition_matrix_microsat(
            kmin=180, kmax=220, motif_size=2, Pgeom=0.2
        ),
    )


def _reference_length_by_sample(
    tree_sequence: tskit.TreeSequence,
    *,
    layout: list[tuple[str, np.ndarray]] | None = None,
) -> dict[str, list[tuple[int, int]]]:
    """Construit un dict {nom_echantillon: [(longueur, nb_sequence), ...]}.

    Args:
        tree_sequence: Un TreeSequence muté du groupe (un locus [M]).
        layout: [(nom_echantillon, np.ndarray[indices_d'individus]), ...] : si fourni, remplace le découpage par échantillon à utilisé par le chemin sériel, où un échantillon n'est plus sa propre échantillon

    Returns:
        Un dict {nom_echantillon: [(longueur, compte), ...]}.
    """
    length_by_sample = {}
    if layout is None:
        layout = compute_population_layout(
            tree_sequence
        )  # liste(tuple(pop,array(indice)))
    # si locus est monomorphe
    if tree_sequence.num_sites == 0:
        for samp_name, sample_ids in layout:
            length_by_sample[samp_name] = [(0, len(sample_ids))]
        return length_by_sample
    else:
        variant = next(tree_sequence.variants())
        tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
        for samp_name, sample_ids in layout:
            length_by_sample[samp_name] = []
            list_taille_sample_ids = tailles[sample_ids].tolist()
            for taille in variant.alleles:
                taille_int = int(taille)
                length_by_sample[samp_name].append(
                    (taille_int, list_taille_sample_ids.count(taille_int))
                )

        return length_by_sample


def test_length_by_sample(microsat_context_te2_xy):
    """Vérifie _length_by_sample sur un locus microsat (Locus_M_A_1_).

    Convertit les codes de génotype en tailles réelles (pb) via
    variant.alleles, puis compte les copies de gène par taille et par
    population -- couvre les deux populations pour éviter de repasser
    silencieusement au bug de découpage corrigé pendant l'écriture de
    cette fonction (indexation par sample_ids plutôt que par un
    intervalle [first_index, last_index] supposé contigu).
    """
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    ts = mutated["Locus_M_A_1_"]
    results = _length_by_sample(ts)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == [
        (201, 5),
        (199, 15),
        (203, 3),
        (195, 13),
        (197, 3),
    ]
    assert results["pop2"] == [
        (201, 4),
        (199, 14),
        (203, 1),
        (195, 16),
        (197, 5),
    ]


def test_length_by_sample_with_same_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    result_with_layout = _length_by_sample(ts, layout=default_layout)
    result_without_layout = _length_by_sample(ts)

    assert list(result_with_layout.keys()) == list(result_without_layout.keys())
    for sample_name in result_with_layout:
        assert (
            result_with_layout[sample_name][0] == result_without_layout[sample_name][0]
        )
        assert (
            result_with_layout[sample_name][1] == result_without_layout[sample_name][1]
        )


def test_length_by_sample_with_different_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    # Create a different layout
    split_layout = []
    for sample_name, node_ids in default_layout:
        half = len(node_ids) // 2
        split_layout.append((f"{sample_name}_a", node_ids[:half]))
        split_layout.append((f"{sample_name}_b", node_ids[half:]))

    result_with_layout = _length_by_sample(ts, layout=split_layout)
    result_without_layout = _length_by_sample(ts)

    assert list(result_with_layout) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]
    assert list(result_without_layout) == ["pop1", "pop2"]

    for sample_name, node_ids in split_layout:
        assert sum(nb_seq for _, nb_seq in result_with_layout[sample_name]) == len(
            node_ids
        )

    # Recollement : couper une population coupe les COMPTES, pas la liste --
    # chaque groupe rend une entrée par allèle distinct, les mêmes allèles
    # pour tous (la boucle interne itère variant.alleles). D'où une somme
    # terme à terme, et non une concaténation comme pour les helpers par
    # individu.
    for sample_name in ("pop1", "pop2"):
        half_a = result_with_layout[f"{sample_name}_a"]
        half_b = result_with_layout[f"{sample_name}_b"]

        # L'alignement des listes est une hypothèse : on la teste au lieu
        # de s'y fier, sinon le zip apparierait silencieusement deux
        # allèles différents.
        assert [taille for taille, _ in half_a] == [taille for taille, _ in half_b]

        rejoined = [
            (taille, count_a + count_b)
            for (taille, count_a), (_, count_b) in zip(half_a, half_b, strict=True)
        ]
        assert rejoined == result_without_layout[sample_name]


def test_count_alleles_per_sample():
    length_per_population = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    results = count_alleles_per_sample(length_per_population)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == 3
    assert results["pop2"] == 4


def test_compute_NAL(microsat_context_te2_xy):
    """Vérifie compute_NAL sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_NAL(mutated.values(), sample_names)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 85 / 10
    assert pytest.approx(results["pop2"]) == 87 / 9


def test_total_genes_copies_per_sample(microsat_context_te2_xy):
    """Vérifie _total_genes_copies_per_sample sur
    toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    length_by_sample = _length_by_sample(mutated["Locus_M_A_1_"])
    results = total_genes_copies_per_sample(length_by_sample)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == 39
    assert results["pop2"] == 40


def test_compute_HET(microsat_context_te2_xy):
    """Vérifie compute_HET sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    tree_sequences1 = [mutated["Locus_M_A_1_"]]

    sample_names = ["pop1", "pop2"]

    results = compute_HET(tree_sequences1, sample_names)

    assert results.keys() == {"pop1", "pop2"}
    assert (
        pytest.approx(results["pop1"])
        == (
            1
            - (5 / 39) ** 2
            - (15 / 39) ** 2
            - (3 / 39) ** 2
            - (13 / 39) ** 2
            - (3 / 39) ** 2
        )
        * 39
        / 38
    )
    assert (
        pytest.approx(results["pop2"])
        == (
            1
            - (4 / 40) ** 2
            - (14 / 40) ** 2
            - (1 / 40) ** 2
            - (16 / 40) ** 2
            - (5 / 40) ** 2
        )
        * 40
        / 39
    )

    tree_sequences2 = list(mutated.values())
    results = compute_HET(tree_sequences2, sample_names)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 8.1632253711 / 10
    assert pytest.approx(results["pop2"]) == 7.0675438596 / 9


# tests pour la stat VAR (Variance of Repeat Number) pour microsatellites


def test_compute_VAR_constants(microsat_context_te2_xy):
    """Vérifie que les constantes de compute_VAR sont correctes."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    length_by_sample = _length_by_sample(mutated["Locus_M_A_1_"])
    raw_sizes, raw_square_sizes, total_counts = _compute_VAR_constants(length_by_sample)

    assert raw_sizes.keys() == {"pop1", "pop2"}
    assert raw_square_sizes.keys() == {"pop1", "pop2"}
    assert total_counts.keys() == {"pop1", "pop2"}

    assert raw_sizes["pop1"] == 201 * 5 + 199 * 15 + 203 * 3 + 195 * 13 + 197 * 3
    assert (
        raw_square_sizes["pop1"]
        == 201**2 * 5 + 199**2 * 15 + 203**2 * 3 + 195**2 * 13 + 197**2 * 3
    )
    assert total_counts["pop1"] == 39

    assert raw_sizes["pop2"] == 201 * 4 + 199 * 14 + 203 * 1 + 195 * 16 + 197 * 5
    assert (
        raw_square_sizes["pop2"]
        == 201**2 * 4 + 199**2 * 14 + 203**2 * 1 + 195**2 * 16 + 197**2 * 5
    )
    assert total_counts["pop2"] == 40


def test_compute_VAR_for_one_sample(microsat_context_te2_xy):
    """Vérifie compute_VAR sur toy_example2_ms_dna_xy pour une seule
    population."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    s, v, n = _compute_VAR_constants(_length_by_sample(mutated["Locus_M_A_1_"]))
    result = _compute_VAR_for_one_sample(s["pop1"], v["pop1"], n["pop1"], 2)

    s1, v1, n1 = (
        201 * 5 + 199 * 15 + 203 * 3 + 195 * 13 + 197 * 3,
        201**2 * 5 + 199**2 * 15 + 203**2 * 3 + 195**2 * 13 + 197**2 * 3,
        39,
    )
    assert pytest.approx(result) == (v1 - s1**2 / n1) / (n1 - 1) / 2**2


def test_compute_VAR(microsat_context_te2_xy):
    """Vérifie compute_VAR sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]
    list_loci = [
        locus for locus in microsat_context_te2_xy.list_loci if locus.ms_or_seq == "M"
    ]
    list_motif_sizes = [locus.motif_size for locus in list_loci]
    results = compute_VAR(mutated.values(), sample_names, list_motif_sizes)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 122.3338141026 / 10
    assert pytest.approx(results["pop2"]) == 131.2087044534 / 9


# tests pour la stat MGW (Mean Gene Width) pour microsatellites


def test_compute_MGW_by_locus(microsat_context_te2_xy):
    """Vérifie _compute_MGW_by_locus sur toy_example2_ms_dna_xy.

    Locus_M_A_1_ : pop1 et pop2 ont chacune leurs 5 allèles présents
    (195/197/199/201/203, tous à compte non nul), donc min=195 et max=203
    pour les deux.
    """
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    result = _compute_MGW_by_locus(_length_by_sample(mutated["Locus_M_A_1_"]), 2)

    assert result.keys() == {"pop1", "pop2"}
    assert result["pop1"] == (5, 1 + (203 - 195) / 2)
    assert result["pop2"] == (5, 1 + (203 - 195) / 2)


def test_compute_MGW(microsat_context_te2_xy):
    """Vérifie compute_MGW sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )
    list_loci = [
        locus for locus in microsat_context_te2_xy.list_loci if locus.ms_or_seq == "M"
    ]
    list_motif_sizes = [locus.motif_size for locus in list_loci]
    sample_names = ["pop1", "pop2"]
    results = compute_MGW(mutated.values(), sample_names, list_motif_sizes)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 0.685483871
    assert pytest.approx(results["pop2"]) == 0.7802690583


## tests pour la stat N2P (Number of Private Alleles) pour microsatellites


def test_compute_N2P_for_one_pair():
    """Vérifie que la fonction _compute_N2P_for_one_pair fonctionne
    correctement."""
    alleles_pop1 = [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)]
    alleles_pop2 = [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)]

    combined_alleles = _compute_N2P_for_one_pair(alleles_pop1, alleles_pop2)

    assert combined_alleles == 5


def test_compute_N2P_for_one_locus():
    """Vérifie que la fonction _compute_N2P_for_one_locus fonctionne
    correctement."""
    length_by_sample = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    combined_alleles_count = _compute_N2P_for_one_locus(
        length_by_sample, ["pop1", "pop2"]
    )

    assert combined_alleles_count == {"1.2": 5.0}


def test_compute_N2P(microsat_context_te2_xy):
    """Vérifie compute_N2P sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_N2P(mutated.values(), sample_names)

    assert results == {"1.2": 11.4}


# tests relatifs à la stat H2P


def _reference_pool_allele_counts_for_two_samples(
    pop_1: list[tuple[int, int]], pop_2: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Fusionne les comptes bruts de deux échantillons, allèle par allèle.

    Ne calcule ni H2P ni une fréquence -- juste n_i·freq_i + n_j·freq_j
    (= compte_i + compte_j, la division par n_i/n_j s'annulant avec la
    multiplication), pour que _compute_H2P_for_one_pair puisse
    réutiliser _compute_HET_for_one_sample dessus (qui fait
    lui-même la division par le total).

    Args:
        pop_1: Liste de tuples (longueur, nb_sequence) pour l'échantillon 1.
        pop_2: Liste de tuples (longueur, nb_sequence) pour l'échantillon 2.

    Returns:
        Liste de tuples (longueur, compte poolé) pour chaque allèle.
    """
    pooled_counts = []
    n1 = sum(count for _, count in pop_1)
    n2 = sum(count for _, count in pop_2)
    n = n1 + n2
    if n <= 1:
        return [(0, 0.0)]
    else:
        for length, count_1 in pop_1:
            count_2 = next((count for len, count in pop_2 if len == length), 0)
            total = count_1 + count_2
            pooled_counts.append((length, total))

    return pooled_counts


def test_pool_allele_counts_vs_reference():
    """Vérifie que la fonction _pool_allele_counts_for_two_samples
    correspond à la référence."""
    alleles_pop1 = [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)]
    alleles_pop2 = [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)]

    pooled_counts = _pool_allele_counts_for_two_samples(alleles_pop1, alleles_pop2)
    reference_pooled_counts = _reference_pool_allele_counts_for_two_samples(
        alleles_pop1, alleles_pop2
    )

    assert pooled_counts == reference_pooled_counts


def test_pool_allele_counts_for_two_samples():
    """Vérifie que la fonction _pool_allele_counts_for_two_samples
    fonctionne correctement."""
    alleles_pop1 = [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)]
    alleles_pop2 = [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)]

    pooled_counts = _pool_allele_counts_for_two_samples(alleles_pop1, alleles_pop2)

    expected_pooled_counts = [
        (201, 0.0),
        (203, 49.0),
        (197, 1.0),
        (193, 27.0),
        (199, 1.0),
        (189, 1.0),
    ]
    assert pooled_counts == expected_pooled_counts


def test_compute_H2P_for_one_pair():
    """Vérifie que la fonction _compute_H2P_for_one_pair fonctionne
    correctement."""
    _length_by_sample = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    h2p_value = _compute_H2P_for_one_pair(_length_by_sample, "pop1", "pop2")

    expected_h2p_value = (
        79
        / 78
        * (
            1
            - (49 / 79) ** 2
            - (1 / 79) ** 2
            - (27 / 79) ** 2
            - (1 / 79) ** 2
            - (1 / 79) ** 2
        )
    )
    assert pytest.approx(h2p_value) == expected_h2p_value

    with pytest.raises(KeyError, match="L'un des échantillons n'est pas"):
        _compute_H2P_for_one_pair(_length_by_sample, "pop1", "pop3")


def test_compute_H2P(microsat_context_te2_xy):
    """Vérifie compute_H2P sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_H2P(mutated.values(), sample_names)

    assert pytest.approx(results["1.2"]) == 0.8068366815


# tests relatifs à la stat V2P


def test_compute_V2P_constants():
    """Vérifie que les constantes de compute_V2P sont correctes."""
    length_by_sample = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    s, v, n = _compute_VAR_constants(length_by_sample)

    s1, v1, n1 = _compute_V2P_constants("pop1", "pop2", s, v, n)

    assert (
        s1
        == 203 * 31
        + 197 * 1
        + 193 * 7
        + 199 * 0
        + 189 * 0
        + 203 * 18
        + 197 * 0
        + 193 * 20
        + 199 * 1
        + 189 * 1
    )
    assert (
        v1
        == 203**2 * 31
        + 197**2 * 1
        + 193**2 * 7
        + 199**2 * 0
        + 189**2 * 0
        + 203**2 * 18
        + 197**2 * 0
        + 193**2 * 20
        + 199**2 * 1
        + 189**2 * 1
    )
    assert n1 == 39 + 40


def test_compute_V2P(microsat_context_te2_xy):
    """Vérifie compute_V2P sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    list_loci = microsat_context_te2_xy.list_loci
    list_motif_sizes = [
        locus.motif_size for locus in list_loci if locus.ms_or_seq == "M"
    ]
    sample_names = ["pop1", "pop2"]

    results = compute_V2P(mutated.values(), sample_names, list_motif_sizes)

    assert pytest.approx(results["1.2"]) == 13.1654476236


# tests relatifs à la stat DAS


def test_compute_identical_pair_for_one_pair():
    """Vérifie que la fonction _compute_identical_pair_for_one_pair fonctionne
    correctement."""
    _length_by_sample = {
        "samp1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "samp2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }
    samp_a, samp_b = "samp1", "samp2"

    identical_count, total_count = _compute_identical_pair_for_one_pair(
        _length_by_sample, samp_a, samp_b
    )

    assert total_count == 39 * 40
    assert identical_count == 31 * 18 + 7 * 20


def test_compute_DAS(microsat_context_te2_xy):
    """Vérifie compute_DAS sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_DAS(mutated.values(), sample_names)

    assert pytest.approx(results["1.2"]) == 0.1870820669


# tests relatifs à la stat DM2


def test_compute_DM2_for_one_locus():
    """Vérifie que la fonction _compute_DM2_for_one_locus fonctionne
    correctement."""
    # test sur un locus où les deux populations sont présentes
    length_by_sample = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    raw_sizes, _, total_counts = _compute_VAR_constants(length_by_sample)

    contribution, new_moy, was_valid = _compute_DM2_for_one_locus(
        "pop1", "pop2", 2, length_by_sample, raw_sizes, total_counts, None
    )

    moy_1 = raw_sizes["pop1"] / total_counts["pop1"]
    moy_2 = raw_sizes["pop2"] / total_counts["pop2"]
    assert contribution == ((moy_2 - moy_1) / 2) ** 2
    assert pytest.approx(new_moy) == (moy_1, moy_2)
    assert was_valid

    # test sur un locus où une des populations est absente.
    length_by_sample_2 = {
        "pop1": [
            (208, 1),
            (212, 4),
            (210, 0),
            (214, 2),
            (206, 0),
            (202, 0),
            (204, 0),
            (198, 6),
            (200, 3),
            (218, 0),
            (216, 2),
            (196, 2),
        ]
    }

    raw_sizes_2, _, total_counts_2 = _compute_VAR_constants(length_by_sample_2)

    contribution_2, new_moy_2, was_valid_2 = _compute_DM2_for_one_locus(
        "pop1",
        "pop2",
        2,
        length_by_sample_2,
        raw_sizes_2,
        total_counts_2,
        (100.0, 105.0),
    )

    assert contribution_2 == ((105.0 - 100.0) / 2) ** 2
    assert new_moy_2 == (100.0, 105.0)
    assert not was_valid_2


def test_compute_DM2(microsat_context_te2_xy):
    """Vérifie compute_DM2 sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    list_loci = microsat_context_te2_xy.list_loci
    list_motif_sizes = [
        locus.motif_size for locus in list_loci if locus.ms_or_seq == "M"
    ]
    sample_names = ["pop1", "pop2"]

    results = compute_DM2(mutated.values(), sample_names, list_motif_sizes)

    assert pytest.approx(results["1.2"]) == 2.4439210429


# test relatifs à la stat FST
def test_length_by_sample_and_individuals(microsat_context_te2_xy):
    """Vérifie _length_by_sample_and_individuals sur un locus microsat
    (Locus_M_A_1_).

    Convertit les codes de génotype en tailles réelles (pb) via
    variant.alleles, puis compte les copies de gène par taille et par
    population -- couvre les deux populations pour éviter de repasser
    silencieusement au bug de découpage corrigé pendant l'écriture de
    cette fonction (indexation par sample_ids plutôt que par un
    intervalle [first_index, last_index] supposé contigu).
    """
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    ts = mutated["Locus_M_A_1_"]
    results = _length_by_sample_and_individuals(ts)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == [
        (195, 195),
        (201, 195),
        (199, 201),
        (195, 201),
        (195, 195),
        (199, 199),
        (199, 199),
        (203, 195),
        (199, 199),
        (203, 195),
        (195, 195),
        (203, 197),
        (199, 199),
        (201, 201),
        (195, 199),
        (199, 197),
        (197, 199),
        (195, 195),
        (199, 199),
        (199, 199),
    ]


def test_length_by_pop_and_individuals_with_same_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    result_with_layout = _length_by_sample_and_individuals(ts, layout=default_layout)
    result_without_layout = _length_by_sample_and_individuals(ts)

    assert list(result_with_layout.keys()) == list(result_without_layout.keys())
    for sample_name in result_with_layout:
        assert result_with_layout[sample_name] == result_without_layout[sample_name]


def test_length_by_pop_and_individuals_with_different_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    # calcul de la ploïdie valable ici car notre fixture n'a que des individus diploïdes
    ploidy = ts.num_samples // ts.num_individuals
    # Create a different layout
    split_layout = []
    for sample_name, node_ids in default_layout:
        half = len(node_ids) // 2
        split_layout.append((f"{sample_name}_a", node_ids[:half]))
        split_layout.append((f"{sample_name}_b", node_ids[half:]))

    result_with_layout = _length_by_sample_and_individuals(ts, layout=split_layout)
    result_without_layout = _length_by_sample_and_individuals(ts)

    assert list(result_with_layout) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]
    assert list(result_without_layout) == ["pop1", "pop2"]

    for sample_name, node_ids in split_layout:
        assert len(result_with_layout[sample_name]) == len(node_ids) // ploidy

    for sample_name in ("pop1", "pop2"):
        half_a = result_with_layout[f"{sample_name}_a"]
        half_b = result_with_layout[f"{sample_name}_b"]

        rejoined = half_a + half_b
        assert rejoined == result_without_layout[sample_name]


def reference_compute_ni_nA_AA_for_one_sample(
    pairs: list[tuple[int, int]], al: int
) -> tuple[int, int, int]:
    """Calcule ni, nA et AA pour un échantillon donné à partir des paires
    d'allèles.

    Args:
        pairs: Liste de tuples (longueur_1, longueur_2) pour chaque individu.
        al: Longueur de l'allèle considéré.

    Returns:
        Un tuple (ni, nA, AA) où
            - ni est le nombre d'individus,
            - nA est la somme par individus du nombre d'éléments de sa paire égaux
            - AA est le nb d'individus dont les deux éléments d ela paire valent al
    """
    ni = len(pairs)
    nA = sum(int(p[0] == al) + int(p[1] == al) for p in pairs)
    AA = sum(1 for p in pairs if p[0] == al and p[1] == al)
    return ni, nA, AA


def test_allele_counts_by_sample():
    """Vérifie que la fonction _allele_counts_by_sample fonctionne correctement."""
    # test sur un échantillon avec 4 individus, dont 3 homozygotes pour l'allèle 203 et 1 hétérozygote pour 203 et 193
    sample_alleles = {
        "sample1": [
            np.int64([203, 203]),
            np.int64([203, 203]),
            np.int64([203, 193]),
            np.int64([203, 203]),
        ]
    }

    allele_counts = _allele_counts_by_sample(sample_alleles)

    assert allele_counts == ({"sample1": 4}, {"sample1": {203: (7, 3), 193: (1, 0)}})

    # test sur un échantillon sans homozygotes
    sample_alleles = {
        "sample2": [
            np.int64([203, 193]),
            np.int64([197, 203]),
            np.int64([203, 193]),
        ]
    }
    allele_counts = _allele_counts_by_sample(sample_alleles)
    assert allele_counts == (
        {"sample2": 3},
        {"sample2": {203: (3, 0), 193: (2, 0), 197: (1, 0)}},
    )


def test_allele_counts_by_sample_vs_oracle():
    """Compare _allele_counts_by_sample à l'oracle sur des entrées aléatoires.

    Peu d'allèles (beaucoup d'homozygotes), trois échantillons de tailles
    différentes, et des allèles absents de certains échantillons (dont un
    absent de tous).
    """
    rng = np.random.default_rng(12345)
    for _ in range(50):
        alleles = [190 + 2 * k for k in range(int(rng.integers(2, 6)))]
        length_by_sample = {
            f"sample{s}": [
                (int(rng.choice(alleles)), int(rng.choice(alleles)))
                for _ in range(int(rng.integers(1, 12)))
            ]
            for s in range(1, 4)
        }
        ni_by_sample, nA_AA_by_sample = _allele_counts_by_sample(length_by_sample)

        all_alleles = {
            al for pairs in length_by_sample.values() for p in pairs for al in p
        }
        all_alleles.add(58)  # allèle présent nulle part
        for name, pairs in length_by_sample.items():
            for al in all_alleles:
                nA, AA = nA_AA_by_sample[name].get(al, (0, 0))
                assert (ni_by_sample[name], nA, AA) == (
                    reference_compute_ni_nA_AA_for_one_sample(pairs, al)
                )


def test_compute_FST_constants_for_two_samples_combined():
    # pairs_1 = [(203, 203), (203, 203), (203, 193), (203, 203)]
    # pairs_2 = [(203, 203), (203, 193), (197, 203), (203, 193)]

    s2G, s2I, s2p = _compute_FST_constants_for_two_samples_combined(4, 7, 3, 4, 5, 1)

    assert s2G == 0.25
    assert s2I == -0.0625
    assert s2p == 0.015625

    # second test avec une allele non_presente
    s2G, s2I, s2p = _compute_FST_constants_for_two_samples_combined(0, 0, 0, 0, 0, 0)
    assert s2G == 0
    assert s2I == 0
    assert s2p == 0


def test_compute_FST_constants_on_all_alleles_for_two_samples():
    """Vérifie que la fonction
    _compute_constants_on_all_alleles_for_two_populations fonctionne
    correctement."""
    _length_by_sample = {
        "pop1": [(203, 203), (203, 203), (203, 193), (203, 203)],
        "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
    }
    ni_by_sample, nA_AA_by_sample = _allele_counts_by_sample(_length_by_sample)

    s1l, s2l, s3l = _compute_FST_constants_on_all_alleles_for_two_samples(
        ni_by_sample, nA_AA_by_sample, "pop1", "pop2"
    )

    assert pytest.approx(s1l) == 0.005208333333333332
    assert pytest.approx(s2l) == -0.078125
    assert pytest.approx(s3l) == 0.421875

    # nouveau test avec des pop non présente
    s1l, s2l, s3l = _compute_FST_constants_on_all_alleles_for_two_samples(
        ni_by_sample, nA_AA_by_sample, "pop3", "pop4"
    )

    assert s1l == 0.0
    assert s2l == 0.0
    assert s3l == 0.0


def test_compute_FST(microsat_context_te2_xy):
    """Vérifie compute_FST sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_FST(mutated.values(), sample_names)

    assert pytest.approx(results["1.2"]) == 0.0119165839


def test_compute_FST_vs_scikit_allel(microsat_context_te1_modified):
    allel = pytest.importorskip("allel")
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te1_modified.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te1_modified,
        seed=42,
    )
    list_loci = microsat_context_te1_modified.list_loci
    g1 = [mutated[locus.name] for locus in list_loci if locus.group == "G1"]
    ours = compute_FST(g1, ["pop1", "pop2", "pop3", "pop4"])["1.2"]
    A = B = C = 0.0
    for ts in g1:
        if ts.num_sites == 0:
            continue
        v = next(ts.variants())
        codes = v.genotypes
        inds = list(ts.individuals())
        G = np.array([[[codes[n] for n in ind.nodes] for ind in inds]])
        pop_of_ind = [ts.node(ind.nodes[0]).population for ind in inds]
        sub = [
            [i for i, p in enumerate(pop_of_ind) if p == 0],
            [i for i, p in enumerate(pop_of_ind) if p == 1],
        ]
        a, b, c = allel.weir_cockerham_fst(
            allel.GenotypeArray(G), subpops=sub, max_allele=int(codes.max())
        )
        A += a.sum()
        B += b.sum()
        C += c.sum()
    assert ours == pytest.approx(A / (A + B + C), abs=1e-9)


# tests relatif à la stat LIK


def test_genotypes_by_sample_and_individuals(microsat_context_te2_xy):
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    ts = mutated["Locus_M_A_1_"]
    results = _genotypes_by_sample_and_individuals(ts)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == [
        (195, 195),
        (201, 195),
        (199, 201),
        (195, 201),
        (195, 195),
        (199, 199),
        (199, 199),
        (203, 195),
        (199, 199),
        (203, 195),
        (195, 195),
        (203, 197),
        (199, 199),
        (201, 201),
        (195, 199),
        (199, 197),
        (197, 199),
        (195, 195),
        (199, 199),
        (199,),
    ]


def test_genotypes_by_pop_and_individuals_with_same_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    result_with_layout = _genotypes_by_sample_and_individuals(ts, layout=default_layout)
    result_without_layout = _genotypes_by_sample_and_individuals(ts)

    assert list(result_with_layout.keys()) == list(result_without_layout.keys())
    for sample_name in result_with_layout:
        assert result_with_layout[sample_name] == result_without_layout[sample_name]


def test_genotypes_by_pop_and_individuals_with_different_layout():
    ts = _mutated_microsat_ts()
    default_layout = compute_population_layout(ts)

    # calcul de la ploïdie valable ici car notre fixture n'a que des individus diploïdes
    ploidy = ts.num_samples // ts.num_individuals
    # Create a different layout
    split_layout = []
    for sample_name, node_ids in default_layout:
        half = len(node_ids) // 2
        split_layout.append((f"{sample_name}_a", node_ids[:half]))
        split_layout.append((f"{sample_name}_b", node_ids[half:]))

    result_with_layout = _genotypes_by_sample_and_individuals(ts, layout=split_layout)
    result_without_layout = _genotypes_by_sample_and_individuals(ts)

    assert list(result_with_layout) == ["pop1_a", "pop1_b", "pop2_a", "pop2_b"]
    assert list(result_without_layout) == ["pop1", "pop2"]

    for sample_name, node_ids in split_layout:
        assert len(result_with_layout[sample_name]) == len(node_ids) // ploidy

    for sample_name in ("pop1", "pop2"):
        half_a = result_with_layout[f"{sample_name}_a"]
        half_b = result_with_layout[f"{sample_name}_b"]

        rejoined = half_a + half_b
        assert rejoined == result_without_layout[sample_name]


# ---------------------------------------------------------------------------
# Lecture des individus par échantillon via `ts.individuals_nodes`.
#
# `_individuals_by_sample`, `_length_by_sample_and_individuals` et
# `_genotypes_by_sample_and_individuals` ont remplacé une boucle
# `for ind in ts.individuals()` (un objet tskit par individu, puis une
# recherche `nodes[0] in nodes_du_layout`). Cette boucle est conservée ici
# comme ORACLE : les fonctions vectorisées doivent donner exactement le même
# résultat, y compris sur ce que les fixtures diploïdes ci-dessus ne voient pas
# (individus haploïdes, mélange de ploïdies, locus monomorphe, layout sériel).
# ---------------------------------------------------------------------------

# {nom du scénario: ((nb_individus, ploïdie), ...)} -- un tuple par échantillon.
_PLOIDY_SCENARIOS = {
    "diploid": ((4, 2), (3, 2)),
    "mixed": ((3, 2), (2, 1)),
    "haploid": ((5, 1), (4, 1)),
    "three_samples": ((2, 2), (3, 1), (2, 2)),
}


def _ts_with_samples(sample_ploidies, mutation_rate, seed=1):
    """Construit un TreeSequence à échantillons de ploïdies variables.

    Chaque échantillon est un `SampleSet` daté de `10 * i` générations, dans UNE
    seule population : c'est la configuration sérielle, où deux échantillons
    partagent une population msprime. Les mutations suivent un modèle de
    microsatellite (allèles = tailles entières).

    Args:
        sample_ploidies: ((nb_individus, ploïdie), ...), un tuple par échantillon.
        mutation_rate: Taux de mutation. 0.0 donne un locus monomorphe.
        seed: Graine de la généalogie et des mutations.

    Returns:
        (ts, counts) où counts est {"s0": nb_individus, "s1": ...}, l'argument
        attendu par `compute_sample_layout`.
    """
    demography = msprime.Demography.isolated_model([200])
    sample_sets = [
        msprime.SampleSet(n, ploidy=ploidy, time=10 * i)
        for i, (n, ploidy) in enumerate(sample_ploidies)
    ]
    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        sequence_length=1,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        model=msprime.TPM(p=0.3, m=0.4, lo=2, hi=25),
        random_seed=seed,
    )
    # Garde-fou : un test « polymorphe » sans mutation serait vide de sens.
    assert (ts.num_sites == 0) == (mutation_rate == 0.0)
    counts = {f"s{i}": n for i, (n, _) in enumerate(sample_ploidies)}
    return ts, counts


def _layout_or_default(ts, counts, use_sample_layout):
    """Layout sériel (par échantillon) ou None (layout par population)."""
    return compute_sample_layout(ts, counts) if use_sample_layout else None


def _reference_individuals_by_sample(ts, layout):
    """Ancienne recherche individu par individu : {nom: [indices d'individus]}."""
    if layout is None:
        layout = compute_population_layout(ts)
    by_sample = {name: [] for name, _ in layout}
    for ind in ts.individuals():
        name = next(pop for pop, nodes in layout if ind.nodes[0] in nodes)
        by_sample[name].append(ind.id)
    return by_sample


def _reference_length_by_sample_and_individuals(ts, layout):
    """Ancienne version de `_length_by_sample_and_individuals` (boucle tskit)."""
    if layout is None:
        layout = compute_population_layout(ts)
    length_by_sample = {name: [] for name, _ in layout}
    if ts.num_sites == 0:
        for ind in ts.individuals():
            name = next(pop for pop, nodes in layout if ind.nodes[0] in nodes)
            length_by_sample[name].append((0, 0))
        return length_by_sample
    variant = next(ts.variants())
    tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
    for ind in ts.individuals():
        nodes = ind.nodes
        name = next(pop for pop, inds in layout if nodes[0] in inds)
        if len(nodes) == 1:
            length_by_sample[name].append(
                (int(tailles[nodes[0]]), int(tailles[nodes[0]]))
            )
        else:
            length_by_sample[name].append(
                (int(tailles[nodes[0]]), int(tailles[nodes[1]]))
            )
    return length_by_sample


def _reference_genotypes_by_sample_and_individuals(ts, layout):
    """Ancienne version de `_genotypes_by_sample_and_individuals` (boucle tskit)."""
    if layout is None:
        layout = compute_population_layout(ts)
    genotype_by_sample = {name: [] for name, _ in layout}
    if ts.num_sites == 0:
        for ind in ts.individuals():
            name = next(pop for pop, nodes in layout if ind.nodes[0] in nodes)
            genotype_by_sample[name].append((0, 0))
        return genotype_by_sample
    variant = next(ts.variants())
    tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
    for ind in ts.individuals():
        nodes = ind.nodes
        name = next(pop for pop, inds in layout if nodes[0] in inds)
        if len(nodes) == 1:
            genotype_by_sample[name].append((int(tailles[nodes[0]]),))
        else:
            genotype_by_sample[name].append(
                (int(tailles[nodes[0]]), int(tailles[nodes[1]]))
            )
    return genotype_by_sample


def test_individuals_by_sample_with_indiv_in_any_sample():
    """Vérifie que `_individuals_by_sample` fonctionne même si un individu est
    assigné à plusieurs échantillons (layout non valide)."""
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS["diploid"], 0.5)
    layout = compute_sample_layout(ts, counts)
    modified_layout = layout[:1]  # s0 et s1 partagent les mêmes individus

    with pytest.raises(ValueError, match="Certains individus ne sont assignés"):
        _individuals_by_sample(ts, layout=modified_layout)


def test_individuals_by_sample_with_indiv_in_two_samples():
    """Vérifie que `_individuals_by_sample` fonctionne même si un individu est
    assigné à deux échantillons."""
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS["diploid"], 0.5)
    layout = compute_sample_layout(ts, counts)
    (name_1, nodes_1), (name_2, nodes_2) = layout

    modified_layout = [
        (name_1, nodes_1),
        (name_2, np.concatenate([nodes_2, nodes_1[:2]])),
    ]
    with pytest.raises(ValueError, match="Certains individus sont déjà assignés"):
        _individuals_by_sample(ts, layout=modified_layout)


def _all_ploidy_cases(test):
    """Applique à un test les trois axes : scénario de ploïdie, polymorphisme, layout."""
    test = pytest.mark.parametrize(
        "use_sample_layout",
        [False, True],
        ids=["population_layout", "sample_layout"],
    )(test)
    test = pytest.mark.parametrize(
        "mutation_rate", [0.0, 0.5], ids=["monomorphic", "polymorphic"]
    )(test)
    return pytest.mark.parametrize("scenario", list(_PLOIDY_SCENARIOS))(test)


@_all_ploidy_cases
def test_length_by_sample_matches_reference(scenario, mutation_rate, use_sample_layout):
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS[scenario], mutation_rate)
    layout = _layout_or_default(ts, counts, use_sample_layout)

    result = _length_by_sample(ts, layout=layout)
    expected = _reference_length_by_sample(ts, layout=layout)

    assert list(result) == list(expected)  # ordre du layout
    for name in expected:
        assert result[name] == expected[name]  # ordre croissant des indices


@_all_ploidy_cases
def test_individuals_by_sample_matches_reference(
    scenario, mutation_rate, use_sample_layout
):
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS[scenario], mutation_rate)
    layout = _layout_or_default(ts, counts, use_sample_layout)

    result = _individuals_by_sample(ts, layout=layout)
    expected = _reference_individuals_by_sample(ts, layout)

    assert list(result) == list(expected)  # ordre du layout
    for name in expected:
        assert result[name].tolist() == expected[name]  # ordre croissant des indices


@_all_ploidy_cases
def test_length_by_sample_and_individuals_matches_reference(
    scenario, mutation_rate, use_sample_layout
):
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS[scenario], mutation_rate)
    layout = _layout_or_default(ts, counts, use_sample_layout)

    result = _length_by_sample_and_individuals(ts, layout=layout)
    expected = _reference_length_by_sample_and_individuals(ts, layout)

    assert list(result) == list(expected)
    assert result == expected
    assert all(type(x) is int for lst in result.values() for t in lst for x in t)


@_all_ploidy_cases
def test_genotypes_by_sample_and_individuals_matches_reference(
    scenario, mutation_rate, use_sample_layout
):
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS[scenario], mutation_rate)
    layout = _layout_or_default(ts, counts, use_sample_layout)

    result = _genotypes_by_sample_and_individuals(ts, layout=layout)
    expected = _reference_genotypes_by_sample_and_individuals(ts, layout)

    assert list(result) == list(expected)
    assert result == expected
    assert all(type(x) is int for lst in result.values() for t in lst for x in t)


def test_individuals_by_sample_returns_one_sorted_array_per_sample():
    """Un tableau d'indices par échantillon, pas une liste qui en contient un."""
    ts, counts = _ts_with_samples(((3, 2), (2, 2)), 0.5)

    result = _individuals_by_sample(ts, layout=compute_sample_layout(ts, counts))

    assert list(result) == ["s0", "s1"]
    assert all(isinstance(indices, np.ndarray) for indices in result.values())
    assert result["s0"].tolist() == [0, 1, 2]
    assert result["s1"].tolist() == [3, 4]


def test_individuals_by_sample_raises_when_an_individual_is_unassigned():
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS["diploid"], 0.5)
    layout = compute_sample_layout(ts, counts)

    # s1 absent du layout : ses individus n'appartiennent à personne.
    with pytest.raises(ValueError, match="individu"):
        _individuals_by_sample(ts, layout=layout[:1])


def test_individuals_by_sample_raises_when_an_individual_is_assigned_twice():
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS["diploid"], 0.5)
    layout = compute_sample_layout(ts, counts)

    # s0 deux fois : ses individus appartiennent à deux entrées du layout.
    with pytest.raises(ValueError, match="individu"):
        _individuals_by_sample(ts, layout=[*layout, ("s0_bis", layout[0][1])])


def test_genotypes_keep_true_ploidy_but_lengths_duplicate_haploids():
    """Sur un locus polymorphe, un haploïde donne (x,) pour les génotypes (LIK
    applique une formule différente selon la ploïdie) et (x, x) pour les longueurs
    (FST traite un haploïde comme un diploïde homozygote)."""
    ts, counts = _ts_with_samples(
        _PLOIDY_SCENARIOS["mixed"], 0.5
    )  # s0 diploïde, s1 haploïde
    layout = compute_sample_layout(ts, counts)

    genotypes = _genotypes_by_sample_and_individuals(ts, layout=layout)
    lengths = _length_by_sample_and_individuals(ts, layout=layout)

    assert {len(t) for t in genotypes["s0"]} == {2}
    assert {len(t) for t in genotypes["s1"]} == {1}
    assert all(first == second for first, second in lengths["s1"])
    assert [(x,) for x, _ in lengths["s1"]] == genotypes["s1"]


def test_genotypes_monomorphic_locus_gives_pairs_even_for_haploid():
    """Comportement historique conservé tel quel : sur un locus sans mutation,
    TOUS les individus reçoivent (0, 0), haploïdes compris (alors que sur un
    locus polymorphe un haploïde reçoit un tuple à 1 élément)."""
    ts, counts = _ts_with_samples(_PLOIDY_SCENARIOS["haploid"], 0.0)

    result = _genotypes_by_sample_and_individuals(
        ts, layout=compute_sample_layout(ts, counts)
    )

    assert result == {"s0": [(0, 0)] * 5, "s1": [(0, 0)] * 4}


def _reference_prepare_triplet_for_admixture(
    mutual_data: list[
        tuple[
            dict[str, list[tuple[int, int]]],
            dict[str, list[tuple[int, int] | tuple[int]]],
        ]
    ],
    focal: str,
    parent1: str,
    parent2: str,
) -> list[tuple[dict[int, float], dict[int, float], list[tuple[int, ...]]]]:
    """Ancienne version conservée comme oracle."""
    prepared_triplet = []
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    for length_by_sample, genotypes_by_sample in mutual_data:
        count_parent1 = {
            length: count
            for length, count in length_by_sample.get(parent1, [])
            if count > 0
        }
        count_parent2 = {
            length: count
            for length, count in length_by_sample.get(parent2, [])
            if count > 0
        }
        total_parent1 = sum(count_parent1.values())
        total_parent2 = sum(count_parent2.values())

        f1 = (
            {length: nb / total_parent1 for length, nb in count_parent1.items()}
            if total_parent1 > 0
            else {}
        )
        f2 = (
            {length: nb / total_parent2 for length, nb in count_parent2.items()}
            if total_parent2 > 0
            else {}
        )
        focal_genotypes = genotypes_by_sample.get(focal, [])

        prepared_triplet.append((f1, f2, focal_genotypes))

    return prepared_triplet


def test_prepared_parental_frequencies():
    """Test the prepared parental frequencies."""
    mutual_data = [
        (
            {"pop1": [(203, 3), (193, 1)], "pop2": [(203, 0), (197, 0)]},
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
            },
        )
    ]
    parent1 = "pop1"
    result1 = _prepared_parental_frequencies(mutual_data, parent1)
    assert result1 == [{203: 0.75, 193: 0.25}]

    parent2 = "pop2"
    result2 = _prepared_parental_frequencies(mutual_data, parent2)
    assert result2 == [{}]


def test_prepare_triplet_for_admixture_vs_reference():
    """Vérifie que `_prepare_triplet_for_admixture` est équivalent à la version de référence."""
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203, 203), (203, 193), (197, 203), (203, 193)],
            },
        )
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"

    result = _prepare_triplet_for_admixture(mutual_data, focal, parent1, parent2)

    expected = _reference_prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )

    assert result == expected


def _vectorize_triplet(prepared_triplet):
    """Référence conservée en tant qu'oracle"""
    # listes vides pour chaque groupe (haploïde : 2 listes, homozygote : 2, hétérozygote : 4)
    f1_haploid, f2_haploid = [], []
    f1_homozygous, f2_homozygous = [], []
    f1_heterozygous_x, f2_heterozygous_x = [], []
    f1_heterozygous_y, f2_heterozygous_y = [], []
    for f1, f2, focal_genotypes in prepared_triplet:
        for genotype in focal_genotypes:
            if len(genotype) == 1:  # haploid
                f1_haploid.append(f1.get(genotype[0], 0))
                f2_haploid.append(f2.get(genotype[0], 0))
            elif (
                len(genotype) == 2 and genotype[0] == genotype[1]
            ):  # homozygous diploid
                f1_homozygous.append(f1.get(genotype[0], 0))
                f2_homozygous.append(f2.get(genotype[0], 0))
            elif (
                len(genotype) == 2 and genotype[0] != genotype[1]
            ):  # heterozygous diploid
                f1_heterozygous_x.append(f1.get(genotype[0], 0))
                f2_heterozygous_x.append(f2.get(genotype[0], 0))
                f1_heterozygous_y.append(f1.get(genotype[1], 0))
                f2_heterozygous_y.append(f2.get(genotype[1], 0))
            else:
                raise ValueError(
                    f"Problème avec le génotype de l'individu : {genotype}"
                )
    # conversion en np.array(dtype=float) et retour
    f1_haploid = np.array(f1_haploid, dtype=float)
    f2_haploid = np.array(f2_haploid, dtype=float)
    f1_homozygous = np.array(f1_homozygous, dtype=float)
    f2_homozygous = np.array(f2_homozygous, dtype=float)
    f1_heterozygous_x = np.array(f1_heterozygous_x, dtype=float)
    f2_heterozygous_x = np.array(f2_heterozygous_x, dtype=float)
    f1_heterozygous_y = np.array(f1_heterozygous_y, dtype=float)
    f2_heterozygous_y = np.array(f2_heterozygous_y, dtype=float)

    focal_genotypes = FocalGenotypeFrequencies(
        haploid_f1=f1_haploid,
        haploid_f2=f2_haploid,
        homozygous_f1=f1_homozygous,
        homozygous_f2=f2_homozygous,
        heterozygous_f1_x=f1_heterozygous_x,
        heterozygous_f2_x=f2_heterozygous_x,
        heterozygous_f1_y=f1_heterozygous_y,
        heterozygous_f2_y=f2_heterozygous_y,
    )

    return focal_genotypes


def _prepare_triplet_for_admixture(
    mutual_data: list[
        tuple[
            dict[str, list[tuple[int, int]]],
            dict[str, list[tuple[int, int] | tuple[int]]],
        ]
    ],
    focal: str,
    parent1: str,
    parent2: str,
    *,
    frequency_cache: dict[str, list[dict[int, float]]] | None = None,
) -> list[tuple[dict[int, float], dict[int, float], list[tuple[int, ...]]]]:
    """Référence gardée comme oracle"""
    prepared_triplet = []
    parent1_frequencies = _parental_frequencies_with_cache(
        mutual_data, parent1, frequency_cache
    )
    parent2_frequencies = _parental_frequencies_with_cache(
        mutual_data, parent2, frequency_cache
    )
    for (f1, f2), (_, genotypes_by_sample) in zip(
        zip(parent1_frequencies, parent2_frequencies, strict=True),
        mutual_data,
        strict=True,
    ):
        focal_genotypes = genotypes_by_sample.get(focal, [])
        prepared_triplet.append((f1, f2, focal_genotypes))

    return prepared_triplet


def test_vectorize_triplet():
    """Vérifie que _vectorize_triplet renvoie bien les bonnes fréqeunces génotypiques"""

    mutual_data1 = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        )
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    prepared_triplet1 = _prepare_triplet_for_admixture(
        mutual_data1, focal, parent1, parent2
    )
    vectorized_triplet1 = _vectorize_triplet(prepared_triplet1)

    assert vectorized_triplet1.haploid_f1.tolist() == [0.75]
    assert vectorized_triplet1.haploid_f2.tolist() == [2 / 3]
    assert vectorized_triplet1.homozygous_f1.tolist() == [0.75]
    assert vectorized_triplet1.homozygous_f2.tolist() == [2 / 3]
    assert vectorized_triplet1.heterozygous_f1_x.tolist() == [0.0, 0.75]
    assert vectorized_triplet1.heterozygous_f1_y.tolist() == [0.75, 0.25]
    assert vectorized_triplet1.heterozygous_f2_x.tolist() == [1 / 3, 2 / 3]
    assert vectorized_triplet1.heterozygous_f2_y.tolist() == [2 / 3, 0]


def _assert_same_focal_frequencies(result, expected):
    for field in dataclasses.fields(result):
        a, b = getattr(result, field.name), getattr(expected, field.name)
        assert a.shape == b.shape and a.dtype == b.dtype, field.name
        assert np.array_equal(a, b), field.name


def test_focal_genotype_frequencies_vs_reference():
    """Vérifie que _focal_genotype_frequencies est équivalent à la version de référence."""
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        ),
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        ),
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )
    freq1 = _prepared_parental_frequencies(mutual_data, parent1)
    freq2 = _prepared_parental_frequencies(mutual_data, parent2)
    freq1_cache = {parent1: freq1}
    freq_cache = {parent1: freq1, parent2: freq2}

    result1 = _focal_genotype_frequencies(
        mutual_data, focal, parent1, parent2, frequency_cache=freq1_cache
    )
    result_empty = _focal_genotype_frequencies(mutual_data, focal, parent1, parent2)
    result_full = _focal_genotype_frequencies(
        mutual_data, focal, parent1, parent2, frequency_cache=freq_cache
    )
    expected = _vectorize_triplet(prepared_triplet)

    _assert_same_focal_frequencies(result1, expected)
    _assert_same_focal_frequencies(result_empty, expected)
    _assert_same_focal_frequencies(result_full, expected)


def test_vectorize_triplet_with_problematic_genotype():
    """Vérifie que _vectorize_triplet renvoie une erreur si on a un génotype problématique"""
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203, 197, 201), (197, 197), (197, 197), (197, 197)],
            },
        )
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )
    with pytest.raises(ValueError, match="Problème avec le génotype de l'individu"):
        _vectorize_triplet(prepared_triplet)


def test_vectorize_triplet_for_focal():
    """Vérifie que _vectorize_triplet renvoie bien les bonnes fréquences génotypiques pour le focal"""
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        )
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    parent_freq1 = _prepared_parental_frequencies(mutual_data, parent1)
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )
    focal_genotypes_per_locus = [
        focal_genotypes for _, _, focal_genotypes in prepared_triplet
    ]
    haploid, homozygous, heterozygous_x, heterozygous_y = _vectorize_parent_for_focal(
        parent_frequencies=parent_freq1,
        focal_genotypes_per_locus=focal_genotypes_per_locus,
    )

    assert haploid.tolist() == [0.75]
    assert homozygous.tolist() == [0.75]
    assert heterozygous_x.tolist() == [0, 0.75]
    assert heterozygous_y.tolist() == [0.75, 0.25]


def test_vectorize_triplet_for_focal_with_problematic_genotype():
    """Vérifie que _vectorize_triplet renvoie une erreur si on a un génotype problématique pour le focal"""
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203, 203, 203), (203, 203), (197, 203), (203, 193)],
            },
        )
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )
    with pytest.raises(ValueError, match="Problème avec le génotype de l'individu"):
        _vectorize_triplet(prepared_triplet)


def test_vectorize_triplet_for_focal_vs_reference():
    mutual_data = [
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        ),
        (
            {
                "pop1": [(203, 3), (193, 1)],
                "pop2": [(203, 2), (197, 1)],
                "pop3": [(203, 0), (197, 0)],
            },
            {
                "pop1": [(203, 203), (203, 193), (203, 203), (203, 203)],
                "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
                "pop3": [(203,), (203, 203), (197, 203), (203, 193)],
            },
        ),
    ]
    focal = "pop3"
    parent1 = "pop1"
    parent2 = "pop2"
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )
    focal_genotypes_per_locus = [
        focal_genotypes for _, _, focal_genotypes in prepared_triplet
    ]
    parent_freq1 = _prepared_parental_frequencies(mutual_data, parent1)
    haploid, homozygous, heterozygous_x, heterozygous_y = _vectorize_parent_for_focal(
        parent_frequencies=parent_freq1,
        focal_genotypes_per_locus=focal_genotypes_per_locus,
    )
    focal_genotypes = _vectorize_triplet(prepared_triplet)
    assert haploid.tolist() == focal_genotypes.haploid_f1.tolist()
    assert homozygous.tolist() == focal_genotypes.homozygous_f1.tolist()
    assert heterozygous_x.tolist() == focal_genotypes.heterozygous_f1_x.tolist()
    assert heterozygous_y.tolist() == focal_genotypes.heterozygous_f1_y.tolist()


def _reference_log_likelihood_admixture(
    data: FocalGenotypeFrequencies,
    a: float,
) -> float:
    """Ancienne version conservée comme oracle"""
    lik = 0.0

    # haploïd
    freq = a * data.haploid_f1 + (1 - a) * data.haploid_f2
    mask = freq > 0
    lik += np.log(freq[mask]).sum()

    # homozygous diploid
    freq = a * data.homozygous_f1 + (1 - a) * data.homozygous_f2
    mask = freq > 0
    lik += np.log(freq[mask] ** 2).sum()

    # heterozygous diploid
    freq_x = a * data.heterozygous_f1_x + (1 - a) * data.heterozygous_f2_x
    freq_y = a * data.heterozygous_f1_y + (1 - a) * data.heterozygous_f2_y
    mask = freq_x * freq_y > 0
    lik += np.log(2 * freq_x[mask] * freq_y[mask]).sum()
    return float(lik)


_ADMIXTURE_CASES = {
    "haploid": (np.array([0.5, 0.5]), np.array([0.25, 0.75])),
    "homozygous": (np.array([0.6, 0.4]), np.array([0.3, 0.7])),
    "heterozygous": (
        np.array([0.2, 0.8]),
        np.array([0.7, 0.3]),
        np.array([0.5, 0.5]),
        np.array([0.4, 0.6]),
    ),
}


_BLOCKS = {
    "haploid": ("haploid_f1", "haploid_f2"),
    "homozygous": ("homozygous_f1", "homozygous_f2"),
    "heterozygous": (
        "heterozygous_f1_x",
        "heterozygous_f1_y",
        "heterozygous_f2_x",
        "heterozygous_f2_y",
    ),
}


@pytest.mark.parametrize("a", [0.0, 0.001, 0.6, 0.999, 1.0])
@pytest.mark.parametrize(
    "empty_blocks", list(itertools.product([False, True], repeat=3))
)
def test_log_likelihood_admixture_vs_reference(empty_blocks, a):
    """Vérifie que _log_likelihood_admixture est équivalent à la version de référence."""
    data = FocalGenotypeFrequencies(
        haploid_f1=np.array([0.5, 0.5]),
        haploid_f2=np.array([0.25, 0.75]),
        homozygous_f1=np.array([0.6, 0.4]),
        homozygous_f2=np.array([0.3, 0.7]),
        heterozygous_f1_x=np.array([0.2, 0.8]),
        heterozygous_f1_y=np.array([0.7, 0.3]),
        heterozygous_f2_x=np.array([0.5, 0.5]),
        heterozygous_f2_y=np.array([0.4, 0.6]),
    )

    for is_empty, fields in zip(empty_blocks, _BLOCKS.values(), strict=True):
        if is_empty:
            data = dataclasses.replace(
                data, **{field: np.array([]) for field in fields}
            )
    assert _log_likelihood_admixture(data, a) == _reference_log_likelihood_admixture(
        data, a
    )

    expected = _reference_log_likelihood_admixture(data, a)
    result = _log_likelihood_admixture(data, a)
    assert result == expected


@pytest.mark.parametrize(
    "field",
    [
        "haploid_f1",
        "haploid_f2",
        "homozygous_f1",
        "homozygous_f2",
        "heterozygous_f1_x",
        "heterozygous_f1_y",
        "heterozygous_f2_x",
        "heterozygous_f2_y",
    ],
)
def test_log_likelihood_admixture_inconsistent_block_raises(field):
    data = FocalGenotypeFrequencies(
        haploid_f1=np.array([0.5, 0.5]),
        haploid_f2=np.array([0.25, 0.75]),
        homozygous_f1=np.array([0.6, 0.4]),
        homozygous_f2=np.array([0.3, 0.7]),
        heterozygous_f1_x=np.array([0.2, 0.8]),
        heterozygous_f1_y=np.array([0.7, 0.3]),
        heterozygous_f2_x=np.array([0.5, 0.5]),
        heterozygous_f2_y=np.array([0.4, 0.6]),
    )
    data = dataclasses.replace(data, **{field: np.array([])})
    with pytest.raises(ValueError, match="broadcast"):
        _reference_log_likelihood_admixture(data, 0.6)
    with pytest.raises(ValueError, match="broadcast"):
        _log_likelihood_admixture(data, 0.6)


def test_compute_num_den_lik_for_one_individual():
    """Vérifie que la fonction _compute_num_den_lik_for_one_individual
    fonctionne correctement."""
    count = {203: 1, 200: 1, 197: 1}
    total_count = 3
    b = 1

    # test haploïde
    genotype_haploid = (203,)
    num, den = _compute_num_den_lik_for_one_individual(
        genotype_haploid, count, total_count, b
    )
    assert num == 2
    assert den == 4

    # test diploïde homozygote
    genotype_diploid_homo = (203, 203)
    num, den = _compute_num_den_lik_for_one_individual(
        genotype_diploid_homo, count, total_count, b
    )
    assert num == 6
    assert den == 20

    # test diploïde heteozygote
    genotype_diploid_hetero = (203, 200)
    num, den = _compute_num_den_lik_for_one_individual(
        genotype_diploid_hetero, count, total_count, b
    )
    assert num == 8
    assert den == 20


def test_compute_LIK_for_one_locus():
    _length_by_sample = {
        "pop1": [(203, 203), (203, 203), (203, 193), (203, 203)],
        "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
    }
    genotypes_by_pop = {
        "pop1": [(203, 203), (203, 203), (203, 193), (203,)],
        "pop2": [(203, 203), (203, 193), (197, 203), (203, 193)],
    }

    lik, present = _compute_LIK_for_one_locus(
        _length_by_sample, genotypes_by_pop, "pop1", "pop2"
    )
    assert pytest.approx(lik) == 1.1175620386325875
    assert present

    # test avec une population absente
    lik, present = _compute_LIK_for_one_locus(
        _length_by_sample, genotypes_by_pop, "pop3", "pop2"
    )
    assert pytest.approx(lik) == 0.0
    assert not present


def test_compute_LIK(microsat_context_te2_xy):
    """Vérifie compute_LIK sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_LIK(mutated.values(), sample_names)

    assert len(results) == 2
    assert (
        results["1.2"] != results["2.1"]
    )  # pas de garanti que ce soit vrai sur tous les datasets
    assert pytest.approx(results["1.2"]) == 1.7822692958


# Test pour la fonction d'entrée principale de calcul des stats microsatellites
def test_compute_all_statistics_microsat(microsat_context_te2_xy):
    """Vérifie compute_all_statistics_microsat sur toy_example2_ms_dna."""
    demography, _ = build_random_demography_for_scenario_index(
        microsat_context_te2_xy.header_text, scenario_index=1, seed=42
    )
    mutated, _ = microsat_mutation_simulation_per_locus(
        demography=demography,
        context=microsat_context_te2_xy,
        seed=42,
    )

    sample_names = ["pop1", "pop2"]

    results = compute_all_statistics_microsat(
        microsat_context_te2_xy.header_text, mutated, sample_names, seed=42
    )

    expected_keys = {"LIK_1_1.2", "MGW_1_1", "FST_1_1.2"}

    assert len(results) == 16
    assert expected_keys.issubset(results.keys())
    assert pytest.approx(results["HET_1_1"]) == 0.8163225371
    assert pytest.approx(results["LIK_1_1.2"]) == 1.7822692958


def _synthetic_aml_inputs(shift):
    """Données synthétiques pour AML : 4 échantillons diploïdes, 2 loci.

    `shift` décale les comptes pour fabriquer deux jeux de données différents.

    Returns:
        (sample_names, length_by_sample_per_locus, genotypes_by_sample_per_locus).
    """
    names = ["pop1", "pop2", "pop3", "pop4"]
    all_genotypes = [
        (200, 200),
        (200, 202),
        (202, 202),
        (202, 204),
        (204, 204),
        (200, 204),
    ]
    length_by_sample_per_locus = []
    genotypes_by_sample_per_locus = []
    for locus in range(2):
        length_by_sample = {}
        genotypes_by_sample = {}
        for k, name in enumerate(names):
            length_by_sample[name] = [
                (200, 2 + k + shift + locus),
                (202, 6 - k + locus),
                (204, 1 + (k * (locus + 1)) % 3 + 2 * shift),
            ]
            genotypes_by_sample[name] = all_genotypes[shift + k % 2 : 5 + locus]
        length_by_sample_per_locus.append(length_by_sample)
        genotypes_by_sample_per_locus.append(genotypes_by_sample)
    return names, length_by_sample_per_locus, genotypes_by_sample_per_locus


def test_aml_parental_frequencies_computed_once_per_parent(monkeypatch):
    """Le cache de `compute_AML_microsat` ne calcule les fréquences qu'une fois par parent.

    4 échantillons, 12 triplets, 24 utilisations d'un parent : sans cache
    effectif, 24 appels ; avec lui, 4 (un par échantillon).
    """
    import bridge.summary_statistics as ss_module

    names, lengths, genotypes = _synthetic_aml_inputs(0)
    calls = []
    original = ss_module._prepared_parental_frequencies

    def spy(mutual_data, parent):
        calls.append(parent)
        return original(mutual_data, parent)

    monkeypatch.setattr(ss_module, "_prepared_parental_frequencies", spy)
    ss_module.compute_AML_microsat(
        [None, None],
        names,
        seed=1,
        length_by_sample_per_locus=lengths,
        genotypes_by_sample_per_locus=genotypes,
    )

    assert sorted(calls) == names


def test_aml_cache_matches_uncached_and_does_not_leak_between_calls(monkeypatch):
    """Avec cache, `compute_AML_microsat` donne le résultat sans cache, sur deux appels successifs.

    Deux jeux de données différents l'un après l'autre : un cache qui
    survivrait à l'appel servirait au second les fréquences du premier.
    """
    import bridge.summary_statistics as ss_module

    results = {}
    for shift in (0, 1):
        names, lengths, genotypes = _synthetic_aml_inputs(shift)
        kwargs = {
            "seed": 1,
            "length_by_sample_per_locus": lengths,
            "genotypes_by_sample_per_locus": genotypes,
        }
        cached = ss_module.compute_AML_microsat([None, None], names, **kwargs)
        with monkeypatch.context() as patch:
            patch.setattr(
                ss_module,
                "_parental_frequencies_with_cache",
                lambda mutual_data, parent, cache: (
                    ss_module._prepared_parental_frequencies(mutual_data, parent)
                ),
            )
            uncached = ss_module.compute_AML_microsat([None, None], names, **kwargs)
        assert cached == uncached
        results[shift] = cached

    assert results[0] != results[1]  # sinon le test ne distingue pas une fuite


def test_aml_parent_vectorization_computed_once_per_focal_parent_pair(monkeypatch):
    """Le cache de vectorisation ne vectorise qu'une fois chaque couple (focal, parent).

    4 échantillons : 12 triplets, chacun utilise deux couples (focal, parent),
    soit 24 utilisations de 12 couples distincts. Sans cache effectif, 24
    appels ; avec lui, 12.
    """
    import bridge.summary_statistics as ss_module

    names, lengths, genotypes = _synthetic_aml_inputs(0)
    calls = []
    original = ss_module._vectorize_parent_for_focal

    def spy(parent_frequencies, focal_genotypes_per_locus):
        calls.append(1)
        return original(parent_frequencies, focal_genotypes_per_locus)

    monkeypatch.setattr(ss_module, "_vectorize_parent_for_focal", spy)
    ss_module.compute_AML_microsat(
        [None, None],
        names,
        seed=1,
        length_by_sample_per_locus=lengths,
        genotypes_by_sample_per_locus=genotypes,
    )

    assert len(calls) == 12


def test_pool_allele_counts_matches_reference_on_random_inputs():
    """`_pool_allele_counts_for_two_samples` donne la même liste que l'ancienne version.

    Listes tirées au hasard : sous-ensembles d'allèles dans un ordre quelconque
    (donc des allèles présents dans un seul des deux échantillons, dont certains
    seulement dans `pop_2`, qui doivent être ignorés), comptes éventuellement
    nuls, listes vides, et cas `n <= 1` qui renvoie `[(0, 0.0)]`.
    """
    rng = np.random.default_rng(2024)
    alleles = [190 + 2 * k for k in range(8)]

    def draw_sample():
        chosen = rng.choice(alleles, size=int(rng.integers(0, 8)), replace=False)
        return [(int(allele), int(rng.integers(0, 6))) for allele in chosen]

    cases = [
        ([], []),
        ([(200, 1)], []),
        ([], [(200, 1)]),
        ([(200, 0), (202, 1)], [(202, 0)]),
        ([(200, 3)], [(204, 5)]),  # allèle de pop_2 absent de pop_1 : ignoré
    ]
    cases += [(draw_sample(), draw_sample()) for _ in range(300)]

    for pop_1, pop_2 in cases:
        result = _pool_allele_counts_for_two_samples(pop_1, pop_2)
        expected = _reference_pool_allele_counts_for_two_samples(pop_1, pop_2)
        assert result == expected, (pop_1, pop_2)
