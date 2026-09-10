"""
Tests des statistiques résumées SNP PoolSeq (bridge/summary_statistics.py).

_prepare_matrices_poolseq : forme/valeurs des matrices (npop, nloci).
compute_all_statistics_poolseq : les 130 stats sont bien présentes, avec au
moins une valeur vérifiée à la main (HWm_1/HWv_1) pour attraper une
régression de formule, pas juste un problème de branchement.
"""

import numpy as np
import pytest
from conftest import (
    OBSERVED_MSS_FILE_TE2,
    OBSERVED_MSS_FILE_TE2_XY,
)

from bridge.ancestry_simulation import (
    dna_mutation_simulation_per_locus,
    microsat_mutation_simulation_per_locus,
)
from bridge.loci_parser import parse_loci_description
from bridge.pipeline import build_random_demography_for_scenario_index
from bridge.summary_statistics import (
    _combined_alleles_for_two_populations,
    _compute_MGW_by_locus,
    _compute_VAR_for_one_population,
    _count_combined_alleles_for_one_locus,
    _genotype_matrix_by_population,
    _length_by_population,
    _prepare_matrices_poolseq,
    _VAR_constants,
    compute_all_statistics_dna,
    compute_all_statistics_microsat,
    compute_all_statistics_poolseq,
    compute_DTA,
    compute_HET,
    compute_HST,
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
    compute_VAR,
    compute_VNS,
    compute_VPD,
    count_alleles_per_population,
    total_genes_copies_per_population,
)


def test_prepare_matrices_poolseq():
    reads_per_locus = [
        {"POP1": (0, 93), "POP2": (0, 100), "POP3": (1, 116), "POP4": (0, 139)},
        {"POP1": (1, 80), "POP2": (0, 90), "POP3": (0, 110), "POP4": (1, 120)},
    ]
    population_names = ["POP1", "POP2", "POP3", "POP4"]
    counts, ns, freq0, freq1 = _prepare_matrices_poolseq(
        reads_per_locus, population_names
    )

    assert counts.shape == (4, 2)
    assert np.array_equal(counts, np.array([[0, 1], [0, 0], [1, 0], [0, 1]]))


def test_compute_all_statistics_poolseq():
    """Vérifie que compute_all_statistics_poolseq renvoie bien un dictionnaire de statistiques
    pour les fichiers POOLSEQ toy_example4 et toy_example5."""
    # Préparer des données fictives pour le test
    reads_per_locus = [
        {"POP1": (0, 93), "POP2": (0, 100), "POP3": (1, 116), "POP4": (0, 139)},
        {"POP1": (1, 80), "POP2": (0, 90), "POP3": (0, 110), "POP4": (1, 120)},
    ]
    population_names = ["POP1", "POP2", "POP3", "POP4"]
    pool_sizes = {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}

    results = compute_all_statistics_poolseq(
        reads_per_locus, population_names, pool_sizes
    )

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


def test_genotype_matrix_by_population(header_text_te2):
    """Vérifie _genotype_matrix_by_population sur un locus <A> (diploïde)
    et un locus <M> (haploïde) du même dataset : la forme retournée doit
    respecter le nombre de sites/samples réels de la TreeSequence, sans
    perte ni doublon de sample entre populations, et le nombre de samples
    par population doit refléter la ploïdie du locus (rapport 2:1 entre
    <A> et <M> -- couvre le bug de ploïdie corrigé le 2026-08-24 dans
    dna_ancestry_parameters_for_heritage)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    def samples_per_population(ts):
        return {
            population.metadata["name"]: ts.samples(population=population.id)
            for population in ts.populations()
        }

    for locus_name in ("Locus_S_A_11_", "Locus_S_M_16_"):
        ts = mutated[locus_name]
        result = _genotype_matrix_by_population(ts)
        expected_samples = samples_per_population(ts)

        assert result.keys() == {"pop1", "pop2"}
        total_samples = 0
        for pop_name, matrix in result.items():
            assert matrix.shape[0] == ts.num_sites
            assert matrix.shape[1] == len(expected_samples[pop_name])
            total_samples += matrix.shape[1]
        assert total_samples == ts.num_samples

    samples_a = _genotype_matrix_by_population(mutated["Locus_S_A_11_"])
    samples_m = _genotype_matrix_by_population(mutated["Locus_S_M_16_"])
    for pop_name in ("pop1", "pop2"):
        assert samples_a[pop_name].shape[1] == 2 * samples_m[pop_name].shape[1]


def test_mean_segregating_sites_per_group(header_text_te2):
    """Vérifie compute_NSS séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre NSS_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NSS(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"pop1": pytest.approx(5.8), "pop2": pytest.approx(5.6)}

    mean_g3 = compute_NSS(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"pop1": pytest.approx(7.2), "pop2": pytest.approx(7.4)}


def test_mean_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    assert compute_NSS([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group(header_text_te2):
    """Vérifie compute_NHA séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre NH_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NHA(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"pop1": pytest.approx(5.2), "pop2": pytest.approx(5.4)}

    mean_g3 = compute_NHA(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"pop1": pytest.approx(6.0), "pop2": pytest.approx(6.0)}


def test_mean_distinct_haplotypes_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = compute_NHA(tree_sequences_empty, population_names)
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_pairwise_differences_per_group(header_text_te2):
    """Vérifie compute_MPD séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre VPD_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_MPD(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {
        "pop1": pytest.approx(1.1997435897435897),
        "pop2": pytest.approx(1.2999999999999998),
    }

    mean_g3 = compute_MPD(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {
        "pop1": pytest.approx(1.5768421052631578),
        "pop2": pytest.approx(1.8084210526315794),
    }


def test_mean_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = compute_MPD(tree_sequences_empty, population_names)
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def test_variance_pairwise_differences_per_group(header_text_te2):
    """Vérifie compute_VPD séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre VPD_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    variance_g2 = compute_VPD(tree_sequences_for_group("G2"), population_names)
    assert variance_g2 == {
        "pop1": pytest.approx(2.0666373720417366),
        "pop2": pytest.approx(1.7601968335472828),
    }

    variance_g3 = compute_VPD(tree_sequences_for_group("G3"), population_names)
    assert variance_g3 == {
        "pop1": pytest.approx(2.6465162907268174),
        "pop2": pytest.approx(2.935839598997494),
    }


def test_variance_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    variance_empty = compute_VPD(tree_sequences_empty, population_names)
    assert variance_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_tajima_d_per_group(header_text_te2):
    """Vérifie compute_DTA séparément sur G2 (<A>) et G3 (<M>)
    de toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque
    chaque `group Gx` du header calcule son propre DTA_i à partir de ses
    seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    dta_g2 = compute_DTA(tree_sequences_for_group("G2"), population_names)
    assert dta_g2 == {
        "pop1": pytest.approx(-0.11361849459656947),
        "pop2": pytest.approx(-0.2621067146690289),
    }

    dta_g3 = compute_DTA(tree_sequences_for_group("G3"), population_names)
    assert dta_g3 == {
        "pop1": pytest.approx(-0.6063248006420837),
        "pop2": pytest.approx(-0.5259556352897572),
    }


def test_mean_tajima_d_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    assert compute_DTA([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_private_segregating_sites_per_group(header_text_te2):
    """Vérifie compute_PSS séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    PSS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    pss_g2 = compute_PSS(tree_sequences_for_group("G2"), population_names)
    assert pss_g2 == {"pop1": pytest.approx(1.2), "pop2": pytest.approx(1.0)}

    pss_g3 = compute_PSS(tree_sequences_for_group("G3"), population_names)
    assert pss_g3 == {"pop1": pytest.approx(3.6), "pop2": pytest.approx(3.8)}


def test_mean_private_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_PSS([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_minor_allele_count_per_group(header_text_te2):
    """Vérifie compute_MNS séparément sur G2 (<A>)
    et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    MNS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mns_g2 = compute_MNS(tree_sequences_for_group("G2"), population_names)
    assert mns_g2 == {
        "pop1": pytest.approx(3.486274509803921),
        "pop2": pytest.approx(4.0),
    }

    mns_g3 = compute_MNS(tree_sequences_for_group("G3"), population_names)
    assert mns_g3 == {
        "pop1": pytest.approx(2.8034188034188032),
        "pop2": pytest.approx(2.7534065934065937),
    }


def test_mean_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_MNS([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_variance_minor_allele_count_per_group(header_text_te2):
    """Vérifie compute_VNS séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    VNS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    vns_g2 = compute_VNS(tree_sequences_for_group("G2"), population_names)
    assert vns_g2 == {
        "pop1": pytest.approx(8.9039600153787),
        "pop2": pytest.approx(18.444444444444446),
    }

    vns_g3 = compute_VNS(tree_sequences_for_group("G3"), population_names)
    assert vns_g3 == {
        "pop1": pytest.approx(1.2056680546424137),
        "pop2": pytest.approx(2.6550416616350683),
    }


def test_variance_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_VNS([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group_pairwize(header_text_te2):
    """Vérifie compute_NH2 séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    NH_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NH2(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"1.2": pytest.approx(6.8)}

    mean_g3 = compute_NH2(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"1.2": pytest.approx(10.0)}


def test_mean_distinct_haplotypes_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_NH2([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_segregating_sites_per_group_pairwize(header_text_te2):
    """Vérifie compute_NS2 séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    NSS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_NS2(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"1.2": pytest.approx(6.8)}

    mean_g3 = compute_NS2(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"1.2": pytest.approx(11.0)}


def test_mean_segregating_sites_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_NS2([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_pairwise_differences_per_group_pairwize(header_text_te2):
    """Vérifie compute_MP2 séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    VPD_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_MP2(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"1.2": pytest.approx(1.2498717948717948)}

    mean_g3 = compute_MP2(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"1.2": pytest.approx(1.6926315789473683)}


def test_mean_pairwise_differences_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_MP2([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_pairwise_differences_between_per_group_pairwize(header_text_te2):
    """Vérifie compute_MPB séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    VPD_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_MPB(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"1.2": pytest.approx(1.28725)}

    mean_g3 = compute_MPB(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"1.2": pytest.approx(1.7454999999999998)}


def test_mean_pairwise_differences_between_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert compute_MPB([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_hst_per_group_pairwize(header_text_te2):
    """Vérifie compute_HST séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    HST_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = compute_HST(tree_sequences_for_group("G2"), population_names)
    assert mean_g2 == {"1.2": pytest.approx(0.029037253935292443)}

    mean_g3 = compute_HST(tree_sequences_for_group("G3"), population_names)
    assert mean_g3 == {"1.2": pytest.approx(0.03028841080070557)}


def test_compute_all_statistics_dna(header_text_te2):
    """Vérifie compute_all_statistics_dna sur toy_example2_ms_dna."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_MSS_FILE_TE2,
        seed=42,
    )

    population_names = ["pop1", "pop2"]

    results = compute_all_statistics_dna(header_text_te2, mutated, population_names)

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
    assert pytest.approx(results["NSS_2_1"]) == 5.8
    assert pytest.approx(results["HST_2_1.2"]) == 0.029037253935292443
    assert pytest.approx(results["NH2_3_1.2"]) == 10.0


# -----------------------------------------------------------------------
# Pour les microsat
# -----------------------------------------------------------------------


def test_length_by_population(header_text_te2_XY):
    """Vérifie _length_by_population sur un locus microsat (Locus_M_A_1_).

    Convertit les codes de génotype en tailles réelles (pb) via
    variant.alleles, puis compte les copies de gène par taille et par
    population -- couvre les deux populations pour éviter de repasser
    silencieusement au bug de découpage corrigé pendant l'écriture de
    cette fonction (indexation par sample_ids plutôt que par un
    intervalle [first_index, last_index] supposé contigu)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    ts = mutated["Locus_M_A_1_"]
    results = _length_by_population(ts)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == [
        (201, 0),
        (203, 31),
        (197, 1),
        (193, 7),
        (199, 0),
        (189, 0),
    ]
    assert results["pop2"] == [
        (201, 0),
        (203, 18),
        (197, 0),
        (193, 20),
        (199, 1),
        (189, 1),
    ]


def test_count_alleles_per_population():
    length_per_population = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    results = count_alleles_per_population(length_per_population)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == 3
    assert results["pop2"] == 4


def test_compute_NAL(header_text_te2_XY):
    """Vérifie compute_NAL sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    population_names = ["pop1", "pop2"]

    results = compute_NAL(mutated.values(), population_names)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == 9.1
    assert pytest.approx(results["pop2"]) == 79 / 9


def test_total_genes_copies_per_population(header_text_te2_XY):
    """Vérifie _total_genes_copies_per_population sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    length_by_population = _length_by_population(mutated["Locus_M_A_1_"])
    results = total_genes_copies_per_population(length_by_population)

    assert results.keys() == {"pop1", "pop2"}
    assert results["pop1"] == 39
    assert results["pop2"] == 40


def test_compute_HET(header_text_te2_XY):
    """Vérifie compute_HET sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    tree_sequences1 = [mutated["Locus_M_A_1_"]]

    population_names = ["pop1", "pop2"]

    results = compute_HET(tree_sequences1, population_names)

    assert results.keys() == {"pop1", "pop2"}
    assert (
        pytest.approx(results["pop1"])
        == (1 - (31 / 39) ** 2 - (7 / 39) ** 2 - (1 / 39) ** 2) * 39 / 38
    )
    assert (
        pytest.approx(results["pop2"])
        == (1 - (20 / 40) ** 2 - (18 / 40) ** 2 - (1 / 40) ** 2 - (1 / 40) ** 2)
        * 40
        / 39
    )

    tree_sequences2 = list(mutated.values())
    results = compute_HET(tree_sequences2, population_names)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 7.929892037786773 / 10
    assert pytest.approx(results["pop2"]) == 7.244601889338731 / 9


def test_VAR_constants(header_text_te2_XY):
    """Vérifie que les constantes de compute_VAR sont correctes."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    length_by_population = _length_by_population(mutated["Locus_M_A_1_"])
    raw_sizes, raw_square_sizes, total_counts = _VAR_constants(length_by_population)

    assert raw_sizes.keys() == {"pop1", "pop2"}
    assert raw_square_sizes.keys() == {"pop1", "pop2"}
    assert total_counts.keys() == {"pop1", "pop2"}

    assert raw_sizes["pop1"] == 203 * 31 + 197 * 1 + 193 * 7 + 199 * 0 + 189 * 0
    assert (
        raw_square_sizes["pop1"]
        == 203**2 * 31 + 197**2 * 1 + 193**2 * 7 + 199**2 * 0 + 189**2 * 0
    )
    assert total_counts["pop1"] == 39

    assert raw_sizes["pop2"] == 203 * 18 + 197 * 0 + 193 * 20 + 199 * 1 + 189 * 1
    assert (
        raw_square_sizes["pop2"]
        == 203**2 * 18 + 197**2 * 0 + 193**2 * 20 + 199**2 * 1 + 189**2 * 1
    )
    assert total_counts["pop2"] == 40


def test_compute_VAR_for_one_population(header_text_te2_XY):
    """Vérifie compute_VAR sur toy_example2_ms_dna_xy pour une seule population."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    s, v, n = _VAR_constants(_length_by_population(mutated["Locus_M_A_1_"]))
    result = _compute_VAR_for_one_population(s["pop1"], v["pop1"], n["pop1"], 2)

    s1, v1, n1 = (
        203 * 31 + 197 * 1 + 193 * 7 + 199 * 0 + 189 * 0,
        203**2 * 31 + 197**2 * 1 + 193**2 * 7 + 199**2 * 0 + 189**2 * 0,
        39,
    )
    assert pytest.approx(result) == (v1 - s1**2 / n1) / (n1 - 1) / 2**2


def test_compute_VAR(header_text_te2_XY):
    """Vérifie compute_VAR sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    population_names = ["pop1", "pop2"]
    list_loci = [
        locus
        for locus in parse_loci_description(header_text_te2_XY)
        if locus.ms_or_seq == "M"
    ]
    list_motif_sizes = [locus.motif_size for locus in list_loci]
    results = compute_VAR(mutated.values(), population_names, list_motif_sizes)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 154.90377867746088 / 10
    assert pytest.approx(results["pop2"]) == 132.11953441295188 / 9


def test_compute_MGW_by_locus(header_text_te2_XY):
    """Vérifie _compute_MGW_by_locus sur toy_example2_ms_dna_xy.

    Locus_M_A_1_ : pop1 a des allèles présents à 203/197/193 (min=193,
    pas 197 -- 189/199/201 sont à compte 0), pop2 à 203/193/199/189
    (min=189, max=203)."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    result = _compute_MGW_by_locus(_length_by_population(mutated["Locus_M_A_1_"]), 2)

    assert result.keys() == {"pop1", "pop2"}
    assert result["pop1"] == (3, 1 + (203 - 193) / 2)
    assert result["pop2"] == (4, 1 + (203 - 189) / 2)


def test_compute_MGW(header_text_te2_XY):
    """Vérifie compute_MGW sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )
    list_loci = [
        locus
        for locus in parse_loci_description(header_text_te2_XY)
        if locus.ms_or_seq == "M"
    ]
    list_motif_sizes = [locus.motif_size for locus in list_loci]
    population_names = ["pop1", "pop2"]
    results = compute_MGW(mutated.values(), population_names, list_motif_sizes)

    assert results.keys() == {"pop1", "pop2"}
    assert pytest.approx(results["pop1"]) == 0.7398373983739838
    assert pytest.approx(results["pop2"]) == 0.7117117117117117


def test_compute_all_statistics_microsat(header_text_te2_XY):
    """Vérifie compute_all_statistics_microsat sur toy_example2_ms_dna."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    population_names = ["pop1", "pop2"]

    with pytest.raises(
        NotImplementedError,
        match="Les statistiques microsat ne sont pas encore implémentées",
    ):
        assert compute_all_statistics_microsat(
            header_text_te2_XY, mutated, population_names
        )


def test_combined_alleles_for_two_populations():
    """Vérifie que la fonction combine_alleles_for_two_populations fonctionne correctement."""
    alleles_pop1 = [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)]
    alleles_pop2 = [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)]

    combined_alleles = _combined_alleles_for_two_populations(alleles_pop1, alleles_pop2)

    assert combined_alleles == 5


def test_count_combined_alleles_for_one_locus():
    """Vérifie que la fonction _count_combined_alleles_for_one_locus fonctionne correctement."""
    length_by_population = {
        "pop1": [(201, 0), (203, 31), (197, 1), (193, 7), (199, 0), (189, 0)],
        "pop2": [(201, 0), (203, 18), (197, 0), (193, 20), (199, 1), (189, 1)],
    }

    combined_alleles_count = _count_combined_alleles_for_one_locus(
        length_by_population, ["pop1", "pop2"]
    )

    assert combined_alleles_count == {"1.2": 5.0}


def test_compute_N2P(header_text_te2_XY):
    """Vérifie compute_N2P sur toy_example2_ms_dna_xy."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2_XY, scenario_index=1, seed=42
    )
    mutated = microsat_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2_XY,
        mss_file_path=OBSERVED_MSS_FILE_TE2_XY,
        seed=42,
    )

    population_names = ["pop1", "pop2"]

    results = compute_N2P(mutated.values(), population_names)

    assert results == {"1.2": 10.8}
