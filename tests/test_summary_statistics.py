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
    OBSERVED_SNP_FILE_TE2,
)

from bridge.ancestry_simulation import dna_mutation_simulation_per_locus
from bridge.loci_parser import parse_loci_description
from bridge.pipeline import build_random_demography_for_scenario_index
from bridge.summary_statistics import (
    _genotype_matrix_by_population,
    _prepare_matrices_poolseq,
    compute_all_statistics_poolseq,
    mean_distinct_haplotypes_per_group,
    mean_distinct_haplotypes_per_group_pairwize,
    mean_minor_allele_count_per_group,
    mean_pairwise_differences_per_group,
    mean_pairwise_differences_per_group_pairwize,
    mean_private_segregating_sites_per_group,
    mean_segregating_sites_per_group,
    mean_segregating_sites_per_group_pairwize,
    mean_tajima_d_per_group,
    variance_minor_allele_count_per_group,
    variance_pairwise_differences_per_group,
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
        mss_file_path=OBSERVED_SNP_FILE_TE2,
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
    """Vérifie mean_segregating_sites_per_group séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre NSS_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_segregating_sites_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {"pop1": pytest.approx(5.8), "pop2": pytest.approx(5.6)}

    mean_g3 = mean_segregating_sites_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {"pop1": pytest.approx(5.6), "pop2": pytest.approx(5.6)}


def test_mean_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    assert mean_segregating_sites_per_group([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group(header_text_te2):
    """Vérifie mean_distinct_haplotypes_per_group séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre NH_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_distinct_haplotypes_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {"pop1": pytest.approx(5.2), "pop2": pytest.approx(5.4)}

    mean_g3 = mean_distinct_haplotypes_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {"pop1": pytest.approx(5.0), "pop2": pytest.approx(3.6)}


def test_mean_distinct_haplotypes_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = mean_distinct_haplotypes_per_group(
        tree_sequences_empty, population_names
    )
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_pairwise_differences_per_group(header_text_te2):
    """Vérifie mean_pairwise_differences_per_group séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre VPD_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_pairwise_differences_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {
        "pop1": pytest.approx(1.1997435897435897),
        "pop2": pytest.approx(1.2999999999999998),
    }

    mean_g3 = mean_pairwise_differences_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {
        "pop1": pytest.approx(1.2252631578947368),
        "pop2": pytest.approx(1.2978947368421054),
    }


def test_mean_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    mean_empty = mean_pairwise_differences_per_group(
        tree_sequences_empty, population_names
    )
    assert mean_empty == {"pop1": 0.0, "pop2": 0.0}


def test_variance_pairwise_differences_per_group(header_text_te2):
    """Vérifie variance_pairwise_differences_per_group séparément sur G2 (<A>) et
    G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes mélangés,
    puisque chaque `group Gx` du header calcule son propre VPD_i à partir
    de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    variance_g2 = variance_pairwise_differences_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert variance_g2 == {
        "pop1": pytest.approx(2.0666373720417366),
        "pop2": pytest.approx(1.7601968335472828),
    }

    variance_g3 = variance_pairwise_differences_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert variance_g3 == {
        "pop1": pytest.approx(1.5865441381230851),
        "pop2": pytest.approx(3.0353940406571986),
    }


def test_variance_pairwise_differences_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    # test si tree_sequence est vide
    tree_sequences_empty = []
    variance_empty = variance_pairwise_differences_per_group(
        tree_sequences_empty, population_names
    )
    assert variance_empty == {"pop1": 0.0, "pop2": 0.0}


def test_mean_tajima_d_per_group(header_text_te2):
    """Vérifie mean_tajima_d_per_group séparément sur G2 (<A>) et G3 (<M>)
    de toy_example2_ms_dna -- jamais les deux groupes mélangés, puisque
    chaque `group Gx` du header calcule son propre DTA_i à partir de ses
    seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    dta_g2 = mean_tajima_d_per_group(tree_sequences_for_group("G2"), population_names)
    assert dta_g2 == {
        "pop1": pytest.approx(-0.11361849459656947),
        "pop2": pytest.approx(-0.2621067146690289),
    }

    dta_g3 = mean_tajima_d_per_group(tree_sequences_for_group("G3"), population_names)
    assert dta_g3 == {
        "pop1": pytest.approx(-0.453016737083286),
        "pop2": pytest.approx(-0.10017554623967702),
    }


def test_mean_tajima_d_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0), comme le `res = 0.0` du C++ avant son `if (nl > 0)`."""
    population_names = ["pop1", "pop2"]
    assert mean_tajima_d_per_group([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_private_segregating_sites_per_group(header_text_te2):
    """Vérifie mean_private_segregating_sites_per_group séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    PSS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    pss_g2 = mean_private_segregating_sites_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert pss_g2 == {"pop1": pytest.approx(1.2), "pop2": pytest.approx(1.0)}

    pss_g3 = mean_private_segregating_sites_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert pss_g3 == {"pop1": pytest.approx(3.8), "pop2": pytest.approx(3.8)}


def test_mean_private_segregating_sites_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert mean_private_segregating_sites_per_group([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_minor_allele_count_per_group(header_text_te2):
    """Vérifie mean_minor_allele_count_per_group séparément sur G2 (<A>)
    et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    MNS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mns_g2 = mean_minor_allele_count_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert mns_g2 == {
        "pop1": pytest.approx(3.486274509803921),
        "pop2": pytest.approx(4.0),
    }

    mns_g3 = mean_minor_allele_count_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert mns_g3 == {
        "pop1": pytest.approx(1.575),
        "pop2": pytest.approx(3.0385964912280703),
    }


def test_mean_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert mean_minor_allele_count_per_group([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_variance_minor_allele_count_per_group(header_text_te2):
    """Vérifie variance_minor_allele_count_per_group séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    VNS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    vns_g2 = variance_minor_allele_count_per_group(
        tree_sequences_for_group("G2"), population_names
    )
    assert vns_g2 == {
        "pop1": pytest.approx(8.9039600153787),
        "pop2": pytest.approx(18.444444444444446),
    }

    vns_g3 = variance_minor_allele_count_per_group(
        tree_sequences_for_group("G3"), population_names
    )
    assert vns_g3 == {
        "pop1": pytest.approx(3.121875),
        "pop2": pytest.approx(5.515358571868267),
    }


def test_variance_minor_allele_count_per_group_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert variance_minor_allele_count_per_group([], population_names) == {
        "pop1": 0.0,
        "pop2": 0.0,
    }


def test_mean_distinct_haplotypes_per_group_pairwize(header_text_te2):
    """Vérifie mean_distinct_haplotypes_per_group_pairwize séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    NH_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_distinct_haplotypes_per_group_pairwize(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {"1.2": pytest.approx(6.8)}

    mean_g3 = mean_distinct_haplotypes_per_group_pairwize(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {"1.2": pytest.approx(6.6)}


def test_mean_distinct_haplotypes_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert mean_distinct_haplotypes_per_group_pairwize([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_segregating_sites_per_group_pairwize(header_text_te2):
    """Vérifie mean_segregating_sites_per_group_pairwize séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    NSS_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_segregating_sites_per_group_pairwize(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {"1.2": pytest.approx(6.8)}

    mean_g3 = mean_segregating_sites_per_group_pairwize(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {"1.2": pytest.approx(9.4)}


def test_mean_segregating_sites_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert mean_segregating_sites_per_group_pairwize([], population_names) == {
        "1.2": 0.0,
    }


def test_mean_pairwise_differences_per_group_pairwize(header_text_te2):
    """Vérifie mean_pairwise_differences_per_group_pairwize séparément sur G2
    (<A>) et G3 (<M>) de toy_example2_ms_dna -- jamais les deux groupes
    mélangés, puisque chaque `group Gx` du header calcule son propre
    VPD_i à partir de ses seuls loci."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text_te2, scenario_index=1, seed=42
    )
    mutated = dna_mutation_simulation_per_locus(
        demography=demography,
        header_text=header_text_te2,
        mss_file_path=OBSERVED_SNP_FILE_TE2,
        seed=42,
    )

    loci_description = parse_loci_description(header_text_te2)
    population_names = ["pop1", "pop2"]

    def tree_sequences_for_group(group_name):
        names = [locus.name for locus in loci_description if locus.group == group_name]
        return [mutated[name] for name in names]

    mean_g2 = mean_pairwise_differences_per_group_pairwize(
        tree_sequences_for_group("G2"), population_names
    )
    assert mean_g2 == {"1.2": pytest.approx(1.2498717948717948)}

    mean_g3 = mean_pairwise_differences_per_group_pairwize(
        tree_sequences_for_group("G3"), population_names
    )
    assert mean_g3 == {"1.2": pytest.approx(1.261578947368421)}


def test_mean_pairwise_differences_per_group_pairwize_empty_defaults_to_zero():
    """Une liste de loci vide ne doit pas faire disparaître de population
    du résultat ni lever d'exception -- chaque population attendue garde
    une valeur (0.0)."""
    population_names = ["pop1", "pop2"]
    assert mean_pairwise_differences_per_group_pairwize([], population_names) == {
        "1.2": 0.0,
    }
