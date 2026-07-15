"""Vérifie ancestry_simulation : construction de l'argument samples pour
msprime.sim_ancestry, simulation de généalogies indépendantes par locus,
et mutation Hudson (exactement une mutation par locus, toujours
polymorphe)."""

import pytest
from conftest import OBSERVED_SNP_FILE, OBSERVED_SNP_FILE_TE5

from bridge.ancestry_simulation import (
    build_male_only_samples_argument,
    build_samples_argument,
    build_sex_stratified_samples_argument,
    simulate_genotypes_for_locus_type,
    simulate_independent_loci,
    simulate_shared_ancestry_loci,
    simulate_snp_genotypes,
)
from bridge.pipeline import build_random_demography_for_scenario_index


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
    samples = build_samples_argument(OBSERVED_SNP_FILE)

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
    samples = build_samples_argument(OBSERVED_SNP_FILE)

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
        build_sex_stratified_samples_argument(OBSERVED_SNP_FILE)  # sexe non renseigné

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
        build_male_only_samples_argument(OBSERVED_SNP_FILE)  # sexe non renseigné

    samples = build_male_only_samples_argument(OBSERVED_SNP_FILE_TE5)
    assert samples == {"pop1": 10, "pop2": 10, "pop3": 10}


def test_simulate_shared_ancestry_loci(header_text):
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)
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
