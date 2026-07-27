"""Vérifie ancestry_simulation : construction de l'argument samples pour
msprime.sim_ancestry, simulation de généalogies indépendantes par locus,
et mutation Hudson (exactement une mutation par locus, toujours
polymorphe)."""

import pytest
from conftest import (
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE4,
    OBSERVED_SNP_FILE_TE5,
)

from bridge.ancestry_simulation import (
    _reindex_reads_by_msprime_name,
    build_male_only_samples_argument,
    build_samples_argument,
    build_sex_stratified_samples_argument,
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
