"""Vérifie pipeline : orchestration de bout en bout (header.txt -> Demography,
point d'entrée -p ./, calcul des statistiques résumées avec filtrage
ALL/HEADER)."""

from pathlib import Path

import msprime
import pytest
from conftest import (
    GENERAL_BINARY_PATH,
    OBSERVED_SNP_FILE_HUMAN,
)

from bridge.header_dataclasses import SnpReplayContext
from bridge.loci_parser import parse_loci_description
from bridge.pipeline import (
    _draw_common_data,
    _simulate_and_compute_statistics,
    _simulate_and_compute_statistics_from_values,
    build_random_demography_for_scenario_index,
    compute_summary_statistics_dna,
    compute_summary_statistics_microsat,
    compute_summary_statistics_mixed,
    compute_summary_statistics_mixed_from_values,
    compute_summary_statistics_snp,
    compute_summary_statistics_snp_from_values,
    read_header_text,
    simulate_particle_genotypes,
    simulate_particle_genotypes_from_values,
)
from bridge.reftable_loop import _group_prior_columns


def test_pipeline_scenario1(header_text):
    """Vérifie que le pipeline complet (header.txt -> Demography) fonctionne de
    bout en bout sur le scénario 1, et que la démographie produite a la
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


def test_simulate_particle_genotypes(snp_context_human):
    """Vérifie le point d'entrée de haut niveau : à partir d'un simple chemin
    de dossier (comme le -p ./ de DIYABC), tout le pipeline doit fonctionner
    sans qu'on ait à lire manuellement header.txt ou le fichier .snp nous-
    mêmes."""
    mutated, values = simulate_particle_genotypes(
        snp_context_human,
        scenario_index=1,
        num_loci=15,
        seed=42,
    )

    mutated_list = list(mutated)
    assert len(mutated_list) == 15
    assert "N1" in values


def test_simulate_particle_genotypes_multi_type(snp_context_te5):
    """Vérifie que simulate_particle_genotypes boucle bien sur TOUS les types
    de locus déclarés dans 'loci description', pas seulement <A> --
    toy_example5 (contrairement à human, <A>-only) déclare 4 types (A/X/Y/M,
    voir reference/toy_example5/headerRF.txt) : num_loci est un.

    compte PAR TYPE (voir
    pipeline._simulate_genotypes_for_all_locus_types), donc on attend
    num_loci * 4 génotypes au total, pas juste num_loci.
    """
    mutated, values = simulate_particle_genotypes(
        snp_context_te5,
        scenario_index=1,
        num_loci=3,
        seed=42,
    )

    mutated_list = list(mutated)
    assert len(mutated_list) == 3 * 4  # 3 loci x 4 types déclarés (A/X/Y/M)
    assert "N1" in values


def test_compute_summary_statistics_snp_multi_type(snp_context_te5):
    """Vérifie que compute_summary_statistics_snp (donc compute_all_statistics_indseq)
    fonctionne aussi sur un dataset multi-type <A>/<X>/<Y>/<M>, pas seulement
    <A> -- 51 statistiques attendues (vs 130 pour human) car toy_example5 n'a
    que 3 populations, pas 4 (moins de paires/triplets)."""
    summary_stats, values = compute_summary_statistics_snp(
        context=snp_context_te5,
        scenario_index=1,
        num_loci=3,
        seed=42,
        stats_filter="ALL",
    )
    assert len(summary_stats) == 51
    assert "N1" in values
    assert not any(v != v for v in summary_stats.values())  # v != v <=> NaN


def test_compute_summary_statistics_poolseq_varies_with_seed(snp_context_te4):
    """Vérifie que compute_summary_statistics_snp simule bien pour PoolSeq (branche
    else de la fonction) au lieu de recopier telles quelles les statistiques de
    l'observé -- régression du bug du 2026-07-23 où l'appel à
    simulate_poolseq_reads_with_mrc_filter avait été supprimé par erreur en
    câblant observed_reads_per_locus, ce qui aurait rendu toutes les particules
    PoolSeq d'un reftable identiques entre elles.

    Deux graines différentes doivent donc tirer des paramètres
    différents ET produire des statistiques différentes.
    """
    stats_seed_1, values_1 = compute_summary_statistics_snp(
        context=snp_context_te4,
        scenario_index=1,
        seed=1,
    )
    stats_seed_2, values_2 = compute_summary_statistics_snp(
        context=snp_context_te4,
        scenario_index=1,
        seed=2,
    )

    assert values_1 != values_2
    assert stats_seed_1 != stats_seed_2


def test_compute_summary_statistics_from_values_poolseq_varies_with_values(
    snp_context_te4,
):
    """Même régression que test_compute_summary_statistics_poolseq_varies_
    with_seed, mais côté compute_summary_statistics_snp_from_values (l'autre
    fonction touchée par le bug du 2026-07-23) : deux jeux de paramètres
    différents (même seed) doivent produire des statistiques différentes."""
    _, values_1 = build_random_demography_for_scenario_index(
        snp_context_te4.header_text, scenario_index=1, seed=1
    )
    _, values_2 = build_random_demography_for_scenario_index(
        snp_context_te4.header_text, scenario_index=1, seed=2
    )
    assert values_1 != values_2  # sinon le test ne prouve rien

    stats_1 = compute_summary_statistics_snp_from_values(
        context=snp_context_te4,
        scenario_index=1,
        values=values_1,
        seed=42,
    )
    stats_2 = compute_summary_statistics_snp_from_values(
        context=snp_context_te4,
        scenario_index=1,
        values=values_2,
        seed=42,
    )

    assert stats_1 != stats_2


def test_read_header_text_prefers_header_txt(tmp_path):
    """Si les deux fichiers sont présents, header.txt doit être lu en priorité
    (config initiale fournie par l'utilisateur), pas headerRF.txt (variante
    produite par un run DIYABC réel)."""
    (tmp_path / "header.txt").write_text("contenu header.txt")
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu header.txt"


def test_read_header_text_falls_back_to_headerRF(tmp_path):
    """Si seul headerRF.txt est présent (ex: reference/Exemple5/), il doit être
    lu en repli."""
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu headerRF.txt"


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie -- "
    "ce test nécessite le binaire 'general' compilé de DIYABC.",
)
def test_compute_summary_statistics_snp_scenario1(tmp_path, snp_context_human):
    """Vérifie que compute_summary_statistics_snp produit bien les 112 statistiques
    résumées attendues (filtre ALL), en déléguant le calcul au vrai binaire C++
    sur des données simulées par notre pipeline."""
    summary_statistics, values = compute_summary_statistics_snp(
        context=snp_context_human,
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


def _replace_group_summary_statistics_section(
    header_text: str, new_section_lines: list[str]
) -> str:
    """Remplace la section 'group summary statistics' de header_text par
    new_section_lines -- même approche que loci_parser.rewrite_loci_count
    pour tester avec un contenu différent sans maintenir un header.txt
    séparé à la main."""
    lines = header_text.splitlines()
    start = next(
        i
        for i, line in enumerate(lines)
        if line.strip().startswith("group summary statistics")
    )
    end = start + 1
    while end < len(lines) and lines[end].strip():
        end += 1
    return "\n".join(lines[:start] + new_section_lines + lines[end:])


def test_compute_summary_statistics_snp_stats_filter_header(
    tmp_path, snp_context_human
):
    """stats_filter='HEADER' ne garde, dans l'ordre de déclaration, que les
    statistiques listées dans 'group summary statistics' -- remplace la section
    obsolète de human/header.txt par un petit sous-ensemble au vocabulaire
    moderne, pour vérifier le filtrage sans dépendre d'un dataset externe."""
    modified_header_text = _replace_group_summary_statistics_section(
        snp_context_human.header_text,
        ["group summary statistics (4)", "group G1 (4)", "ML1p 1 2", "HWm 1 2"],
    )
    (tmp_path / "header.txt").write_text(modified_header_text)
    (tmp_path / OBSERVED_SNP_FILE_HUMAN.name).symlink_to(OBSERVED_SNP_FILE_HUMAN)

    snp_context_human_modified = SnpReplayContext(
        header_text=modified_header_text,
        snp_path=snp_context_human.snp_path,
        snp_file_type=snp_context_human.snp_file_type,
        loci_description=snp_context_human.loci_description,
        counts_per_sample=snp_context_human.counts_per_sample,
        sex_ratio=snp_context_human.sex_ratio,
        maf_ratio=snp_context_human.maf_ratio,
        mrc_ratio=snp_context_human.mrc_ratio,
        reads_observed=snp_context_human.reads_observed,
        sexes_per_sample=snp_context_human.sexes_per_sample,
    )

    summary_stats, values = compute_summary_statistics_snp(
        context=snp_context_human_modified,
        scenario_index=1,
        num_loci=10,
        seed=1,
        stats_filter="HEADER",
    )

    assert list(summary_stats.keys()) == ["ML1p_1", "ML1p_2", "HWm_1", "HWm_2"]
    assert "N1" in values


def test_compute_summary_statistics_stats_filter_header_raises_on_unknown_names(
    snp_context_human,
):
    """stats_filter='HEADER' sur le vrai human/header.txt (vocabulaire obsolète
    HP0/HM1/...) doit lever une ValueError explicite plutôt que de produire
    silencieusement un reftable vide ou incomplet."""
    with pytest.raises(ValueError, match="non calculées"):
        compute_summary_statistics_snp(
            context=snp_context_human,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="HEADER",
        )


def test_compute_summary_statistics_unknown_stats_filter_raises(snp_context_human):
    with pytest.raises(NotImplementedError, match="stats_filter"):
        compute_summary_statistics_snp(
            context=snp_context_human,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="BOGUS",
        )


# -------------------------------------------------------------
# Rejeu d'un jeu IndSeq SÉRIEL (1 population, plusieurs dates)
# -------------------------------------------------------------

_SERIAL_INDSEQ_HEADER = """toy.snp
1 parameters and 0 summary statistics

2 scenarios: 6 5
scenario 1 [0.5] (1)
Npast
0 sample 1
50 sample 1
200 sample 1
500 sample 1
scenario 2 [0.5] (1)
Npast
0 sample 1
50 sample 1
200 sample 1
500 sample 1

historical parameters priors (1,0)
Npast N UN[10,50000,0,0]

loci description (1)
5 <A> G1 from 1

group summary statistics (0)
"""


@pytest.mark.parametrize("maf_ratio", [0.0, 0.05])
def test_simulate_particle_genotypes_from_values_serial_indseq(maf_ratio):
    """Le rejeu d'un jeu SÉRIEL (une population échantillonnée à 4 dates, comme
    `human_seriel`) doit construire ses `SampleSet` depuis le scénario, comme le
    fait le jumeau de tirage `simulate_particle_genotypes`.

    Régression : cette construction avait disparu de la branche de rejeu au
    commit 392d193 (25/09) et rien ne l'a signalé pendant dix jours, faute de
    test rejouant un jeu IndSeq sériel. Sans elle, msprime reçoit les noms
    pop1..pop4 des blocs POP du .snp alors que la démographie n'a qu'une
    population, et lève `KeyError: Population with name 'pop2' not found`.
    Les deux valeurs de MAF couvrent les deux chemins (`hudson` et filtré).
    """
    context = SnpReplayContext(
        header_text=_SERIAL_INDSEQ_HEADER,
        snp_path=Path("toy.snp"),
        snp_file_type="IND",
        loci_description=parse_loci_description(_SERIAL_INDSEQ_HEADER),
        counts_per_sample={"POP1": 6, "POP2": 4, "POP3": 5, "POP4": 3},
        sex_ratio=0.5,
        maf_ratio=maf_ratio,
        mrc_ratio=None,
        reads_observed=None,
        sexes_per_sample={},
    )

    genotypes = list(
        simulate_particle_genotypes_from_values(context, 2, {"Npast": 1000.0}, seed=7)
    )

    assert len(genotypes) == 5
    # Un échantillon par date, chacun avec 2 copies de gènes par individu.
    assert {name: len(g) for name, g in genotypes[0].items()} == {
        "pop1": 12,
        "pop2": 8,
        "pop3": 10,
        "pop4": 6,
    }


# -------------------------------------------------------------
# Tests pour la partie DNA
# -------------------------------------------------------------


def test_compute_summary_statistics_dna(dna_context_te2):
    stats, _, _ = compute_summary_statistics_dna(
        context=dna_context_te2,
        scenario_index=1,
        seed=42,
    )

    assert len(stats) == 42
    assert stats["NSS_2_1"] == pytest.approx(9.8)
    assert stats["HST_2_1.2"] == pytest.approx(0.012362823348065015)
    assert stats["NH2_3_1.2"] == pytest.approx(5.0)


# -------------------------------------------------------------
# Tests pour la partie Microsat
# -------------------------------------------------------------


def test_compute_summary_statistics_microsat(microsat_context_te2_xy):

    summary_stats, _, _ = compute_summary_statistics_microsat(
        context=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    assert len(summary_stats) == 16
    assert summary_stats["FST_1_1.2"] == pytest.approx(0.007081532724996798)
    assert summary_stats["LIK_1_1.2"] == pytest.approx(1.8115736105446298)
    assert summary_stats["LIK_1_2.1"] == pytest.approx(1.7855891126163377)


# -------------------------------------------------------------
# Tests pour la partie Microsat + DNA
# -------------------------------------------------------------


def test_draw_common_data(microsat_context_te2_xy):
    """Vérifie que _draw_common_data extrait correctement les données
    communes nécessaires à la simulation et au calcul des statistiques pour
    les scénarios ADN et microsat."""
    demography, sample_sets, values = _draw_common_data(
        context=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    assert isinstance(demography, msprime.Demography)
    assert isinstance(sample_sets, list)
    assert all(isinstance(s, msprime.SampleSet) for s in sample_sets)
    assert "N1" in values


def test_simulate_and_compute_statistics(microsat_context_te2_xy):
    """Vérifie que _simulate_and_compute_statistics calcule correctement les statistiques
    résumées pour les scénarios ADN et microsat."""
    demography, sample_sets, values = _draw_common_data(
        context=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    summary_stats, _ = _simulate_and_compute_statistics(
        context=microsat_context_te2_xy,
        demography=demography,
        sample_sets=sample_sets,
        type_of_data="microsat",
        seed=42,
    )

    assert isinstance(summary_stats, dict)
    assert "FST_1_1.2" in summary_stats
    assert "LIK_1_1.2" in summary_stats
    assert "LIK_1_2.1" in summary_stats


def test_simulate_and_compute_statistics_with_wrong_type(microsat_context_te2_xy):
    """Vérifie que _simulate_and_compute_statistics lève une ValueError si le type de données
    est incorrect."""
    demography, sample_sets, values = _draw_common_data(
        context=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    with pytest.raises(ValueError, match="Type de données"):
        _simulate_and_compute_statistics(
            context=microsat_context_te2_xy,
            demography=demography,
            sample_sets=sample_sets,
            type_of_data="invalid_type",
            seed=42,
        )


def test_compute_summary_statistics_mixed(microsat_context_te2_xy, dna_context_te2_xy):
    """Vérifie que le calcul des statistiques résumées mixtes (ADN + microsat) fonctionne correctement et que les résultats sont cohérents avec les calculs séparés pour chaque type de données."""
    summary_stats_microsat, _, _ = compute_summary_statistics_microsat(
        context=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    summary_stats_dna, _, _ = compute_summary_statistics_dna(
        context=dna_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    summary_stats_mixed, _, _ = compute_summary_statistics_mixed(
        context_dna=dna_context_te2_xy,
        context_microsat=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
    )

    assert summary_stats_mixed == {**summary_stats_dna, **summary_stats_microsat}


def test_compute_summary_statistics_mixed_with_stats_filter_Header(
    microsat_context_te2_xy,
    dna_context_te2_xy,
):
    """Vérifie que le statistics filter 'HEADER' fonctionne correctement pour le calcul des statistiques résumées mixtes (ADN + microsat)."""
    summary_stats_mixed, _, _ = compute_summary_statistics_mixed(
        context_dna=dna_context_te2_xy,
        context_microsat=microsat_context_te2_xy,
        scenario_index=1,
        seed=42,
        stats_filter="HEADER",
    )

    assert "HST_3_1.2" in summary_stats_mixed
    assert "LIK_1_1.2" in summary_stats_mixed
    assert len(summary_stats_mixed) == 58


def test_compute_summary_statistics_mixed_values_with_different_header(
    microsat_context_te2_xy,
    dna_context_te2,
):
    """Vérifie que l'erreur est bien levé par compute_summary_statistics_mixed si les deux contextes n'ont pas le même header.txt"""

    with pytest.raises(ValueError, match="doivent avoir le même header.txt"):
        compute_summary_statistics_mixed(
            context_dna=dna_context_te2,
            context_microsat=microsat_context_te2_xy,
            scenario_index=1,
            seed=42,
        )


def test_compute_summary_statistics_mixed_from_values(
    microsat_context_te2_xy, dna_context_te2_xy
):
    """Vérifie que compute_summary_statistics_mixed_from_values fonctionne correctement et que les résultats sont cohérents avec les calculs séparés pour chaque type de données."""

    demography, sample_sets, values = _draw_common_data(
        microsat_context_te2_xy, 1, seed=42
    )
    _, nested = _simulate_and_compute_statistics(
        microsat_context_te2_xy,
        demography,
        sample_sets,
        type_of_data="microsat",
        seed=42,
    )
    group_priors_values = {
        c: nested[g][p]
        for c, g, p in _group_prior_columns(microsat_context_te2_xy.header_text)
    }

    summary_stats_microsat = _simulate_and_compute_statistics_from_values(
        context=microsat_context_te2_xy,
        demography=demography,
        sample_sets=sample_sets,
        type_of_data="microsat",
        group_priors_values=group_priors_values,
        seed=42,
    )

    summary_stats_dna = _simulate_and_compute_statistics_from_values(
        context=dna_context_te2_xy,
        demography=demography,
        sample_sets=sample_sets,
        type_of_data="dna",
        group_priors_values=group_priors_values,
        seed=42,
    )

    summary_stats_mixed = compute_summary_statistics_mixed_from_values(
        context_dna=dna_context_te2_xy,
        context_microsat=microsat_context_te2_xy,
        scenario_index=1,
        values=values,
        group_priors_values=group_priors_values,
        seed=42,
    )

    assert summary_stats_mixed == {**summary_stats_dna, **summary_stats_microsat}
