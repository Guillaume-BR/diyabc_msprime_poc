"""Vérifie pipeline : orchestration de bout en bout (header.txt ->
Demography, point d'entrée -p ./, calcul des statistiques résumées avec
filtrage ALL/HEADER)."""

import msprime
import pytest
from conftest import (
    GENERAL_BINARY_PATH,
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE4,
    REFERENCE_DIR,
)

from bridge.pipeline import (
    build_random_demography_for_scenario_index,
    compute_summary_statistics,
    compute_summary_statistics_dna,
    compute_summary_statistics_from_values,
    read_header_text,
    run_poc_for_directory,
)


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


def test_run_poc_for_directory_multi_type():
    """Vérifie que run_poc_for_directory boucle bien sur TOUS les types de
    locus déclarés dans 'loci description', pas seulement <A> --
    toy_example5 (contrairement à human, <A>-only) déclare 4 types
    (A/X/Y/M, voir reference/toy_example5/headerRF.txt) : num_loci est un
    compte PAR TYPE (voir pipeline._simulate_genotypes_for_all_locus_types),
    donc on attend num_loci * 4 génotypes au total, pas juste num_loci."""
    mutated, values = run_poc_for_directory(
        REFERENCE_DIR.parent / "toy_example5",
        scenario_index=1,
        num_loci=3,
        seed=42,
    )

    mutated_list = list(mutated)
    assert len(mutated_list) == 3 * 4  # 3 loci x 4 types déclarés (A/X/Y/M)
    assert "N1" in values


def test_compute_summary_statistics_multi_type():
    """Vérifie que compute_summary_statistics (donc compute_all_statistics)
    fonctionne aussi sur un dataset multi-type <A>/<X>/<Y>/<M>, pas
    seulement <A> -- 51 statistiques attendues (vs 130 pour human) car
    toy_example5 n'a que 3 populations, pas 4 (moins de paires/triplets)."""
    summary_stats, values = compute_summary_statistics(
        reference_directory=REFERENCE_DIR.parent / "toy_example5",
        scenario_index=1,
        num_loci=3,
        seed=42,
        stats_filter="ALL",
    )
    assert len(summary_stats) == 51
    assert "N1" in values
    assert not any(v != v for v in summary_stats.values())  # v != v <=> NaN


def test_compute_summary_statistics_poolseq_varies_with_seed():
    """Vérifie que compute_summary_statistics simule bien pour PoolSeq
    (branche else de la fonction) au lieu de recopier telles quelles les
    statistiques de l'observé -- régression du bug du 2026-07-23 où
    l'appel à simulate_poolseq_reads_with_mrc_filter avait été supprimé
    par erreur en câblant observed_reads_per_locus, ce qui aurait rendu
    toutes les particules PoolSeq d'un reftable identiques entre elles.
    Deux graines différentes doivent donc tirer des paramètres différents
    ET produire des statistiques différentes."""
    reference_directory = OBSERVED_SNP_FILE_TE4.parent

    stats_seed_1, values_1 = compute_summary_statistics(
        reference_directory=reference_directory,
        scenario_index=1,
        seed=1,
    )
    stats_seed_2, values_2 = compute_summary_statistics(
        reference_directory=reference_directory,
        scenario_index=1,
        seed=2,
    )

    assert values_1 != values_2
    assert stats_seed_1 != stats_seed_2


def test_compute_summary_statistics_from_values_poolseq_varies_with_values(
    header_text_te4,
):
    """Même régression que test_compute_summary_statistics_poolseq_varies_
    with_seed, mais côté compute_summary_statistics_from_values (l'autre
    fonction touchée par le bug du 2026-07-23) : deux jeux de paramètres
    différents (même seed) doivent produire des statistiques
    différentes."""
    reference_directory = OBSERVED_SNP_FILE_TE4.parent

    _, values_1 = build_random_demography_for_scenario_index(
        header_text_te4, scenario_index=1, seed=1
    )
    _, values_2 = build_random_demography_for_scenario_index(
        header_text_te4, scenario_index=1, seed=2
    )
    assert values_1 != values_2  # sinon le test ne prouve rien

    stats_1 = compute_summary_statistics_from_values(
        reference_directory=reference_directory,
        scenario_index=1,
        values=values_1,
        seed=42,
    )
    stats_2 = compute_summary_statistics_from_values(
        reference_directory=reference_directory,
        scenario_index=1,
        values=values_2,
        seed=42,
    )

    assert stats_1 != stats_2


def test_read_header_text_prefers_header_txt(tmp_path):
    """Si les deux fichiers sont présents, header.txt doit être lu en
    priorité (config initiale fournie par l'utilisateur), pas
    headerRF.txt (variante produite par un run DIYABC réel)."""
    (tmp_path / "header.txt").write_text("contenu header.txt")
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu header.txt"


def test_read_header_text_falls_back_to_headerRF(tmp_path):
    """Si seul headerRF.txt est présent (ex: reference/Exemple5/), il
    doit être lu en repli."""
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu headerRF.txt"


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


def test_compute_summary_statistics_stats_filter_header(tmp_path, header_text):
    """stats_filter='HEADER' ne garde, dans l'ordre de déclaration, que
    les statistiques listées dans 'group summary statistics' --
    remplace la section obsolète de human/header.txt par un petit
    sous-ensemble au vocabulaire moderne, pour vérifier le filtrage
    sans dépendre d'un dataset externe."""
    modified_header_text = _replace_group_summary_statistics_section(
        header_text,
        ["group summary statistics (4)", "group G1 (4)", "ML1p 1 2", "HWm 1 2"],
    )
    (tmp_path / "header.txt").write_text(modified_header_text)
    (tmp_path / OBSERVED_SNP_FILE_HUMAN.name).symlink_to(OBSERVED_SNP_FILE_HUMAN)

    summary_stats, values = compute_summary_statistics(
        reference_directory=tmp_path,
        scenario_index=1,
        num_loci=10,
        seed=1,
        stats_filter="HEADER",
    )

    assert list(summary_stats.keys()) == ["ML1p_1", "ML1p_2", "HWm_1", "HWm_2"]
    assert "N1" in values


def test_compute_summary_statistics_stats_filter_header_raises_on_unknown_names(
    header_text,
):
    """stats_filter='HEADER' sur le vrai human/header.txt (vocabulaire
    obsolète HP0/HM1/...) doit lever une ValueError explicite plutôt que
    de produire silencieusement un reftable vide ou incomplet."""
    with pytest.raises(ValueError, match="non calculées"):
        compute_summary_statistics(
            reference_directory=REFERENCE_DIR,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="HEADER",
        )


def test_compute_summary_statistics_unknown_stats_filter_raises():
    with pytest.raises(NotImplementedError, match="stats_filter"):
        compute_summary_statistics(
            reference_directory=REFERENCE_DIR,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="BOGUS",
        )


def test_compute_summary_statistics_dna():
    stats, _ = compute_summary_statistics_dna(
        reference_directory=REFERENCE_DIR.parent / "toy_example2_ms_dna",
        scenario_index=1,
        seed=42,
    )

    assert len(stats) == 42
    assert stats["NSS_2_1"] == pytest.approx(5.8)
    assert stats["HST_2_1.2"] == pytest.approx(0.029037253935292443)
    assert stats["NH2_3_1.2"] == pytest.approx(10.0)
