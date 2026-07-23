"""Vérifie observed_data : comptage d'individus par population, mapping
indice-de-scénario -> nom réel de population, sex-ratio et sexe par
individu (ces deux derniers nécessaires pour les futurs loci <X>/<Y>/<M>,
voir notes/exploration.md)."""

import pytest
from conftest import (
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE3,
    OBSERVED_SNP_FILE_TE4,
    OBSERVED_SNP_FILE_TE5,
)

from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    coalescence_coefficient,
    count_samples_per_population,
    detect_snp_file_type,
    individual_sexes_per_population,
    observed_reads,
    parse_maf_ratio,
    parse_mrc_ratio,
    parse_sex_ratio,
    population_index_to_name,
)


def test_detect_snp_file_type():
    """Vérifie que le type de fichier .snp : INDSEQ (individus par ligne) ou POOLSEQ (pools par ligne)
    est correctement détecté pour les fichiers de référence human et toy_example4."""
    assert detect_snp_file_type(OBSERVED_SNP_FILE_HUMAN) == "IND"
    assert detect_snp_file_type(OBSERVED_SNP_FILE_TE4) == "POOL"


def test_count_samples_per_population():
    """Vérifie que le comptage retrouve bien les 4 populations à 30
    individus chacune, annoncées en commentaire dans le fichier."""
    counts_human = count_samples_per_population(OBSERVED_SNP_FILE_HUMAN)
    counts_te4 = count_samples_per_population(OBSERVED_SNP_FILE_TE4)

    assert set(counts_human.keys()) == {"ASW", "YRI", "CHB", "GBR"}
    assert all(n == 30 for n in counts_human.values())
    assert counts_te4 == {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}


def test_population_index_to_name():
    """Vérifie le mapping indice de scénario (1-indexed) -> nom réel de
    population, dans l'ordre d'apparition du fichier .snp."""
    mapping = population_index_to_name(OBSERVED_SNP_FILE_HUMAN)

    assert mapping == {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"}


def test_parse_sex_ratio():
    """Vérifie que le parsing du fichier snp renvoie bien la bonne proportion"""
    sex_ratio_human = parse_sex_ratio(OBSERVED_SNP_FILE_HUMAN)
    sex_ratio_te5 = parse_sex_ratio(OBSERVED_SNP_FILE_TE5)
    assert sex_ratio_human == 0.5
    assert sex_ratio_te5 == 0.428571 / (1 + 0.428571)


def test_parse_maf_ratio():
    """Vérifie que le parsing du fichier snp renvoie bien la bonne proportion"""
    maf_ratio_human = parse_maf_ratio(OBSERVED_SNP_FILE_HUMAN)
    maf_ratio_te5 = parse_maf_ratio(OBSERVED_SNP_FILE_TE5)
    maf_ratio_te3 = parse_maf_ratio(OBSERVED_SNP_FILE_TE3)
    assert maf_ratio_human == 0.0
    assert maf_ratio_te5 == 0.0
    assert maf_ratio_te3 == 0.05


def test_parse_mrc_ratio():
    """Vérifie que le parsing du fichier snp renvoie bien la bonne proportion"""
    mrc_ratio_human = parse_mrc_ratio(OBSERVED_SNP_FILE_HUMAN)
    mrc_ratio_te5 = parse_mrc_ratio(OBSERVED_SNP_FILE_TE5)
    mrc_ratio_te3 = parse_mrc_ratio(OBSERVED_SNP_FILE_TE3)
    mrc_ratio_te4 = parse_mrc_ratio(OBSERVED_SNP_FILE_TE4)
    assert mrc_ratio_human == 1
    assert mrc_ratio_te5 == 1
    assert mrc_ratio_te3 == 1
    assert mrc_ratio_te4 == 5


def test_observed_reads(header_text_te4):
    """Vérifie que le parsing du fichier snp renvoie bien le nombre de reads
    observés par population, pour les fichiers POOLSEQ toy_example4 et
    toy_example5, APRÈS purge des loci sous le seuil MRC (<MRC=5> ici) --
    reproduit DataC::purgelocMRCPOOLSEQ, qui élimine ces loci au
    chargement, avant toute utilisation. Le tout premier locus BRUT du
    fichier ({"POP1": (0, 93), "POP2": (0, 100), "POP3": (1, 116),
    "POP4": (0, 139)}) a min(somme reads1, somme reads2) = min(1, 447) =
    1 < 5 -- il doit donc être absent du résultat ; le premier élément
    retourné est le premier locus du fichier qui passe réellement le
    seuil."""
    max_loci = parse_loci_description(header_text_te4).total_loci["A"]
    observed_reads_te4 = observed_reads(OBSERVED_SNP_FILE_TE4, num_loci=max_loci)
    observed_reads_te4_total = observed_reads(OBSERVED_SNP_FILE_TE4, num_loci=None)
    assert len(observed_reads_te4) == 100
    assert (
        len(observed_reads_te4_total) == 14388
    )  # 10x plus de loci bruts dans le fichier
    assert observed_reads_te4[0] == {
        "POP1": (52, 112),
        "POP2": (29, 95),
        "POP3": (45, 106),
        "POP4": (18, 93),
    }
    with pytest.raises(ValueError, match="format POOLSEQ"):
        observed_reads(OBSERVED_SNP_FILE_HUMAN)  # fichier INDSEQ


def test_individual_sexes_per_population():
    """Sur human (dataset <A>-only), le sexe n'est jamais renseigné : les
    120 individus doivent tous ressortir en "9" (sexe inconnu, cf.
    data.cpp:702-704). Sur toy_example5 (qui a des loci <X>/<Y>/<M>), le
    sexe est réellement renseigné : on doit retrouver 10 M et 10 F par
    population, cohérent avec les 20 individus par population comptés par
    test_count_samples_per_population."""
    sexes_human = individual_sexes_per_population(OBSERVED_SNP_FILE_HUMAN)

    assert set(sexes_human.keys()) == {"ASW", "YRI", "CHB", "GBR"}
    assert all(sexes == ["9"] * 30 for sexes in sexes_human.values())

    sexes_te5 = individual_sexes_per_population(OBSERVED_SNP_FILE_TE5)

    assert set(sexes_te5.keys()) == {"P1", "P2", "P3"}
    for sexes in sexes_te5.values():
        assert len(sexes) == 20
        assert sexes.count("M") == 10
        assert sexes.count("F") == 10


def test_coalescence_coefficient():
    """Vérifie que le coefficient de coalescence est bien calculé pour
    human (dataset <A>-only) et toy_example5 (qui a des loci <X>/<Y>/<M>).
    """
    sexe_ratio_human = parse_sex_ratio(OBSERVED_SNP_FILE_HUMAN)
    assert coalescence_coefficient("A", sexe_ratio_human) == 16 * sexe_ratio_human * (
        1 - sexe_ratio_human
    )

    sexe_ratio_te5 = parse_sex_ratio(OBSERVED_SNP_FILE_TE5)
    assert coalescence_coefficient("X", sexe_ratio_te5) == 18 * sexe_ratio_te5 * (
        1 - sexe_ratio_te5
    ) / (1 + sexe_ratio_te5)

    with pytest.raises(
        NotImplementedError,
        match="Type de locus inconnu pour le calcul du coefficient de coalescence : 'Z'",
    ):
        coalescence_coefficient("Z", sexe_ratio_te5)  # type de locus inconnu
