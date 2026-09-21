"""
Fixtures et constantes partagées par tous les modules de tests, qui
portent tous sur le même dataset de référence (human/header.txt +
human_snp_all22chr_maf5.snp -- scénario 1 décortiqué à la main avec le
mentor, voir notes/exploration.md) plus le dataset toy_example5 (qui a de
vrais loci <X>/<Y>/<M> et un sex-ratio non trivial, pour tester
observed_data.py).
"""

import os
from pathlib import Path

import pytest

from bridge.header_dataclasses import (
    DnaReplayContext,
    MicrosatReplayContext,
    SnpReplayContext,
)
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    allele_bounds_per_locus,
    base_frequency_by_locus,
    count_samples_per_population,
    detect_snp_file_type,
    individual_sexes_per_population,
    observed_count_population,
    observed_microsatellites,
    observed_reads,
    observed_sequences,
    parse_maf_ratio,
    parse_mrc_ratio,
    parse_sex_ratio,
)

REFERENCE_DIR = Path(__file__).parent.parent / "reference"
GENERAL_BINARY_PATH = os.environ.get("DIYABC_GENERAL_PATH")

OBSERVED_SNP_FILE_HUMAN = REFERENCE_DIR / "human" / "human_snp_all22chr_maf5.snp"

OBSERVED_MSS_FILE_TE1 = (
    REFERENCE_DIR
    / "toy_example1_ms"
    / "pseudo_observed_DATASET_toy_example1_microsatellites_one_pop_multiple_samples_over_time_001.mss"
)

OBSERVED_MSS_FILE_TE2 = (
    REFERENCE_DIR
    / "toy_example2_ms_dna"
    / "pseudo_observed_DATASET_toy_example2_microsatellites_DNAsequences_ancestral_admixture_unsampled_pops_001.mss"
)

OBSERVED_MSS_FILE_TE2_XY = (
    REFERENCE_DIR
    / "toy_example2_ms_dna_XY"
    / "pseudo_observed_DATASET_toy_example2_microsatellites_DNAsequences_ancestral_admixture_unsampled_pops_001.mss"
)

OBSERVED_SNP_FILE_TE3 = (
    REFERENCE_DIR
    / "toy_example3"
    / "pseudo_observed_DATASET_SNP_INDSEQ_4pops_Scenario3_MER.snp"
)

OBSERVED_SNP_FILE_TE3_SCENARIO1 = (
    REFERENCE_DIR
    / "toy_example3_scenario1"
    / "pseudo_observed_DATASET_SNP_INDSEQ_4pops_Scenario3_MER.snp"
)

OBSERVED_SNP_FILE_TE4 = (
    REFERENCE_DIR
    / "toy_example4"
    / "pseudo_observed_DATASET_SNP_POOLSEQ_4pops_Scenario3_MER.snp"
)

OBSERVED_SNP_FILE_TE5 = (
    REFERENCE_DIR / "toy_example5" / "simu_dataset_test_divergence_admixture_001.snp"
)


@pytest.fixture
def header_text() -> str:
    return (REFERENCE_DIR / "human" / "header.txt").read_text()


@pytest.fixture
def header_text_te5() -> str:
    return (OBSERVED_SNP_FILE_TE5.parent / "headerRF.txt").read_text()


@pytest.fixture
def header_text_te4() -> str:
    return (OBSERVED_SNP_FILE_TE4.parent / "headerRF.txt").read_text()


@pytest.fixture
def header_text_te3_scenario1() -> str:
    """toy_example3, scénario 1 isolé (split+admixture, 8 priors) --
    seul dataset de ce projet avec un vrai filtre MAF actif (<MAF=0.05>,
    contrairement à human/toy_example5 qui sont <MAF=hudson>)."""
    return (OBSERVED_SNP_FILE_TE3_SCENARIO1.parent / "headerRF.txt").read_text()


@pytest.fixture
def header_text_te1() -> str:
    path_te1 = REFERENCE_DIR / "toy_example1_ms" / "headerRF.txt"
    return path_te1.read_text()


@pytest.fixture
def header_text_te1_modified() -> str:
    path_te1 = REFERENCE_DIR / "toy_example1_ms_modified" / "headerRF.txt"
    return path_te1.read_text()


@pytest.fixture
def header_text_te2() -> str:
    path_te2 = REFERENCE_DIR / "toy_example2_ms_dna" / "headerRF.txt"
    return path_te2.read_text()


@pytest.fixture
def header_text_te2_XY() -> str:
    path_te2_XY = REFERENCE_DIR / "toy_example2_ms_dna_XY" / "headerRF.txt"
    return path_te2_XY.read_text()


@pytest.fixture
def snp_context_human(header_text) -> SnpReplayContext:
    snp_path = OBSERVED_SNP_FILE_HUMAN
    snp_file_type = detect_snp_file_type(snp_path)
    loci_description = parse_loci_description(header_text)
    count_samples = count_samples_per_population(snp_path)
    sex_ratio = parse_sex_ratio(snp_path)
    maf_ratio = parse_maf_ratio(snp_path)
    mrc_ratio = parse_mrc_ratio(snp_path)
    reads_observed = None
    sexes_per_population = individual_sexes_per_population(snp_path)
    return SnpReplayContext(
        header_text=header_text,
        snp_path=snp_path,
        snp_file_type=snp_file_type,
        loci_description=loci_description,
        count_samples=count_samples,
        sex_ratio=sex_ratio,
        maf_ratio=maf_ratio,
        mrc_ratio=mrc_ratio,
        reads_observed=reads_observed,
        sexes_per_population=sexes_per_population,
    )


@pytest.fixture
def snp_context_te4(header_text_te4) -> SnpReplayContext:
    snp_path = OBSERVED_SNP_FILE_TE4
    snp_file_type = detect_snp_file_type(snp_path)
    loci_description = parse_loci_description(header_text_te4)
    count_samples = count_samples_per_population(snp_path)
    sex_ratio = parse_sex_ratio(snp_path)
    maf_ratio = parse_maf_ratio(snp_path)
    mrc_ratio = parse_mrc_ratio(snp_path)
    reads_observed = observed_reads(snp_path)
    sexes_per_population = (
        individual_sexes_per_population(snp_path) if snp_file_type == "IND" else {}
    )
    return SnpReplayContext(
        header_text=header_text_te4,
        snp_path=snp_path,
        snp_file_type=snp_file_type,
        loci_description=loci_description,
        count_samples=count_samples,
        sex_ratio=sex_ratio,
        maf_ratio=maf_ratio,
        mrc_ratio=mrc_ratio,
        reads_observed=reads_observed,
        sexes_per_population=sexes_per_population,
    )


@pytest.fixture
def snp_context_te5(header_text_te5) -> SnpReplayContext:
    snp_path = OBSERVED_SNP_FILE_TE5
    snp_file_type = detect_snp_file_type(snp_path)
    loci_description = parse_loci_description(header_text_te5)
    count_samples = count_samples_per_population(snp_path)
    sex_ratio = parse_sex_ratio(snp_path)
    maf_ratio = parse_maf_ratio(snp_path)
    mrc_ratio = parse_mrc_ratio(snp_path)
    reads_observed = None
    sexes_per_population = individual_sexes_per_population(snp_path)
    return SnpReplayContext(
        header_text=header_text_te5,
        snp_path=snp_path,
        snp_file_type=snp_file_type,
        loci_description=loci_description,
        count_samples=count_samples,
        sex_ratio=sex_ratio,
        maf_ratio=maf_ratio,
        mrc_ratio=mrc_ratio,
        reads_observed=reads_observed,
        sexes_per_population=sexes_per_population,
    )


@pytest.fixture
def dna_context_te2(header_text_te2) -> DnaReplayContext:
    list_loci = parse_loci_description(header_text_te2)
    dna_observed = observed_sequences(OBSERVED_MSS_FILE_TE2, list_loci)
    return DnaReplayContext(
        header_text=header_text_te2,
        mss_path=OBSERVED_MSS_FILE_TE2,
        list_loci=list_loci,
        dna_observed=dna_observed,
        frequencies_per_locus=base_frequency_by_locus(dna_observed),
        samples_default=observed_count_population(OBSERVED_MSS_FILE_TE2),
        sex_ratio=parse_sex_ratio(OBSERVED_MSS_FILE_TE2),
    )


@pytest.fixture
def dna_context_te2_xy(header_text_te2_XY) -> DnaReplayContext:
    list_loci = parse_loci_description(header_text_te2_XY)
    dna_observed = observed_sequences(OBSERVED_MSS_FILE_TE2_XY, list_loci)
    return DnaReplayContext(
        header_text=header_text_te2_XY,
        mss_path=OBSERVED_MSS_FILE_TE2_XY,
        list_loci=list_loci,
        dna_observed=dna_observed,
        frequencies_per_locus=base_frequency_by_locus(dna_observed),
        samples_default=observed_count_population(OBSERVED_MSS_FILE_TE2_XY),
        sex_ratio=parse_sex_ratio(OBSERVED_MSS_FILE_TE2_XY),
    )


@pytest.fixture
def microsat_context_te2_xy(header_text_te2_XY) -> MicrosatReplayContext:
    list_loci = parse_loci_description(header_text_te2_XY)
    microsat_observed = observed_microsatellites(OBSERVED_MSS_FILE_TE2_XY, list_loci)
    return MicrosatReplayContext(
        header_text=header_text_te2_XY,
        mss_path=OBSERVED_MSS_FILE_TE2_XY,
        list_loci=list_loci,
        microsat_observed=microsat_observed,
        bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
        samples_default=observed_count_population(OBSERVED_MSS_FILE_TE2_XY),
        sex_ratio=parse_sex_ratio(OBSERVED_MSS_FILE_TE2_XY),
    )


@pytest.fixture
def microsat_context_te1_modified(header_text_te1_modified) -> MicrosatReplayContext:
    list_loci = parse_loci_description(header_text_te1_modified)
    microsat_observed = observed_microsatellites(OBSERVED_MSS_FILE_TE1, list_loci)
    return MicrosatReplayContext(
        header_text=header_text_te1_modified,
        mss_path=OBSERVED_MSS_FILE_TE1,
        list_loci=list_loci,
        microsat_observed=microsat_observed,
        bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
        samples_default=observed_count_population(OBSERVED_MSS_FILE_TE1),
        sex_ratio=parse_sex_ratio(OBSERVED_MSS_FILE_TE1),
    )
