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

REFERENCE_DIR = Path(__file__).parent.parent / "reference" / "human"
GENERAL_BINARY_PATH = os.environ.get("DIYABC_GENERAL_PATH")

OBSERVED_SNP_FILE_HUMAN = REFERENCE_DIR / "human_snp_all22chr_maf5.snp"
OBSERVED_MSS_FILE_TE2 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example2_ms_dna"
    / "pseudo_observed_DATASET_toy_example2_microsatellites_DNAsequences_ancestral_admixture_unsampled_pops_001.mss"
)
OBSERVED_SNP_FILE_TE5 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example5"
    / "simu_dataset_test_divergence_admixture_001.snp"
)
OBSERVED_SNP_FILE_TE3 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example3"
    / "pseudo_observed_DATASET_SNP_INDSEQ_4pops_Scenario3_MER.snp"
)
OBSERVED_SNP_FILE_TE3_SCENARIO1 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example3_scenario1"
    / "pseudo_observed_DATASET_SNP_INDSEQ_4pops_Scenario3_MER.snp"
)

OBSERVED_SNP_FILE_TE4 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example4"
    / "pseudo_observed_DATASET_SNP_POOLSEQ_4pops_Scenario3_MER.snp"
)

OBSERVED_MSS_FILE_TE2_XY = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example2_ms_dna_XY"
    / "pseudo_observed_DATASET_toy_example2_microsatellites_DNAsequences_ancestral_admixture_unsampled_pops_001.mss"
)


@pytest.fixture
def header_text() -> str:
    return (REFERENCE_DIR / "header.txt").read_text()


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
    path_te1 = REFERENCE_DIR.parent / "toy_example1_ms" / "headerRF.txt"
    return path_te1.read_text()


@pytest.fixture
def header_text_te2() -> str:
    path_te2 = REFERENCE_DIR.parent / "toy_example2_ms_dna" / "headerRF.txt"
    return path_te2.read_text()


@pytest.fixture
def header_text_te2_XY() -> str:
    path_te2_XY = REFERENCE_DIR.parent / "toy_example2_ms_dna_XY" / "headerRF.txt"
    return path_te2_XY.read_text()
