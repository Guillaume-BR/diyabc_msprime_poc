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

OBSERVED_SNP_FILE = REFERENCE_DIR / "human_snp_all22chr_maf5.snp"
OBSERVED_SNP_FILE_TE5 = (
    Path(__file__).parent.parent
    / "reference"
    / "toy_example5"
    / "simu_dataset_test_divergence_admixture_001.snp"
)


@pytest.fixture
def header_text() -> str:
    return (REFERENCE_DIR / "header.txt").read_text()


@pytest.fixture
def header_text_te5() -> str:
    return (OBSERVED_SNP_FILE_TE5.parent / "headerRF.txt").read_text()
