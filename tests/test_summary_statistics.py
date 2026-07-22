"""
Tests des statistiques résumées SNP PoolSeq (bridge/summary_statistics.py).

_prepare_matrices_poolseq : forme/valeurs des matrices (npop, nloci).
compute_all_statistics_poolseq : les 130 stats sont bien présentes, avec au
moins une valeur vérifiée à la main (HWm_1/HWv_1) pour attraper une
régression de formule, pas juste un problème de branchement.
"""

import numpy as np

from bridge.summary_statistics import (
    _prepare_matrices_poolseq,
    compute_all_statistics_poolseq,
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
