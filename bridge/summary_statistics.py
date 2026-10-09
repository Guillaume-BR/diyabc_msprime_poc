"""
Implémentation Python des statistiques résumées calculées par DIYABC
(HeaderC::calstatobs / ParticleC::docalstat, src-JMC-C++/sumstat.cpp) --
SNP (IndSeq et PoolSeq) et séquences ADN.

PROTOCOLE DE VALIDATION : pour chaque formule implémentée ici, on vérifie
qu'elle produit les mêmes valeurs que le vrai binaire `general` sur les
MÊMES données en entrée -- comparaison exacte (à la précision float32
près), pas statistique.

Quatre familles, quatre formats d'entrée, quatre points d'entrée
(compute_all_statistics_<type> ci-dessous) :
  - SNP IndSeq : liste de dicts {nom_echantillon: [génotypes haploïdes
    0/1]}, un dict par locus -- la forme produite par
    ancestry_simulation.simulate_snp_genotypes.
  - SNP PoolSeq : liste de dicts {nom_echantillon: (nreads_dérivé,
    nreads_total)}, un dict par locus -- la forme produite par
    ancestry_simulation.simulate_poolseq_reads_with_mrc_filter.
  - Séquences ADN : dict {nom_locus: tskit.TreeSequence mutée} -- la
    forme produite par ancestry_simulation.dna_mutation_simulation_
    per_locus, tranchée par échantillon via _genotype_matrix_by_sample
    (format complètement différent des deux précédents, pas de liste de
    génotypes 0/1).
  - Microsats : dict {nom_locus: tskit.TreeSequence mutée} -- la forme
    produite par ancestry_simulation.microsat_mutation_simulation_per_locus,
    tranchée par échantillon via _genotype_matrix_by_sample.

Organisation : une fonction par famille de statistiques, suivant
exactement la nomenclature de sumstat.cpp (cal_snfl, cal_snhw, cal_snhb,
cal_snfsti...) pour faciliter la traçabilité entre le code Python et sa
source C++ de référence.
"""

import random
from itertools import combinations, permutations

import numpy as np
import tskit

from bridge.ancestry_simulation import compute_population_layout
from bridge.configuration import (
    _GROUP_STAT_SEED_OFFSET,
    _LIKELIHOOD_SEED_OFFSET,
    _locus_seed,
)
from bridge.header_dataclasses import FocalGenotypeFrequencies
from bridge.loci_parser import parse_loci_description

# ---------------------------------------------------------------------------
# Utilitaires scalaires (conservés pour la traçabilité et les tests unitaires)
# ---------------------------------------------------------------------------


def _allele_freq(haploid_genotypes: list[int]) -> float:
    """Calcule la fréquence de l'allèle dérivé (1) dans un échantillon.

    Équivalent de locuslist[loc].freq[pop][1] dans le code C++.

    Args:
        haploid_genotypes: Les génotypes (0/1) d'un échantillon, un
            locus.

    Returns:
        La fréquence de l'allèle dérivé (nan si échantillon vide).
    """
    n = len(haploid_genotypes)
    if n == 0:
        return float("nan")
    return sum(haploid_genotypes) / n


def _q1(haploid_genotypes: list[int]) -> float:
    """Calcule la probabilité d'identité par état intra-échantillon.

    Tirage SANS remise -- formule exacte de sumstat.cpp::q1 (cas SNP,
    bias=False) : `q1 = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1))`, où
    y1, y2 = comptes d'allèles 0 et 1 (= freq * n).

    Args:
        haploid_genotypes: Les génotypes (0/1) d'un échantillon, un
            locus.

    Returns:
        q1 (nan si échantillon de taille <= 1).
    """
    n = len(haploid_genotypes)
    if n <= 1:
        return float("nan")
    y2 = sum(haploid_genotypes)
    y1 = n - y2
    return (y1 * (y1 - 1) + y2 * (y2 - 1)) / (n * (n - 1))


def _q2(haploid_genotypes_a: list[int], haploid_genotypes_b: list[int]) -> float:
    """Calcule la probabilité d'identité par état inter-échantillons.

    Formule exacte de sumstat.cpp::q2 (cas SNP) :
    `q2 = (y11*y21 + y12*y22) / (n1*n2)`.

    Args:
        haploid_genotypes_a: Les génotypes (0/1) de l'échantillon A, un
            locus.
        haploid_genotypes_b: Les génotypes (0/1) de l'échantillon B,
            même locus.

    Returns:
        q2 (nan si l'une des deux échantillons est vide).
    """
    n1 = len(haploid_genotypes_a)
    n2 = len(haploid_genotypes_b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    y12 = sum(haploid_genotypes_a)
    y11 = n1 - y12
    y22 = sum(haploid_genotypes_b)
    y21 = n2 - y22
    return (y11 * y21 + y12 * y22) / (n1 * n2)


# ---------------------------------------------------------------------------
# Infrastructure numpy partagée
# ---------------------------------------------------------------------------


def _prepare_matrices(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construit les matrices (n_sample, n_loci) de comptes et fréquences.

    Appelé UNE SEULE FOIS dans compute_all_statistics_indseq et transmis via
    _mats à toutes les familles de statistiques -- évite de reconstruire
    les matrices (n_sample x n_loci) une fois par famille.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon, dans l'ordre voulu
            pour les lignes des matrices.

    Returns:
        Le tuple (counts, ns, freq0, freq1) :
        counts[i, l] = nb d'allèles dérivés (1) dans pop i au locus l ;
        ns[i, l] = nb total de lignées dans pop i au locus l ;
        freq1 = counts / ns, freq0 = 1 - freq1.
    """
    counts = np.array(
        [[sum(lg[p]) for lg in genotypes_per_locus] for p in sample_names],
        dtype=float,
    )
    ns = np.array(
        [[len(lg[p]) for lg in genotypes_per_locus] for p in sample_names],
        dtype=float,
    )
    freq1 = counts / ns
    freq0 = 1.0 - freq1
    return counts, ns, freq0, freq1


def _prepare_matrices_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]], sample_names: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construit les matrices (n_sample, n_loci) de comptes et fréquences pour POOLSEQ.

    Équivalent PoolSeq de _prepare_matrices -- même Returns, mais
    `counts`/`ns` viennent des lectures observées (reads), pas de
    génotypes individuels.

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon, dans l'ordre voulu
            pour les lignes des matrices.

    Returns:
        Le tuple (counts, ns, freq0, freq1) :
        counts[i, l] = nb de lectures dérivées observées dans pop i au
        locus l ; ns[i, l] = nb total de lectures observées dans pop i
        au locus l ; freq1 = counts / ns, freq0 = 1 - freq1.
    """
    counts = np.array(
        [
            [reads_per_locus[loc][p][0] for loc in range(len(reads_per_locus))]
            for p in sample_names
        ],
        dtype=float,
    )
    ns = np.array(
        [
            [reads_per_locus[loc][p][1] for loc in range(len(reads_per_locus))]
            for p in sample_names
        ],
        dtype=float,
    )
    freq1 = counts / ns
    freq0 = 1.0 - freq1
    return counts, ns, freq0, freq1


def _forward_fill(
    values: np.ndarray, valid: np.ndarray, fill: float = 0.0
) -> np.ndarray:
    """Propage la dernière valeur valide vue aux positions non valides
    (forward-fill vectorisé).

    Reproduit le comportement de la variable C++ non ré-initialisée
    x_prev dans cal_snfstd et cal_snnei. Implémentation : searchsorted
    sur les indices valides, O(n log n) tout-numpy, sans boucle Python.

    Args:
        values: Les valeurs, certaines invalides.
        valid: Masque booléen, même longueur que values.
        fill: Valeur utilisée pour les positions initiales sans aucune
            valeur valide précédente.

    Returns:
        Le tableau avec les positions invalides remplacées par la
        dernière valeur valide vue (ou `fill`).
    """
    n = len(values)
    if not valid.any():
        return np.full(n, fill)
    valid_idx = np.where(valid)[0]  # positions valides
    last = (
        np.searchsorted(valid_idx, np.arange(n), side="right") - 1
    )  # dernier valide <= i
    has_prev = last >= 0
    return np.where(has_prev, values[valid_idx[np.maximum(last, 0)]], fill)


def _half_sorted_by_pairs(v: list[int]) -> bool:
    """Filtre HALF de DIYABC (cal_snaml / cal_snf3r / cal_snf4r).

    Args:
        v: Une permutation d'indices.

    Returns:
        True si v est "half-sorted by pairs".
    """
    n = len(v)
    for i in range(n - 1, 0, -2):
        if v[i - 1] > v[i]:
            return False
        if i - 2 > 0 and v[i - 3] > v[i - 1]:
            return False
    return True


def _half_arrangements(n: int, r: int) -> list[list[int]]:
    """Calcule les arrangements HALF de r éléments parmi n.

    Ordre reproduit empiriquement depuis DIYABC (cal_snaml, cal_snf3r,
    cal_snf4r).

    Args:
        n: Le nombre d'éléments parmi lesquels arranger.
        r: La taille de chaque arrangement.

    Returns:
        La liste des arrangements HALF, chacun une liste de r indices.
    """
    result = []
    for combo in sorted(combinations(range(n), r), reverse=True):
        for perm in sorted(set(permutations(combo))):
            if _half_sorted_by_pairs(list(perm)):
                result.append(list(perm))
    return result


# ---------------------------------------------------------------------------
# ML1 : proportion de loci monomorphes par échantillon (cal_snfl, n_sample=1)
# ---------------------------------------------------------------------------


def compute_ML1(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML1p_i : proportion de loci monomorphes dans l'échantillon i.

    Un locus est monomorphe si sum==0 (fixé ancestral) ou sum==n (fixé
    dérivé).

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Matrices (counts, ns, freq0, freq1) déjà calculées par
            _prepare_matrices (voir compute_all_statistics_indseq). Si None,
            calculées ici.

    Returns:
        Un dict {"ML1p_i": valeur}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    mono = (counts == 0) | (counts == ns)  # (n_sample, n_loci) booléen
    return {f"ML1p_{i + 1}": float(mono[i].mean()) for i in range(len(sample_names))}


# ---------------------------------------------------------------------------
# ML2 : proportion de loci fixés identiquement sur les paires (cal_snfl, n_sample=2)
# ---------------------------------------------------------------------------


def compute_ML2(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML2p_i.j : proportion de loci fixés au même allèle dans la paire (i, j).

    Référence : cal_snfl(n_sample=2) -- freq_a == freq_b ∈ {0, 1}.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"ML2p_i.j": valeur}, une entrée par paire d'échantillons.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n = len(sample_names)
    results = {}
    for i in range(n):
        for j in range(i + 1, n):
            both_zero = (counts[i] == 0) & (counts[j] == 0)
            both_one = (counts[i] == ns[i]) & (counts[j] == ns[j])
            results[f"ML2p_{i + 1}.{j + 1}"] = float((both_zero | both_one).mean())
    return results


# ---------------------------------------------------------------------------
# ML3 : proportion de loci fixés identiquement sur les triplets (cal_snfl, n_sample=3)
# ---------------------------------------------------------------------------


def compute_ML3(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML3p_i.j.k : même logique que ML2, sur les triplets d'échantillons.

    Référence : cal_snfl(n_sample=3).

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"ML3p_i.j.k": valeur}, une entrée par triplet d'échantillons.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n = len(sample_names)
    results = {}
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                all_zero = (counts[i] == 0) & (counts[j] == 0) & (counts[k] == 0)
                all_one = (
                    (counts[i] == ns[i]) & (counts[j] == ns[j]) & (counts[k] == ns[k])
                )
                results[f"ML3p_{i + 1}.{j + 1}.{k + 1}"] = float(
                    (all_zero | all_one).mean()
                )
    return results


# ---------------------------------------------------------------------------
# HW / HB : hétérozygotie intra- et inter-échantillon (cal_snhw, cal_snhb)
# ---------------------------------------------------------------------------


def compute_HW_HB(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """HWm_i/HWv_i (intra-pop) et HBm_i.j/HBv_i.j (inter-pop).

    HW = 1 - q1, HB = 1 - q2 (formules de sumstat.cpp vectorisées) :
    q1[i,l] = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1)) [sans remise] ;
    q2[i,j,l] = (y1_i*y1_j + y2_i*y2_j) / (n_i*n_j). HWv et HBv
    utilisent ddof=1 (validé contre le C++).

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"HWm_i": ..., "HWv_i": ..., "HBm_i.j": ...,
        "HBv_i.j": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    y1, y2 = ns - counts, counts

    hw = 1.0 - (y1 * (y1 - 1) + y2 * (y2 - 1)) / (ns * (ns - 1))  # (n_sample, n_loci)

    n_sample = len(sample_names)
    results = {}
    for i in range(n_sample):
        results[f"HWm_{i + 1}"] = float(hw[i].mean())
        results[f"HWv_{i + 1}"] = float(hw[i].var(ddof=1))

    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            hb = 1.0 - (y1[i] * y1[j] + y2[i] * y2[j]) / (ns[i] * ns[j])
            key = f"{i + 1}.{j + 1}"
            results[f"HBm_{key}"] = float(hb.mean())
            results[f"HBv_{key}"] = float(hb.var(ddof=1))

    return results


def compute_HW_HB_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_HW_HB.

    Correction de biais de lecture Q1 nécessaire côté PoolSeq (la
    profondeur de séquençage `pool_sizes` n'est pas le vrai nombre de
    copies de gène échantillonnées) -- voir la formule `q1` ci-dessous,
    absente du chemin IndSeq.

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool},
            voir observed_data._parse_pool_header_line.
        _mats: Matrices (counts, ns, freq0, freq1) déjà calculées par
            _prepare_matrices_poolseq. Si None, calculées ici.

    Returns:
        Un dict {"HWm_i": ..., "HWv_i": ..., "HBm_i.j": ...,
        "HBv_i.j": ...}.
    """
    results = {}
    (
        counts,
        ns,
        _,
        _,
    ) = _mats or _prepare_matrices_poolseq(reads_per_locus, sample_names)
    n_sample = len(sample_names)
    ## Calcul de HWm et HWv pour chaque échantillon
    for i in range(n_sample):
        r1 = counts[i]
        c1 = ns[i]
        r2 = c1 - r1
        s1 = r1 * (r1 - 1)
        s2 = r2 * (r2 - 1)
        np_i = pool_sizes[sample_names[i]]
        q1 = ((np_i / (c1 * (c1 - 1))) * (s1 + s2) - 1) / (np_i - 1)
        hw = 1 - q1
        results[f"HWm_{i + 1}"] = float(hw.mean())
        results[f"HWv_{i + 1}"] = float(hw.var(ddof=1))

    # Calcul de HBm et HBv pour chaque paire d'échantillons
    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            r11 = counts[i]
            c1 = ns[i]
            r12 = c1 - r11
            r21 = counts[j]
            c2 = ns[j]
            r22 = c2 - r21
            q2 = (r11 * r21 + r12 * r22) / (c1 * c2)
            hb = 1 - q2

            key = f"{i + 1}.{j + 1}"
            results[f"HBm_{key}"] = float(hb.mean())
            results[f"HBv_{key}"] = float(hb.var(ddof=1))

    return results


# ---------------------------------------------------------------------------
# FST1 : FST échantillon-spécifique (cal_snfsti)
# ---------------------------------------------------------------------------


def compute_FST1(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST1m_i = 1 - HWm_i / HBmoy_global, FST1v_i = HWv_i / HBmoy_global².

    HBmoy_global = moyenne de TOUS les HBm (toutes paires confondues) --
    confirmé dans cal_snfsti (sumstat.cpp), pas seulement les paires de
    samp_i. FST1v est une propagation d'erreur analytique, pas une
    variance empirique.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"FST1m_i": ..., "FST1v_i": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    y1, y2 = ns - counts, counts

    hw = 1.0 - (y1 * (y1 - 1) + y2 * (y2 - 1)) / (ns * (ns - 1))  # (n_sample, n_loci)

    n_sample = len(sample_names)
    all_hbm = []
    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            hb = 1.0 - (y1[i] * y1[j] + y2[i] * y2[j]) / (ns[i] * ns[j])
            all_hbm.append(float(hb.mean()))

    hbmoy = float(np.mean(all_hbm)) if all_hbm else float("nan")

    results = {}
    for i in range(n_sample):
        hwm = float(hw[i].mean())
        hwv = float(hw[i].var(ddof=1))
        if hbmoy != 0:
            results[f"FST1m_{i + 1}"] = 1.0 - hwm / hbmoy
            results[f"FST1v_{i + 1}"] = hwv / (hbmoy**2)
        else:
            results[f"FST1m_{i + 1}"] = float("nan")
            results[f"FST1v_{i + 1}"] = float("nan")

    return results


def compute_FST1_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST1.

    cal_snfsti n'a pas de branche type==3 : elle combine juste des
    HW/HB déjà calculés. Duplique donc ici la formule q1/q2 poolseq de
    compute_HW_HB_poolseq, exactement comme compute_FST1 duplique déjà
    la formule q1/q2 IndSeq plutôt que d'appeler compute_HW_HB -- même
    style que l'existant.

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.
        _mats: Voir compute_HW_HB_poolseq.

    Returns:
        Un dict {"FST1m_i": ..., "FST1v_i": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(reads_per_locus, sample_names)
    n_sample = len(sample_names)

    hw_by_sample = []
    for i in range(n_sample):
        np_i = pool_sizes[sample_names[i]]
        r1, c1 = counts[i], ns[i]
        r2 = c1 - r1
        s1, s2 = r1 * (r1 - 1), r2 * (r2 - 1)
        q1 = ((np_i / (c1 * (c1 - 1))) * (s1 + s2) - 1) / (np_i - 1)
        hw_by_sample.append(1.0 - q1)

    all_hbm = []
    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            r11, c1 = counts[i], ns[i]
            r12 = c1 - r11
            r21, c2 = counts[j], ns[j]
            r22 = c2 - r21
            q2 = (r11 * r21 + r12 * r22) / (c1 * c2)
            hb = 1.0 - q2
            all_hbm.append(float(hb.mean()))

    hbmoy = float(np.mean(all_hbm)) if all_hbm else float("nan")

    results = {}
    for i in range(n_sample):
        hwm = float(hw_by_sample[i].mean())
        hwv = float(hw_by_sample[i].var(ddof=1))
        if hbmoy != 0:
            results[f"FST1m_{i + 1}"] = 1.0 - hwm / hbmoy
            results[f"FST1v_{i + 1}"] = hwv / (hbmoy**2)
        else:
            results[f"FST1m_{i + 1}"] = float("nan")
            results[f"FST1v_{i + 1}"] = float("nan")

    return results


# ---------------------------------------------------------------------------
# FST2/3/4 via Weir & Cockerham vectorisé (cal_snfstd)
# ---------------------------------------------------------------------------


def _fst_wc(loci, samples, _counts=None, _ns=None):
    """Calcule le FST de Weir & Cockerham, vectorisé sur tous les loci.

    Formule identique à cal_snfstd, toutes les opérations par-locus
    faites en numpy sur des vecteurs de longueur n_loci.

    Args:
        loci: Liste de dicts {nom_echantillon: [génotype, ...]}, un
            dict par locus.
        samples: Les noms d'échantillon à inclure dans ce calcul.
        _counts: Matrice (len(samples), n_loci) de comptes d'allèles
            dérivés, déjà calculée -- passée comme slice de la matrice
            globale depuis compute_FST2/3/4 pour éviter de reconstruire
            les comptes locus par locus pour chaque sous-ensemble. Si
            None, calculée ici.
        _ns: Matrice (len(samples), n_loci) de tailles d'échantillon,
            même principe que _counts.

    Returns:
        Le tuple (FSTm, FSTv).
    """
    n_loci = len(loci)
    n_sample = len(samples)

    if _counts is not None and _ns is not None:
        counts, ns = _counts, _ns
    else:
        counts = np.array([[sum(lg[p]) for lg in loci] for p in samples], dtype=float)
        ns = np.array([[len(lg[p]) for lg in loci] for p in samples], dtype=float)

    p1 = counts / ns
    p0 = 1.0 - p1

    S_1 = ns.sum(axis=0)
    S_2 = (ns**2).sum(axis=0)
    n_d = float(n_sample)

    pi0 = (ns * p0).sum(axis=0) / S_1
    pi1 = (ns * p1).sum(axis=0) / S_1

    SSI = (ns * p0 * (1 - p0) + ns * p1 * (1 - p1)).sum(axis=0)
    SSP = (ns * (p0 - pi0) ** 2 + ns * (p1 - pi1) ** 2).sum(axis=0)

    n_c = (S_1 - S_2 / S_1) / (n_d - 1.0)
    MSI = SSI / (S_1 - n_d)
    MSP = SSP / (n_d - 1.0)
    num = MSP - MSI
    den = MSP + (n_c - 1.0) * MSI

    valid = np.abs(den) > 0
    ratio = np.where(valid, num / np.where(valid, den, 1.0), 0.0)
    xs = _forward_fill(ratio, valid, fill=0.0)

    numt = num.sum()
    dent = den.sum()
    fstm = numt / dent if abs(dent) > 0 else 0.0

    sw2diff = n_loci * (n_loci - 1)
    mean = xs.mean()
    fstv = ((xs - mean) ** 2).sum() * n_loci / sw2diff if sw2diff > 0 else 0.0

    return float(fstm), float(fstv)


def compute_FST2(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST2m_i.j / FST2v_i.j : Weir & Cockerham par paire.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"FST2m_i.j": ..., "FST2v_i.j": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}
    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            key = f"{i + 1}.{j + 1}"
            m, v = _fst_wc(
                genotypes_per_locus,
                [sample_names[i], sample_names[j]],
                _counts=counts[[i, j]],
                _ns=ns[[i, j]],
            )
            results[f"FST2m_{key}"] = m
            results[f"FST2v_{key}"] = v
    return results


def _fst_wc_poolseq(samples, pool_sizes, _counts, _ns):
    """Variante PoolSeq de _fst_wc (cal_snfstd, branche grouplist[gr].type==3).

    Calcul en DEUX passes (contrairement à l'IndSeq) : pi1/pi2 (moyennes
    pondérées par la profondeur de lecture) doivent être connues avant de
    calculer SSP. C_1/C_1_star mélangent la profondeur de lecture (`n`,
    variable par locus) et la VRAIE taille du pool (`c`, constante par
    échantillon) -- c'est ce mélange qui constitue la correction propre à
    PoolSeq (le terme de variance intra-pool supplémentaire).

    L'agrégation finale (ratio de sommes num/den + variance via
    _forward_fill) est IDENTIQUE à _fst_wc -- confirmé par l'exploration
    C++, même code d'agrégation pour les deux types d'échantillon.

    Args:
        samples: Les noms d'échantillon à inclure dans ce calcul.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.
        _counts: Matrice (len(samples), n_loci) de lectures dérivées
            (nreads1), slice de la matrice globale -- même contrat que
            _fst_wc.
        _ns: Matrice (len(samples), n_loci) de profondeur de lecture
            (nreads_total), même principe que _counts.

    Returns:
        Le tuple (FSTm, FSTv).
    """
    x1, n = _counts, _ns  # (n_sample, n_loci) : reads allèle1, profondeur de lecture
    x2 = n - x1
    n_loci = n.shape[1]
    c = np.array([pool_sizes[p] for p in samples], dtype=float).reshape(-1, 1)

    # --- Passe 1 ---
    term = n / c + (c - 1) / c  # (n_sample, n_loci)
    C_1 = term.sum(axis=0)  # (n_loci,)
    C_1_star = (n * term).sum(axis=0)  # (n_loci,)
    R_1 = n.sum(axis=0)
    R_2 = (n * n).sum(axis=0)
    SSI = (x1 - x1 * x1 / n + x2 - x2 * x2 / n).sum(axis=0)

    pi1 = x1.sum(axis=0) / R_1
    pi2 = x2.sum(axis=0) / R_1
    C_1_star = C_1_star / R_1

    # --- Passe 2 (a besoin de pi1/pi2 de la passe 1) ---
    r1 = x1 / n - pi1
    r2 = x2 / n - pi2
    SSP = (n * (r1 * r1 + r2 * r2)).sum(axis=0)

    n_c = (R_1 - R_2 / R_1) / (C_1 - C_1_star)
    MSI = SSI / (R_1 - C_1)
    MSP = SSP / (C_1 - C_1_star)
    num = MSP - MSI
    den = MSP + (n_c - 1.0) * MSI

    # --- Agrégation : identique à _fst_wc ---
    valid = np.abs(den) > 0
    ratio = np.where(valid, num / np.where(valid, den, 1.0), 0.0)
    xs = _forward_fill(ratio, valid, fill=0.0)

    numt = num.sum()
    dent = den.sum()
    fstm = numt / dent if abs(dent) > 0 else 0.0

    sw2diff = n_loci * (n_loci - 1)
    mean = xs.mean()
    fstv = ((xs - mean) ** 2).sum() * n_loci / sw2diff if sw2diff > 0 else 0.0

    return float(fstm), float(fstv)


def compute_FST2_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST2 : FST2m_i.j / FST2v_i.j par paire.

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.
        _mats: Voir compute_HW_HB_poolseq.

    Returns:
        Un dict {"FST2m_i.j": ..., "FST2v_i.j": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(reads_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}
    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            key = f"{i + 1}.{j + 1}"
            samples = [sample_names[i], sample_names[j]]
            m, v = _fst_wc_poolseq(
                samples, pool_sizes, _counts=counts[[i, j]], _ns=ns[[i, j]]
            )
            results[f"FST2m_{key}"] = m
            results[f"FST2v_{key}"] = v
    return results


def compute_FST3_FST4_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST3_FST4_FSTG : FST3/FST4 sur
    triplets/quadruplets (COMB).

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.
        _mats: Voir compute_HW_HB_poolseq.

    Returns:
        Un dict {"FST3m_i.j.k": ..., "FST3v_i.j.k": ...} et, s'il y a
        au moins 4 échantillons, {"FST4m_i.j.k.l": ...,
        "FST4v_i.j.k.l": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(reads_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    for combo in combinations(range(n_sample), 3):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        samples = [sample_names[i] for i in idx]
        m, v = _fst_wc_poolseq(samples, pool_sizes, _counts=counts[idx], _ns=ns[idx])
        results[f"FST3m_{key}"] = m
        results[f"FST3v_{key}"] = v

    for combo in combinations(range(n_sample), 4):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        samples = [sample_names[i] for i in idx]
        m, v = _fst_wc_poolseq(samples, pool_sizes, _counts=counts[idx], _ns=ns[idx])
        results[f"FST4m_{key}"] = m
        results[f"FST4v_{key}"] = v

    return results


def compute_FST3_FST4_FSTG(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST3/FST4 : Weir & Cockerham sur triplets et quadruplets (COMB).

    Ne calcule PAS `FSTG` malgré son nom, seulement FST3/FST4 --
    `FSTG` (FST global, tous échantillons combinés d'un coup, même
    `cal_snfstd` que FST2/3/4 avec n_sample=0, voir statdefs.cpp:184) n'est
    implémenté nulle part dans ce module. Non bloquant en pratique :
    aucun header.txt de ce dépôt ne le déclare dans sa section 'group
    summary statistics', donc `stats_filter="HEADER"` ne le demande
    jamais.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"FST3m_i.j.k": ..., "FST3v_i.j.k": ...} et, s'il y a
        au moins 4 échantillons, {"FST4m_i.j.k.l": ...,
        "FST4v_i.j.k.l": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    for combo in combinations(range(n_sample), 3):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        m, v = _fst_wc(
            genotypes_per_locus,
            [sample_names[i] for i in idx],
            _counts=counts[idx],
            _ns=ns[idx],
        )
        results[f"FST3m_{key}"] = m
        results[f"FST3v_{key}"] = v

    for combo in combinations(range(n_sample), 4):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        m, v = _fst_wc(
            genotypes_per_locus,
            [sample_names[i] for i in idx],
            _counts=counts[idx],
            _ns=ns[idx],
        )
        results[f"FST4m_{key}"] = m
        results[f"FST4v_{key}"] = v

    return results


# ---------------------------------------------------------------------------
# NEI : distance de Nei (1972) vectorisée (cal_snnei)
# ---------------------------------------------------------------------------


def compute_NEI(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire, vectorisée.

    `NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)`. x_prev
    persiste si denom==0 (comportement C++ non réinitialisé) --
    reproduit via _forward_fill.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"NEIm_i.j": ..., "NEIv_i.j": ...}.
    """
    counts, ns, freq0, freq1 = _mats or _prepare_matrices(
        genotypes_per_locus, sample_names
    )
    n_loci = len(genotypes_per_locus)
    f, g = freq0, freq1
    norm = np.sqrt(f * f + g * g)  # (n_sample, n_loci)
    n_sample = len(sample_names)
    results = {}

    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            denom = norm[i] * norm[j]
            valid = denom > 0
            nei = np.where(
                valid,
                1.0 - (f[i] * f[j] + g[i] * g[j]) / np.where(valid, denom, 1.0),
                0.0,
            )
            xs = _forward_fill(nei, valid, fill=0.0)

            key = f"{i + 1}.{j + 1}"
            sw2diff = n_loci * (n_loci - 1)
            mean = xs.mean()
            results[f"NEIm_{key}"] = float(mean)
            results[f"NEIv_{key}"] = (
                float(((xs - mean) ** 2).sum() * n_loci / sw2diff)
                if sw2diff > 0
                else 0.0
            )

    return results


# ---------------------------------------------------------------------------
# AML : admixture maximum likelihood sur triplets HALF (cal_snaml)
# ---------------------------------------------------------------------------


def compute_AML(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """AMLm / AMLv : coefficient d'admixture ML sur triplets HALF.

    `aml = (f3-f2)/(f1-f2)` clampé à [0,1]. Les loci non informatifs
    (f1==f2, w=0) sont exclus de la moyenne -- équivalent au Welford
    pondéré avec w ∈ {0,1}, ce qui réduit à mean/var(ddof=1) sur le
    sous-ensemble informatif.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"AMLm_h.p1.p2": ..., "AMLv_h.p1.p2": ...}, une entrée
        par arrangement HALF de 3 échantillons.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    for t in _half_arrangements(n_sample, 3):
        h, p1, p2 = t[0], t[1], t[2]
        key = f"{h + 1}.{p1 + 1}.{p2 + 1}"

        f1 = freq0[p1]  # freq allèle 0 dans parent 1  (n_loci,)
        f2 = freq0[p2]  # freq allèle 0 dans parent 2
        f3 = freq0[h]  # freq allèle 0 dans l'hybride

        diff = f1 - f2
        informative = diff != 0  # w=1 si parents diffèrent

        aml_raw = np.where(
            informative, (f3 - f2) / np.where(informative, diff, 1.0), 0.5
        )
        x_inf = np.clip(aml_raw, 0.0, 1.0)[informative]

        results[f"AMLm_{key}"] = float(x_inf.mean()) if len(x_inf) > 0 else 0.0
        results[f"AMLv_{key}"] = float(x_inf.var(ddof=1)) if len(x_inf) > 1 else 0.0

    return results


# ---------------------------------------------------------------------------
# F3 / F4 : statistiques de Patterson vectorisées (cal_snf3r, cal_snf4r)
# ---------------------------------------------------------------------------


def compute_F3(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """F3m/F3v sur triplets HALF.

    `F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)` (pop0=hybride,
    pop1/2=parents). Tous les loci ont w=1 → mean/var(ddof=1)
    directement.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"F3m_i0.i1.i2": ..., "F3v_i0.i1.i2": ...}, une entrée
        par arrangement HALF de 3 échantillons.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    # --- F3 ---
    for t in _half_arrangements(n_sample, 3):
        i0, i1, i2 = t[0], t[1], t[2]
        key = f"{i0 + 1}.{i1 + 1}.{i2 + 1}"

        np_ = ns[i0]  # nb lignées dans l'hybride  (n_loci,)
        f1 = freq0[i0]  # freq allèle 0 dans l'hybride
        f2 = freq0[i1]  # freq allèle 0 dans parent 1
        f3 = freq0[i2]  # freq allèle 0 dans parent 2

        alpha = np.where(np_ > 1, f1 * (1 - f1) / np.where(np_ > 1, np_ - 1, 1.0), 0.0)
        x_vals = (f1 - f2) * (f1 - f3) - alpha

        results[f"F3m_{key}"] = float(x_vals.mean())
        results[f"F3v_{key}"] = float(x_vals.var(ddof=1))

    return results


def compute_F3_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_F3 (cal_snf3r, branche grouplist[gr].type==3).

    `alpha = ((np*a1p*(a1p-1))/(c1p*(c1p-1)) - a1p/c1p) / (np-1)`
    (pop0=hybride), `F3 = alpha + betaBC - betaAB - betaAC`, avec
    `beta_XY = (aXp*aYp)/(cXp*cYp)`.

    `np` = taille du pool (VRAIE, pas la profondeur de lecture) de la
    échantillon hybride -- vient de pool_sizes, pas de `ns`/`_mats`
    (contrairement à ns, qui est la profondeur de lecture, variable par
    locus). Agrégation identique à compute_F3 (mean/var(ddof=1) simples
    sur les loci, pas de ratio de sommes).

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.
        _mats: Voir compute_HW_HB_poolseq.

    Returns:
        Un dict {"F3m_i0.i1.i2": ..., "F3v_i0.i1.i2": ...}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(reads_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    for t in _half_arrangements(n_sample, 3):
        i0, i1, i2 = t[0], t[1], t[2]
        key = f"{i0 + 1}.{i1 + 1}.{i2 + 1}"

        np_i0 = pool_sizes[sample_names[i0]]

        a1p, c1p = counts[i0], ns[i0]
        a2p, c2p = counts[i1], ns[i1]
        a3p, c3p = counts[i2], ns[i2]

        alpha = ((np_i0 * a1p * (a1p - 1)) / (c1p * (c1p - 1)) - a1p / c1p) / (
            np_i0 - 1
        )
        betaAB = (a1p * a2p) / (c1p * c2p)
        betaAC = (a1p * a3p) / (c1p * c3p)
        betaBC = (a2p * a3p) / (c2p * c3p)
        x_vals = alpha + betaBC - betaAB - betaAC

        results[f"F3m_{key}"] = float(x_vals.mean())
        results[f"F3v_{key}"] = float(x_vals.var(ddof=1))

    return results


def compute_F4(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """F4m/F4v sur quadruplets HALF.

    `F4 = (a-b)*(c-d)`. Tous les loci ont w=1 → mean/var(ddof=1)
    directement.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.
        _mats: Voir compute_ML1.

    Returns:
        Un dict {"F4m_ia.ib.ic.id": ..., "F4v_ia.ib.ic.id": ...}, une
        entrée par arrangement HALF de 4 échantillons.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(genotypes_per_locus, sample_names)
    n_sample = len(sample_names)
    results = {}

    # --- F4 ---
    for t in _half_arrangements(n_sample, 4):
        ia, ib, ic, id_ = t[0], t[1], t[2], t[3]
        key = f"{ia + 1}.{ib + 1}.{ic + 1}.{id_ + 1}"

        a = freq0[ia]
        b = freq0[ib]
        c = freq0[ic]
        d = freq0[id_]
        x_vals = (a - b) * (c - d)

        results[f"F4m_{key}"] = float(x_vals.mean())
        results[f"F4v_{key}"] = float(x_vals.var(ddof=1))

    return results


# ---------------------------------------------------------------------------
# Statistiques pour les séquences ADN
# ---------------------------------------------------------------------------


def _genotype_matrix_by_sample(
    tree_sequence: tskit.TreeSequence,
    *,
    layout: list[tuple[str, np.ndarray]] | None = None,
) -> dict[str, np.ndarray]:
    """Découpe la matrice de génotypes d'un locus ADN par échantillon.

    genotype_matrix() n'est appelé qu'UNE FOIS pour toute la
    TreeSequence, puis tranché par échantillon via fancy indexing (pas
    de reconstruction par sample).

    Args:
        tree_sequence: La TreeSequence mutée d'un locus.
        layout: [(nom_echantillon, np.ndarray[indices_d'individus]), ...] : si fourni, remplace le découpage par échantillon à utilisé par le chemin sériel, où un échantillon n'est plus sa propre échantillon

    Returns:
        Un dict {nom_echantillon: matrice (n_sites, n_samples_pop)} --
        convention native de tskit (genotype_matrix() est déjà (sites,
        samples)), pas de transposition.
    """
    genotype_matrix = tree_sequence.genotype_matrix()
    if layout is None:
        layout = compute_population_layout(tree_sequence)
    return {
        samp_name: genotype_matrix[:, sample_ids] for samp_name, sample_ids in layout
    }


# ---------------------------------------------------------------------------
# NSS : nombre de sites ségrégeants par échantillon
# ---------------------------------------------------------------------------


def _segregating_sites_mask(matrix: np.ndarray) -> np.ndarray:
    """Calcule le masque des sites ségrégeants d'une matrice de génotypes.

    Factorisé hors de _count_segregating_sites pour être réutilisé par
    PSS (_private_segregating_sites_per_locus), qui a besoin du masque
    par site, pas seulement du compte agrégé.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Masque booléen (longueur n_sites) : True si le site n'est pas
        identique chez tous les échantillons.
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    return np.any(matrix != matrix[:, [0]], axis=1)


def _count_segregating_sites(matrix: np.ndarray) -> int:
    """Compte le nombre de sites polymorphes d'une matrice de génotypes.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Le nombre de sites ségrégeants.

    Raises:
        ValueError: Si matrix.shape[1] == 0 (échantillon sans échantillon).
    """
    return int(np.sum(_segregating_sites_mask(matrix)))


def compute_NSS(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule NSS_i (cal_nss1p) : pour chaque échantillon, la moyenne du
    nombre de sites ségrégeants sur tous les loci du groupe passé en
    argument (un groupe = les TreeSequences des loci séquence d'un même
    `group Gx` du header, ex. les 5 loci <A> de G2).

    `sample_names` fixe explicitement les clés du dict retourné (comme
    compute_ML1/compute_HW_HB) -- chaque échantillon attendu a toujours une
    valeur (0.0 par défaut, comme le `res = 0.0` du C++), même si
    `tree_sequences` est vide, plutôt que d'être silencieusement absente du
    résultat.

    Suppose que toutes les échantillons de `sample_names` sont présentes
    sur tous les loci du groupe (divise par `len(tree_sequences)`, pas par
    un décompte par échantillon comme le `nl` du C++ -- lève un KeyError si
    ce n'est pas le cas plutôt que d'exclure silencieusement ce locus,
    contrairement au C++) -- vérifié vrai sur toy_example2_ms_dna, pas
    garanti en général.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    num_loci = len(tree_sequences)
    mean_segregating_sites = {samp_name: 0.0 for samp_name in sample_names}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)

    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            matrix = genotype_matrices[samp_name]
            mean_segregating_sites[samp_name] += _count_segregating_sites(matrix)

    if num_loci > 0:
        for samp_name in sample_names:
            mean_segregating_sites[samp_name] /= num_loci

    return mean_segregating_sites


# ---------------------------------------------------------------------------
# NDH : nombre d'haplotypes distincts par échantillon
# ---------------------------------------------------------------------------


def _count_distinct_haplotypes(matrix: np.ndarray) -> int:
    """Compte le nombre d'haplotypes distincts d'une matrice de génotypes.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Le nombre d'haplotypes distincts (colonnes uniques).

    Raises:
        ValueError: Si matrix.shape[1] == 0 (échantillon sans échantillon).
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    # np.unique(axis=1) déduplique les colonnes (les haplotypes) --
    # le résultat a la forme (n_sites, n_haplotypes_distincts).
    distinct_haplotypes = np.unique(matrix, axis=1)
    return distinct_haplotypes.shape[1]


def compute_NHA(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule le nombre moyen d'haplotypes distincts par échantillon sur un
    groupe de loci.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    num_loci = len(tree_sequences)
    mean_distinct_haplotypes = {samp_name: 0.0 for samp_name in sample_names}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)

    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            matrix = genotype_matrices[samp_name]
            mean_distinct_haplotypes[samp_name] += _count_distinct_haplotypes(matrix)

    if num_loci > 0:
        for samp_name in sample_names:
            mean_distinct_haplotypes[samp_name] /= num_loci

    return mean_distinct_haplotypes


# ---------------------------------------------------------------------------
# MDP/VDP : moyenne et variance de différences de paires (pairwise differences) par échantillon
# ---------------------------------------------------------------------------


def _pairwise_hamming_distances(matrix: np.ndarray) -> np.ndarray:
    """Calcule les distances de Hamming par paire d'échantillons.

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Un vecteur 1D de longueur C(n_samples, 2) -- une valeur par
        paire (i, j) avec i < j, pas la matrice carrée
        (n_samples, n_samples) complète.

    Raises:
        ValueError: Si matrix.shape[1] == 0 (échantillon sans échantillon).
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")

    matrix_hammming = (matrix[:, :, None] != matrix[:, None, :]).sum(axis=0)
    triangle_sup = np.triu_indices(matrix_hammming.shape[0], k=1)
    return matrix_hammming[triangle_sup]


def compute_MPD(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule MPD_i (cal_mpd1p) : pour chaque échantillon, la moyenne du
    nombre de différences par paire (distance de Hamming) sur tous les
    loci du groupe passé en argument (un groupe = les TreeSequences des
    loci séquence d'un même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Contrairement à compute_NSS/compute_NHA
    (qui divisent par `len(tree_sequences)`), le
    dénominateur ici est un compteur PAR ECHANTILLON (comme le `nl` de
    cal_mpd1p) : un locus où un échantillon compte moins de deux séquences
    donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
    vide), auquel cas ce locus est exclu du calcul pour cet échantillon
    -- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
    laisser un `nan` (moyenne d'un vecteur vide) empoisonner le résultat.

    Suppose quand même que toutes les échantillons de `sample_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `compute_*` de ce module.

    Reproduit depuis le 02/10 un BUG de DIYABC : un locus sans aucune
        mutation est exclu du dénominateur, alors qu'il devrait y compter
        avec une contribution nulle. Voir le bloc de commentaires dans le
        corps de la fonction, et "G2 `MPD`/`VPD` residual" dans CLAUDE.md.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    mean_pairwise_differences = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        # ---------------------------------------------------------------
        # BUG DIYABC REPRODUIT -- voir "G2 `MPD`/`VPD` residual" dans
        # CLAUDE.md (Open work) et notes/exploration.md, 02/10 (suite).
        #
        # Sur un locus où la simulation n'a tiré AUCUNE mutation, DIYABC
        # produit des séquences VIDES : `init_dnaseq` (particuleC.cpp:1777)
        # part de `string dna = ""` et sa branche `else` ne remplit la
        # séquence que si `dnatrue` -- ce qui est faux sur le chemin
        # reftable. Or la chaîne vide EST son marqueur de donnée manquante
        # (`#define SEQMISSING ""`, particuleC.hpp:17). `cal_mpdpl` ne retient
        # donc aucune paire (`ndd = 0`) et `cal_mpd1p` exclut le locus de son
        # dénominateur (`if (nd > 0)`). `cal_nss1p`, lui, le compte quand même,
        # parce que son `OK` vient de `samplesize()`, qui lit les données
        # OBSERVÉES et non les séquences simulées.
        # Résultat : MPD est moyenné sur les seuls loci polymorphes, NSS sur
        # tous. Quand theta -> 0, MPD -> 1/a_79 = 0.202 au lieu de 0.
        #
        # VERSION CORRECTE, à restaurer en supprimant les deux lignes
        # `if ts.num_sites == 0: continue` ci-dessous : un locus sans
        # mutation est un locus réellement monomorphe ; il doit compter au
        # dénominateur avec une contribution nulle, ce que la théorie exige
        # (E[pi] = theta, qui tend vers 0 avec theta) et ce que ce code
        # faisait avant le 02/10.
        # ---------------------------------------------------------------
        if ts.num_sites == 0:
            continue
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            matrix = genotype_matrices[samp_name]
            pairwise_distances = _pairwise_hamming_distances(matrix)
            if len(pairwise_distances) > 0:
                mean_pairwise_differences[samp_name] += pairwise_distances.mean()
                valid_loci_count[samp_name] += 1

    for samp_name in sample_names:
        if valid_loci_count[samp_name] > 0:
            mean_pairwise_differences[samp_name] /= valid_loci_count[samp_name]

    return mean_pairwise_differences


def compute_VPD(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule VPD_i (cal_vpd1p) : pour chaque échantillon, la variance du
    nombre de différences par paire (distance de Hamming) sur tous les
    loci du groupe passé en argument (un groupe = les TreeSequences des
    loci séquence d'un même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Contrairement à compute_NSS/compute_NHA
    (qui divisent par `len(tree_sequences)`), le
    dénominateur ici est un compteur PAR ECHANTILLON (comme le `nl` de
    cal_vpd1p) : un locus où un échantillon compte moins de deux séquences
    donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
    vide), auquel cas ce locus est exclu du calcul pour cet échantillon
    -- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
    laisser un `nan` (variance d'un vecteur vide) empoisonner le résultat.

    Suppose quand même que toutes les échantillons de `sample_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `compute_*` de ce module.

    Reproduit depuis le 02/10 un BUG de DIYABC : un locus sans aucune
        mutation est exclu du dénominateur, alors qu'il devrait y compter
        avec une contribution nulle. Voir le bloc de commentaires dans le
        corps de la fonction, et "G2 `MPD`/`VPD` residual" dans CLAUDE.md.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_variance}.
    """
    variance_pairwise_differences = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        # ---------------------------------------------------------------
        # BUG DIYABC REPRODUIT -- voir "G2 `MPD`/`VPD` residual" dans
        # CLAUDE.md (Open work) et notes/exploration.md, 02/10 (suite).
        #
        # Sur un locus où la simulation n'a tiré AUCUNE mutation, DIYABC
        # produit des séquences VIDES : `init_dnaseq` (particuleC.cpp:1777)
        # part de `string dna = ""` et sa branche `else` ne remplit la
        # séquence que si `dnatrue` -- ce qui est faux sur le chemin
        # reftable. Or la chaîne vide EST son marqueur de donnée manquante
        # (`#define SEQMISSING ""`, particuleC.hpp:17). `cal_vpd1p` ne retient
        # donc aucune paire (`ndd = 0`) et `cal_vpd1p` exclut le locus de son
        # dénominateur (`if (nd > 1)`, un cran plus strict encore). `cal_nss1p`, lui, le compte quand même,
        # parce que son `OK` vient de `samplesize()`, qui lit les données
        # OBSERVÉES et non les séquences simulées.
        # Résultat : VPD est moyennée sur les seuls loci polymorphes, donc
        # gonflée du même facteur que MPD.
        #
        # VERSION CORRECTE, à restaurer en supprimant les deux lignes
        # `if ts.num_sites == 0: continue` ci-dessous : un locus sans
        # mutation est un locus réellement monomorphe ; il doit compter au
        # dénominateur avec une contribution nulle, ce que la théorie exige
        # (E[pi] = theta, qui tend vers 0 avec theta) et ce que ce code
        # faisait avant le 02/10.
        # ---------------------------------------------------------------
        if ts.num_sites == 0:
            continue
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            matrix = genotype_matrices[samp_name]
            pairwise_distances = _pairwise_hamming_distances(matrix)
            if len(pairwise_distances) > 1:
                variance_pairwise_differences[samp_name] += pairwise_distances.var(
                    ddof=1
                )
                valid_loci_count[samp_name] += 1

    for samp_name in sample_names:
        if valid_loci_count[samp_name] > 0:
            variance_pairwise_differences[samp_name] /= valid_loci_count[samp_name]

    return variance_pairwise_differences


# ---------------------------------------------------------------------------
# DTA : distance de Tajima par échantillon
# ---------------------------------------------------------------------------


def _tajima_constants(n_samples: int) -> tuple[float, float, float]:
    """Calcule les constantes a1, e1, e2 pour le D de Tajima (cal_dta1pl,
    lignes 1566-1575) à partir du nombre d'échantillons (n_samples) -- ne
    dépend que de n_samples, jamais des données elles-mêmes.

    Args:
        n_samples: Nombre d'échantillons (>= 2).

    Returns:
        Le tuple (a1, e1, e2) -- a1 est aussi réutilisé directement dans
        la formule finale du D (S / a1), e1/e2 sont les coefficients de
        la variance sous neutralité. b1, b2, c1, c2, a2 sont des étapes
        intermédiaires purement internes, jamais réutilisées ailleurs.
    """
    a1 = sum(1.0 / i for i in range(1, n_samples))
    a2 = sum(1.0 / (i * i) for i in range(1, n_samples))
    b1 = (n_samples + 1) / (n_samples - 1) / 3.0
    b2 = 2 * ((n_samples**2 + n_samples) + 3.0) / 9.0 / (n_samples**2 - n_samples)
    c1 = b1 - 1.0 / a1
    c2 = b2 - ((n_samples + 2) / a1 / n_samples) + (a2 / a1 / a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return a1, e1, e2


def _tajima_d_per_locus(matrix: np.ndarray) -> float | None:
    """D de Tajima (cal_dta1pl) pour UN échantillon sur UN locus.

    `D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))` -- pi = MPD (moyenne des
    différences par paire), S = NSS (nombre de sites ségrégeants),
    a1/e1/e2 = _tajima_constants(n_samples).

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Le D de Tajima, ou `None` si `n_samples < 2` (pas assez
        d'échantillons pour calculer ne serait-ce qu'une paire -- ce
        locus doit être EXCLU de la moyenne du groupe, `OKK = false`
        côté C++). Retourne `0.0` (pas `None`) si `n_samples >= 2` mais
        qu'aucun site n'est ségrégeant (S == 0, le dénominateur de la
        formule est nul) -- ce locus reste INCLUS dans la moyenne du
        groupe avec une valeur de 0.0, contrairement au cas précédent :
        le C++ ne repasse jamais `OKK` à false dans ce cas (lignes
        1579/1594-1598) -- distinction délibérée, pas une
        simplification de notre part.
    """
    n_samples = matrix.shape[1]
    if n_samples < 2:
        return None

    S = _count_segregating_sites(matrix)
    a1, e1, e2 = _tajima_constants(n_samples)
    denominator = e1 * S + e2 * S * (S - 1.0)
    if denominator <= 0:
        return 0.0

    pi = _pairwise_hamming_distances(matrix).mean()
    return (pi - S / a1) / np.sqrt(denominator)


def compute_DTA(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule DTA_i (cal_dta1p) : pour chaque échantillon, la moyenne du
    D de Tajima (_tajima_d_per_locus) sur tous les loci du groupe passé
    en argument (un groupe = les TreeSequences des loci séquence d'un
    même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Comme compute_MPD/variance_pairwise_
    differences_per_group, le dénominateur est un compteur PAR
    SAMPLE (le `nl` de cal_dta1p) : un locus où `_tajima_d_per_locus`
    retourne `None` (moins de 2 échantillons) est exclu -- ni ajouté à la
    somme, ni compté. Un locus où `_tajima_d_per_locus` retourne `0.0`
    (0 site ségrégeant) reste, lui, INCLUS dans le compte (voir
    _tajima_d_per_locus pour la distinction).

    Suppose quand même que toutes les échantillons de `sample_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `compute_*` de ce module.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Les layouts à utiliser pour chaque locus, ou None si aucun n'est fourni.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    mean_tajima_d = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            matrix = genotype_matrices[samp_name]
            tajima_d = _tajima_d_per_locus(matrix)
            if tajima_d is not None:
                mean_tajima_d[samp_name] += tajima_d
                valid_loci_count[samp_name] += 1

    for samp_name in sample_names:
        if valid_loci_count[samp_name] > 0:
            mean_tajima_d[samp_name] /= valid_loci_count[samp_name]

    return mean_tajima_d


# ---------------------------------------------------------------------------
# PSS : sites ségrégeants "privés" par échantillon (cal_pss1p)
# ---------------------------------------------------------------------------


def _private_segregating_sites_per_locus(
    genotype_matrices: dict[str, np.ndarray], target_sample: str
) -> int:
    """Compte les sites ségrégeants "privés" de `target_sample` sur UN locus
    (cal_pss1p).

    Un site ségrégeant privé est ségrégeant dans `target_sample` mais
    NULLE PART ailleurs, parmi TOUTES les échantillons de
    `genotype_matrices` (pas seulement ceux d'un même groupe -- le
    C++ compare à `this->nsample`, le nombre total d'échantillons du
    dataset).

    Toutes les matrices de `genotype_matrices` viennent du même
    `genotype_matrix()` (juste tranchées par colonnes, voir
    _genotype_matrix_by_sample) -- la ligne `i` désigne donc le MÊME
    site physique pour toutes les échantillons. Contrairement au C++, qui
    compare des listes d'indices de sites variables de longueurs
    différentes par une recherche imbriquée (`ssa[sample][j] ==
    ssa[sa][k]`), on peut donc comparer les masques booléens position par
    position directement -- pas de recherche d'égalité nécessaire.

    Args:
        genotype_matrices: Dict {nom_echantillon: matrice} pour TOUTES
            les échantillons du dataset (voir _genotype_matrix_by_sample).
        target_sample: L'échantillon pour lequel compter les sites
            privés.

    Returns:
        Le nombre de sites ségrégeants privés de target_sample.
    """
    target_mask = _segregating_sites_mask(genotype_matrices[target_sample])
    other_masks = [
        _segregating_sites_mask(matrix)
        for samp_name, matrix in genotype_matrices.items()
        if samp_name != target_sample
    ]
    segregating_elsewhere = (
        np.logical_or.reduce(other_masks) if other_masks else np.zeros_like(target_mask)
    )
    return int(np.sum(target_mask & ~segregating_elsewhere))


def compute_PSS(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule PSS_i (cal_pss1p) : pour chaque échantillon, la moyenne du
    nombre de sites ségrégeants privés (_private_segregating_sites_per_ locus)
    sur tous les loci du groupe passé en argument.

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. `sample_names` DOIT couvrir
    TOUTES les échantillons du dataset (pas seulement ceux d'un groupe),
    puisque `_private_segregating_sites_per_locus` compare `target_sample` à
    toutes les autres échantillons présents dans `genotype_matrices`.

    Contrairement à compute_NSS et aux autres
    fonctions `compute_*` de ce module, le dénominateur ici est
    `len(tree_sequences)` SANS AUCUNE exclusion (`nl` s'incrémente sans
    condition dans cal_pss1p, ligne 1624 -- pas de garde-fou du tout,
    même pas le `samplesize > 0` de NSS/NHA).

    Suppose que toutes les échantillons de `sample_names` sont
    présentes sur tous les loci du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus (tous ceux du
            dataset).
        layouts: Les layouts à utiliser pour chaque locus, ou None si aucun n'est fourni.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    mean_pss = {samp_name: 0.0 for samp_name in sample_names}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            mean_pss[samp_name] += _private_segregating_sites_per_locus(
                genotype_matrices, samp_name
            )

    if num_loci > 0:
        for samp_name in sample_names:
            mean_pss[samp_name] /= num_loci

    return mean_pss


# ---------------------------------------------------------------------------
# MNS/VNS : moyenne/variance du compte de l'allèle minoritaire (afs, cal_mns1p/cal_vns1p)
# ---------------------------------------------------------------------------


def _minor_allele_counts_at_segregating_sites(matrix: np.ndarray) -> np.ndarray:
    """Compte l'allèle minoritaire à chaque site ségrégeant.

    Pour chaque site ségrégeant d'une matrice de génotypes, compte le(s)
    allèle(s) le(s) moins fréquent(s) parmi les échantillons (afs,
    `nf[jj]` après tri croissant et saut des zéros --
    équivalent à `min(comptes des valeurs effectivement présentes)`,
    généralisation multi-allélique du "minor allele count" puisqu'un
    site ADN peut avoir jusqu'à 4 états A/C/G/T, pas seulement 2 comme
    un SNP).

    Args:
        matrix: Matrice de génotypes (n_sites, n_samples).

    Returns:
        Un vecteur 1D de longueur = nombre de sites SÉGRÉGEANTS (pas
        n_sites) -- un site fixé (une seule valeur parmi les
        échantillons) n'est pas inclus, comme `afs` qui ne pousse rien
        dans `t_afs` quand `jj >= 3` (une seule base présente).

    Raises:
        ValueError: Si matrix.shape[1] == 0 (échantillon sans échantillon).
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    minor_counts = []
    for site in matrix:
        _, counts = np.unique(site, return_counts=True)
        if len(counts) > 1:
            minor_counts.append(counts.min())
    return np.array(minor_counts, dtype=float)


def compute_MNS(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule MNS_i (cal_mns1p) : pour chaque échantillon, la moyenne,
    sur les loci du groupe, de la moyenne (par locus) des comptes
    d'allèle minoritaire aux sites ségrégeants
    (_minor_allele_counts_at_segregating_sites).

    Un locus sans site ségrégeant contribue 0.0 (la boucle C++ sur
    `t_afs` ne s'exécute simplement pas -- même effet qu'un vecteur
    vide ici). Comme PSS (et contrairement à MPD/VPD/DTA), AUCUNE
    exclusion de locus : `nl` = `len(tree_sequences)` sans condition
    (cal_mns1p, `nl++` inconditionnel, ligne 1705) -- pas besoin de
    compteur par échantillon.

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    mean_mns = {samp_name: 0.0 for samp_name in sample_names}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            minor_counts = _minor_allele_counts_at_segregating_sites(
                genotype_matrices[samp_name]
            )
            if len(minor_counts) > 0:
                mean_mns[samp_name] += minor_counts.mean()

    if num_loci > 0:
        for samp_name in sample_names:
            mean_mns[samp_name] /= num_loci

    return mean_mns


def compute_VNS(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule VNS_i (cal_vns1p) : pour chaque échantillon, la moyenne,
    sur les loci du groupe, de la variance (par locus) des comptes
    d'allèle minoritaire aux sites ségrégeants.

    ATTENTION : variance BIAISÉE (ddof=0, division par n, PAS n-1) --
    contrairement à compute_VPD (VPD) qui
    utilise ddof=1. Vérifié explicitement contre cal_vns1p, ligne 1731 :
    `v = (sx2 - sx*sx/a) / a`, pas `/ (a-1)`.

    Un locus avec moins de 2 sites ségrégeants contribue 0.0 (`v = 0.0`
    explicite, ligne 1734 -- pas assez de points pour une variance).
    Comme MNS/PSS, AUCUNE exclusion de locus au niveau du groupe (`nl`
    s'incrémente sans condition, ligne 1736) : le `if (v > 0.0) res +=
    v` du C++ (ligne 1738) n'exclut rien du dénominateur, il évite
    seulement d'ajouter un `v` négatif -- ce qui ne peut mathématiquement
    pas arriver pour une vraie variance (`v >= 0` toujours), donc ce
    garde-fou n'a aucun effet observable et n'est pas reproduit ici.

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {nom_echantillon: valeur_moyenne}.
    """
    variance_vns = {samp_name: 0.0 for samp_name in sample_names}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for samp_name in sample_names:
            minor_counts = _minor_allele_counts_at_segregating_sites(
                genotype_matrices[samp_name]
            )
            if len(minor_counts) > 1:
                variance_vns[samp_name] += minor_counts.var()

    if num_loci > 0:
        for samp_name in sample_names:
            variance_vns[samp_name] /= num_loci

    return variance_vns


# --------------------------------------------------------------------------
# NH2 : Nombre d'ahplotypes distincts pour un regroupement de deux échantillons
# --------------------------------------------------------------------------


def compute_NH2(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule NH2_ij (cal_nh2p) : pour chaque paire d'échantillons, la
    moyenne du nombre d'haplotypes distincts sur tous les loci du groupe
    passé en argument (un groupe = les TreeSequences des loci séquence
    d'un même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Un dict {"i.j": valeur_moyenne}, une entrée par paire d'échantillons.
    """
    num_loci = len(tree_sequences)

    pairs = [
        (ia, ib)
        for ia in range(len(sample_names))
        for ib in range(ia + 1, len(sample_names))
    ]
    mean_distinct_haplotypes = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)

    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for ia, ib in pairs:
            samp_a = sample_names[ia]
            samp_b = sample_names[ib]
            combined_matrix = np.hstack(
                (genotype_matrices[samp_a], genotype_matrices[samp_b])
            )
            key = f"{ia + 1}.{ib + 1}"
            mean_distinct_haplotypes[key] += _count_distinct_haplotypes(combined_matrix)

    if num_loci > 0:
        for key in mean_distinct_haplotypes:
            mean_distinct_haplotypes[key] /= num_loci

    return mean_distinct_haplotypes


# ---------------------------------------------------------------------------
# NS2 : Nombre de sites ségrégeants pour un regroupement de deux échantillons
# ---------------------------------------------------------------------------


def compute_NS2(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule NS2_ij (cal_ns2p) : pour chaque paire d'échantillons, la
    moyenne du nombre de sites ségrégeants sur tous les loci du groupe
    passé en argument (un groupe = les TreeSequences des loci séquence
    d'un même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en

    Returns:
        Un dict {"i.j": valeur_moyenne}, une entrée par paire d'échantillons.
    """
    pairs = [
        (ia, ib)
        for ia in range(len(sample_names))
        for ib in range(ia + 1, len(sample_names))
    ]
    mean_segregating_sites = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for ia, ib in pairs:
            samp_a = sample_names[ia]
            samp_b = sample_names[ib]
            combined_matrix = np.hstack(
                (genotype_matrices[samp_a], genotype_matrices[samp_b])
            )
            key = f"{ia + 1}.{ib + 1}"
            mean_segregating_sites[key] += _count_segregating_sites(combined_matrix)

    if num_loci > 0:
        for key in mean_segregating_sites:
            mean_segregating_sites[key] /= num_loci

    return mean_segregating_sites


# ---------------------------------------------------------------------------
# MP2 : Moyenne de différences par paire pour un regroupement de deux échantillons
# ---------------------------------------------------------------------------


def _mean_pairwise_differences_within_per_locus(
    genotype_matrices: dict[str, np.ndarray], samp_a: str, samp_b: str
) -> float:
    """Calcule MP2 "within" pour un locus : ratio poolé des sommes, PAS la
    moyenne des deux MPD.

    Additionne les distances de Hamming intra-échantillon de samp_a ET
    samp_b (jamais entre les deux), puis divise par le nombre total de
    paires des deux côtés -- équivalent à la moyenne des deux MPD
    seulement si samp_a et samp_b ont la même taille d'échantillon (voir
    compute_MP2).

    Args:
        genotype_matrices: Dict {nom_echantillon: matrice}, au moins
            samp_a et samp_b.
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.

    Returns:
        Le ratio poolé (somme des différences / somme des paires),
        0.0 si aucune des deux échantillons n'a de paire.
    """
    distances_a = _pairwise_hamming_distances(genotype_matrices[samp_a])
    distances_b = _pairwise_hamming_distances(genotype_matrices[samp_b])
    total_pairs = len(distances_a) + len(distances_b)
    total_differences = distances_a.sum() + distances_b.sum()
    if total_pairs > 0:
        return total_differences / total_pairs
    return float(total_differences)  # 0.0 en pratique, comme le C++


def compute_MP2(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule MP2_ij (cal_mp2p) : pour chaque paire d'échantillons, la
    moyenne du nombre de différences par paire (distance de Hamming)
    sur tous les loci du groupe passé en argument (un groupe = les
    TreeSequences des loci séquence d'un même `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en

    Returns:
        Un dict {"i.j": valeur_moyenne}, une entrée par paire d'échantillons.
    """
    pairs = [
        (ia, ib)
        for ia in range(len(sample_names))
        for ib in range(ia + 1, len(sample_names))
    ]
    mean_pairwise_differences = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for ia, ib in pairs:
            samp_a = sample_names[ia]
            samp_b = sample_names[ib]
            key = f"{ia + 1}.{ib + 1}"
            mean_pairwise_differences[key] += (
                _mean_pairwise_differences_within_per_locus(
                    genotype_matrices, samp_a, samp_b
                )
            )

    if num_loci > 0:
        for key in mean_pairwise_differences:
            mean_pairwise_differences[key] /= num_loci

    return mean_pairwise_differences


# ---------------------------------------------------------------------------
# MPB : Moyenne de différences par paire entre deux échantillons
# ---------------------------------------------------------------------------


def _pairwise_hamming_distances_between(
    matrix_a: np.ndarray, matrix_b: np.ndarray
) -> np.ndarray:
    """Calcule les distances de Hamming par paire entre deux matrices de
    génotypes.

    Contrairement à _pairwise_hamming_distances (matrice unique, une
    triangulaire à extraire), ici TOUTE paire (p dans matrix_a, q dans
    matrix_b) est valide -- jamais d'auto-comparaison, donc pas de
    triangle à extraire.

    Args:
        matrix_a: Matrice de génotypes (n_sites, n_samples_a).
        matrix_b: Matrice de génotypes (n_sites, n_samples_b), même
            n_sites que matrix_a.

    Returns:
        La matrice 2D (n_samples_a, n_samples_b) des distances de
        Hamming, pas un vecteur aplati.

    Raises:
        ValueError: Si matrix_a et matrix_b n'ont pas le même nombre de
            sites, ou si l'une des deux est vide (0 échantillon).
    """
    if matrix_a.shape[0] != matrix_b.shape[0]:
        raise ValueError("Les matrices doivent avoir le même nombre de sites.")
    if matrix_a.shape[1] == 0 or matrix_b.shape[1] == 0:
        raise ValueError("Une des matrices de génotypes est vide.")

    # Calculer les distances de Hamming entre toutes les paires d'échantillons
    distances = (matrix_a[:, :, None] != matrix_b[:, None, :]).sum(axis=0)
    return distances


def _mean_pairwise_differences_between_per_locus(
    genotype_matrices: dict[str, np.ndarray], samp_a: str, samp_b: str
) -> float:
    """Calcule la moyenne des différences par paire (MPB) entre deux
    échantillons pour un locus.

    Args:
        genotype_matrices: Dict {nom_echantillon: matrice}, au moins
            samp_a et samp_b.
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.

    Returns:
        La moyenne des distances de Hamming entre chaque échantillon de
        samp_a et chaque échantillon de samp_b.

    Raises:
        KeyError: Si samp_a ou samp_b n'est pas dans genotype_matrices.
    """
    if samp_a not in genotype_matrices or samp_b not in genotype_matrices:
        raise KeyError(
            "L'un des échantillons n'est pas présent dans les matrices de génotypes."
        )
    return _pairwise_hamming_distances_between(
        genotype_matrices[samp_a], genotype_matrices[samp_b]
    ).mean()


def compute_MPB(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule MPB_ij (cal_mpb2p) : pour chaque paire d'échantillons, la
    moyenne du nombre de différences par paire (distance de Hamming)
    entre les deux échantillons sur tous les loci du groupe passé en
    argument (un groupe = les TreeSequences des loci séquence d'un même
    `group Gx` du header).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque échantillon attendu a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    échantillons de `sample_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en
            fournit pas.

    Returns:
        Un dict {"i.j": valeur_moyenne}, une entrée par paire d'échantillons.
    """
    pairs = [
        (ia, ib)
        for ia in range(len(sample_names))
        for ib in range(ia + 1, len(sample_names))
    ]
    mean_pairwise_differences_between = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    num_loci = len(tree_sequences)
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for ia, ib in pairs:
            samp_a = sample_names[ia]
            samp_b = sample_names[ib]
            key = f"{ia + 1}.{ib + 1}"
            distances_between = _pairwise_hamming_distances_between(
                genotype_matrices[samp_a], genotype_matrices[samp_b]
            )
            mean_pairwise_differences_between[key] += distances_between.mean()

    if num_loci > 0:
        for key in mean_pairwise_differences_between:
            mean_pairwise_differences_between[key] /= num_loci

    return mean_pairwise_differences_between


# ---------------------------------------------------------------------------
# HST : Comparaison de la diversité entre deux échantillons à la diversité au
# sein de ces échantillons
# ---------------------------------------------------------------------------


def compute_HST(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule HST_ij (cal_fst2p) : pour chaque paire d'échantillons,
    une mesure de différenciation type FST à partir des loci du groupe
    passé en argument (un groupe = les TreeSequences des loci séquence
    d'un même `group Gx` du header).

    HST = num/den, où num = somme sur les loci de (Hb_locus - Hw_locus)
    et den = somme sur les loci de Hb_locus (Hb = MPB par-locus, Hw =
    MPW par-locus) -- un RATIO DE SOMMES accumulées sur tout le groupe,
    PAS une moyenne de valeurs par-locus divisée par num_loci (contraire
    à compute_MP2/compute_MPB) -- même schéma d'agrégation que _fst_wc
    (FST2/FST3/FST4, côté SNP).

    `sample_names` fixe explicitement les clés du dict retourné --
    chaque paire attendue a toujours une valeur (0.0 par défaut, y
    compris si `den == 0` pour cette paire), même si `tree_sequences`
    est vide. Suppose que toutes les échantillons de `sample_names`
    sont présentes sur tous les loci du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
            chacune).
        sample_names: Le nom des échantillons attendus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en
            fournit pas.

    Returns:
        Un dict {"i.j": valeur}, une entrée par paire d'échantillons.
    """
    pairs = [
        (ia, ib)
        for ia in range(len(sample_names))
        for ib in range(ia + 1, len(sample_names))
    ]
    mean_hst = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    num = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    den = {f"{ia + 1}.{ib + 1}": 0.0 for ia, ib in pairs}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        genotype_matrices = _genotype_matrix_by_sample(ts, layout=layout)
        for ia, ib in pairs:
            samp_a = sample_names[ia]
            samp_b = sample_names[ib]
            key = f"{ia + 1}.{ib + 1}"
            mpb = _mean_pairwise_differences_between_per_locus(
                genotype_matrices, samp_a, samp_b
            )
            mpw = _mean_pairwise_differences_within_per_locus(
                genotype_matrices, samp_a, samp_b
            )
            num[key] += mpb - mpw
            den[key] += mpb

    for key in mean_hst:
        if den[key] > 0:
            mean_hst[key] = num[key] / den[key]

    return mean_hst


# ---------------------------------------------------------------------------
# Statistiques pour les Microsatellites (microsat)
# ---------------------------------------------------------------------------

# Helper nécessaire pour le reste des stattistiques


def _length_by_sample(
    tree_sequence: tskit.TreeSequence,
    *,
    layout: list[tuple[str, np.ndarray]] | None = None,
) -> dict[str, list[tuple[int, int]]]:
    """Construit un dict {nom_echantillon: [(longueur, nb_sequence), ...]}.

    Args:
        tree_sequence: Un TreeSequence muté du groupe (un locus [M]).
        layout: [(nom_echantillon, np.ndarray[indices_d'individus]), ...] : si fourni, remplace le découpage par échantillon à utilisé par le chemin sériel, où un échantillon n'est plus sa propre échantillon

    Returns:
        Un dict {nom_echantillon: [(longueur, compte), ...]}.
    """
    length_by_sample = {}
    if layout is None:
        layout = compute_population_layout(
            tree_sequence
        )  # liste(tuple(pop,array(indice)))
    # si locus est monomorphe
    if tree_sequence.num_sites == 0:
        for samp_name, sample_ids in layout:
            length_by_sample[samp_name] = [(0, len(sample_ids))]
        return length_by_sample
    else:
        variant = next(tree_sequence.variants())
        tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
        for samp_name, sample_ids in layout:
            length_by_sample[samp_name] = []
            for taille in variant.alleles:
                taille_int = int(taille)
                length_by_sample[samp_name].append(
                    (taille_int, list(tailles[sample_ids]).count(taille_int))
                )

        return length_by_sample


# NAL : mean number of alleles across loci


def count_alleles_per_sample(
    length_by_sample: dict[str, list[tuple[int, int]]],
) -> dict[str, int]:
    """Compte le nombre de tuples ayant un compte > 0.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.

    Returns:
        Dict {nom_echantillon: nombre d'allèles distincts}.
    """
    allele_counts = {}
    for samp_name in length_by_sample:
        allele_counts[samp_name] = len(
            [result for result in length_by_sample[samp_name] if result[1] > 0]
        )
    return allele_counts


def compute_NAL(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule NAL_i : pour chaque échantillon, la moyenne du nombre d'allèles distincts
    sur tous les loci du groupe passé en argument (un groupe = les TreeSequences des loci séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dicts {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.

    Returns:
        Dict {nom_echantillon: NAL}.
    """

    allele_counts = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)

    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]

    for _, length_by_sample in zip(
        tree_sequences, length_by_sample_per_locus, strict=True
    ):
        counts = count_alleles_per_sample(length_by_sample)
        for samp_name in counts:
            allele_counts[samp_name] += counts[samp_name]
            valid_loci_count[samp_name] += 1

    # Calcul de la moyenne pour chaque échantillon
    for samp_name in sample_names:
        allele_counts[samp_name] /= (
            valid_loci_count[samp_name] if valid_loci_count[samp_name] > 0 else 1
        )

    return allele_counts


# HET : mean gene diversity across loci


def total_genes_copies_per_sample(
    length_by_sample: dict[str, list[tuple[int, int]]],
) -> dict[str, int]:
    """Compte le nombre total d'allèles distincts pour chaque échantillon.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.

    Returns:
        Dict {nom_echantillon: nombre total d'allèles distincts}.
    """
    total_counts = {}
    for samp_name in length_by_sample:
        total_counts[samp_name] = sum(
            [result[1] for result in length_by_sample[samp_name]]
        )
    return total_counts


def _compute_HET_for_one_sample(
    total_count: int, _lengths_counts: list[tuple[int, int]]
) -> float:
    """Calcule HET pour un échantillon donné à partir du nombre total d'allèles
    et des comptes par longueur.

    Args:
        total_count: Nombre total d'allèles distincts pour l'échantillon.
        _lengths_counts: Liste de tuples (longueur, nb_sequence) pour l'échantillon.

    Returns:
        La diversité génétique HET pour l'échantillon.
    """
    if total_count <= 1:
        return 0.0
    return (
        1 - sum((count / total_count) ** 2 for _, count in _lengths_counts if count > 0)
    ) * (total_count / (total_count - 1))


def compute_HET(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule HET_i : pour chaque échantillon, la moyenne de la diversité génétique
    sur tous les loci du groupe passé en argument (un groupe = les TreeSequences des loci
    séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dicts {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.

    Returns:
        Dict {nom_echantillon: HET}.
    """
    gene_diversity = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, length_by_sample in zip(
        tree_sequences, length_by_sample_per_locus, strict=True
    ):
        total_counts = total_genes_copies_per_sample(length_by_sample)
        for samp_name in total_counts:
            if total_counts[samp_name] > 1:
                gene_diversity[samp_name] += _compute_HET_for_one_sample(
                    total_counts[samp_name], length_by_sample[samp_name]
                )
                valid_loci_count[samp_name] += 1
    # Calcul de la moyenne pour chaque échantillon
    for samp_name in sample_names:
        gene_diversity[samp_name] /= (
            valid_loci_count[samp_name] if valid_loci_count[samp_name] > 0 else 1
        )

    return gene_diversity


# VAR : mean allele size variance across loci


def _compute_VAR_constants(
    length_by_sample: dict[str, list[tuple[int, int]]],
) -> tuple[dict[str, float], dict[str, int], dict[str, int]]:
    """Calcule s = somme des tailles brutes (en pb, pas les comptes par valeur distincte), v = somme des tailles brutes au carré
    et n = nombre total d'allèles.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.
        motif_size: Taille du motif pour le locus.

    Returns:
        Tuple de trois dictionnaires : {nom_echantillon: s}, {nom_echantillon: v}, {nom_echantillon: n}.
    """
    raw_sizes = {
        samp_name: sum(
            length * count for length, count in length_by_sample[samp_name] if count > 0
        )
        for samp_name in length_by_sample
    }
    raw_square_sizes = {
        samp_name: sum(
            length**2 * count
            for length, count in length_by_sample[samp_name]
            if count > 0
        )
        for samp_name in length_by_sample
    }
    total_counts = total_genes_copies_per_sample(length_by_sample)
    return raw_sizes, raw_square_sizes, total_counts


def _compute_VAR_for_one_sample(
    raw_sizes: float, raw_square_sizes: float, total_count: int, motif_size: int
) -> float:
    """Calcule VAR pour un échantillon donné à partir des sommes brutes et du
    nombre total d'allèles.

    Args:
        raw_sizes: Somme des tailles brutes (en pb) pour l'échantillon.
        raw_square_sizes: Somme des tailles brutes au carré pour l'échantillon.
        total_count: Nombre total d'allèles distincts pour l'échantillon.
        motif_size: Taille du motif pour le locus.

    Returns:
        La variance de la taille des allèles VAR pour l'échantillon.
    """
    if total_count <= 1:
        return 0.0
    return (
        (raw_square_sizes - (raw_sizes**2) / total_count)
        / (total_count - 1)
        / (motif_size**2)
    )


def compute_VAR(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    list_motif_sizes: list[int],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule VAR_i : pour chaque échantillon, la moyenne de la variance de la taille des allèles
    sur tous les loci du groupe passé en argument (un groupe = les TreeSequences des loci séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        list_motif_sizes: Liste des tailles de motifs pour chaque locus.
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Dict {nom_echantillon: VAR}.
    """

    allele_size_variance = {samp_name: 0.0 for samp_name in sample_names}
    valid_loci_count = {samp_name: 0 for samp_name in sample_names}
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, motif_size, length_by_sample in zip(
        tree_sequences, list_motif_sizes, length_by_sample_per_locus, strict=True
    ):
        raw_sizes, raw_square_sizes, total_counts = _compute_VAR_constants(
            length_by_sample
        )
        for samp_name in length_by_sample:
            if total_counts[samp_name] > 1:
                allele_size_variance[samp_name] += _compute_VAR_for_one_sample(
                    raw_sizes[samp_name],
                    raw_square_sizes[samp_name],
                    total_counts[samp_name],
                    motif_size,
                )
                valid_loci_count[samp_name] += 1

    # Calcul de la moyenne pour chaque échantillon
    for samp_name in sample_names:
        allele_size_variance[samp_name] /= (
            valid_loci_count[samp_name] if valid_loci_count[samp_name] > 0 else 1
        )

    return allele_size_variance


# MGW : mean M index across loci


def _compute_MGW_by_locus(
    length_by_sample: dict[str, list[tuple[int, int]]], motif_size: int
) -> dict[str, tuple[float, float]]:
    """Calcule (num, den) de MGW pour un locus, par échantillon.

    cal_mgw1p accumule num et den séparément sur tous les loci du groupe
    et ne divise qu'une seule fois à la fin (num_total/den_total) --
    PAS une moyenne des ratios num/den par locus. Cette brique retourne
    donc les deux termes séparés plutôt qu'un ratio déjà calculé, pour
    que compute_MGW puisse les accumuler correctement.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.
        motif_size: Taille du motif pour ce locus.

    Returns:
        Dict {nom_echantillon: (num, den)}.
    """
    result = {}
    for samp_name in length_by_sample:
        lengths_present = [
            length for length, count in length_by_sample[samp_name] if count > 0
        ]
        if not lengths_present:
            continue
        num_alleles = len(lengths_present)
        den = 1 + (max(lengths_present) - min(lengths_present)) / motif_size
        result[samp_name] = (num_alleles, den)
    return result


def compute_MGW(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    list_motif_sizes: list[int],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule MGW_i : pour chaque échantillon, la moyenne de l'indice M
    sur tous les loci du groupe passé en argument (un groupe = les
    TreeSequences des loci séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        list_motif_sizes: Liste des tailles de motifs pour chaque locus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.

    Returns:
        Dict {nom_echantillon: MGW}.
    """
    num_sum = {samp_name: 0.0 for samp_name in sample_names}
    den_sum = {samp_name: 0.0 for samp_name in sample_names}

    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, motif_size, length_by_sample in zip(
        tree_sequences, list_motif_sizes, length_by_sample_per_locus, strict=True
    ):
        for samp_name, (num, den) in _compute_MGW_by_locus(
            length_by_sample, motif_size
        ).items():
            num_sum[samp_name] += num
            den_sum[samp_name] += den

    # cal_mgw1p : un seul ratio num_total/den_total, pas une moyenne de
    # ratios par locus -- voir docstring de _compute_MGW_by_locus.
    mgw_values = {
        samp_name: num_sum[samp_name] / den_sum[samp_name]
        if den_sum[samp_name] > 0
        else 0.0
        for samp_name in sample_names
    }

    return mgw_values


# N2P - mean number of alleles across loci (two samples)


def _compute_N2P_for_one_pair(
    samp_a: list[tuple[int, int]], samp_b: list[tuple[int, int]]
) -> float:
    """Calcule N2P_ij pour une paire d'échantillons à partir des listes de
    tuples (longueur, nb_sequence).

    Args:
        samp_a: Liste de tuples (longueur, nb_sequence) pour la première échantillon.
        samp_b: Liste de tuples (longueur, nb_sequence) pour la seconde échantillon.

    Returns:
        La moyenne du nombre d'allèles distincts entre les deux échantillons.
    """
    alleles_a = {length for length, count in samp_a if count > 0}
    alleles_b = {length for length, count in samp_b if count > 0}
    combined_alleles = alleles_a.union(alleles_b)
    return len(combined_alleles)


def _compute_N2P_for_one_locus(
    length_by_sample: dict[str, list[tuple[int, int]]], sample_names: list[str]
) -> dict[str, float]:
    """Calcule le nombre total d'allèles distincts pour un locus donné.

    Args:
        length_by_sample: Dictionnaire {nom_echantillon: [tuples (longueur, nb_sequence)]}.

    Returns:
        Le nombre total d'allèles distincts.
    """
    pairs = pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]

    combined_alleles = {f"{i + 1}.{j + 1}": 0 for i, j in pairs}
    for i, j in pairs:
        key = f"{i + 1}.{j + 1}"
        if sample_names[i] in length_by_sample and sample_names[j] in length_by_sample:
            lengths_a = length_by_sample[sample_names[i]]
            lengths_b = length_by_sample[sample_names[j]]
            combined_alleles[key] += _compute_N2P_for_one_pair(lengths_a, lengths_b)
        elif sample_names[i] in length_by_sample:
            lengths_a = length_by_sample[sample_names[i]]
            combined_alleles[key] += count_alleles_per_sample(length_by_sample)[
                sample_names[i]
            ]
        elif sample_names[j] in length_by_sample:
            lengths_b = length_by_sample[sample_names[j]]
            combined_alleles[key] += count_alleles_per_sample(length_by_sample)[
                sample_names[j]
            ]
    return combined_alleles


def compute_N2P(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotype_by_sample_per_locus: list[dict[str, list[tuple[int] | tuple[int, int]]]]
    | None = None,
) -> dict[str, float]:
    """Calcule N2P_i_j : pour chaque paire d'échantillons, la moyenne du nombre d'allèles distincts
    sur tous les loci du groupe passé en argument (un groupe = les TreeSequences des loci séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.
        genotype_by_sample_per_locus: non utilisé mais pour conserver la compatibilité avec le dispatch.

    Returns:
        Dict {"i.j": N2P}.
    """
    valid_loci = {}
    all_combined_alleles = {}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, length_by_sample in zip(
        tree_sequences, length_by_sample_per_locus, strict=True
    ):
        combined_alleles = _compute_N2P_for_one_locus(length_by_sample, sample_names)
        for key in combined_alleles:
            valid_loci.setdefault(key, 0)
            valid_loci[key] += 1
            all_combined_alleles.setdefault(key, 0.0)
            all_combined_alleles[key] += combined_alleles[key]

    for key in all_combined_alleles:
        all_combined_alleles[key] /= valid_loci[key]

    return all_combined_alleles


# H2P - mean gene diversity across loci (two samples)


def _pool_allele_counts_for_two_samples(
    pop_1: list[tuple[int, int]], pop_2: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Fusionne les comptes bruts de deux échantillons, allèle par allèle.

    Ne calcule ni H2P ni une fréquence -- juste n_i·freq_i + n_j·freq_j
    (= compte_i + compte_j, la division par n_i/n_j s'annulant avec la
    multiplication), pour que _compute_H2P_for_one_pair puisse
    réutiliser _compute_HET_for_one_sample dessus (qui fait
    lui-même la division par le total).

    Args:
        pop_1: Liste de tuples (longueur, nb_sequence) pour l'échantillon 1.
        pop_2: Liste de tuples (longueur, nb_sequence) pour l'échantillon 2.

    Returns:
        Liste de tuples (longueur, compte poolé) pour chaque allèle.
    """
    pooled_counts = []
    n1 = sum(count for _, count in pop_1)
    n2 = sum(count for _, count in pop_2)
    n = n1 + n2
    if n <= 1:
        return [(0, 0.0)]
    else:
        for length, count_1 in pop_1:
            count_2 = next((count for len, count in pop_2 if len == length), 0)
            total = count_1 + count_2
            pooled_counts.append((length, total))

    return pooled_counts


def _compute_H2P_for_one_pair(
    length_by_sample: dict[str, list[tuple[int, int]]], samp_a: str, samp_b: str
) -> float:
    """Calcule H2P_ij pour une paire d'échantillons à partir des listes de
    tuples (longueur, nb_sequence).

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.

    Returns:
        La diversité génétique H2P pour la paire d'échantillons.
    """
    if samp_a not in length_by_sample or samp_b not in length_by_sample:
        raise KeyError(
            "L'un des échantillons n'est pas présent dans les matrices de génotypes."
        )
    pooled_counts = _pool_allele_counts_for_two_samples(
        length_by_sample[samp_a], length_by_sample[samp_b]
    )
    total_count = sum(count for _, count in pooled_counts)
    if total_count <= 1:
        return 0.0
    return (
        1 - sum((count / total_count) ** 2 for _, count in pooled_counts if count > 0)
    ) * (total_count / (total_count - 1))


def compute_H2P(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotype_by_sample_per_locus: list[dict[str, list[tuple[int] | tuple[int, int]]]]
    | None = None,
) -> dict[str, float]:
    """Calcule la diversité génétique H2P pour toutes les paires
    d'échantillons.

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.
        genotype_by_sample_per_locus: non utilisé mais pour conserver la compatibilité avec le dispatch.

    Returns:
        Dict {i.j: H2P} pour chaque paire d'échantillons.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]

    H2P_values = {f"{i + 1}.{j + 1}": [] for i, j in pairs}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]

    for _, length_by_sample in zip(
        tree_sequences, length_by_sample_per_locus, strict=True
    ):
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            if (
                sample_names[i] in length_by_sample
                and sample_names[j] in length_by_sample
            ):
                H2P_value = _compute_H2P_for_one_pair(
                    length_by_sample, sample_names[i], sample_names[j]
                )
                H2P_values[key].append(H2P_value)
            else:
                total_counts = total_genes_copies_per_sample(length_by_sample)
                if sample_names[i] in length_by_sample:
                    _lengths_counts = length_by_sample[sample_names[i]]
                    H2P_value = _compute_HET_for_one_sample(
                        total_counts[sample_names[i]], _lengths_counts
                    )
                    H2P_values[key].append(H2P_value)
                elif sample_names[j] in length_by_sample:
                    _lengths_counts = length_by_sample[sample_names[j]]
                    H2P_value = _compute_HET_for_one_sample(
                        total_counts[sample_names[j]], _lengths_counts
                    )
                    H2P_values[key].append(H2P_value)
    for key in H2P_values:
        H2P_values[key] = (
            sum(H2P_values[key]) / len(H2P_values[key])
            if len(H2P_values[key]) > 0
            else 0.0
        )
    return H2P_values


# V2P : mean allele size variance across loci (two samples)


def _compute_V2P_constants(
    population1: str,
    population2: str,
    raw_sizes: dict[str, float],
    raw_square_sizes: dict[str, int],
    total_counts: dict[str, int],
) -> tuple[float, float, int]:
    """Calcule les constantes nécessaires pour V2P pour une paire
    d'échantillons.

    Args:
        population1: Nom du premier échantillon.
        population2: Nom du second échantillon.
        raw_sizes: Dict {nom_echantillon: somme des tailles brutes}.
        raw_square_sizes: Dict {nom_echantillon: somme des tailles brutes au carré}.
        total_counts: Dict {nom_echantillon: nombre total d'allèles distincts}.

    Returns:
        Tuple (raw_size_sum, raw_square_size_sum, total_count_sum) pour la paire
    """
    raw_size_sum = raw_sizes.get(population1, 0.0) + raw_sizes.get(population2, 0.0)
    raw_square_size_sum = raw_square_sizes.get(population1, 0) + raw_square_sizes.get(
        population2, 0
    )
    total_count_sum = total_counts.get(population1, 0) + total_counts.get(
        population2, 0
    )
    return raw_size_sum, raw_square_size_sum, total_count_sum


def compute_V2P(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    list_motif_sizes: list[int],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule la variance de la taille des allèles entre deux échantillons.

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        list_motif_sizes: Liste des tailles de motifs pour chaque locus.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.

    Returns:
        Dict {"i.j": V2P} pour chaque paire d'échantillons.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]

    V2P_values = {f"{i + 1}.{j + 1}": [] for i, j in pairs}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, motif_size, length_by_sample in zip(
        tree_sequences, list_motif_sizes, length_by_sample_per_locus, strict=True
    ):
        raw_sizes, raw_square_sizes, total_counts = _compute_VAR_constants(
            length_by_sample
        )
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            samp_a = sample_names[i]
            samp_b = sample_names[j]
            raw_size_sum, raw_square_size_sum, total_count_sum = _compute_V2P_constants(
                samp_a, samp_b, raw_sizes, raw_square_sizes, total_counts
            )
            if total_count_sum > 1:
                V2P_value = _compute_VAR_for_one_sample(
                    raw_size_sum, raw_square_size_sum, total_count_sum, motif_size
                )
                V2P_values[key].append(V2P_value)

    for key in V2P_values:
        V2P_values[key] = (
            sum(V2P_values[key]) / len(V2P_values[key])
            if len(V2P_values[key]) > 0
            else 0.0
        )

    return V2P_values


# [DAS] - shared allele distance between two samples (Chakraborty and Jin 1993)


def _compute_identical_pair_for_one_pair(
    length_by_sample: dict[str, list[tuple[int, int]]], samp_a: str, samp_b: str
) -> tuple[int, int]:
    """Calcule le nombre de paires d'allèles identiques et le nombre total de
    paires possibles pour une paire d'échantillons.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}.
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.

    Returns:
        Tuple (identical_count, total_count) pour la paire d'échantillons.
    """
    counts_a = {
        length: count for length, count in length_by_sample.get(samp_a, []) if count > 0
    }
    counts_b = {
        length: count for length, count in length_by_sample.get(samp_b, []) if count > 0
    }
    identical_count = sum(
        counts_a.get(length, 0) * counts_b.get(length, 0)
        for length in set(counts_a) | set(counts_b)
    )
    total_count = sum(counts_a.values()) * sum(
        counts_b.values()
    )  # nombre total de paires possibles
    return identical_count, total_count


def compute_DAS(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotype_by_sample_per_locus: list[dict[str, list[tuple[int] | tuple[int, int]]]]
    | None = None,
) -> dict[str, float]:
    """Calcule la distance d'allèle partagée (DAS) entre toutes les paires
    d'échantillons.

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Un layout par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.
        genotype_by_sample_per_locus: non utilisé mais pour conserver la compatibilité avec le dispatch.
    Returns:
        Dict {i.j: DAS} pour chaque paire d'échantillons.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]
    identical_sum = {f"{i + 1}.{j + 1}": 0 for i, j in pairs}
    total_sum = {f"{i + 1}.{j + 1}": 0 for i, j in pairs}

    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, length_by_sample in zip(
        tree_sequences, length_by_sample_per_locus, strict=True
    ):
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            identical_count, total_count = _compute_identical_pair_for_one_pair(
                length_by_sample, sample_names[i], sample_names[j]
            )
            identical_sum[key] += identical_count
            total_sum[key] += total_count

    return {
        key: identical_sum[key] / total_sum[key] if total_sum[key] > 0 else 0.0
        for key in identical_sum
    }


# [DM2] - distance between two samples (Goldstein et al. 1995)


def _compute_DM2_for_one_locus(
    samp_a: str,
    samp_b: str,
    motif_size: int,
    length_by_sample: dict[str, list[tuple[int, int]]],
    raw_sizes: dict[str, float],
    total_counts: dict[str, int],
    previous_moy: tuple[float, float] | None,
) -> tuple[float, tuple[float, float] | None, bool]:
    """Calcule la contribution de DM2 à UN locus, pour une paire
    d'échantillons.

    Reproduit fidèlement un bug de cal_dmu2p (sumstat.cpp) : dans le
    C++, le buffer moy[] est alloué UNE SEULE FOIS avant la boucle sur
    les loci, et n'est réécrit que si les deux échantillons ont des
    échantillons à ce locus (sasize*sasize1 > 0) -- sinon il garde les
    valeurs (moy_a, moy_b) du DERNIER locus valide. Mais la ligne
    d'accumulation (dmu2 += sqr((moy[1]-moy[0])/motif_size)) est en
    dehors du if qui protège le calcul de moy[] -- elle s'exécute donc
    à CHAQUE locus, y compris avec des valeurs de moy[] périmées (d'un
    autre locus), divisées par le motif_size du locus COURANT. `nl`,
    lui, n'est incrémenté que quand les deux échantillons sont
    présentes -- numérateur et dénominateur sont donc désynchronisés
    dès qu'un tel locus existe.

    Ce n'est pas un choix statistique voulu (aucune justification
    biologique à réutiliser le delta-mu d'un locus différent) mais un
    bug de portée de variable en C++, de la même famille que le bug
    mutsit/sitefix déjà documenté pour sample_site_rates (header.cpp)
    -- confirmé avec l'utilisateur le 2026-09-11. Reproduit ici pour
    coller bit à bit à la sortie réelle de DIYABC ; la version "voulue"
    (sauter proprement le locus invalide) est gardée en commentaire
    ci-dessous, au cas où une future décision serait de corriger plutôt
    que reproduire ce comportement.

    Args:
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.
        motif_size: Taille du motif pour CE locus.
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}
            pour CE locus.
        raw_sizes: Dict {nom_echantillon: somme des tailles brutes} pour
            CE locus (sortie de _compute_VAR_constants).
        total_counts: Dict {nom_echantillon: nombre total de copies de
            gène} pour CE locus (sortie de _compute_VAR_constants).
        previous_moy: Le (moy_a, moy_b) du dernier locus valide
            rencontré pour cette paire, ou None si aucun locus valide
            n'a encore été rencontré (tout premier locus du groupe --
            cas non observé sur nos datasets, où le premier locus a
            toujours les deux échantillons présents).

    Returns:
        Tuple (contribution, new_moy, was_valid) :
            contribution: le terme à ajouter à la somme dmu2 pour ce
                locus (0.0 si previous_moy vaut encore None).
            new_moy: (moy_a, moy_b) mis à jour -- recalculé si les deux
                échantillons sont présentes à ce locus, sinon identique
                à previous_moy (le bug reproduit).
            was_valid: True si les deux échantillons étaient présentes à
                ce locus (donc si ce locus doit compter dans nl).
    """
    both_present = samp_a in length_by_sample and samp_b in length_by_sample

    if both_present:
        moy_a = raw_sizes[samp_a] / total_counts[samp_a]
        moy_b = raw_sizes[samp_b] / total_counts[samp_b]
        new_moy = (moy_a, moy_b)
    else:
        # BUG REPRODUIT (cal_dmu2p) : moy[] n'est pas recalculé, on
        # réutilise les valeurs périmées du dernier locus valide.
        # Version "voulue" (non utilisée ici) :
        #     return 0.0, previous_moy, False
        new_moy = previous_moy

    if new_moy is None:
        return 0.0, new_moy, both_present

    moy_a, moy_b = new_moy
    contribution = ((moy_b - moy_a) / motif_size) ** 2
    return contribution, new_moy, both_present


def compute_DM2(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    list_motif_sizes: list[int],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
) -> dict[str, float]:
    """Calcule DM2_i_j : (delta mu)^2 de Goldstein et al.

    (1995), pour chaque paire d'échantillons, sur tous les loci du groupe passé en
    argument.

    Reproduit fidèlement le bug d'accumulation de cal_dmu2p -- voir le
    docstring de _compute_DM2_for_one_locus. L'ORDRE d'itération sur
    les loci compte ici, contrairement aux autres stats microsat de ce
    fichier : chaque locus peut réutiliser l'état (moy) du précédent,
    donc `tree_sequences`/`list_motif_sizes` doivent être dans l'ordre
    réel du groupe (celui du header), pas un ordre arbitraire.

    Args:
        tree_sequences: Liste de TreeSequences, DANS L'ORDRE du groupe.
        sample_names: Liste des noms d'échantillon.
        list_motif_sizes: Liste des tailles de motifs, un par locus,
            dans le même ordre que tree_sequences.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dictionnaires {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus.

    Returns:
        Dict {"i.j": DM2}.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]
    dmu2_sum = {f"{i + 1}.{j + 1}": 0.0 for i, j in pairs}
    valid_loci_count = {f"{i + 1}.{j + 1}": 0 for i, j in pairs}
    previous_moy: dict[str, tuple[float, float] | None] = {
        f"{i + 1}.{j + 1}": None for i, j in pairs
    }
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, motif_size, length_by_sample in zip(
        tree_sequences, list_motif_sizes, length_by_sample_per_locus, strict=True
    ):
        raw_sizes, _, total_counts = _compute_VAR_constants(length_by_sample)
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            samp_a, samp_b = sample_names[i], sample_names[j]
            contribution, new_moy, was_valid = _compute_DM2_for_one_locus(
                samp_a,
                samp_b,
                motif_size,
                length_by_sample,
                raw_sizes,
                total_counts,
                previous_moy[key],
            )
            dmu2_sum[key] += contribution
            previous_moy[key] = new_moy
            if was_valid:
                valid_loci_count[key] += 1

    return {
        key: dmu2_sum[key] / valid_loci_count[key] if valid_loci_count[key] > 0 else 0.0
        for key in dmu2_sum
    }


# FST : between two samples (Weir and Cockerham 1984)


def _length_by_sample_and_individuals(
    tree_sequence: tskit.TreeSequence,
    *,
    layout: list[tuple[str, np.ndarray]] | None = None,
) -> dict[str, list[tuple[int, int]]]:
    """Calcule la longueur des séquences pour chaque individupar échantillon.
    La ploidie de l'individu est détectée par le nombre de noeud dans l'arbre
    via tree_sequence.individuals(). On retournera à chaque fois un tuple
    (longueur_1,longueur_2) pour chaque individu et longueur_1 sera répétée si
    l'individu est haploïde.

    Args:
        tree_sequence: Un objet TreeSequence de tskit.
        layout: [(nom_echantillon, np.ndarray[indices_d'individus]), ...] : si fourni, remplace le découpage par échantillon à utilisé par le chemin sériel, où un échantillon n'est plus sa propre échantillon

    Returns:
        Dict {nom_echantillon: [(longueur_1, longueur_2), ...]}.
    """
    if layout is None:
        layout = compute_population_layout(tree_sequence)
    length_by_sample = {pop: [] for pop, _ in layout}
    if tree_sequence.num_sites == 0:
        for ind in tree_sequence.individuals():
            nodes = ind.nodes
            échantillon = next(pop for pop, inds in layout if nodes[0] in inds)
            length_by_sample[échantillon].append((0, 0))
        return length_by_sample

    variant = next(tree_sequence.variants())
    tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
    for ind in tree_sequence.individuals():
        nodes = ind.nodes
        échantillon = next(pop for pop, inds in layout if nodes[0] in inds)
        if len(nodes) == 1:
            length_by_sample[échantillon].append(
                (int(tailles[nodes[0]]), int(tailles[nodes[0]]))
            )
        else:
            length_by_sample[échantillon].append(
                (int(tailles[nodes[0]]), int(tailles[nodes[1]]))
            )

    return length_by_sample


def _compute_ni_nA_AA_for_one_sample(
    pairs: list[tuple[int, int]], al: int
) -> tuple[int, int, int]:
    """Calcule ni, nA et AA pour un échantillon donné à partir des paires
    d'allèles.

    Args:
        pairs: Liste de tuples (longueur_1, longueur_2) pour chaque individu.
        al: Longueur de l'allèle considéré.

    Returns:
        Un tuple (ni, nA, AA) où
            - ni est le nombre d'individus,
            - nA est la somme par individus du nombre d'éléments de sa paire égaux
            - AA est le nb d'individus dont les deux éléments d ela paire valent al
    """
    ni = len(pairs)
    nA = sum(int(p[0] == al) + int(p[1] == al) for p in pairs)
    AA = sum(1 for p in pairs if p[0] == al and p[1] == al)
    return ni, nA, AA


def _compute_FST_constants_for_two_samples_combined(
    pairs_1: list[tuple[int, int]], pairs_2: list[tuple[int, int]], al: int
) -> tuple[int, int, int]:
    """Calcule les constantes nécessaires pour FST pour une paire
    d'échantillons combinés.

    Args:
        pairs_1: Liste de tuples (longueur_1, longueur_2) pour la première échantillon.
        pairs_2: Liste de tuples (longueur_1, longueur_2) pour la seconde échantillon.
        al: Longueur de l'allèle considéré.

    Returns:
        Un tuple (s2G,s2I, s2P)
    """
    ni_1, nA_1, AA_1 = _compute_ni_nA_AA_for_one_sample(pairs_1, al)
    ni_2, nA_2, AA_2 = _compute_ni_nA_AA_for_one_sample(pairs_2, al)

    sni = ni_1 + ni_2
    sni2 = ni_1**2 + ni_2**2
    sniA = nA_1 + nA_2
    sniAA = AA_1 + AA_2
    s2A = nA_1**2 / (2 * ni_1) + nA_2**2 / (2 * ni_2) if ni_1 > 0 and ni_2 > 0 else 0.0

    nc = sni - (sni2 / sni) if sni > 0 else 0.0

    if (sni * nc) > 0:
        MSG = (0.5 * sniA - sniAA) / sni
        MSI = (0.5 * sniA + sniAA - s2A) / (sni - 2.0)
        MSP = s2A - 0.5 * sniA**2 / sni
        s2G = MSG
        s2I = 0.5 * (MSI - MSG)
        s2P = (MSP - MSI) / (2.0 * nc)
        return s2G, s2I, s2P
    else:
        return 0.0, 0.0, 0.0


def _compute_FST_constants_on_all_alleles_for_two_samples(
    length_by_sample: dict[str, list[tuple[int, int]]], samp_a: str, samp_b: str
) -> tuple[float, float, float]:
    """Calcule les constantes nécessaires pour FST pour une paire
    d'échantillons sur tous les allèles.

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur_1, longueur_2), ...]}.
        samp_a: Nom du premier échantillon.
        samp_b: Nom du second échantillon.

    Returns:
        Un tuple (s2G_total, s2I_total, s2P_total) pour la paire d'échantillons.
    """
    pairs_1 = length_by_sample.get(samp_a, [])
    pairs_2 = length_by_sample.get(samp_b, [])

    unique_alleles = set(length for pair in pairs_1 + pairs_2 for length in pair)

    s1l = 0.0
    s2l = 0.0
    s3l = 0.0

    for al in unique_alleles:
        s2G, s2I, s2P = _compute_FST_constants_for_two_samples_combined(
            pairs_1, pairs_2, al
        )
        s1l += s2P
        s2l += s2P + s2I
        s3l += s2P + s2I + s2G

    return s1l, s2l, s3l


def compute_FST(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotype_by_sample_per_locus: list[dict[str, list[tuple[int] | tuple[int, int]]]]
    | None = None,
) -> dict[str, float]:
    """Calcule FST_i_j : pour chaque paire d'échantillons, la moyenne de FST sur tous les loci du groupe passé en argument (un groupe = les TreeSequences des loci séquence d'un même `group Gx` du header).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: non utilisé pour l'instant ; accepté pour que le dispatch appelle toutes les statistiques de ce dict de la même façon.
        genotype_by_sample_per_locus: non utilisé pour l'instant ; accepté pour que le dispatch appelle toutes les statistiques de ce dict de la même façon.
    Returns:
        Dict {"i.j": FST}.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(i + 1, len(sample_names))
    ]
    s1 = {f"{i + 1}.{j + 1}": 0.0 for i, j in pairs}
    # s2_sum = {f"{i + 1}.{j + 1}": 0.0 for i, j in pairs}
    s3 = {f"{i + 1}.{j + 1}": 0.0 for i, j in pairs}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    for ts, layout in zip(tree_sequences, layouts, strict=True):
        length_by_sample = _length_by_sample_and_individuals(ts, layout=layout)
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            samp_a, samp_b = sample_names[i], sample_names[j]
            # calcul de nc
            ni_1 = len(length_by_sample.get(samp_a, []))
            ni_2 = len(length_by_sample.get(samp_b, []))
            sni = ni_1 + ni_2
            sni2 = ni_1**2 + ni_2**2
            nc = sni - sni2 / sni if sni > 0 else 0.0
            # mise à jour des constantes
            s1l, _, s3l = _compute_FST_constants_on_all_alleles_for_two_samples(
                length_by_sample, samp_a, samp_b
            )
            s1[key] += s1l * nc
            s3[key] += s3l * nc
    return {key: s1[key] / s3[key] if s3[key] > 0 else 0.0 for key in s1}


# [LIK] - mean index of classification (two samples) (Rannala and Moutain 1997; Pascual et al. 2007)


def _genotypes_by_sample_and_individuals(
    tree_sequence: tskit.TreeSequence,
    *,
    layout: list[tuple[str, np.ndarray]] | None = None,
) -> dict[str, list[tuple[int, int] | tuple[int]]]:
    """Calcule les tailles d'allèles de chaque individu, par échantillon.

    Contrairement à _length_by_sample_and_individuals (qui duplique la
    valeur des individus haploïdes en une paire), on garde ici la
    VRAIE ploïdie -- un tuple à 1 élément pour un individu haploïde, à
    2 éléments pour un diploïde -- puisque cal_lik2p applique une
    formule différente selon le cas, pas une formule unique qui se
    prête à la duplication.

    Args:
        tree_sequence: Un objet TreeSequence de tskit.
        layout: [(nom_echantillon, np.ndarray[indices_d'individus]), ...] : si fourni, remplace le découpage par échantillon à utilisé par le chemin sériel, où un échantillon n'est plus sa propre échantillon

    Returns:
        Dict {nom_echantillon: [(longueur_1, longueur_2) ou (longueur_1,), ...]}.
    """
    if layout is None:
        layout = compute_population_layout(tree_sequence)
    genotype_by_sample = {pop: [] for pop, _ in layout}
    if tree_sequence.num_sites == 0:
        for ind in tree_sequence.individuals():
            nodes = ind.nodes
            échantillon = next(pop for pop, inds in layout if nodes[0] in inds)
            genotype_by_sample[échantillon].append((0, 0))
        return genotype_by_sample
    else:
        variant = next(tree_sequence.variants())
        tailles = np.array([int(a) for a in variant.alleles])[variant.genotypes]
        for ind in tree_sequence.individuals():
            nodes = ind.nodes
            échantillon = next(pop for pop, inds in layout if nodes[0] in inds)
            if len(nodes) == 1:
                genotype_by_sample[échantillon].append((tailles[nodes[0]],))
            else:
                genotype_by_sample[échantillon].append(
                    (tailles[nodes[0]], tailles[nodes[1]])
                )

        return genotype_by_sample


def _compute_num_den_lik_for_one_individual(
    genotype: tuple[int] | tuple[int, int],
    count_j: dict[int, int],
    total_count_j: int,
    b: float,
) -> tuple[float, float]:
    """Calcule (num_lik, den_lik) de LIK pour un individu, selon sa ploïdie.

    Reproduit cal_lik2p (sumstat.cpp) : trois formules distinctes selon
    que l'individu est haploïde, diploïde homozygote ou diploïde
    hétérozygote -- pas une formule unique applicable aux trois cas
    (contrairement à FST, voir _genotypes_by_sample_and_individuals).

    Args:
        genotype: Tailles d'allèles de l'individu -- un tuple à 1
            élément s'il est haploïde, à 2 s'il est diploïde.
        count_j: Dict {taille: compte} pour l'échantillon de référence
            (samp_j, celle dont on utilise les fréquences).
        total_count_j: Nombre total de copies de gène dans samp_j.
        b: Pseudo-compte (1/nal, nal = nombre d'allèles distincts dans
            le dataset à ce locus).

    Returns:
        Tuple (num_lik, den_lik) pour cet individu.
    """

    if len(genotype) == 1:  # haploid
        num_lik = count_j.get(genotype[0], 0) + b
        den_lik = total_count_j + 1
    elif len(genotype) == 2 and genotype[0] == genotype[1]:  # homozygous diploid
        num_lik = (1 + b + count_j.get(genotype[0], 0)) * (
            count_j.get(genotype[0], 0) + b
        )
        den_lik = (total_count_j + 2) * (total_count_j + 1)
    elif len(genotype) == 2 and genotype[0] != genotype[1]:  # heterozygous diploid
        num_lik = (
            2 * (b + count_j.get(genotype[0], 0)) * (b + count_j.get(genotype[1], 0))
        )
        den_lik = (total_count_j + 2) * (total_count_j + 1)
    else:
        raise ValueError(f"Problème avec le génotype de l'individu : {genotype}")

    return num_lik, den_lik


def _compute_LIK_for_one_locus(
    length_by_sample: dict[str, list[tuple[int, int]]],
    genotypes_by_sample: dict[str, list[tuple[int, ...]]],
    samp_i: str,
    samp_j: str,
) -> tuple[float, bool]:
    """Calcule la contribution de LIK_i_j à UN locus (sens i -> j uniquement).

    Teste les individus de samp_i (échantillon testé) contre les
    fréquences de samp_j (échantillon de référence) -- ASYMÉTRIQUE,
    contrairement à toutes les autres stats microsat à deux
    échantillons : LIK_i_j et LIK_j_i utilisent des rôles inversés et ne
    sont pas censées être égales. `nal`/`b` sont calculés sur l'union
    des allèles de samp_i et samp_j présents à ce locus (comme il n'y a
    que 2 échantillons dans ce projet, ça correspond exactement à
    `nal` du C++, qui pool en théorie sur TOUTES les échantillons du
    dataset).

    Args:
        length_by_sample: Dict {nom_echantillon: [(longueur, nb_sequence), ...]}
            pour CE locus (sortie de _length_by_sample_and_individuals).
        genotypes_by_sample: Dict {nom_echantillon: [génotype par individu, ...]}
            pour CE locus (sortie de _genotypes_by_sample_and_individuals).
        samp_i: Population testée.
        samp_j: Population de référence (dont on utilise les fréquences).

    Returns:
        Tuple (contribution, was_valid) : la contribution de ce locus
        (déjà multipliée par a = 1/nombre d'individus de samp_i), et
        True si samp_i ET samp_j ont des données à ce locus (sinon ce
        locus doit être exclu, contribution=0.0).
    """
    if samp_i not in length_by_sample or samp_j not in length_by_sample:
        return 0.0, False

    count_j = {length: count for length, count in length_by_sample[samp_j]}
    total_count_j = sum(count_j.values())

    # nal du C++ (cal_lik2p) : un allèle n'est compté que si sa fréquence
    # SOMMÉE SUR TOUS LES ÉCHANTILLONS est non nulle. `length_by_sample` porte
    # déjà toutes les échantillons, et ses clés viennent de variant.alleles --
    # donc elles incluent des états créés par une mutation puis écrasés, que
    # plus aucun échantillon ne porte. Les compter gonfle nal et écrase b.
    total_counts = {}
    for rows in length_by_sample.values():
        for length, count in rows:
            total_counts[length] = total_counts.get(length, 0) + count
    nb_allele = sum(1 for count in total_counts.values() if count > 0)
    b = 1 / nb_allele if nb_allele > 0 else 0.0

    a = (
        1 / len(genotypes_by_sample[samp_i])
        if len(genotypes_by_sample[samp_i]) > 0
        else 0.0
    )

    likelihood = 0

    for genotype in genotypes_by_sample[samp_i]:
        num_lik, den_lik = _compute_num_den_lik_for_one_individual(
            genotype, count_j, total_count_j, b
        )
        likelihood -= np.log10(num_lik / den_lik)

    return likelihood * a, True


def compute_LIK(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotypes_by_sample_per_locus: list[dict[str, list[tuple[int, ...]]]] | None = None,
) -> dict[str, float]:
    """Calcule LIK_i_j : indice de vraisemblance d'assignation (Rannala &
    Mountain 1997 ; Pascual et al. 2007), pour chaque paire ORDONNÉE de
    échantillons, moyenné sur tous les loci du groupe.

    Contrairement à NAL/HET/VAR/N2P/H2P/V2P/DAS/DM2/FST (paires non
    ordonnées i<j), LIK est asymétrique : `pairs` couvre toutes les
    paires ordonnées i != j (donc "i.j" ET "j.i" pour chaque
    combinaison), pas seulement i<j. Moyenne par locus (comme
    NAL/HET/VAR), pas un ratio de sommes (comme MGW/DAS/FST).

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dicts {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.
        genotypes_by_sample_per_locus: Liste de dicts {nom_echantillon: [(génome, ...), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.

    Returns:
        Dict {"i.j": LIK}, une entrée par paire ORDONNÉE d'échantillons.
    """
    pairs = [
        (i, j)
        for i in range(len(sample_names))
        for j in range(len(sample_names))
        if i != j
    ]
    likelihood_sum = {f"{i + 1}.{j + 1}": 0.0 for i, j in pairs}
    valid_loci_count = {f"{i + 1}.{j + 1}": 0 for i, j in pairs}
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    if genotypes_by_sample_per_locus is None:
        genotypes_by_sample_per_locus = [
            _genotypes_by_sample_and_individuals(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, length_by_sample, genotypes_by_sample in zip(
        tree_sequences,
        length_by_sample_per_locus,
        genotypes_by_sample_per_locus,
        strict=True,
    ):
        for i, j in pairs:
            key = f"{i + 1}.{j + 1}"
            samp_i, samp_j = sample_names[i], sample_names[j]
            likelihood, both_present = _compute_LIK_for_one_locus(
                length_by_sample, genotypes_by_sample, samp_i, samp_j
            )
            likelihood_sum[key] += likelihood
            if both_present:
                valid_loci_count[key] += 1
    return {
        key: likelihood_sum[key] / valid_loci_count[key]
        if valid_loci_count[key] > 0
        else 0.0
        for key in likelihood_sum
    }


# Stat AML - Maximum likelihood coefficient of admixture (Choisy et al. 2004)


def _extract_mutual_data(
    tree_sequences: list[tskit.TreeSequence],
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotypes_by_sample_per_locus: list[dict[str, list[tuple[int, ...]]]] | None = None,
) -> list[tuple[dict[str, list[tuple[int, int]]], dict[str, list[tuple[int, ...]]]]]:
    """Prépare les données nécessaires pour le calcul de la log-vraisemblance d'admixture.

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dicts {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.
        genotypes_by_sample_per_locus: Liste de dicts {nom_echantillon: [(génome, ...), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.
    Returns:
        Une liste de tuples pour chaque locus, contenant :
        - dict des longueurs par échantillon,
        - dict des génotypes par échantillon.
    """
    mutual_data = []
    if layouts is None:
        layouts = [None] * len(tree_sequences)
    if length_by_sample_per_locus is None:
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    if genotypes_by_sample_per_locus is None:
        genotypes_by_sample_per_locus = [
            _genotypes_by_sample_and_individuals(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts, strict=True)
        ]
    for _, length_by_sample, genotypes_by_sample in zip(
        tree_sequences,
        length_by_sample_per_locus,
        genotypes_by_sample_per_locus,
        strict=True,
    ):
        mutual_data.append((length_by_sample, genotypes_by_sample))

    return mutual_data


def _prepare_triplet_for_admixture(
    mutual_data: list[
        tuple[
            dict[str, list[tuple[int, int]]],
            dict[str, list[tuple[int, int] | tuple[int]]],
        ]
    ],
    focal: str,
    parent1: str,
    parent2: str,
) -> list[tuple[dict[int, float], dict[int, float], list[tuple[int, ...]]]]:
    """Prépare les loci pour le calcul de la log-vraisemblance d'admixture.
    Pour chaque locus, on calcule les fréquences des allèles dans les
    échantillons parentaux et on récupère les génotypes des individus de
    l'échantillon focal.

    Args:
        mutual_data: Liste de tuples contenant les données mutuelles pour chaque locus.
        focal: Nom de l'échantillon focal.
        parent1: Nom du premier échantillon parental.
        parent2: Nom de la deuxième échantillon parental.


    Returns:
        Une liste de tuples pour chaque locus, contenant :
        - dict des fréquences des allèles dans parent1,
        - dict des fréquences des allèles dans parent2,
        - liste des génotypes des individus de l'échantillon focal.
    """
    prepared_triplet = []
    # Un layout par locus, ou None partout si l'appelant n'en fournit pas.
    for length_by_sample, genotypes_by_sample in mutual_data:
        count_parent1 = {
            length: count
            for length, count in length_by_sample.get(parent1, [])
            if count > 0
        }
        count_parent2 = {
            length: count
            for length, count in length_by_sample.get(parent2, [])
            if count > 0
        }
        total_parent1 = sum(count_parent1.values())
        total_parent2 = sum(count_parent2.values())

        f1 = (
            {length: nb / total_parent1 for length, nb in count_parent1.items()}
            if total_parent1 > 0
            else {}
        )
        f2 = (
            {length: nb / total_parent2 for length, nb in count_parent2.items()}
            if total_parent2 > 0
            else {}
        )
        focal_genotypes = genotypes_by_sample.get(focal, [])

        prepared_triplet.append((f1, f2, focal_genotypes))

    return prepared_triplet


def _vectorize_triplet(prepared_triplet):
    # listes vides pour chaque groupe (haploÃ¯de : 2 listes, homozygote : 2, hétérozygote : 4)
    f1_haploid, f2_haploid = [], []
    f1_homozygous, f2_homozygous = [], []
    f1_heterozygous_x, f2_heterozygous_x = [], []
    f1_heterozygous_y, f2_heterozygous_y = [], []
    for f1, f2, focal_genotypes in prepared_triplet:
        for genotype in focal_genotypes:
            if len(genotype) == 1:  # haploid
                f1_haploid.append(f1.get(genotype[0], 0))
                f2_haploid.append(f2.get(genotype[0], 0))
            elif (
                len(genotype) == 2 and genotype[0] == genotype[1]
            ):  # homozygous diploid
                f1_homozygous.append(f1.get(genotype[0], 0))
                f2_homozygous.append(f2.get(genotype[0], 0))
            elif (
                len(genotype) == 2 and genotype[0] != genotype[1]
            ):  # heterozygous diploid
                f1_heterozygous_x.append(f1.get(genotype[0], 0))
                f2_heterozygous_x.append(f2.get(genotype[0], 0))
                f1_heterozygous_y.append(f1.get(genotype[1], 0))
                f2_heterozygous_y.append(f2.get(genotype[1], 0))
            else:
                raise ValueError(
                    f"Problème avec le génotype de l'individu : {genotype}"
                )
    # conversion en np.array(dtype=float) et retour
    f1_haploid = np.array(f1_haploid, dtype=float)
    f2_haploid = np.array(f2_haploid, dtype=float)
    f1_homozygous = np.array(f1_homozygous, dtype=float)
    f2_homozygous = np.array(f2_homozygous, dtype=float)
    f1_heterozygous_x = np.array(f1_heterozygous_x, dtype=float)
    f2_heterozygous_x = np.array(f2_heterozygous_x, dtype=float)
    f1_heterozygous_y = np.array(f1_heterozygous_y, dtype=float)
    f2_heterozygous_y = np.array(f2_heterozygous_y, dtype=float)

    focal_genotypes = FocalGenotypeFrequencies(
        haploid_f1=f1_haploid,
        haploid_f2=f2_haploid,
        homozygous_f1=f1_homozygous,
        homozygous_f2=f2_homozygous,
        heterozygous_f1_x=f1_heterozygous_x,
        heterozygous_f2_x=f2_heterozygous_x,
        heterozygous_f1_y=f1_heterozygous_y,
        heterozygous_f2_y=f2_heterozygous_y,
    )

    return focal_genotypes


def _log_likelihood_admixture(
    data: FocalGenotypeFrequencies,
    a: float,
) -> float:
    """Un seul float, la log-vraisemblance totale (somme sur tous les loci du
    groupe, somme sur tous les individus de focal) : équivalent de li[rep] du
    C++ pour UN a donné, pas encore le couple (li0, delta).

    Args:
        data: dataclass FocalGenotypeFrequencies contenant les fréquences des allèles pour les individus de l'échantillon focal.
        a: Valeur du coefficient d'admixture.

    Returns:
        La valeur de la statistique.
    """
    lik = 0.0

    # haploïd
    freq = a * data.haploid_f1 + (1 - a) * data.haploid_f2
    mask = freq > 0
    lik += np.log(freq[mask]).sum()

    # homozygous diploid
    freq = a * data.homozygous_f1 + (1 - a) * data.homozygous_f2
    mask = freq > 0
    lik += np.log(freq[mask] ** 2).sum()

    # heterozygous diploid
    freq_x = a * data.heterozygous_f1_x + (1 - a) * data.heterozygous_f2_x
    freq_y = a * data.heterozygous_f1_y + (1 - a) * data.heterozygous_f2_y
    mask = freq_x * freq_y > 0
    lik += np.log(2 * freq_x[mask] * freq_y[mask]).sum()

    return float(lik)


def _pente_lik(
    data: FocalGenotypeFrequencies,
    i0: int,
) -> tuple[float, float]:
    """
    Sortie : (li0, delta) où li0 = _log_likelihood_admixture_mixture(..., a=0.001*i0) et delta = li(a=0.001*(i0+1)) - li0  :
    exactement le (li[0], li[1]-li[0]) que pente_lik retourne en C++.
    Args:
        data: Dataclass FocalGenotypeFrequencies contenant les fréquences des allèles pour les individus de l'échantillon focal.
        i0: Indice du point d'évaluation pour le calcul de la pente.
    Returns:
        La valeur de la statistique et la pente.
    """
    lik_a = _log_likelihood_admixture(data, a=0.001 * i0)
    lik_a_plus = _log_likelihood_admixture(data, a=0.001 * (i0 + 1))
    delta = lik_a_plus - lik_a

    return lik_a, delta


def _compute_AML_one_triplet(
    mutual_data: list[
        tuple[dict[str, list[tuple[int, int]]], dict[str, list[tuple[int, ...]]]]
    ],
    focal: str,
    parent1: str,
    parent2: str,
    seed: int,
) -> float:
    """Calcule le coefficient d'admixture maximum de vraisemblance (AML) pour
    un triplet d'échantillons.

    Args:
        mutual_data: Liste de tuples contenant les données mutuelles pour chaque locus.
        focal: Nom de l'échantillon focal.
        parent1: Nom du premier échantillon parental.
        parent2: Nom de la deuxième échantillon parental.
        seed: Seed pour la génération aléatoire (pour les cas où AML ne peut pas être calculé).

    Returns:
        La statistique AML pour le triplet donné.
    """
    prepared_triplet = _prepare_triplet_for_admixture(
        mutual_data, focal, parent1, parent2
    )

    data = _vectorize_triplet(prepared_triplet)

    i1, i2 = 1, 998
    lik1, p1 = _pente_lik(data, i1)
    lik2, p2 = _pente_lik(data, i2)

    if abs(lik1) + abs(lik2) < 1e-10:
        rng = random.Random(seed)
        return rng.uniform(0, 1)
    elif p1 < 0 and p2 < 0:
        return 0.0
    elif p1 > 0 and p2 > 0:
        return 1.0
    else:
        while i2 - i1 > 1:
            i3 = (i1 + i2) // 2
            lik3, p3 = _pente_lik(data, i3)
            if p1 * p3 < 0:
                i2 = i3
                p2 = p3
                lik2 = lik3
            else:
                i1 = i3
                p1 = p3
                lik1 = lik3

    if lik1 > lik2:
        return 0.001 * i1
    else:
        return 0.001 * i2


def compute_AML_microsat(
    tree_sequences: list[tskit.TreeSequence],
    sample_names: list[str],
    seed: int = 0,
    *,
    layouts: list[list[tuple[str, np.ndarray]]] | None = None,
    length_by_sample_per_locus: list[dict[str, list[tuple[int, int]]]] | None = None,
    genotypes_by_sample_per_locus: list[dict[str, list[tuple[int, ...]]]] | None = None,
) -> dict[str, float]:
    """Calcule le coefficient d'admixture maximum de vraisemblance (AML) pour
    chaque triplet d'échantillons.

    Args:
        tree_sequences: Liste de TreeSequences (un arbre par locus).
        sample_names: Liste des noms d'échantillon.
        seed: Seed pour la génération aléatoire (pour les cas où AML ne peut pas être calculé).
        layouts: Liste des layouts, un par locus, ou None partout si l'appelant n'en fournit pas.
        length_by_sample_per_locus: Liste de dicts {nom_echantillon: [(longueur, nb_sequence), ...]} pour chaque locus, ou None si l'appelant n'en fournit pas.

    Returns:
        Dict {"i.j.k": AML} pour chaque triplet d'échantillons.
    """
    n_sample = len(sample_names)
    results = {}
    mutual_data = _extract_mutual_data(
        tree_sequences,
        layouts=layouts,
        length_by_sample_per_locus=length_by_sample_per_locus,
        genotypes_by_sample_per_locus=genotypes_by_sample_per_locus,
    )
    for i, t in enumerate(_half_arrangements(n_sample, 3)):
        h, p1, p2 = t[0], t[1], t[2]
        key = f"{h + 1}.{p1 + 1}.{p2 + 1}"
        focal, parent1, parent2 = (
            sample_names[h],
            sample_names[p1],
            sample_names[p2],
        )
        results[key] = _compute_AML_one_triplet(
            mutual_data,
            focal,
            parent1,
            parent2,
            _locus_seed(seed, _LIKELIHOOD_SEED_OFFSET, i),
        )
    return results


# ---------------------------------------------------------------------------
# Point d'entrée principal
# ---------------------------------------------------------------------------


def compute_all_statistics_indseq(
    genotypes_per_locus: list[dict[str, list[int]]],
    sample_names: list[str],
) -> dict[str, float]:
    """Calcule les 130 statistiques résumées SNP (IndSeq).

    Les matrices (n_sample x n_loci) de comptes et fréquences sont construites
    une seule fois (_prepare_matrices) et transmises à toutes les familles
    de statistiques via _mats.

    Args:
        genotypes_per_locus: Liste de dicts {nom_echantillon:
            [génotype, ...]}, un dict par locus.
        sample_names: Les noms d'échantillon.

    Returns:
        Un dict {nom_stat: valeur} -- même format que parse_statobs().
    """
    mats = _prepare_matrices(genotypes_per_locus, sample_names)
    results = {}
    results.update(compute_ML1(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_ML2(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_ML3(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_HW_HB(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_FST1(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_FST2(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_NEI(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_AML(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_F3(genotypes_per_locus, sample_names, _mats=mats))
    results.update(compute_F4(genotypes_per_locus, sample_names, _mats=mats))
    results.update(
        compute_FST3_FST4_FSTG(genotypes_per_locus, sample_names, _mats=mats)
    )
    return results


def compute_all_statistics_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    sample_names: list[str],
    pool_sizes: dict[str, int],
) -> dict[str, float]:
    """Calcule les statistiques résumées SNP pour PoolSeq.

    Les matrices (n_sample x n_loci) de comptes et tailles d'échantillon sont
    construites une seule fois (_prepare_matrices_poolseq) et transmises
    à toutes les familles de statistiques via _mats.

    Args:
        reads_per_locus: Liste de dicts {nom_echantillon: (nreads_dérivé,
            nreads_total)}, un dict par locus.
        sample_names: Les noms d'échantillon.
        pool_sizes: Dict {nom_echantillon: taille_haploïde du pool}.

    Returns:
        Un dict {nom_stat: valeur} -- même format que parse_statobs().
    """
    mats = _prepare_matrices_poolseq(reads_per_locus, sample_names)

    results = {}
    results.update(compute_ML1(None, sample_names, _mats=mats))
    results.update(compute_ML2(None, sample_names, _mats=mats))
    results.update(compute_ML3(None, sample_names, _mats=mats))
    results.update(
        compute_HW_HB_poolseq(
            reads_per_locus, sample_names, pool_sizes=pool_sizes, _mats=mats
        )
    )
    results.update(compute_NEI(reads_per_locus, sample_names, _mats=mats))
    results.update(compute_F4(reads_per_locus, sample_names, _mats=mats))
    results.update(
        compute_F3_poolseq(reads_per_locus, sample_names, pool_sizes, _mats=mats)
    )
    results.update(
        compute_FST2_poolseq(reads_per_locus, sample_names, pool_sizes, _mats=mats)
    )
    results.update(
        compute_FST3_FST4_poolseq(reads_per_locus, sample_names, pool_sizes, _mats=mats)
    )
    results.update(
        compute_FST1_poolseq(reads_per_locus, sample_names, pool_sizes, _mats=mats)
    )
    results.update(compute_AML(reads_per_locus, sample_names, _mats=mats))
    return results


_DNA_PER_SAMPLE_STATS = {
    "NHA": compute_NHA,
    "NSS": compute_NSS,
    "MPD": compute_MPD,
    "VPD": compute_VPD,
    "DTA": compute_DTA,
    "PSS": compute_PSS,
    "MNS": compute_MNS,
    "VNS": compute_VNS,
}

_DNA_PAIRWISE_STATS = {
    "NH2": compute_NH2,
    "NS2": compute_NS2,
    "MP2": compute_MP2,
    "MPB": compute_MPB,
    "HST": compute_HST,
}


def compute_all_statistics_dna(
    header_text: str,
    tree_sequences_by_locus: dict[str, tskit.TreeSequence],
    sample_names: list[str],
    *,
    seed: int = 0,
    layouts_by_locus: dict[str, list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule les 13 statistiques résumées ADN pour chaque `group Gx` séquence
    (`[S]`) du header, et retourne un dict {nom_colonne: valeur} utilisant les
    VRAIS noms de colonnes DIYABC (`STAT_<groupe>_<pop-ou- paire>`, ex.
    `NSS_2_1`, `NH2_3_1.2`) -- vérifié caractère pour caractère contre la
    sortie réelle de `diyabc` sur `toy_example2_ms_dna`
    (`STAT_<groupe>_<suffixe>` quand il y a plusieurs groupes, `STAT_<suffixe>`
    seul sinon -- même convention que
    `stats_group_parser.parse_requested_statistic_names`).

    Args:
        header_text: contenu de header.txt/headerRF.txt (pour
            `parse_loci_description`, qui donne le groupe de chaque
            locus).
        tree_sequences_by_locus: {nom_locus: TreeSequence mutée} --
            la sortie de `dna_mutation_simulation_per_locus`. Les loci
            microsat (`ms_or_seq == "M"`) présents dans le header sont
            ignorés ici (pas de code de simulation microsat).
        sample_names: toutes les échantillons du dataset, dans
            l'ordre "pop1".."popN" (leur position dans cette liste,
            pas leur nom, détermine l'indice numérique utilisé dans
            les noms de colonnes).

    Returns:
        Un dict {nom_colonne_diyabc: valeur}.
    """
    loci_by_group: dict[str, list[str]] = {}
    list_loci = parse_loci_description(header_text)
    for locus in list_loci:
        if locus.ms_or_seq != "S":
            continue
        loci_by_group.setdefault(locus.group, []).append(locus.name)

    all_groups = {locus.group for locus in list_loci}
    multi_group = len(all_groups) > 1

    results = {}
    for group_label, locus_names in loci_by_group.items():
        group_number = group_label[1:]  # "G2" -> "2", comme stats_group_parser.py
        tree_sequences = [tree_sequences_by_locus[name] for name in locus_names]

        layouts = (
            None
            if layouts_by_locus is None
            else [layouts_by_locus[name] for name in locus_names]
        )
        for stat_name, stat_fn in _DNA_PER_SAMPLE_STATS.items():
            for samp_name, value in stat_fn(
                tree_sequences, sample_names, layouts=layouts
            ).items():
                samp_index = sample_names.index(samp_name) + 1
                key = (
                    f"{stat_name}_{group_number}_{samp_index}"
                    if multi_group
                    else f"{stat_name}_{samp_index}"
                )
                results[key] = value

        for stat_name, stat_fn in _DNA_PAIRWISE_STATS.items():
            for pair_key, value in stat_fn(
                tree_sequences, sample_names, layouts=layouts
            ).items():
                key = (
                    f"{stat_name}_{group_number}_{pair_key}"
                    if multi_group
                    else f"{stat_name}_{pair_key}"
                )
                results[key] = value

    return results


_MICROSAT_PER_SAMPLE_WITHOUT_MOTIF_SIZE = {
    "NAL": compute_NAL,
    "HET": compute_HET,
}

_MICROSAT_PAIRWISE_WITHOUT_MOTIF_SIZE = {
    "N2P": compute_N2P,
    "H2P": compute_H2P,
    "DAS": compute_DAS,
    "FST": compute_FST,
    "LIK": compute_LIK,
}

_MICROSAT_PER_SAMPLE_WITH_MOTIF_SIZE = {
    "VAR": compute_VAR,
    "MGW": compute_MGW,
}

_MICROSAT_PAIRWISE_WITH_MOTIF_SIZE = {
    "V2P": compute_V2P,
    "DM2": compute_DM2,
}

_MICROSAT_TRIPLET_STATS = {"AML": compute_AML_microsat}


def compute_all_statistics_microsat(
    header_text: str,
    tree_sequences_by_locus: dict[str, tskit.TreeSequence],
    sample_names: list[str],
    *,
    seed: int,
    layouts_by_locus: dict[str, list[tuple[str, np.ndarray]]] | None = None,
) -> dict[str, float]:
    """Calcule les statistiques résumées microsat pour chaque `group Gx`
    microsat (`[M]`) du header, et retourne un dict {nom_colonne: valeur}
    utilisant les VRAIS noms de colonnes DIYABC -- vérifié caractère pour
    caractère contre la sortie réelle de `diyabc` sur `toy_example2_ms_dna`
    (`STAT_<groupe>_<suffixe>` quand il y a plusieurs groupes, `STAT_<suffixe>`
    seul sinon -- même convention que
    `stats_group_parser.parse_requested_statistic_names`).

    Args:
        header_text: contenu de header.txt/headerRF.txt (pour
            `parse_loci_description`, qui donne le groupe de chaque
            locus).
        tree_sequences_by_locus: {nom_locus: TreeSequence mutée} --
            la sortie de `microsat_mutation_simulation_per_locus`. Les loci ADN (`ms_or_seq == "S"`) présents dans le header sont ignorés ici (pas de code de simulation ADN).
        sample_names: toutes les échantillons du dataset, dans
            l'ordre "pop1".."popN" (leur position dans cette liste,
            pas leur nom, détermine l'indice numérique utilisé dans
            les noms de colonnes).
    Returns:
        Un dict {nom_colonne_diyabc: valeur}.
    """
    loci_by_group: dict[str, list[str]] = {}
    motif_sizes_by_locus: dict[str, int] = {}
    list_loci = parse_loci_description(header_text)
    for locus in list_loci:
        if locus.ms_or_seq != "M":
            continue
        loci_by_group.setdefault(locus.group, []).append(locus.name)
        motif_sizes_by_locus[locus.name] = locus.motif_size

    all_groups = {locus.group for locus in list_loci}
    multi_group = len(all_groups) > 1

    results = {}

    for group_label, locus_names in loci_by_group.items():
        group_number = group_label[1:]  # "G2" -> "2", comme stats_group_parser.py
        tree_sequences = [tree_sequences_by_locus[name] for name in locus_names]
        layouts = (
            None
            if layouts_by_locus is None
            else [layouts_by_locus[name] for name in locus_names]
        )
        layouts_for_lengths = (
            layouts if layouts is not None else [None] * len(tree_sequences)
        )
        length_by_sample_per_locus = [
            _length_by_sample(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts_for_lengths, strict=True)
        ]

        genotypes_by_sample_per_locus = [
            _genotypes_by_sample_and_individuals(ts, layout=layout)
            for ts, layout in zip(tree_sequences, layouts_for_lengths, strict=True)
        ]
        for stat_name, stat_fn in _MICROSAT_PER_SAMPLE_WITHOUT_MOTIF_SIZE.items():
            for samp_name, value in stat_fn(
                tree_sequences,
                sample_names,
                layouts=layouts,
                length_by_sample_per_locus=length_by_sample_per_locus,
                genotypes_by_sample_per_locus=genotypes_by_sample_per_locus,
            ).items():
                samp_index = sample_names.index(samp_name) + 1
                key = (
                    f"{stat_name}_{group_number}_{samp_index}"
                    if multi_group
                    else f"{stat_name}_{samp_index}"
                )
                results[key] = value
        for stat_name, stat_fn in _MICROSAT_PAIRWISE_WITHOUT_MOTIF_SIZE.items():
            for stat_index, value in stat_fn(
                tree_sequences,
                sample_names,
                layouts=layouts,
                length_by_sample_per_locus=length_by_sample_per_locus,
            ).items():
                key = (
                    f"{stat_name}_{group_number}_{stat_index}"
                    if multi_group
                    else f"{stat_name}_{stat_index}"
                )
                results[key] = value

        motif_sizes = [motif_sizes_by_locus[name] for name in locus_names]
        for stat_name, stat_fn in _MICROSAT_PER_SAMPLE_WITH_MOTIF_SIZE.items():
            for samp_name, value in stat_fn(
                tree_sequences,
                sample_names,
                motif_sizes,
                layouts=layouts,
                length_by_sample_per_locus=length_by_sample_per_locus,
            ).items():
                samp_index = sample_names.index(samp_name) + 1
                key = (
                    f"{stat_name}_{group_number}_{samp_index}"
                    if multi_group
                    else f"{stat_name}_{samp_index}"
                )
                results[key] = value
        for stat_name, stat_fn in _MICROSAT_PAIRWISE_WITH_MOTIF_SIZE.items():
            for stat_index, value in stat_fn(
                tree_sequences,
                sample_names,
                motif_sizes,
                layouts=layouts,
                length_by_sample_per_locus=length_by_sample_per_locus,
            ).items():
                key = (
                    f"{stat_name}_{group_number}_{stat_index}"
                    if multi_group
                    else f"{stat_name}_{stat_index}"
                )
                results[key] = value

        for stat_name, stat_fn in _MICROSAT_TRIPLET_STATS.items():
            for stat_index, value in stat_fn(
                tree_sequences,
                sample_names,
                _locus_seed(seed, _GROUP_STAT_SEED_OFFSET, int(group_number)),
                layouts=layouts,
                length_by_sample_per_locus=length_by_sample_per_locus,
                genotypes_by_sample_per_locus=genotypes_by_sample_per_locus,
            ).items():  # pour éviter d'avoir la même graine pour différents groupes
                key = (
                    f"{stat_name}_{group_number}_{stat_index}"
                    if multi_group
                    else f"{stat_name}_{stat_index}"
                )
                results[key] = value
    return results
