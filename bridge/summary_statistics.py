"""
Implémentation Python des statistiques résumées SNP calculées par DIYABC
(HeaderC::calstatobs / ParticleC::docalstat, src-JMC-C++/sumstat.cpp).

PROTOCOLE DE VALIDATION : pour chaque formule implémentée ici, on vérifie
qu'elle produit les mêmes valeurs que le vrai binaire `general` sur les
MÊMES données en entrée -- comparaison exacte (à la précision float32
près), pas statistique.

Structure d'entrée attendue : liste de dicts {nom_population: [génotypes
haploïdes 0/1]}, un dict par locus -- la forme produite par
ancestry_simulation.simulate_snp_genotypes.

Organisation : une fonction par famille de statistiques, suivant
exactement la nomenclature de sumstat.cpp (cal_snfl, cal_snhw, cal_snhb,
cal_snfsti...) pour faciliter la traçabilité entre le code Python et sa
source C++ de référence.
"""

from itertools import combinations, permutations

import numpy as np
import tskit

from bridge.ancestry_simulation import compute_population_layout

# ---------------------------------------------------------------------------
# Utilitaires scalaires (conservés pour la traçabilité et les tests unitaires)
# ---------------------------------------------------------------------------


def _allele_freq(haploid_genotypes: list[int]) -> float:
    """Fréquence de l'allèle dérivé (1) dans une population -- équivalent
    de locuslist[loc].freq[pop][1] dans le code C++."""
    n = len(haploid_genotypes)
    if n == 0:
        return float("nan")
    return sum(haploid_genotypes) / n


def _q1(haploid_genotypes: list[int]) -> float:
    """Probabilité d'identité par état intra-population, tirage SANS
    remise -- formule exacte de sumstat.cpp::q1 (cas SNP, bias=False) :
        q1 = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1))
    où y1, y2 = comptes d'allèles 0 et 1 (= freq * n).
    """
    n = len(haploid_genotypes)
    if n <= 1:
        return float("nan")
    y2 = sum(haploid_genotypes)
    y1 = n - y2
    return (y1 * (y1 - 1) + y2 * (y2 - 1)) / (n * (n - 1))


def _q2(haploid_genotypes_a: list[int], haploid_genotypes_b: list[int]) -> float:
    """Probabilité d'identité par état inter-populations -- formule exacte
    de sumstat.cpp::q2 (cas SNP) :
        q2 = (y11*y21 + y12*y22) / (n1*n2)
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
    population_names: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construit les matrices (npop, nloci) de comptes et fréquences.

    Retourne (counts, ns, freq0, freq1) :
      counts[i, l] = nb d'allèles dérivés (1) dans pop i au locus l
      ns[i, l]     = nb total de lignées dans pop i au locus l
      freq1 = counts / ns,  freq0 = 1 - freq1

    Appelé UNE SEULE FOIS dans compute_all_statistics et transmis via _mats
    à toutes les familles de statistiques -- évite de reconstruire les
    matrices (npop × nloci) une fois par famille.
    """
    counts = np.array(
        [[sum(lg[p]) for lg in genotypes_per_locus] for p in population_names],
        dtype=float,
    )
    ns = np.array(
        [[len(lg[p]) for lg in genotypes_per_locus] for p in population_names],
        dtype=float,
    )
    freq1 = counts / ns
    freq0 = 1.0 - freq1
    return counts, ns, freq0, freq1


def _prepare_matrices_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Construit les matrices (npop, nloci) de comptes et tailles d'échantillon pour POOLSEQ.

    Retourne (counts, ns) :
      counts[i, l] = nb d'allèles dérivés (1) observés dans pop i au locus l
      ns[i, l]     = nb total de reads observés dans pop i au locus l
    """
    counts = np.array(
        [
            [reads_per_locus[loc][p][0] for loc in range(len(reads_per_locus))]
            for p in population_names
        ],
        dtype=float,
    )
    ns = np.array(
        [
            [reads_per_locus[loc][p][1] for loc in range(len(reads_per_locus))]
            for p in population_names
        ],
        dtype=float,
    )
    freq1 = counts / ns
    freq0 = 1.0 - freq1
    return counts, ns, freq0, freq1


def _forward_fill(
    values: np.ndarray, valid: np.ndarray, fill: float = 0.0
) -> np.ndarray:
    """Propagation 'forward-fill' vectorisée : aux positions où valid est False,
    propage la dernière valeur valide vue (ou `fill` si aucune encore).

    Reproduit le comportement de la variable C++ non ré-initialisée x_prev
    dans cal_snfstd et cal_snnei.  Implémentation : searchsorted sur les
    indices valides, O(n log n) tout-numpy, sans boucle Python.
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


def _halfsortedbypairs(v: list[int]) -> bool:
    """Filtre HALF de DIYABC (cal_snaml / cal_snf3r / cal_snf4r) :
    vrai si la permutation v est 'half-sorted by pairs'."""
    n = len(v)
    for i in range(n - 1, 0, -2):
        if v[i - 1] > v[i]:
            return False
        if i - 2 > 0 and v[i - 3] > v[i - 1]:
            return False
    return True


def _half_arrangements(n: int, r: int) -> list[list[int]]:
    """Arrangements HALF de r éléments parmi n -- ordre reproduit
    empiriquement depuis DIYABC (cal_snaml, cal_snf3r, cal_snf4r)."""
    result = []
    for combo in sorted(combinations(range(n), r), reverse=True):
        for perm in sorted(set(permutations(combo))):
            if _halfsortedbypairs(list(perm)):
                result.append(list(perm))
    return result


# ---------------------------------------------------------------------------
# ML1 : proportion de loci monomorphes par population (cal_snfl, npop=1)
# ---------------------------------------------------------------------------


def compute_ML1(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML1p_i : proportion de loci monomorphes dans la population i.
    Un locus est monomorphe si sum==0 (fixé dérivé) ou sum==n (fixé ancestral).
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    mono = (counts == 0) | (counts == ns)  # (npop, nloci) booléen
    return {
        f"ML1p_{i + 1}": float(mono[i].mean()) for i in range(len(population_names))
    }


# ---------------------------------------------------------------------------
# ML2 : proportion de loci fixés identiquement sur les paires (cal_snfl, npop=2)
# ---------------------------------------------------------------------------


def compute_ML2(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML2p_i.j : proportion de loci fixés au même allèle dans la paire (i, j).
    Référence : cal_snfl(npop=2) -- freq_a == freq_b ∈ {0, 1}.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    n = len(population_names)
    results = {}
    for i in range(n):
        for j in range(i + 1, n):
            both_zero = (counts[i] == 0) & (counts[j] == 0)
            both_one = (counts[i] == ns[i]) & (counts[j] == ns[j])
            results[f"ML2p_{i + 1}.{j + 1}"] = float((both_zero | both_one).mean())
    return results


# ---------------------------------------------------------------------------
# ML3 : proportion de loci fixés identiquement sur les triplets (cal_snfl, npop=3)
# ---------------------------------------------------------------------------


def compute_ML3(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """ML3p_i.j.k : même logique que ML2, sur les triplets de populations.
    Référence : cal_snfl(npop=3).
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    n = len(population_names)
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
# HW / HB : hétérozygotie intra- et inter-population (cal_snhw, cal_snhb)
# ---------------------------------------------------------------------------


def compute_HW_HB(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """HWm_i/HWv_i (intra-pop) et HBm_i.j/HBv_i.j (inter-pop).

    HW = 1 - q1,  HB = 1 - q2  (formules de sumstat.cpp vectorisées).
    q1[i,l] = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1))    [sans remise]
    q2[i,j,l] = (y1_i*y1_j + y2_i*y2_j) / (n_i*n_j)
    HWv et HBv utilisent ddof=1 (validé contre le C++).
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    y1, y2 = ns - counts, counts

    hw = 1.0 - (y1 * (y1 - 1) + y2 * (y2 - 1)) / (ns * (ns - 1))  # (npop, nloci)

    npop = len(population_names)
    results = {}
    for i in range(npop):
        results[f"HWm_{i + 1}"] = float(hw[i].mean())
        results[f"HWv_{i + 1}"] = float(hw[i].var(ddof=1))

    for i in range(npop):
        for j in range(i + 1, npop):
            hb = 1.0 - (y1[i] * y1[j] + y2[i] * y2[j]) / (ns[i] * ns[j])
            key = f"{i + 1}.{j + 1}"
            results[f"HBm_{key}"] = float(hb.mean())
            results[f"HBv_{key}"] = float(hb.var(ddof=1))

    return results


def compute_HW_HB_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    population_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    results = {}
    (
        counts,
        ns,
        _,
        _,
    ) = _mats or _prepare_matrices_poolseq(reads_per_locus, population_names)
    npop = len(population_names)
    ## Calcul de HWm et HWv pour chaque population
    for i in range(npop):
        r1 = counts[i]
        c1 = ns[i]
        r2 = c1 - r1
        s1 = r1 * (r1 - 1)
        s2 = r2 * (r2 - 1)
        np_i = pool_sizes[population_names[i]]
        q1 = ((np_i / (c1 * (c1 - 1))) * (s1 + s2) - 1) / (np_i - 1)
        hw = 1 - q1
        results[f"HWm_{i + 1}"] = float(hw.mean())
        results[f"HWv_{i + 1}"] = float(hw.var(ddof=1))

    # Calcul de HBm et HBv pour chaque paire de populations
    for i in range(npop):
        for j in range(i + 1, npop):
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
# FST1 : FST population-spécifique (cal_snfsti)
# ---------------------------------------------------------------------------


def compute_FST1(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST1m_i = 1 - HWm_i / HBmoy_global,  FST1v_i = HWv_i / HBmoy_global².

    HBmoy_global = moyenne de TOUS les HBm (toutes paires confondues) --
    confirmé dans cal_snfsti (sumstat.cpp), pas seulement les paires de pop_i.
    FST1v est une propagation d'erreur analytique, pas une variance empirique.
    """
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    y1, y2 = ns - counts, counts

    hw = 1.0 - (y1 * (y1 - 1) + y2 * (y2 - 1)) / (ns * (ns - 1))  # (npop, nloci)

    npop = len(population_names)
    all_hbm = []
    for i in range(npop):
        for j in range(i + 1, npop):
            hb = 1.0 - (y1[i] * y1[j] + y2[i] * y2[j]) / (ns[i] * ns[j])
            all_hbm.append(float(hb.mean()))

    hbmoy = float(np.mean(all_hbm)) if all_hbm else float("nan")

    results = {}
    for i in range(npop):
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
    population_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST1 (cal_snfsti n'a pas de branche
    type==3 : elle combine juste des HW/HB déjà calculés). Duplique donc
    ici la formule q1/q2 poolseq de compute_HW_HB_poolseq, exactement
    comme compute_FST1 duplique déjà la formule q1/q2 IndSeq plutôt que
    d'appeler compute_HW_HB -- même style que l'existant.
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(
        reads_per_locus, population_names
    )
    npop = len(population_names)

    hw_by_pop = []
    for i in range(npop):
        np_i = pool_sizes[population_names[i]]
        r1, c1 = counts[i], ns[i]
        r2 = c1 - r1
        s1, s2 = r1 * (r1 - 1), r2 * (r2 - 1)
        q1 = ((np_i / (c1 * (c1 - 1))) * (s1 + s2) - 1) / (np_i - 1)
        hw_by_pop.append(1.0 - q1)

    all_hbm = []
    for i in range(npop):
        for j in range(i + 1, npop):
            r11, c1 = counts[i], ns[i]
            r12 = c1 - r11
            r21, c2 = counts[j], ns[j]
            r22 = c2 - r21
            q2 = (r11 * r21 + r12 * r22) / (c1 * c2)
            hb = 1.0 - q2
            all_hbm.append(float(hb.mean()))

    hbmoy = float(np.mean(all_hbm)) if all_hbm else float("nan")

    results = {}
    for i in range(npop):
        hwm = float(hw_by_pop[i].mean())
        hwv = float(hw_by_pop[i].var(ddof=1))
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


def _fst_wc(loci, pops, _counts=None, _ns=None):
    """Weir & Cockerham vectorisé sur tous les loci.
    Retourne (FSTm, FSTv). Formule identique à cal_snfstd, toutes les
    opérations par-locus faites en numpy sur des vecteurs de longueur nloci.

    _counts, _ns : matrices (len(pops), nloci) pré-calculées -- passées
    comme slices de la matrice globale depuis compute_FST2/3/4 pour éviter
    de reconstruire les comptes locus par locus pour chaque sous-ensemble.
    """
    nloci = len(loci)
    npop = len(pops)

    if _counts is not None and _ns is not None:
        counts, ns = _counts, _ns
    else:
        counts = np.array([[sum(lg[p]) for lg in loci] for p in pops], dtype=float)
        ns = np.array([[len(lg[p]) for lg in loci] for p in pops], dtype=float)

    p1 = counts / ns
    p0 = 1.0 - p1

    S_1 = ns.sum(axis=0)
    S_2 = (ns**2).sum(axis=0)
    n_d = float(npop)

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

    sw2diff = nloci * (nloci - 1)
    mean = xs.mean()
    fstv = ((xs - mean) ** 2).sum() * nloci / sw2diff if sw2diff > 0 else 0.0

    return float(fstm), float(fstv)


def compute_FST2(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST2m_i.j / FST2v_i.j : Weir & Cockerham par paire."""
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    npop = len(population_names)
    results = {}
    for i in range(npop):
        for j in range(i + 1, npop):
            key = f"{i + 1}.{j + 1}"
            m, v = _fst_wc(
                genotypes_per_locus,
                [population_names[i], population_names[j]],
                _counts=counts[[i, j]],
                _ns=ns[[i, j]],
            )
            results[f"FST2m_{key}"] = m
            results[f"FST2v_{key}"] = v
    return results


def _fst_wc_poolseq(pops, pool_sizes, _counts, _ns):
    """Variante PoolSeq de _fst_wc (cal_snfstd, branche grouplist[gr].type==3).

    Calcul en DEUX passes (contrairement à l'IndSeq) : pi1/pi2 (moyennes
    pondérées par la profondeur de lecture) doivent être connues avant de
    calculer SSP. C_1/C_1_star mélangent la profondeur de lecture (`n`,
    variable par locus) et la VRAIE taille du pool (`c`, constante par
    population) -- c'est ce mélange qui constitue la correction propre à
    PoolSeq (le terme de variance intra-pool supplémentaire).

    L'agrégation finale (ratio de sommes num/den + variance via
    _forward_fill) est IDENTIQUE à _fst_wc -- confirmé par l'exploration
    C++, même code d'agrégation pour les deux types de population.

    _counts, _ns : matrices (len(pops), nloci) pré-calculées, slices de
    la matrice globale (nreads1, nreads_total) -- même contrat que
    _fst_wc.
    """
    x1, n = _counts, _ns  # (npop, nloci) : reads allèle1, profondeur de lecture
    x2 = n - x1
    nloci = n.shape[1]
    c = np.array([pool_sizes[p] for p in pops], dtype=float).reshape(-1, 1)

    # --- Passe 1 ---
    term = n / c + (c - 1) / c  # (npop, nloci)
    C_1 = term.sum(axis=0)  # (nloci,)
    C_1_star = (n * term).sum(axis=0)  # (nloci,)
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

    sw2diff = nloci * (nloci - 1)
    mean = xs.mean()
    fstv = ((xs - mean) ** 2).sum() * nloci / sw2diff if sw2diff > 0 else 0.0

    return float(fstm), float(fstv)


def compute_FST2_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    population_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST2 : FST2m_i.j / FST2v_i.j par paire."""
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(
        reads_per_locus, population_names
    )
    npop = len(population_names)
    results = {}
    for i in range(npop):
        for j in range(i + 1, npop):
            key = f"{i + 1}.{j + 1}"
            pops = [population_names[i], population_names[j]]
            m, v = _fst_wc_poolseq(
                pops, pool_sizes, _counts=counts[[i, j]], _ns=ns[[i, j]]
            )
            results[f"FST2m_{key}"] = m
            results[f"FST2v_{key}"] = v
    return results


def compute_FST3_FST4_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    population_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_FST3_FST4_FSTG : FST3/FST4 sur
    triplets/quadruplets (COMB)."""
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(
        reads_per_locus, population_names
    )
    npop = len(population_names)
    results = {}

    for combo in combinations(range(npop), 3):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        pops = [population_names[i] for i in idx]
        m, v = _fst_wc_poolseq(pops, pool_sizes, _counts=counts[idx], _ns=ns[idx])
        results[f"FST3m_{key}"] = m
        results[f"FST3v_{key}"] = v

    for combo in combinations(range(npop), 4):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        pops = [population_names[i] for i in idx]
        m, v = _fst_wc_poolseq(pops, pool_sizes, _counts=counts[idx], _ns=ns[idx])
        results[f"FST4m_{key}"] = m
        results[f"FST4v_{key}"] = v

    return results


def compute_FST3_FST4_FSTG(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """FST3/FST4 : Weir & Cockerham sur triplets et quadruplets (COMB)."""
    counts, ns, _, _ = _mats or _prepare_matrices(genotypes_per_locus, population_names)
    npop = len(population_names)
    results = {}

    for combo in combinations(range(npop), 3):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        m, v = _fst_wc(
            genotypes_per_locus,
            [population_names[i] for i in idx],
            _counts=counts[idx],
            _ns=ns[idx],
        )
        results[f"FST3m_{key}"] = m
        results[f"FST3v_{key}"] = v

    for combo in combinations(range(npop), 4):
        idx = list(combo)
        key = ".".join(str(i + 1) for i in idx)
        m, v = _fst_wc(
            genotypes_per_locus,
            [population_names[i] for i in idx],
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
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire, vectorisée.
    NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)
    x_prev persiste si denom==0 (comportement C++ non réinitialisé) --
    reproduit via _forward_fill.
    """
    counts, ns, freq0, freq1 = _mats or _prepare_matrices(
        genotypes_per_locus, population_names
    )
    nloci = len(genotypes_per_locus)
    f, g = freq0, freq1
    norm = np.sqrt(f * f + g * g)  # (npop, nloci)

    results = {}
    npop = len(population_names)
    for i in range(npop):
        for j in range(i + 1, npop):
            denom = norm[i] * norm[j]
            valid = denom > 0
            nei = np.where(
                valid,
                1.0 - (f[i] * f[j] + g[i] * g[j]) / np.where(valid, denom, 1.0),
                0.0,
            )
            xs = _forward_fill(nei, valid, fill=0.0)

            key = f"{i + 1}.{j + 1}"
            sw2diff = nloci * (nloci - 1)
            mean = xs.mean()
            results[f"NEIm_{key}"] = float(mean)
            results[f"NEIv_{key}"] = (
                float(((xs - mean) ** 2).sum() * nloci / sw2diff)
                if sw2diff > 0
                else 0.0
            )

    return results


# ---------------------------------------------------------------------------
# AML : admixture maximum likelihood sur triplets HALF (cal_snaml)
# ---------------------------------------------------------------------------


def compute_AML(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """AMLm / AMLv : coefficient d'admixture ML sur triplets HALF.
    aml = (f3-f2)/(f1-f2) clampé à [0,1].
    Les loci non informatifs (f1==f2, w=0) sont exclus de la moyenne --
    équivalent au Welford pondéré avec w ∈ {0,1}, ce qui réduit à
    mean/var(ddof=1) sur le sous-ensemble informatif.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(
        genotypes_per_locus, population_names
    )
    npop = len(population_names)
    results = {}

    for t in _half_arrangements(npop, 3):
        h, p1, p2 = t[0], t[1], t[2]
        key = f"{h + 1}.{p1 + 1}.{p2 + 1}"

        f1 = freq0[p1]  # freq allèle 0 dans parent 1  (nloci,)
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
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """F3m/F3v sur triplets HALF,
    F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)   [pop0=hybride, pop1/2=parents]
    Tous les loci ont w=1 → mean/var(ddof=1) directement.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(
        genotypes_per_locus, population_names
    )
    npop = len(population_names)
    results = {}

    # --- F3 ---
    for t in _half_arrangements(npop, 3):
        i0, i1, i2 = t[0], t[1], t[2]
        key = f"{i0 + 1}.{i1 + 1}.{i2 + 1}"

        np_ = ns[i0]  # nb lignées dans l'hybride  (nloci,)
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
    population_names: list[str],
    pool_sizes: dict[str, int],
    _mats=None,
) -> dict[str, float]:
    """Variante PoolSeq de compute_F3 (cal_snf3r, branche grouplist[gr].type==3).

    alpha = ((np*a1p*(a1p-1))/(c1p*(c1p-1)) - a1p/c1p) / (np-1)  [pop0=hybride]
    F3 = alpha + betaBC - betaAB - betaAC, avec beta_XY = (aXp*aYp)/(cXp*cYp)

    np = taille du pool (VRAIE, pas la profondeur de lecture) de la
    population hybride -- vient de pool_sizes, pas de `ns`/`_mats`
    (contrairement à ns, qui est la profondeur de lecture, variable par
    locus). Agrégation identique à compute_F3 (mean/var(ddof=1) simples
    sur les loci, pas de ratio de sommes).
    """
    counts, ns, _, _ = _mats or _prepare_matrices_poolseq(
        reads_per_locus, population_names
    )
    npop = len(population_names)
    results = {}

    for t in _half_arrangements(npop, 3):
        i0, i1, i2 = t[0], t[1], t[2]
        key = f"{i0 + 1}.{i1 + 1}.{i2 + 1}"

        np_i0 = pool_sizes[population_names[i0]]

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
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """
    F4m/F4v sur quadruplets HALF.
    F4 = (a-b)*(c-d)
    Tous les loci ont w=1 → mean/var(ddof=1) directement.
    """
    counts, ns, freq0, _ = _mats or _prepare_matrices(
        genotypes_per_locus, population_names
    )
    npop = len(population_names)
    results = {}

    # --- F4 ---
    for t in _half_arrangements(npop, 4):
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


def _genotype_matrix_by_population(
    tree_sequence: tskit.TreeSequence,
) -> dict[str, np.ndarray]:
    """Crée une sous-matrice de génotypes par population.

    Retourne {nom_pop: matrice (n_sites, n_samples_pop)} -- convention
    native de tskit (genotype_matrix() est déjà (sites, samples)), pas
    de transposition. genotype_matrix() n'est appelé qu'UNE FOIS pour
    toute la TreeSequence, puis tranché par population via fancy
    indexing (pas de reconstruction par sample).
    """
    genotype_matrix = tree_sequence.genotype_matrix()
    layout = compute_population_layout(tree_sequence)
    return {pop_name: genotype_matrix[:, sample_ids] for pop_name, sample_ids in layout}


# ---------------------------------------------------------------------------
# NSS : nombre de sites ségrégeants par population
# ---------------------------------------------------------------------------


def _segregating_sites_mask(matrix: np.ndarray) -> np.ndarray:
    """Masque booléen (longueur n_sites) : True si le site n'est pas
    identique chez tous les échantillons de la matrice (n_sites,
    n_samples). Factorisé hors de _count_segregating_sites pour être
    réutilisé par PSS (_private_segregating_sites_per_locus), qui a
    besoin du masque par site, pas seulement du compte agrégé."""
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    return np.any(matrix != matrix[:, [0]], axis=1)


def _count_segregating_sites(matrix: np.ndarray) -> int:
    """Compte le nombre de sites polymorphes dans une matrice de génotypes
    (n_sites, n_samples)."""
    return int(np.sum(_segregating_sites_mask(matrix)))


def mean_segregating_sites_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule NSS_i (cal_nss1p) : pour chaque population, la moyenne du
    nombre de sites ségrégeants sur tous les loci du groupe passé en
    argument (un groupe = les TreeSequences des loci séquence d'un même
    `group Gx` du header, ex. les 5 loci <A> de G2).

    `population_names` fixe explicitement les clés du dict retourné (comme
    compute_ML1/compute_HW_HB) -- chaque population attendue a toujours une
    valeur (0.0 par défaut, comme le `res = 0.0` du C++), même si
    `tree_sequences` est vide, plutôt que d'être silencieusement absente du
    résultat.

    Suppose que toutes les populations de `population_names` sont présentes
    sur tous les loci du groupe (divise par `len(tree_sequences)`, pas par
    un décompte par population comme le `nl` du C++ -- lève un KeyError si
    ce n'est pas le cas plutôt que d'exclure silencieusement ce locus,
    contrairement au C++) -- vérifié vrai sur toy_example2_ms_dna, pas
    garanti en général.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    num_loci = len(tree_sequences)
    mean_segregating_sites = {pop_name: 0.0 for pop_name in population_names}
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            matrix = genotype_matrices[pop_name]
            mean_segregating_sites[pop_name] += _count_segregating_sites(matrix)

    if num_loci > 0:
        for pop_name in population_names:
            mean_segregating_sites[pop_name] /= num_loci

    return mean_segregating_sites


# ---------------------------------------------------------------------------
# NDH : nombre d'haplotypes distincts par population
# ---------------------------------------------------------------------------


def _count_distinct_haplotypes(matrix: np.ndarray) -> int:
    """Compte le nombre d'haplotypes distincts dans une matrice de génotypes
    (n_sites, n_samples)."""
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    # np.unique(axis=1) déduplique les colonnes (les haplotypes) --
    # le résultat a la forme (n_sites, n_haplotypes_distincts).
    distinct_haplotypes = np.unique(matrix, axis=1)
    return distinct_haplotypes.shape[1]


def mean_distinct_haplotypes_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule le nombre moyen d'haplotypes distincts par population sur un groupe de loci.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    num_loci = len(tree_sequences)
    mean_distinct_haplotypes = {pop_name: 0.0 for pop_name in population_names}
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            matrix = genotype_matrices[pop_name]
            mean_distinct_haplotypes[pop_name] += _count_distinct_haplotypes(matrix)

    if num_loci > 0:
        for pop_name in population_names:
            mean_distinct_haplotypes[pop_name] /= num_loci

    return mean_distinct_haplotypes


# ---------------------------------------------------------------------------
# MDP/VDP : moyenne et variance de différences de paires (pairwise differences) par population
# ---------------------------------------------------------------------------


def _pairwise_hamming_distances(matrix: np.ndarray) -> np.ndarray:
    """Calcule les distances de Hamming par paire d'échantillons, pour une
    matrice de génotypes (n_sites, n_samples). Retourne un vecteur 1D de
    longueur C(n_samples, 2) -- une valeur par paire (i, j) avec i < j,
    pas la matrice carrée (n_samples, n_samples) complète."""

    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")

    matrix_hammming = (matrix[:, :, None] != matrix[:, None, :]).sum(axis=0)
    triangle_sup = np.triu_indices(matrix_hammming.shape[0], k=1)
    return matrix_hammming[triangle_sup]


def mean_pairwise_differences_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule MPD_i (cal_mpd1p) : pour chaque population, la moyenne du
    nombre de différences par paire (distance de Hamming) sur tous les
    loci du groupe passé en argument (un groupe = les TreeSequences des
    loci séquence d'un même `group Gx` du header).

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Contrairement à mean_segregating_sites_per_group/mean_distinct_
    haplotypes_per_group (qui divisent par `len(tree_sequences)`), le
    dénominateur ici est un compteur PAR POPULATION (comme le `nl` de
    cal_mpd1p) : un locus où une population a moins de 2 échantillons
    donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
    vide), auquel cas ce locus est exclu du calcul pour cette population
    -- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
    laisser un `nan` (moyenne d'un vecteur vide) empoisonner le résultat.

    Suppose quand même que toutes les populations de `population_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `mean_*_per_group` de ce module.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    mean_pairwise_differences = {pop_name: 0.0 for pop_name in population_names}
    valid_loci_count = {pop_name: 0 for pop_name in population_names}
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            matrix = genotype_matrices[pop_name]
            pairwise_distances = _pairwise_hamming_distances(matrix)
            if len(pairwise_distances) > 0:
                mean_pairwise_differences[pop_name] += pairwise_distances.mean()
                valid_loci_count[pop_name] += 1

    for pop_name in population_names:
        if valid_loci_count[pop_name] > 0:
            mean_pairwise_differences[pop_name] /= valid_loci_count[pop_name]

    return mean_pairwise_differences


def variance_pairwise_differences_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule VPD_i (cal_vpd1p) : pour chaque population, la variance du
    nombre de différences par paire (distance de Hamming) sur tous les
    loci du groupe passé en argument (un groupe = les TreeSequences des
    loci séquence d'un même `group Gx` du header).

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Contrairement à mean_segregating_sites_per_group/mean_distinct_
    haplotypes_per_group (qui divisent par `len(tree_sequences)`), le
    dénominateur ici est un compteur PAR POPULATION (comme le `nl` de
    cal_vpd1p) : un locus où une population a moins de 2 échantillons
    donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
    vide), auquel cas ce locus est exclu du calcul pour cette population
    -- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
    laisser un `nan` (variance d'un vecteur vide) empoisonner le résultat.

    Suppose quand même que toutes les populations de `population_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `mean_*_per_group` de ce module.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_variance}
    """
    variance_pairwise_differences = {pop_name: 0.0 for pop_name in population_names}
    valid_loci_count = {pop_name: 0 for pop_name in population_names}
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            matrix = genotype_matrices[pop_name]
            pairwise_distances = _pairwise_hamming_distances(matrix)
            if len(pairwise_distances) > 1:
                variance_pairwise_differences[pop_name] += pairwise_distances.var(
                    ddof=1
                )
                valid_loci_count[pop_name] += 1

    for pop_name in population_names:
        if valid_loci_count[pop_name] > 0:
            variance_pairwise_differences[pop_name] /= valid_loci_count[pop_name]

    return variance_pairwise_differences


# ---------------------------------------------------------------------------
# DTA : distance de Tajima par population
# ---------------------------------------------------------------------------


def _tajima_constants(n_samples: int) -> tuple[float, float, float]:
    """Calcule les constantes a1, e1, e2 pour le D de Tajima (cal_dta1pl,
    lignes 1566-1575) à partir du nombre d'échantillons (n_samples) --
    ne dépend que de n_samples, jamais des données elles-mêmes.

    Args:
        n_samples (int): Nombre d'échantillons (>= 2).

    Returns:
        tuple: (a1, e1, e2) -- a1 est aussi réutilisé directement dans la
        formule finale du D (S / a1), e1/e2 sont les coefficients de la
        variance sous neutralité. b1, b2, c1, c2, a2 sont des étapes
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
    """D de Tajima (cal_dta1pl) pour UNE population sur UN locus, à
    partir de sa matrice de génotypes (n_sites, n_samples).

    D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))  -- pi = MPD (moyenne des
    différences par paire), S = NSS (nombre de sites ségrégeants),
    a1/e1/e2 = _tajima_constants(n_samples).

    Retourne `None` si `n_samples < 2` (pas assez d'échantillons pour
    calculer ne serait-ce qu'une paire -- ce locus doit être EXCLU de la
    moyenne du groupe, `OKK = false` côté C++). Retourne `0.0` si
    `n_samples >= 2` mais qu'aucun site n'est ségrégeant (S == 0, le
    dénominateur de la formule est nul) -- ce locus reste INCLUS dans la
    moyenne du groupe avec une valeur de 0.0, contrairement au cas
    précédent : le C++ ne repasse jamais `OKK` à false dans ce cas
    (lignes 1579/1594-1598) -- distinction délibérée, pas une
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


def mean_tajima_d_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule DTA_i (cal_dta1p) : pour chaque population, la moyenne du
    D de Tajima (_tajima_d_per_locus) sur tous les loci du groupe passé
    en argument (un groupe = les TreeSequences des loci séquence d'un
    même `group Gx` du header).

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut,
    comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

    Comme mean_pairwise_differences_per_group/variance_pairwise_
    differences_per_group, le dénominateur est un compteur PAR
    POPULATION (le `nl` de cal_dta1p) : un locus où `_tajima_d_per_locus`
    retourne `None` (moins de 2 échantillons) est exclu -- ni ajouté à la
    somme, ni compté. Un locus où `_tajima_d_per_locus` retourne `0.0`
    (0 site ségrégeant) reste, lui, INCLUS dans le compte (voir
    _tajima_d_per_locus pour la distinction).

    Suppose quand même que toutes les populations de `population_names`
    sont présentes (au moins 1 échantillon) sur tous les loci du groupe
    -- lève un KeyError si ce n'est pas le cas, comme les autres
    fonctions `mean_*_per_group` de ce module.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    mean_tajima_d = {pop_name: 0.0 for pop_name in population_names}
    valid_loci_count = {pop_name: 0 for pop_name in population_names}
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            matrix = genotype_matrices[pop_name]
            tajima_d = _tajima_d_per_locus(matrix)
            if tajima_d is not None:
                mean_tajima_d[pop_name] += tajima_d
                valid_loci_count[pop_name] += 1

    for pop_name in population_names:
        if valid_loci_count[pop_name] > 0:
            mean_tajima_d[pop_name] /= valid_loci_count[pop_name]

    return mean_tajima_d


# ---------------------------------------------------------------------------
# PSS : sites ségrégeants "privés" par population (cal_pss1p)
# ---------------------------------------------------------------------------


def _private_segregating_sites_per_locus(
    genotype_matrices: dict[str, np.ndarray], target_pop: str
) -> int:
    """Nombre de sites ségrégeants "privés" de `target_pop` sur UN locus
    (cal_pss1p) : sites ségrégeants dans `target_pop` mais NULLE PART
    ailleurs, parmi TOUTES les populations de `genotype_matrices` (pas
    seulement celles d'un même groupe -- le C++ compare à `this->nsample`,
    le nombre total de populations du dataset).

    Toutes les matrices de `genotype_matrices` viennent du même
    `genotype_matrix()` (juste tranchées par colonnes, voir
    _genotype_matrix_by_population) -- la ligne `i` désigne donc le MÊME
    site physique pour toutes les populations. Contrairement au C++, qui
    compare des listes d'indices de sites variables de longueurs
    différentes par une recherche imbriquée (`ssa[sample][j] ==
    ssa[sa][k]`), on peut donc comparer les masques booléens position par
    position directement -- pas de recherche d'égalité nécessaire.
    """
    target_mask = _segregating_sites_mask(genotype_matrices[target_pop])
    other_masks = [
        _segregating_sites_mask(matrix)
        for pop_name, matrix in genotype_matrices.items()
        if pop_name != target_pop
    ]
    segregating_elsewhere = (
        np.logical_or.reduce(other_masks) if other_masks else np.zeros_like(target_mask)
    )
    return int(np.sum(target_mask & ~segregating_elsewhere))


def mean_private_segregating_sites_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule PSS_i (cal_pss1p) : pour chaque population, la moyenne du
    nombre de sites ségrégeants privés (_private_segregating_sites_per_
    locus) sur tous les loci du groupe passé en argument.

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. `population_names` DOIT couvrir
    TOUTES les populations du dataset (pas seulement celles d'un groupe),
    puisque `_private_segregating_sites_per_locus` compare `target_pop` à
    toutes les autres populations présentes dans `genotype_matrices`.

    Contrairement à mean_segregating_sites_per_group et aux autres
    fonctions `mean_*_per_group` de ce module, le dénominateur ici est
    `len(tree_sequences)` SANS AUCUNE exclusion (`nl` s'incrémente sans
    condition dans cal_pss1p, ligne 1624 -- pas de garde-fou du tout,
    même pas le `samplesize > 0` de NSS/NHA).

    Suppose que toutes les populations de `population_names` sont
    présentes sur tous les loci du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues (toutes celles du dataset).

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    mean_pss = {pop_name: 0.0 for pop_name in population_names}
    num_loci = len(tree_sequences)
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            mean_pss[pop_name] += _private_segregating_sites_per_locus(
                genotype_matrices, pop_name
            )

    if num_loci > 0:
        for pop_name in population_names:
            mean_pss[pop_name] /= num_loci

    return mean_pss


# ---------------------------------------------------------------------------
# MNS/VNS : moyenne/variance du compte de l'allèle minoritaire (afs, cal_mns1p/cal_vns1p)
# ---------------------------------------------------------------------------


def _minor_allele_counts_at_segregating_sites(matrix: np.ndarray) -> np.ndarray:
    """Pour chaque site ségrégeant d'une matrice de génotypes (n_sites,
    n_samples), calcule le compte du/des allèle(s) le(s) moins
    fréquent(s) parmi les échantillons (afs, `nf[jj]` après tri croissant
    et saut des zéros -- équivalent à `min(comptes des valeurs
    effectivement présentes)`, généralisation multi-allélique du "minor
    allele count" puisqu'un site ADN peut avoir jusqu'à 4 états A/C/G/T,
    pas seulement 2 comme un SNP).

    Retourne un vecteur 1D de longueur = nombre de sites SÉGRÉGEANTS
    (pas n_sites) -- un site fixé (une seule valeur parmi les
    échantillons) n'est pas inclus, comme `afs` qui ne pousse rien dans
    `t_afs` quand `jj >= 3` (une seule base présente).
    """
    if matrix.shape[1] == 0:
        raise ValueError("La matrice de génotypes est vide.")
    minor_counts = []
    for site in matrix:
        _, counts = np.unique(site, return_counts=True)
        if len(counts) > 1:
            minor_counts.append(counts.min())
    return np.array(minor_counts, dtype=float)


def mean_minor_allele_count_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule MNS_i (cal_mns1p) : pour chaque population, la moyenne,
    sur les loci du groupe, de la moyenne (par locus) des comptes
    d'allèle minoritaire aux sites ségrégeants
    (_minor_allele_counts_at_segregating_sites).

    Un locus sans site ségrégeant contribue 0.0 (la boucle C++ sur
    `t_afs` ne s'exécute simplement pas -- même effet qu'un vecteur
    vide ici). Comme PSS (et contrairement à MPD/VPD/DTA), AUCUNE
    exclusion de locus : `nl` = `len(tree_sequences)` sans condition
    (cal_mns1p, `nl++` inconditionnel, ligne 1705) -- pas besoin de
    compteur par population.

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    populations de `population_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    mean_mns = {pop_name: 0.0 for pop_name in population_names}
    num_loci = len(tree_sequences)
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            minor_counts = _minor_allele_counts_at_segregating_sites(
                genotype_matrices[pop_name]
            )
            if len(minor_counts) > 0:
                mean_mns[pop_name] += minor_counts.mean()

    if num_loci > 0:
        for pop_name in population_names:
            mean_mns[pop_name] /= num_loci

    return mean_mns


def variance_minor_allele_count_per_group(
    tree_sequences: list[tskit.TreeSequence],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule VNS_i (cal_vns1p) : pour chaque population, la moyenne,
    sur les loci du groupe, de la variance (par locus) des comptes
    d'allèle minoritaire aux sites ségrégeants.

    ATTENTION : variance BIAISÉE (ddof=0, division par n, PAS n-1) --
    contrairement à variance_pairwise_differences_per_group (VPD) qui
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

    `population_names` fixe explicitement les clés du dict retourné --
    chaque population attendue a toujours une valeur (0.0 par défaut),
    même si `tree_sequences` est vide. Suppose que toutes les
    populations de `population_names` sont présentes sur tous les loci
    du groupe -- lève un KeyError sinon.

    Args:
        tree_sequences (list[tskit.TreeSequence]): Liste des TreeSequences à analyser.
        population_names (list[str]): Populations attendues.

    Returns:
        dict: {pop_name: valeur_moyenne}
    """
    variance_vns = {pop_name: 0.0 for pop_name in population_names}
    num_loci = len(tree_sequences)
    for ts in tree_sequences:
        genotype_matrices = _genotype_matrix_by_population(ts)
        for pop_name in population_names:
            minor_counts = _minor_allele_counts_at_segregating_sites(
                genotype_matrices[pop_name]
            )
            if len(minor_counts) > 1:
                variance_vns[pop_name] += minor_counts.var()

    if num_loci > 0:
        for pop_name in population_names:
            variance_vns[pop_name] /= num_loci

    return variance_vns


# ---------------------------------------------------------------------------
# Point d'entrée principal
# ---------------------------------------------------------------------------


def compute_all_statistics(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule les 130 statistiques résumées SNP et retourne un dict
    {nom_stat: valeur} -- même format que parse_statobs().

    Les matrices (npop × nloci) de comptes et fréquences sont construites
    une seule fois (_prepare_matrices) et transmises à toutes les familles
    de statistiques via _mats.
    """
    mats = _prepare_matrices(genotypes_per_locus, population_names)
    results = {}
    results.update(compute_ML1(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_ML2(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_ML3(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_HW_HB(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_FST1(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_FST2(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_NEI(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_AML(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_F3(genotypes_per_locus, population_names, _mats=mats))
    results.update(compute_F4(genotypes_per_locus, population_names, _mats=mats))
    results.update(
        compute_FST3_FST4_FSTG(genotypes_per_locus, population_names, _mats=mats)
    )
    return results


def compute_all_statistics_poolseq(
    reads_per_locus: list[dict[str, tuple[int, int]]],
    population_names: list[str],
    pool_sizes: dict[str, int],
) -> dict[str, float]:
    """Calcule les statistiques résumées SNP pour POOLSEQ et retourne un dict
    {nom_stat: valeur} -- même format que parse_statobs().

    Les matrices (npop × nloci) de comptes et tailles d'échantillon sont construites
    une seule fois (_prepare_matrices_poolseq) et transmises à toutes les familles
    de statistiques via _mats.
    """
    mats = _prepare_matrices_poolseq(reads_per_locus, population_names)

    results = {}
    results.update(compute_ML1(None, population_names, _mats=mats))
    results.update(compute_ML2(None, population_names, _mats=mats))
    results.update(compute_ML3(None, population_names, _mats=mats))
    results.update(
        compute_HW_HB_poolseq(
            reads_per_locus, population_names, pool_sizes=pool_sizes, _mats=mats
        )
    )
    results.update(compute_NEI(reads_per_locus, population_names, _mats=mats))
    results.update(compute_F4(reads_per_locus, population_names, _mats=mats))
    results.update(
        compute_F3_poolseq(reads_per_locus, population_names, pool_sizes, _mats=mats)
    )
    results.update(
        compute_FST2_poolseq(reads_per_locus, population_names, pool_sizes, _mats=mats)
    )
    results.update(
        compute_FST3_FST4_poolseq(
            reads_per_locus, population_names, pool_sizes, _mats=mats
        )
    )
    results.update(
        compute_FST1_poolseq(reads_per_locus, population_names, pool_sizes, _mats=mats)
    )
    results.update(compute_AML(reads_per_locus, population_names, _mats=mats))
    return results
