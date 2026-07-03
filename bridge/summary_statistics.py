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
    Un locus est monomorphe si sum==0 (fixé ancestral) ou sum==n (fixé dérivé).
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


def compute_F3_F4(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
    _mats=None,
) -> dict[str, float]:
    """F3m/F3v sur triplets HALF, F4m/F4v sur quadruplets HALF.

    F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)   [pop0=hybride, pop1/2=parents]
    F4 = (a-b)*(c-d)

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
    results.update(compute_F3_F4(genotypes_per_locus, population_names, _mats=mats))
    results.update(
        compute_FST3_FST4_FSTG(genotypes_per_locus, population_names, _mats=mats)
    )
    return results
