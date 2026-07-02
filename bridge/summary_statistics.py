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

from typing import Iterator
import numpy as np


# ---------------------------------------------------------------------------
# Utilitaires internes
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
    y2 = sum(haploid_genotypes)  # count of allele 1
    y1 = n - y2  # count of allele 0
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
    y12 = sum(haploid_genotypes_a)  # allele 1 in pop a
    y11 = n1 - y12  # allele 0 in pop a
    y22 = sum(haploid_genotypes_b)  # allele 1 in pop b
    y21 = n2 - y22  # allele 0 in pop b
    return (y11 * y21 + y12 * y22) / (n1 * n2)


# ---------------------------------------------------------------------------
# ML1 : proportion de loci monomorphes, par population
# Référence : sumstat.cpp::cal_snfl (npop=1)
# ---------------------------------------------------------------------------


def compute_ML1(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule ML1p_i pour chaque population i : proportion de loci
    monomorphes dans cette population.

    Un locus est monomorphe si freq(allele 0) == 0 OU == 1, c'est-à-dire
    si sum(genos) == 0 (que des ancestraux) ou sum(genos) == n (que des
    dérivés) -- traduit de la condition "freq[pop][0] == 0.0 or == 1.0"
    de cal_snfl dans sumstat.cpp.
    """
    monomorphic_counts = {pop: 0 for pop in population_names}
    total_loci = len(genotypes_per_locus)

    for locus_genotypes in genotypes_per_locus:
        for pop in population_names:
            genos = locus_genotypes[pop]
            n = len(genos)
            s = sum(genos)
            if s == 0 or s == n:
                monomorphic_counts[pop] += 1

    return {
        f"ML1p_{i + 1}": monomorphic_counts[pop] / total_loci
        for i, pop in enumerate(population_names)
    }


def compute_ML2(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """ML2p_i.j : proportion de loci monomorphes ET identiques
    (même allèle fixé) dans la paire (pop_i, pop_j).
    Référence : cal_snfl(npop=2) -- un locus est fixé sur la paire si
    freq0_pop_i == freq0_pop_j ET vaut 0.0 ou 1.0."""
    results = {}
    total = len(genotypes_per_locus)
    for i, pa in enumerate(population_names):
        for j, pb in enumerate(population_names):
            if j <= i:
                continue
            count = 0
            for lg in genotypes_per_locus:
                na, nb = len(lg[pa]), len(lg[pb])
                sa, sb = sum(lg[pa]), sum(lg[pb])
                freq_a = sa / na if na else float("nan")
                freq_b = sb / nb if nb else float("nan")
                fixed_a = freq_a == 0.0 or freq_a == 1.0
                if fixed_a and freq_a == freq_b:
                    count += 1
            results[f"ML2p_{i + 1}.{j + 1}"] = count / total
    return results


def compute_ML3(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """ML3p_i.j.k : même logique que ML2, sur les triplets de populations.
    Référence : cal_snfl(npop=3)."""
    results = {}
    total = len(genotypes_per_locus)
    n = len(population_names)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                pa, pb, pc = (
                    population_names[i],
                    population_names[j],
                    population_names[k],
                )
                count = 0
                for lg in genotypes_per_locus:
                    na, nb, nc = len(lg[pa]), len(lg[pb]), len(lg[pc])
                    freq_a = sum(lg[pa]) / na if na else float("nan")
                    freq_b = sum(lg[pb]) / nb if nb else float("nan")
                    freq_c = sum(lg[pc]) / nc if nc else float("nan")
                    if (freq_a == 0.0 or freq_a == 1.0) and freq_a == freq_b == freq_c:
                        count += 1
                results[f"ML3p_{i + 1}.{j + 1}.{k + 1}"] = count / total
    return results


# ---------------------------------------------------------------------------
# HW : hétérozygotie intra-population (within)
# HB : hétérozygotie inter-population (between, par paire)
# Référence : sumstat.cpp::cal_snhw, cal_snhb
# ---------------------------------------------------------------------------


def compute_HW_HB(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule HWm_i (moyenne) et HWv_i (variance) pour chaque population,
    et HBm_i.j (moyenne) et HBv_i.j (variance) pour chaque paire.

    HW = 1 - q1 (hétérozygotie intra-pop)
    HB = 1 - q2 (hétérozygotie inter-pop)

    Référence exacte : cal_snhw et cal_snhb dans sumstat.cpp.
    """
    npop = len(population_names)
    hw_per_locus = {pop: [] for pop in population_names}
    hb_per_locus = {
        (pa, pb): []
        for i, pa in enumerate(population_names)
        for pb in population_names[i + 1 :]
    }

    for locus_genotypes in genotypes_per_locus:
        for pop in population_names:
            hw_per_locus[pop].append(1 - _q1(locus_genotypes[pop]))
        for i, pa in enumerate(population_names):
            for pb in population_names[i + 1 :]:
                hb_per_locus[(pa, pb)].append(
                    1 - _q2(locus_genotypes[pa], locus_genotypes[pb])
                )

    results = {}
    for i, pop in enumerate(population_names):
        vals = hw_per_locus[pop]
        results[f"HWm_{i + 1}"] = float(np.mean(vals))
        results[f"HWv_{i + 1}"] = float(np.var(vals, ddof=1))

    for i, pa in enumerate(population_names):
        for j, pb in enumerate(population_names):
            if j <= i:
                continue
            vals = hb_per_locus[(pa, pb)]
            key = f"{i + 1}.{j + 1}"
            results[f"HBm_{key}"] = float(np.mean(vals))
            results[f"HBv_{key}"] = float(np.var(vals, ddof=1))

    return results


# ---------------------------------------------------------------------------
# FST1 : FST population-spécifique (Weir & Goudet 2017 / Hivert et al. 2018)
# Référence : sumstat.cpp::cal_snfsti
# FST1 = 1 - HW / HBmoy, où HBmoy est la moyenne des HB impliquant cette pop
# ---------------------------------------------------------------------------


def compute_FST1(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule FST1m_i et FST1v_i pour chaque population i.

    Formule exacte lue dans cal_snfsti (sumstat.cpp) :
        FST1m_i = 1 - HWm_i / HBmoy_global
        FST1v_i = HWv_i / (HBmoy_global²)

    où HBmoy_global = moyenne de TOUS les HBm (toutes les paires de
    populations, pas seulement celles impliquant pop_i) -- confirmé par
    le code C++ qui additionne tous les HBm sans filtre de population
    (boucle sur grouplist[gr].sumstat, condition curstat == Hbstat).

    FST1v n'est PAS une variance empirique locus par locus -- c'est une
    formule analytique de propagation d'erreur : HWv_i / HBmoy².
    Découverte en relisant précisément cal_snfsti :
        stsnp.mx2 = Hwv / (Hbmoy * Hbmoy)
    """
    npop = len(population_names)
    pairs = [
        (pa, pb)
        for i, pa in enumerate(population_names)
        for pb in population_names[i + 1 :]
    ]

    hw_per_locus = {pop: [] for pop in population_names}
    hb_per_locus = {pair: [] for pair in pairs}

    for locus_genotypes in genotypes_per_locus:
        for pop in population_names:
            hw_per_locus[pop].append(1 - _q1(locus_genotypes[pop]))
        for pa, pb in pairs:
            hb_per_locus[(pa, pb)].append(
                1 - _q2(locus_genotypes[pa], locus_genotypes[pb])
            )

    # HBmoy_global = moyenne de TOUS les HBm (toutes paires)
    all_hbm = [float(np.mean(hb_per_locus[pair])) for pair in pairs]
    hbmoy_global = float(np.mean(all_hbm))

    results = {}
    for i, pop in enumerate(population_names):
        hwm = float(np.mean(hw_per_locus[pop]))
        hwv = float(np.var(hw_per_locus[pop], ddof=1))

        fst1m = (1 - hwm / hbmoy_global) if hbmoy_global != 0 else float("nan")
        fst1v = (hwv / (hbmoy_global**2)) if hbmoy_global != 0 else float("nan")

        results[f"FST1m_{i + 1}"] = fst1m
        results[f"FST1v_{i + 1}"] = fst1v

    return results


#def _fst_wc(loci, pops):
#    """Weir & Cockerham sur un ensemble de populations -- formule générale
#    de cal_snfstd (sumstat.cpp), npop quelconque.
#    Retourne (FST2m, FST2v) via ratio sum(num)/sum(den) et Welford."""
#    x_prev = 0.0
#    xs, numt, dent = [], 0.0, 0.0
#    for lg in loci:
#        S_1 = S_2 = SSI = SSP = 0.0
#        pi_hat = [0.0, 0.0]
#        samples_data = []
#        n0 = 0
#        for pop in pops:
#            g = lg[pop]
#            n = float(len(g))
#            if n > 0:
#                s = float(sum(g))
#                p = [1 - s / n, s / n]
#                S_1 += n
#                S_2 += n * n
#                for k in range(2):
#                    pi_hat[k] += n * p[k]
#                    SSI += n * p[k] * (1 - p[k])
#                samples_data.append((n, p))
#            else:
#                n0 += 1
#        if not samples_data:
#            xs.append(x_prev)
#            continue
#        for k in range(2):
#            pi_hat[k] /= S_1
#        for n, p in samples_data:
#            for k in range(2):
#                r = p[k] - pi_hat[k]
#                SSP += n * (r * r)
#        n_d = float(len(samples_data))
#        n_c = (S_1 - S_2 / S_1) / (n_d - 1.0)
#        MSI = SSI / (S_1 - n_d)
#        MSP = SSP / (n_d - 1.0)
#        num = MSP - MSI
#        den = MSP + (n_c - 1.0) * MSI
#        if abs(den) > 0:
#            x_prev = num / den
#        xs.append(x_prev)
#        numt += num
#        dent += den
#    fstm = numt / dent if abs(dent) > 0 else 0.0
#    n = len(xs)
#    sw2diff = n * n - n
#    mx2 = sum((x - sum(xs) / n) ** 2 for x in xs)
#    fstv = mx2 * n / sw2diff if sw2diff > 0 else 0.0
#    return fstm, fstv

def _fst_wc(loci, pops):
    """Weir & Cockerham vectorisé sur tous les loci.
    Retourne (FSTm, FSTv). Formule identique à cal_snfstd, mais toutes
    les opérations par-locus sont faites en numpy sur des vecteurs de
    longueur n_loci au lieu d'une boucle Python."""
    nloci = len(loci)
    npop = len(pops)

    # Matrice des comptes d'allèle 1 : shape (npop, nloci)
    counts = np.array([[sum(lg[p]) for lg in loci] for p in pops], dtype=float)
    ns = np.array([[len(lg[p]) for lg in loci] for p in pops], dtype=float)

    p1 = counts / ns          # freq allèle 1, shape (npop, nloci)
    p0 = 1.0 - p1             # freq allèle 0

    S_1 = ns.sum(axis=0)      # (nloci,)
    S_2 = (ns**2).sum(axis=0)
    n_d = float(npop)

    # pi_hat pour chaque allèle : moyenne pondérée sur les pops
    pi0 = (ns * p0).sum(axis=0) / S_1
    pi1 = (ns * p1).sum(axis=0) / S_1

    SSI = (ns * p0 * (1 - p0) + ns * p1 * (1 - p1)).sum(axis=0)
    SSP = (ns * (p0 - pi0)**2 + ns * (p1 - pi1)**2).sum(axis=0)

    n_c = (S_1 - S_2 / S_1) / (n_d - 1.0)
    MSI = SSI / (S_1 - n_d)
    MSP = SSP / (n_d - 1.0)
    num = MSP - MSI
    den = MSP + (n_c - 1.0) * MSI

    # x = num/den quand den != 0, sinon persiste la valeur précédente
    x = np.zeros(nloci)
    valid = np.abs(den) > 0
    ratio = np.where(valid, num / np.where(valid, den, 1.0), 0.0)
    # Propagation "forward-fill" de x_prev quand den==0
    x_prev = 0.0
    xs = np.empty(nloci)
    for i in range(nloci):
        if valid[i]:
            x_prev = ratio[i]
        xs[i] = x_prev

    numt = num.sum()
    dent = den.sum()
    fstm = numt / dent if abs(dent) > 0 else 0.0

    n = nloci
    sw2diff = n * n - n
    mean = xs.mean()
    mx2 = ((xs - mean) ** 2).sum()
    fstv = mx2 * n / sw2diff if sw2diff > 0 else 0.0

    return float(fstm), float(fstv)


def compute_FST2(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """FST2 via _fst_wc -- même code que FST3/FST4."""
    from itertools import combinations

    results = {}
    for i, pa in enumerate(population_names):
        for j, pb in enumerate(population_names[i + 1 :], i + 1):
            key = f"{i + 1}.{j + 1}"
            m, v = _fst_wc(genotypes_per_locus, [pa, pb])
            results[f"FST2m_{key}"] = m
            results[f"FST2v_{key}"] = v
    return results


def compute_FST3_FST4_FSTG(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """FST3/FST4/FSTG : Weir & Cockerham généralisé à 3, 4, toutes pops.
    Référence : cal_snfstd(npop=3/4/0) dans sumstat.cpp.
    Ordre COMB pour FST3/FST4/FSTG."""
    from itertools import combinations

    results = {}
    npop = len(population_names)

    # FST3 -- triplets COMB
    for combo in combinations(range(npop), 3):
        pops = [population_names[i] for i in combo]
        key = ".".join(str(i + 1) for i in combo)
        m, v = _fst_wc(genotypes_per_locus, pops)
        results[f"FST3m_{key}"] = m
        results[f"FST3v_{key}"] = v

    # FST4 -- quadruplets COMB
    for combo in combinations(range(npop), 4):
        pops = [population_names[i] for i in combo]
        key = ".".join(str(i + 1) for i in combo)
        m, v = _fst_wc(genotypes_per_locus, pops)
        results[f"FST4m_{key}"] = m
        results[f"FST4v_{key}"] = v

    # FSTG -- toutes populations
    #m, v = _fst_wc(genotypes_per_locus, population_names)
    #results["FSTGm"] = m
    #results["FSTGv"] = v

    return results


# ----------------------------------------------------------------------------
# NEI
# ----------------------------------------------------------------------------
#def compute_NEI(
#    genotypes_per_locus: list[dict[str, list[int]]],
#    population_names: list[str],
#) -> dict[str, float]:
#    """NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire.
#    Référence : cal_snnei dans sumstat.cpp.
#    NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)
#    x_prev persiste si n==0 (comportement C++ non réinitialisé)."""
#    import math
#
#    results = {}
#    pairs = [
#        (pa, pb)
#        for i, pa in enumerate(population_names)
#        for pb in population_names[i + 1 :]
#    ]
#
#    for pa, pb in pairs:
#        x_prev = 0.0
#        xs = []
#        for lg in genotypes_per_locus:
#            ga, gb = lg[pa], lg[pb]
#            na, nb = len(ga), len(gb)
#            if na > 0 and nb > 0:
#                fi = 1 - sum(ga) / na  # freq allele 0 in pa
#                gi = sum(ga) / na  # freq allele 1 in pa
#                fj = 1 - sum(gb) / nb  # freq allele 0 in pb
#
#                gj = sum(gb) / nb  # freq allele 1 in pb
#                denom = math.sqrt(fi * fi + gi * gi) * math.sqrt(fj * fj + gj * gj)
#                if denom > 0:
#                    x_prev = 1 - (fi * fj + gi * gj) / denom
#            xs.append(x_prev)
#
#        i = population_names.index(pa) + 1
#        j = population_names.index(pb) + 1
#        key = f"{i}.{j}"
#        n = len(xs)
#        sw2diff = n * n - n
#        mx2 = sum((x - sum(xs) / n) ** 2 for x in xs)
#        results[f"NEIm_{key}"] = sum(xs) / n
#        results[f"NEIv_{key}"] = mx2 * n / sw2diff if sw2diff > 0 else 0.0
#    return results

def compute_NEI(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire, vectorisé.
    NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)
    x_prev persiste si n==0."""
    pops = population_names
    nloci = len(genotypes_per_locus)

    # Fréquences allèle 0 (f) et allèle 1 (g) pour chaque pop : (npop, nloci)
    counts = np.array([[sum(lg[p]) for lg in genotypes_per_locus] for p in pops], dtype=float)
    ns = np.array([[len(lg[p]) for lg in genotypes_per_locus] for p in pops], dtype=float)
    g = counts / ns          # freq allèle 1
    f = 1.0 - g              # freq allèle 0
    norm = np.sqrt(f*f + g*g)  # (npop, nloci)

    results = {}
    for i in range(len(pops)):
        for j in range(i+1, len(pops)):
            denom = norm[i] * norm[j]
            valid = denom > 0
            nei = np.where(valid, 1.0 - (f[i]*f[j] + g[i]*g[j]) / np.where(valid, denom, 1.0), 0.0)

            # forward-fill quand invalide
            xs = np.empty(nloci)
            x_prev = 0.0
            for k in range(nloci):
                if valid[k]:
                    x_prev = nei[k]
                xs[k] = x_prev

            key = f"{i+1}.{j+1}"
            n = nloci
            sw2diff = n*n - n
            mean = xs.mean()
            mx2 = ((xs - mean)**2).sum()
            results[f"NEIm_{key}"] = float(mean)
            results[f"NEIv_{key}"] = float(mx2 * n / sw2diff) if sw2diff > 0 else 0.0
    return results
# ----------------------------------------------------------------------------
# AML : admixture maximum likelihood sur triplets
# ----------------------------------------------------------------------------


def compute_AML(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """AMLm et AMLv : coefficient d'admixture ML sur triplets.
    Référence : cal_snaml dans sumstat.cpp.
    aml = (f3 - f2) / (f1 - f2), clampé à [0,1] si hors bornes.
    w=0 si f1==f2 (locus non informatif, exclu de la moyenne pondérée).
    Ordre des triplets : HALF (halfsortedbypairs), reproduit empiriquement.
    samp[0]=hybride, samp[1]=parent1, samp[2]=parent2."""
    from itertools import permutations, combinations

    def halfsortedbypairs(v):
        n = len(v)
        for i in range(n - 1, 0, -2):
            if not (v[i - 1] <= v[i]):
                return False
            if (i - 2) > 0 and not (v[i - 3] <= v[i - 1]):
                return False
        return True

    n = len(population_names)
    triplets = []
    for combo in sorted(combinations(range(n), 3), reverse=True):
        for perm in sorted(set(permutations(combo))):
            if halfsortedbypairs(list(perm)):
                triplets.append(list(perm))

    results = {}
    for t in triplets:
        hybrid = population_names[t[0]]
        p1 = population_names[t[1]]
        p2 = population_names[t[2]]
        key = f"{t[0] + 1}.{t[1] + 1}.{t[2] + 1}"

        x_prev = 0.0
        sw = sw2 = mx = mx2 = 0.0
        for lg in genotypes_per_locus:
            n1 = len(lg[p1])
            n2 = len(lg[p2])
            if n1 > 0 and n2 > 0:
                f1 = 1 - sum(lg[p1]) / n1  # freq allele 0 in parent1
                f2 = 1 - sum(lg[p2]) / n2  # freq allele 0 in parent2
                n3 = len(lg[hybrid])
                f3 = 1 - sum(lg[hybrid]) / n3 if n3 > 0 else 0.0
                w = 1.0
                if f1 != f2:
                    aml = (f3 - f2) / (f1 - f2)
                    x_prev = max(0.0, min(1.0, aml)) if (aml < 0 or aml > 1) else aml
                else:
                    x_prev = 0.5
                    w = 0.0
            # Welford pondéré
            if w > 0:
                sw += w
                sw2 += w * w
                mo = mx
                mx += (w / sw) * (x_prev - mo)
                mx2 += w * (x_prev - mo) * (x_prev - mx)

        sw2diff = sw * sw - sw2
        results[f"AMLm_{key}"] = float(mx)
        results[f"AMLv_{key}"] = float(mx2 * sw / sw2diff) if sw2diff > 1e-9 else 0.0

    return results


# ----------------------------------------------------------------------------
# F3-F4 : Patterson statistics sur triplets et quadruplets
# ----------------------------------------------------------------------------


def compute_F3_F4(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """F3m/F3v : statistique de Patterson f3 sur triplets (HALF).
    F4m/F4v : statistique de Patterson f4 sur quadruplets (HALF).
    Référence : cal_snf3r et cal_snf4r dans sumstat.cpp (branche SNP).
    F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)  où np = nb lignées pop hybride
    F4 = (a-b)*(c-d)
    """
    from itertools import permutations, combinations

    def halfsortedbypairs(v):
        n = len(v)
        for i in range(n - 1, 0, -2):
            if not (v[i - 1] <= v[i]):
                return False
            if (i - 2) > 0 and not (v[i - 3] <= v[i - 1]):
                return False
        return True

    def get_half_arrangements(n, r):
        result = []
        for combo in sorted(combinations(range(n), r), reverse=True):
            for perm in sorted(set(permutations(combo))):
                if halfsortedbypairs(list(perm)):
                    result.append(list(perm))
        return result

    npop = len(population_names)
    results = {}

    # --- F3 ---
    def welford_stats(xs_w):
        """Accumule Welford pondéré, retourne (mean, var) via cal_varL."""
        sw = sw2 = mx = mx2 = 0.0
        for x, w in xs_w:
            if w > 0:
                sw += w
                sw2 += w * w
                mo = mx
                mx += (w / sw) * (x - mo)
                mx2 += w * (x - mo) * (x - mx)
        sw2diff = sw * sw - sw2
        var = mx2 * sw / sw2diff if sw2diff > 1e-9 else 0.0
        return mx, var

    for t in get_half_arrangements(npop, 3):
        pop0 = population_names[t[0]]  # hybride (sample)
        pop1 = population_names[t[1]]  # parent1 (sample1)
        pop2 = population_names[t[2]]  # parent2 (sample2)
        key = f"{t[0] + 1}.{t[1] + 1}.{t[2] + 1}"
        xs_w = []
        x_prev = 0.0
        for lg in genotypes_per_locus:
            n1, n2 = len(lg[pop1]), len(lg[pop2])
            if n1 > 0 and n2 > 0:
                np_ = float(len(lg[pop0]))
                f1 = 1 - sum(lg[pop0]) / np_ if np_ > 0 else 0.0
                f2 = 1 - sum(lg[pop1]) / n1
                f3 = 1 - sum(lg[pop2]) / n2
                alpha = f1 * (1 - f1) / (np_ - 1) if np_ > 1 else 0.0
                x_prev = (f1 - f2) * (f1 - f3) - alpha
            xs_w.append((x_prev, 1.0))
        mx, vx = welford_stats(xs_w)
        results[f"F3m_{key}"] = mx
        results[f"F3v_{key}"] = vx

    # --- F4 ---
    for t in get_half_arrangements(npop, 4):
        pa = population_names[t[0]]
        pb = population_names[t[1]]
        pc = population_names[t[2]]
        pd = population_names[t[3]]
        key = f"{t[0] + 1}.{t[1] + 1}.{t[2] + 1}.{t[3] + 1}"
        xs_w = []
        x_prev = 0.0
        for lg in genotypes_per_locus:
            n1, n2, n3 = len(lg[pb]), len(lg[pc]), len(lg[pd])
            if n1 > 0 and n2 > 0 and n3 > 0:
                a = 1 - sum(lg[pa]) / len(lg[pa]) if lg[pa] else 0.0
                b = 1 - sum(lg[pb]) / n1
                c = 1 - sum(lg[pc]) / n2
                d = 1 - sum(lg[pd]) / n3
                x_prev = (a - b) * (c - d)
            xs_w.append((x_prev, 1.0))
        mx, vx = welford_stats(xs_w)
        results[f"F4m_{key}"] = mx
        results[f"F4v_{key}"] = vx

    return results


# ------------------------------------------------------------------
# Point d'entrée principal
# ---------------------------------------------------------------------------


def compute_all_statistics(
    genotypes_per_locus: list[dict[str, list[int]]],
    population_names: list[str],
) -> dict[str, float]:
    """Calcule toutes les statistiques implémentées et retourne un dict
    unifié {nom_stat: valeur} -- même format que parse_statobs().

    LIMITES ACTUELLES (à étendre au fil des validations) :
    - ML1 seulement (pas ML2/ML3)
    - HW, HB (moyenne et variance)
    - FST1 seulement (pas FST2/FST3/FST4/FSTG)
    - NEI, AML, F3, F4 : NON ENCORE IMPLÉMENTÉS
    """
    results = {}
    results.update(compute_ML1(genotypes_per_locus, population_names))
    results.update(compute_ML2(genotypes_per_locus, population_names))
    results.update(compute_ML3(genotypes_per_locus, population_names))
    results.update(compute_HW_HB(genotypes_per_locus, population_names))
    results.update(compute_FST1(genotypes_per_locus, population_names))
    results.update(compute_FST2(genotypes_per_locus, population_names))
    results.update(compute_NEI(genotypes_per_locus, population_names))
    results.update(compute_AML(genotypes_per_locus, population_names))
    results.update(compute_F3_F4(genotypes_per_locus, population_names))
    results.update(compute_FST3_FST4_FSTG(genotypes_per_locus, population_names))
    return results
