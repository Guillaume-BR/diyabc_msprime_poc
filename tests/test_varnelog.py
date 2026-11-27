"""Teste les foncitons présentes dans varnelog.py"""

import math

import msprime
import numpy as np
import pytest
from scipy import stats

from bridge.demography_builder import rescale_demography
from bridge.varnelog import (
    _check_varnelog_parameters,
    _compute_horizon_Dstar,
    add_varnelog,
    varnelog_segments,
    verhulst_N,
)


def test_verhulst_N():
    s1 = 1
    s2 = 10
    r = 3
    # Cas 1 : NF < NS croissance
    NF = 20
    NS = 50

    assert verhulst_N(s1, NF, NS, r) < verhulst_N(s2, NF, NS, r)
    assert verhulst_N(0, NF, NS, r) == pytest.approx(NF)
    assert verhulst_N(1000, NF, NS, r) == pytest.approx(NS)
    assert verhulst_N(s1, NF, NS, r=3) == pytest.approx(46.5, abs=1e-1)

    # Cas 2 : NF > NS décroissance
    NF = 50
    NS = 20
    assert verhulst_N(s1, NF, NS, r) > verhulst_N(s2, NF, NS, r)
    # Cas 3 : NF = NS constant
    NF = 20
    NS = 20
    assert verhulst_N(s1, NF, NS, r) == verhulst_N(s2, NF, NS, r)


def test_compute_horizon_Dstar():
    # Cas abs(A) > 1e-12
    NF = 10
    NS = 100
    r = 3
    Dstar = _compute_horizon_Dstar(NF, NS, r)
    assert verhulst_N(Dstar, NF, NS, r) == pytest.approx(NS / 1.001)

    # Cas abs(A) <= 1e-12
    NF = 5
    NS = 5
    r = 3
    assert _compute_horizon_Dstar(NF, NS, r) == 1


# Tests concernant _check_varnelog_parameters


def test_check_varnelog_parameters():
    # Test ok
    NF = 10
    NS = 10
    r = 0.2
    assert _check_varnelog_parameters(NF, NS, r) is None

    # r < 0
    NF = 10
    NS = 10
    r = -0.2
    with pytest.raises(ValueError, match="Le taux de croissance"):
        _check_varnelog_parameters(NF, NS, r)

    # NF < 1
    NF = 0
    NS = 10
    r = 0.2
    with pytest.raises(ValueError, match="La taille des populations"):
        _check_varnelog_parameters(NF, NS, r)

    # NS < 1
    NF = 10
    NS = -10
    r = 0.2
    with pytest.raises(ValueError, match="La taille des populations"):
        _check_varnelog_parameters(NF, NS, r)


# Tests concernant varnelog_segments

varnelog_parameters = [
    (0, 60, 70, 20, 5000, 0.2, 20),  #'normal'
    (0, 60, 70, 20, 5000, 0.01, 20),  # little r
    (0, 60, 70, 20, 5000, 10, 20),  # huge r
    (0, 60, 60, 20, 5000, 0.2, 20),  # t1=tF
    (0, 40, 50, 20, 5000, 0.2, 20),  # t1<tF
    (0, 60, 70, 20000, 5000, 0.2, 20),  # NF > NS
    (0, 60, 70, 5000, 5000, 0.2, 20),  # NF = NS
    (0, 60, 60, 20000, 5000, 0.2, 20),  # t1=tF et NF>NS
]


@pytest.mark.parametrize("t0, t1, tF, NF, NS, r, K", varnelog_parameters)
def test_varnelog_segments_are_contiguous(t0, t1, tF, NF, NS, r, K):
    results = varnelog_segments(t0, t1, tF, NF, NS, r, K)

    for k in range(len(results) - 1):
        assert results[k][1] == results[k + 1][0]
        assert results[k][0] < results[k + 1][0]
        assert results[k][1] < results[k + 1][1]
    assert results[0][0] == t0


@pytest.mark.parametrize("t0, t1, tF, NF, NS, r, K", varnelog_parameters)
def test_varnelog_segments_isolate_founding_generation(t0, t1, tF, NF, NS, r, K):
    results = varnelog_segments(t0, t1, tF, NF, NS, r, K)
    assert results[-1][1] == t1
    if t1 == tF:
        assert results[-1] == (tF - 1, tF, NF)


def test_varnelog_segments_with_wrong_fondation_time():
    with pytest.raises(ValueError, match="dépasse la fondation"):
        varnelog_segments(t0=0, t1=40, tF=30, NF=20, NS=5000, r=0.2, K=20)


@pytest.mark.parametrize("t0, t1, tF, NF, NS, r, K", varnelog_parameters)
def test_varnelog_segments_count_is_bounded(t0, t1, tF, NF, NS, r, K):
    results = varnelog_segments(t0, t1, tF, NF, NS, r, K)
    assert len(results) <= K + 2


@pytest.mark.parametrize("t0, t1, tF, NF, NS, r, K", varnelog_parameters)
def test_varnelog_segments_harmonic_mean(t0, t1, tF, NF, NS, r, K):
    results = varnelog_segments(t0, t1, tF, NF, NS, r, K)
    for a, b, N in results:
        s1 = tF - b
        s2 = tF - a
        inv = sum(1 / verhulst_N(k, NF, NS, r) for k in range(s1, s2))
        moyenne = (s2 - s1) / inv
        assert pytest.approx(moyenne, abs=0.5) == N


def test_varnelog_segments_with_negative_growth_rate():
    with pytest.raises(ValueError, match="Le taux de croissance r=-2"):
        varnelog_segments(t0=2, t1=3, tF=3, NF=3, NS=3, r=-2, K=5)


def test_varnelog_segments_with_negative_NF():
    with pytest.raises(ValueError, match="La taille des populations fondatrice"):
        varnelog_segments(t0=2, t1=3, tF=3, NF=-3, NS=3, r=2, K=5)


def test_varnelog_segments_with_negative_NS():
    with pytest.raises(ValueError, match="La taille des populations fondatrice"):
        varnelog_segments(t0=2, t1=3, tF=3, NF=3, NS=-3, r=2, K=5)


def test_add_varnelog_with_inversed_time():
    demography = msprime.Demography()
    with pytest.raises(ValueError, match="VARNELOG : il faut tL < tF"):
        add_varnelog(
            demography=demography, pop="pop3", tL=4, tF=3, NF=3, NS=3, r=2, K=3
        )


# Tests concernant add_varnelog


def _build_demography(tL, tF, NF, NS, r, mode):
    demography = msprime.Demography()
    demography.add_population(name="pop1", initial_size=NS)
    demography.add_population(name="pop2", initial_size=100)
    add_varnelog(demography, "pop1", tL, tF, NF, NS, r, mode=mode, K=20)
    demography.add_population_split(time=tF, derived=["pop1"], ancestral="pop2")
    demography.sort_events()
    return demography


def test_add_varnelog_exact_mode():
    tL = 0
    tF = 30
    NF = 20
    NS = 5000
    r = 0.2
    demography = _build_demography(tL=tL, tF=tF, NF=NF, NS=NS, r=r, mode="exact")
    epochs = demography.debug().epochs
    epochs_before_tF = [epoch for epoch in epochs if epoch.start_time < tF]
    assert len(epochs_before_tF) == tF - tL
    for epoch in epochs_before_tF:
        pop1 = epoch.populations[demography["pop1"].id]
        start_size = pop1.start_size
        assert start_size == pytest.approx(
            verhulst_N(tF - 1 - epoch.start_time, NF, NS, r)
        )

    for epoch in epochs:
        pop2 = epoch.populations[demography["pop2"].id]
        assert pop2.start_size == pop2.end_size


def test_add_varnelog_exact_mode_negative_r():
    tL = 0
    tF = 30
    NF = 20
    NS = 5000
    r = -0.2
    with pytest.raises(ValueError, match="Le taux de croissance r=-0.2"):
        _build_demography(tL, tF, NF, NS, r, mode="exact")


def test_add_varnelog_exact_mode_wrong_NF():
    tL = 0
    tF = 30
    NF = 0.5
    NS = 5000
    r = 0.2
    with pytest.raises(ValueError, match="La taille des populations fondatrice"):
        _build_demography(tL, tF, NF, NS, r, mode="exact")


def test_add_varnelog_harmonic_mode():
    tL = 0
    tF = 30
    NF = 20
    NS = 5000
    r = 0.2
    K = 20
    demography = _build_demography(tL=tL, tF=tF, NF=NF, NS=NS, r=r, mode="harmonic")
    epochs = demography.debug().epochs
    epochs_before_tF = [epoch for epoch in epochs if epoch.start_time < tF]
    assert 1 <= len(epochs_before_tF) <= K + 2
    for epoch in epochs_before_tF:
        pop1 = epoch.populations[demography["pop1"].id]
        s1 = int(tF - epoch.end_time)
        s2 = int(tF - epoch.start_time)
        inv = sum(1 / verhulst_N(k, NF, NS, r) for k in range(s1, s2))
        moyenne = (s2 - s1) / inv
        start_size = pop1.start_size
        assert start_size == pytest.approx(moyenne, abs=0.5)

    for epoch in epochs:
        pop2 = epoch.populations[demography["pop2"].id]
        assert pop2.start_size == pop2.end_size


def test_add_varnelog_after_tF():
    tL = 0
    tF = 30
    NF = 20
    NS = 5000
    r = 0.2
    for mode in ["exact", "harmonic"]:
        demography = _build_demography(tL=tL, tF=tF, NF=NF, NS=NS, r=r, mode=mode)
        epochs = demography.debug().epochs
        epochs_after_tF = [epoch for epoch in epochs if epoch.end_time >= tF]
        for epoch in epochs_after_tF:
            pop1 = epoch.populations[demography["pop1"].id]
            assert pop1.start_size == pop1.end_size

        for epoch in epochs:
            pop2 = epoch.populations[demography["pop2"].id]
            assert pop2.start_size == pop2.end_size


def test_add_varnelog_with_positive_tL():
    tL = 10
    tF = 50
    NF = 20
    NS = 5000
    r = 0.2
    for mode in ["exact", "harmonic"]:
        demography = _build_demography(tL=tL, tF=tF, NF=NF, NS=NS, r=r, mode=mode)
        epochs = demography.debug().epochs
        epochs_before_tL = [epoch for epoch in epochs if epoch.start_time < tL]

        for epoch in epochs_before_tL:
            pop1 = epoch.populations[demography["pop1"].id]
            assert pop1.start_size == NS

        for epoch in epochs:
            pop2 = epoch.populations[demography["pop2"].id]
            assert pop2.start_size == pop2.end_size


def test_add_varnelog_exact_long_horizon():
    tL = 10
    tF = 100
    NF = 20
    NS = 5000
    r = 0.2

    demography = _build_demography(tL=tL, tF=tF, NF=NF, NS=NS, r=r, mode="exact")
    epochs = demography.debug().epochs
    pop1_id = demography["pop1"].id

    kmax = math.ceil(_compute_horizon_Dstar(NF, NS, r))
    # Entre tL et tF - kmax (exclu) : époque ajoutée à NS.
    epochs_at_NS = [e for e in epochs if tL <= e.start_time < tF - kmax]
    assert len(epochs_at_NS) >= 1
    for epoch in epochs_at_NS:
        assert epoch.populations[pop1_id].start_size == NS

    # À tF - kmax : première génération de la phase logistique, k = kmax - 1.
    junction = [e for e in epochs if e.start_time == tF - kmax]
    assert len(junction) == 1
    size = junction[0].populations[pop1_id].start_size
    assert size == pytest.approx(verhulst_N(kmax - 1, NF, NS, r))
    # Au-delà de l'horizon (k = kmax >= Dstar), N est à moins de 0,1 % de NS :
    # c'est ce qui autorise à remplacer la suite par NS.
    assert verhulst_N(kmax, NF, NS, r) == pytest.approx(NS, rel=1e-3)


def _fraction_colaesced_with_varnelog(tL, tF, NF, NS, r, mode, ploidy, seed, factor):
    n_replicates = 20000
    demography = _build_demography(tL, tF, NF, NS, r, mode)
    if factor != 1.0:
        demography = rescale_demography(demography, factor)
    n_individuals = 1 if ploidy == 2 else 2  # toujours 2 lignées échantillonnées
    replicates = msprime.sim_ancestry(
        samples={"pop1": n_individuals},
        demography=demography,
        ploidy=ploidy,
        sequence_length=1,
        num_replicates=n_replicates,
        random_seed=seed,
    )
    tmrca = np.array([ts.first().tmrca(0, 1) for ts in replicates])
    before_foundation = tmrca[tmrca <= tF]
    return len(before_foundation), n_replicates


params = [(0, "exact"), (0, "harmonic"), (10, "exact"), (10, "harmonic")]


@pytest.mark.parametrize("tL,mode", params)
def test_varnelog_effect(tL, mode):
    tF = 100
    NF = 20
    NS = 5000
    r = 0.2
    n_coalesced, n = _fraction_colaesced_with_varnelog(
        tL, tF, NF, NS, r, mode=mode, ploidy=2, seed=42, factor=1
    )

    if tL > 0:
        somme = sum(1 / (2 * verhulst_N(k, NF, NS, r)) for k in range(tF - tL)) + tL / (
            2 * NS
        )
    else:
        somme = sum(1 / (2 * verhulst_N(k, NF, NS, r)) for k in range(tF))

    p = 1 - math.exp(-somme)
    result = stats.binomtest(n_coalesced, n, p)

    assert result.pvalue > 1e-3, (
        f"{n_coalesced}/{n} coalescences à t={tF}, attendu {p:.3f}"
    )


def test_add_varnelog_with_null_growth_rate():
    demography = msprime.Demography()
    with pytest.raises(ValueError, match="Le taux de croissance r=0"):
        add_varnelog(
            demography=demography, pop="pop3", tL=3, tF=4, NF=3, NS=3, r=0, K=3
        )


def test_add_varnelog_with_wrong_mode():
    demography = msprime.Demography()
    with pytest.raises(ValueError, match="mode = 'harmonic' ou 'exact'"):
        add_varnelog(
            demography=demography,
            pop="pop3",
            tL=3,
            tF=4,
            NF=3,
            NS=3,
            r=0.2,
            K=3,
            mode="coucou",
        )
