"""Teste les foncitons présentes dans varnelog.py"""

import msprime
import pytest

from bridge.varnelog import (
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
    add_varnelog(demography, "pop1", tL, tF, NF, NS, r, mode=mode)
    demography.add_population_split(time=tF, derived=["pop1"], ancestral="pop2")
    return demography


def test_add_varnelog_exact_mode():
    # demography = _build_demography(tL=0, tF=30, NF=20, NS=5000, r=0.2)
    # epochs = demography.debug().epochs
    pass


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
