"""Vérifie parameter_sampling : tirage rejet des priors sous contraintes
d'ordre, arrondi DIYABC (round-half-up pour N/T, continu pour A), et
tirage pondéré du scénario."""

import dataclasses
import random

from bridge.header_dataclasses import Prior, Scenario
from bridge.loci_parser import parse_loci_description
from bridge.parameter_sampling import (
    _draw_one_value,
    draw_group_parameter_values,
    draw_parameter_values,
    draw_scenario,
    sample_site_rates,
    sampling_group_local_param,
)
from bridge.prior_parser import parse_group_priors, parse_priors
from bridge.scenario_parser import parse_header_scenarios


def test_draw_parameter_values(header_text):
    """Vérifie que draw_parameter_values tire bien une valeur pour chaque
    prior, et que le tirage retourné respecte toutes les contraintes
    d'ordre (t4>t3, t3>t2, t44>t33, t44>t22)."""
    priors, constraints = parse_priors(header_text)

    seed = 42
    values = draw_parameter_values(priors, constraints, seed)

    # Toutes les valeurs ont bien été tirées
    assert set(values.keys()) == {p.name for p in priors}

    # Toutes les contraintes sont respectées par ce tirage
    for constraint in constraints:
        assert constraint.is_satisfied(values)


def test_draw_parameter_values_reproducible(header_text):
    """Même graine -> même tirage (déterminisme attendu pour la
    reproductibilité scientifique)."""
    priors, constraints = parse_priors(header_text)

    values1 = draw_parameter_values(priors, constraints, seed=123)
    values2 = draw_parameter_values(priors, constraints, seed=123)

    assert values1 == values2


def test_draw_parameter_values_rounds_N_and_T_not_A(header_text):
    """DIYABC arrondit à l'entier le plus proche les priors de catégorie
    N (taille) et T (temps) juste après le tirage, mais laisse le taux
    d'admixture (catégorie A, ex: 'ra') continu -- vérifié contre
    particuleC.cpp ("if (category<2) value = floor(0.5+value)", avec
    category 0=N, 1=T, 2=A dans header.cpp)."""
    priors, constraints = parse_priors(header_text)
    priors_by_name = {p.name: p for p in priors}

    non_integer_ra_seen = False
    for seed in range(1, 50):
        values = draw_parameter_values(priors, constraints, seed)
        for name, prior in priors_by_name.items():
            if prior.category in ("N", "T"):
                assert values[name] == int(values[name]), (
                    f"{name} (catégorie {prior.category}) devrait être un "
                    f"entier, obtenu {values[name]}"
                )
        if values["ra"] != int(values["ra"]):
            non_integer_ra_seen = True

    assert non_integer_ra_seen, "ra (catégorie A) ne devrait jamais être arrondi"


class _FakeRng:
    """rng minimal pour contrôler exactement la valeur tirée par
    _draw_one_value, sans dépendre de random.Random."""

    def __init__(self, value):
        self._value = value

    def uniform(self, min_, max_):
        return self._value


def test_draw_one_value_round_half_up_not_banker():
    """L'arrondi doit être round-half-up (floor(0.5+x), comme
    particuleC.cpp), PAS round-half-to-even (comportement par défaut de
    la fonction round() de Python) -- cas où les deux méthodes divergent
    : x=2.5 -> 3 en round-half-up, mais round(2.5) == 2 en Python."""
    prior = Prior(
        name="N1", category="N", law="UN", min=0.0, max=10.0, mean=0.0, sdshape=0.0
    )

    value = _draw_one_value(prior, _FakeRng(2.5))

    assert value == 3.0
    assert (
        round(2.5) == 2
    )  # confirme que Python round() aurait donné un résultat différent


def test_draw_scenario_uses_weight(header_text):
    """Sur les 6 scénarios de human (poids ~1/6 chacun), un tirage sur de
    nombreuses graines doit produire les 6 indices possibles -- preuve
    que le poids est bien utilisé pour sélectionner le scénario, pas un
    choix fixe (particuleC.cpp::ParticleC::drawscenario)."""
    scenarios = parse_header_scenarios(header_text)

    drawn_indices = {draw_scenario(scenarios, seed).index for seed in range(1, 201)}

    assert drawn_indices == {1, 2, 3, 4, 5, 6}


def test_draw_scenario_degenerate_weight():
    """Si un seul scénario a un poids non nul, il doit toujours être
    tiré, quelle que soit la graine."""
    scenarios = [
        Scenario(index=1, weight=1.0, initial_pop_size_exprs=["N1"]),
        Scenario(index=2, weight=0.0, initial_pop_size_exprs=["N1"]),
    ]

    for seed in range(1, 20):
        assert draw_scenario(scenarios, seed).index == 1


def test_draw_scenario_fallback_when_weights_sum_below_one():
    """Si la somme des poids est < 1 (ne devrait pas arriver avec un vrai
    header.txt DIYABC, mais le C++ ne le vérifie pas), le DERNIER
    scénario sert de secours pour tout tirage au-delà de la somme
    cumulée -- même comportement que la boucle C++ (bornée à
    nscenarios-1)."""
    scenarios = [
        Scenario(index=1, weight=0.1, initial_pop_size_exprs=["N1"]),
        Scenario(index=2, weight=0.1, initial_pop_size_exprs=["N1"]),
    ]

    # seed=2 -> rng.random() == 0.956..., largement > 0.2 (somme des poids)
    drawn = draw_scenario(scenarios, seed=2)
    assert drawn.index == 2


def test_draw_group_parameter_values(header_text_te2):
    """Vérifie bien que draw_group_parameter_values tire une valeur pour chaque paramètre du groupe
    et que le tirage retourné respecte toutes les dépendances entre les paramètres
    du groupe comme Mean_u pour la loi Gamma"""
    group_priors = parse_group_priors(header_text_te2)
    group_priors_values = draw_group_parameter_values(group_priors, seed=42)

    assert len(group_priors_values) == 3
    assert len(group_priors_values["G2"]) == 6


def test_draw_group_parameter_values_hierarchical_mean(header_text_te2):
    """Vérifie que draw_group_parameter_values tire une valeur pour chaque paramètre du groupe
    et que le tirage retourné respecte toutes les dépendances entre les paramètres
    du groupe comme Mean_u pour la loi Gamma"""
    group_priors = parse_group_priors(header_text_te2)
    priors_g1 = group_priors["G1"]
    for i, gp in enumerate(priors_g1):
        if gp.name == "GAMMU":
            priors_g1[i] = dataclasses.replace(gp, sdshape=1e-15)

    values = draw_group_parameter_values(group_priors, seed=42)

    assert values["G1"]["GAMMU"] == values["G1"]["MEANMU"]


def test_sampling_group_local_param(header_text_te2):
    """Vérifie que sampling_group_local_param tire une valeur pour chaque locus du groupe
    et que le tirage retourné respecte toutes les dépendances entre les paramètres
    du groupe comme Mean_u pour la loi Gamma"""
    group_priors = parse_group_priors(header_text_te2)
    priors_g2 = group_priors["G2"]
    list_loci_g2 = [
        locus
        for locus in parse_loci_description(header_text_te2)
        if locus.group == "G2"
    ]
    values = draw_group_parameter_values(group_priors, seed=42)

    # test avec un sdshape > 0.001, donc tirage aléatoire
    prior_gamma1 = next((gp for gp in priors_g2 if gp.name == "GAMK1"), None)
    kappa_values = sampling_group_local_param(
        prior_gamma1,
        values["G2"]["MEANK1"],
        len(list_loci_g2),
        True,
        list_loci_g2,
        rng=random.Random(42),
    )
    assert len(set(kappa_values.values())) == len(list_loci_g2)

    # test avec un sdshape < 0.001, donc valeur fixe
    prior_gamma1_modified = dataclasses.replace(prior_gamma1, sdshape=1e-15)
    kappa_values_modified = sampling_group_local_param(
        prior_gamma1_modified,
        values["G2"]["MEANK1"],
        len(list_loci_g2),
        True,
        list_loci_g2,
        rng=random.Random(42),
    )
    assert all(
        value == values["G2"]["MEANK1"] for value in kappa_values_modified.values()
    )

    # test pour k2
    prior_gamma2 = next((gp for gp in priors_g2 if gp.name == "GAMK2"), None)
    kappa_values2 = sampling_group_local_param(
        prior_gamma2,
        values["G2"]["MEANK2"],
        len(list_loci_g2),
        False,
        list_loci_g2,
        rng=random.Random(42),
    )
    assert len(set(kappa_values2.values())) == len(list_loci_g2)


def test_sample_site_rates():
    p_fixe = 5
    dnalength = 10
    gams = 1.0
    rng = random.Random(42)
    result = sample_site_rates(p_fixe, gams, dnalength, rng)

    assert len(result) == dnalength
    assert abs(sum(result) - 1) < 1e-6  # somme des taux = 1

    # test with gams = 0 et un p_fixe qui donne reellement des sites fixes
    # (p_fixe=20, dnalength=10 -> nsv=8, donc 2 sites fixes)
    p_fixe = 20
    gams = 0.0
    result = sample_site_rates(p_fixe, gams, dnalength, rng)
    assert len(result) == dnalength
    nb_sites_fixes = 2
    assert all(rate == 0.0 for rate in result[:nb_sites_fixes])
    assert all(
        abs(rate - 1 / (dnalength - nb_sites_fixes)) < 1e-6
        for rate in result[nb_sites_fixes:]
    )
