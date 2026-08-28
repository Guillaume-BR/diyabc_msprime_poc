"""Vérifie prior_parser : extraction des priors et contraintes d'ordre
depuis header.txt, et la règle de filtrage des priors quasi-constants."""

import pytest

from bridge.header_dataclasses import GroupPrior, OrderConstraint, Prior
from bridge.prior_parser import (
    get_parameter_used_by_model,
    is_constant_prior,
    parse_group_priors,
    parse_priors,
)


def test_priors_and_constraints(header_text):
    priors, constraints = parse_priors(header_text)

    assert len(priors) == 21
    assert len(constraints) == 4

    priors_by_name = {p.name: p for p in priors}
    assert priors_by_name["N1"].category == "N"
    assert priors_by_name["N1"].law == "UN"
    assert (priors_by_name["N1"].min, priors_by_name["N1"].max) == (1000.0, 100000.0)

    assert priors_by_name["t1"].category == "T"
    assert (priors_by_name["t1"].min, priors_by_name["t1"].max) == (1.0, 30.0)

    assert OrderConstraint(param1="t4", operator=">", param2="t3") in constraints
    assert OrderConstraint(param1="t3", operator=">", param2="t2") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t33") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t22") in constraints


def test_parse_priors_no_draw_until(header_text_te1):
    """Vérifie le parsing de la section 'historical parameters priors' de
    toy_example1, qui n'a pas de section 'DRAW UNTIL'."""
    priors, constraints = parse_priors(header_text_te1)

    assert len(priors) == 3
    assert constraints == []


def test_parse_group_priors(header_text_te2):
    group_priors = parse_group_priors(header_text_te2)
    assert len(group_priors) == 3
    assert group_priors["G1"][0].name == "MEANMU"
    assert group_priors["G1"][0].law == "UN"
    assert group_priors["G3"][-1].p_fixe == 10
    assert group_priors["G3"][-1].gams == 2
    assert group_priors["G3"][-1].model is True


def test_is_constant_prior():
    """Vérifie la règle de filtrage des priors quasi-constants."""
    # Cas normal : large intervalle, jamais constant
    normal_prior = Prior(
        name="N1",
        category="N",
        law="UN",
        min=1000.0,
        max=100000.0,
        mean=0.0,
        sdshape=0.0,
    )
    assert is_constant_prior(normal_prior) is False

    # Cas dégénéré : min == max, clairement constant
    constant_prior = Prior(
        name="X", category="N", law="UN", min=100.0, max=100.0, mean=0.0, sdshape=0.0
    )
    assert is_constant_prior(constant_prior) is True

    # Cas limite : différence infime, sous le seuil
    near_constant_prior = Prior(
        name="Y",
        category="N",
        law="UN",
        min=100.0,
        max=100.00001,
        mean=0.0,
        sdshape=0.0,
    )
    assert is_constant_prior(near_constant_prior) is True

    # Cas limite inverse : différence juste au-dessus du seuil
    barely_variable_prior = Prior(
        name="Z", category="N", law="UN", min=100.0, max=100.1, mean=0.0, sdshape=0.0
    )
    assert is_constant_prior(barely_variable_prior) is False


def test_get_parameter_used_by_model(header_text_te2):
    """Vérifie que la fonction get_parameter_used_by_model retourne bien
    le bon tuple de booléen pour le nom de modèle donné.
    """
    group_priors = parse_group_priors(header_text_te2)
    last_prior = group_priors["G3"][-1]
    assert get_parameter_used_by_model(last_prior) == (True, False)

    with pytest.raises(NotImplementedError, match="Le modèle mutationnel"):
        get_parameter_used_by_model(group_priors["G1"][0])

    prior_without_name_model = GroupPrior(
        group="G1",
        ms_or_seq=None,
        name=None,
        law=None,
        min=None,
        max=None,
        mean=None,
        sdshape=None,
        model=True,
        name_model=None,
        p_fixe=5,
        gams=6,
    )
    with pytest.raises(NotImplementedError, match="Le nom du modèle"):
        get_parameter_used_by_model(prior_without_name_model)
