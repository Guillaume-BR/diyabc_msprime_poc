"""Vérifie prior_parser : extraction des priors et contraintes d'ordre depuis
header.txt, et la règle de filtrage des priors quasi-constants."""

import pytest

from bridge.header_dataclasses import DirichletPrior, GroupPrior, OrderConstraint, Prior
from bridge.prior_parser import (
    get_parameter_used_by_model,
    is_constant_prior,
    parse_group_priors,
    parse_priors,
)


def test_parse_priors_and_constraints(header_text):
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


# def test_parse_priors_with_invalid_bounds(header_text):
#    # min > max
#    header_text_invalid1 = header_text.replace(
#        "N1 N UN[1000.0,100000.0,0.0,0.0]", "N1 N UN[100000.0,1000.0,0.0,0.0]"
#    )
#    with pytest.raises(ValueError, match="borne min > max"):
#        parse_priors(header_text_invalid1)

# min == max
# header_text_invalid2 = header_text.replace(
#    "N1 N UN[1000.0,100000.0,0.0,0.0]", "N1 N UN[50000.0,50000.0,0.0,0.0]"
# )
# with pytest.raises(ValueError, match="Domaine vide ou réduit à un point"):
#    parse_priors(header_text_invalid2)


def test_parse_priors_with_dirichlet_two_parameters(header_text):
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra A DI[0.05,0.95,0.0,0.0]"
    )
    priors, _ = parse_priors(header_text)

    dirichlet_priors = next(p for p in priors if isinstance(p, DirichletPrior))

    assert dirichlet_priors.names == ["ra"]
    assert dirichlet_priors.category == "A"
    assert dirichlet_priors.min == 0.05
    assert dirichlet_priors.max == 0.95


def test_parse_priors_with_dirichlet_three_parameters(header_text):
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    # insertion de ra,rb A DI[0.05,0.95,0.0,0.0]
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )
    priors, _ = parse_priors(header_text)

    dirichlet_priors = next(p for p in priors if isinstance(p, DirichletPrior))

    assert dirichlet_priors.names == ["ra", "rb"]
    assert dirichlet_priors.category == "A"
    assert dirichlet_priors.min == 0.05
    assert dirichlet_priors.max == 0.95


def test_parse_priors_with_invalid_dirichlet(header_text):

    # Borne min
    header_text_invalid1 = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra A DI[0.6,0.95,0.0,0.0]"
    )
    with pytest.raises(ValueError, match="Domaine vide ou réduit à un point"):
        parse_priors(header_text_invalid1)

    header_text_invalid1_bis = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra A DI[0.5,0.95,0.0,0.0]"
    )
    with pytest.raises(ValueError, match="Domaine vide ou réduit à un point"):
        parse_priors(header_text_invalid1_bis)

    # borne max
    header_text_invalid2 = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra A DI[0.05,0.25,0.0,0.0]"
    )
    with pytest.raises(
        ValueError, match="Les parts ne peuvent sommer à 1 avec la borne max"
    ):
        parse_priors(header_text_invalid2)

    header_text_invalid2_bis = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra A DI[0.05,0.5,0.0,0.0]"
    )
    with pytest.raises(
        ValueError, match="Les parts ne peuvent sommer à 1 avec la borne max"
    ):
        parse_priors(header_text_invalid2_bis)

    # Mauvais type
    header_text_invalid3 = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra T DI[0.05,0.25,0.0,0.0]"
    )
    with pytest.raises(ValueError, match="Le prior DI ne supporte"):
        parse_priors(header_text_invalid3)

    # Prior Dirichlet avec deux fois le même nom
    header_text_invalid4 = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,ra A DI[0.05,0.95,0.0,0.0]"
    )
    with pytest.raises(
        ValueError,
        match="Un même nom apparaît plusieurs fois dans",
    ):
        parse_priors(header_text_invalid4)

    # header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    line_index = header_text.find("ra A UN[0.05,0.95,0.0,0.0]")
    header_text_invalid5 = (
        header_text[:line_index]
        + "ra,rb A DI[0.05,0.95,0.0,0.0] \n"
        + header_text[line_index:]
    )
    with pytest.raises(
        ValueError, match="Les noms de priors classiques et de Dirichlet DI"
    ):
        parse_priors(header_text_invalid5)


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


def test_parse_group_priors_model_error(header_text_te2):
    """Vérifie que parse_group_priors lève une exception si un paramètre de
    modèle est présent mais que son nom n'est pas indiqué."""
    # On modifie le header pour supprimer le nom du modèle
    header_text_modified = header_text_te2.replace("MODEL", "MODLE")
    with pytest.raises(ValueError, match="Ligne de modèle de la section"):
        parse_group_priors(header_text_modified)

    header_text_modified2 = header_text_te2.replace("MODEL K2P 10 2.00", "MODEL")
    assert header_text_modified2 != header_text_te2
    with pytest.raises(ValueError, match="Ligne de modèle de la section group priors"):
        parse_group_priors(header_text_modified2)

    header_text_modified2 = header_text_te2.replace("MODEL K2P 10 2.00", "MODEL K2P 10")
    assert header_text_modified2 != header_text_te2
    with pytest.raises(ValueError, match="malformée"):
        parse_group_priors(header_text_modified2)

    header_text_modified3 = header_text_te2.replace(
        "MODEL K2P 10 2.00", "MODEL K2P x 2.00"
    )
    assert header_text_modified3 != header_text_te2
    with pytest.raises(ValueError, match="malformée"):
        parse_group_priors(header_text_modified3)

    # on place "MODEL" avant les priors de groupe, ce qui est incorrect
    lines = header_text_te2.splitlines()
    i_group = lines.index("group G2 [S]")
    i_model = next(
        i for i in range(i_group, len(lines)) if lines[i].startswith("MODEL")
    )
    model_line = lines.pop(i_model)  # on retire la ligne MODEL de sa place
    lines.insert(i_group + 1, model_line)  # et on la met juste après "group G2 [S]"
    header_model_first = "\n".join(lines)
    with pytest.raises(
        ValueError,
        match="La section 'MODEL' doit être placée après les sections 'group'",
    ):
        parse_group_priors(header_model_first)


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
    """Vérifie que la fonction get_parameter_used_by_model retourne bien le bon
    tuple de booléen pour le nom de modèle donné."""
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
