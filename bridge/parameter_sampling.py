"""
Tirage des valeurs numériques des priors, avec retirage si les contraintes
d'ordre (ex: "t4>t3") ne sont pas respectées -- équivalent du mécanisme
"DRAW UNTIL" observé dans header.txt.

Les priors de catégorie N (taille) et T (temps) sont arrondis à
l'entier le plus proche juste après le tirage, comme DIYABC
(particuleC.cpp, voir _draw_one_value) -- seul le taux d'admixture (A)
reste continu.

Tirage des paramètres de group priors, avec gestion des dépendances entre les priors
d'un même groupe (ex: MEANMU et GAMMU pour la loi GA)
"""

import dataclasses
import math
import random

from bridge.scenario_types import (
    GroupPrior,
    LociDescriptionDetailed,
    OrderConstraint,
    Prior,
    Scenario,
)

# Décalage de graine pour draw_group_parameter_values, distinct de celui utilisé
# par draw_parameter_values (qui n'en a pas besoin, une seule instance de
# random.Random par appel) -- évite que les deux tirages soient corrélés si un
# appelant utilise la même seed de base pour les deux, même bug que celui déjà
# rencontré et corrigé le 2026-07-16 (corrélation scénario/premier prior en
# multi-scénario). Valeur choisie loin des autres offsets déjà utilisés dans le
# projet (_MAF_REJECTION_SEED_OFFSET=2_000_000, _MRC_REJECTION_SEED_OFFSET=
# 3_000_000, _BINOMIAL_SEED_OFFSET=4_000_000 dans ancestry_simulation.py ;
# _SCENARIO_DRAW_SEED_OFFSET=50_000_000 dans reftable_loop.py).
_GROUP_PRIOR_SEED_OFFSET = 10_000_000


class ConstraintsNotSatisfiedError(Exception):
    """Levée si aucun tirage valide n'a été trouvé en max_attempts essais."""


def draw_scenario(scenarios: list[Scenario], seed: int) -> Scenario:
    """Tire un scénario parmi `scenarios`, pondéré par son `weight`
    (le "prior_proba" de particuleC.cpp::ParticleC::drawscenario).

    Reproduit exactement l'algorithme C++ : tirage d'un nombre uniforme
    `ra` dans [0,1), puis balayage de la somme cumulée des poids jusqu'à
    ce qu'elle atteigne ou dépasse `ra` -- même logique d'inversion de
    CDF que _draw_single_mutation_edge_child (ancestry_simulation.py),
    vérifiée boundary-compatible avec la boucle C++ ("while ra > sp").

    Ne normalise PAS les poids (comme le C++, qui ne le fait pas non
    plus) : si leur somme est < 1, le DERNIER scénario de la liste sert
    de secours pour tout `ra` au-delà de la somme cumulée -- même
    comportement que la boucle C++, bornée à nscenarios-1.
    """
    rng = random.Random(seed)
    ra = rng.random()

    cumulative = 0.0
    for scenario in scenarios:
        cumulative += scenario.weight
        if ra <= cumulative:
            return scenario

    return scenarios[-1]


def _draw_one_value(prior: Prior, rng: random.Random) -> float:
    """Tire une valeur pour un prior donné, selon sa loi et ses bornes.
    Lève NotImplementedError si la loi n'est pas encore implémentée.
    """
    if prior.bounds[0] == prior.bounds[1]:
        # bornes identiques : pas de tirage, valeur fixée
        value = float(prior.bounds[0])
    elif prior.law == "UN":  # uniform
        min_, max_, *_ = prior.bounds
        value = rng.uniform(min_, max_)
    elif prior.law == "LU":  # log-uniform
        min_, max_, *_ = prior.bounds
        value = math.exp(rng.uniform(math.log(min_), math.log(max_)))
    elif prior.law == "NO":  # normal
        min_, max_, mean, sdshape = prior.bounds
        while True:
            value = rng.gauss(mean, sdshape)
            if min_ <= value <= max_:
                break
    elif prior.law == "LN":  # log-gaussian
        min_, max_, mean, sdshape = prior.bounds
        while True:
            value = math.exp(rng.gauss(math.log(mean), math.log(sdshape)))
            if min_ <= value <= max_:
                break
    elif prior.law == "GA":  # gamma
        min_, max_, mean, sdshape = prior.bounds
        if mean < 1e-12:
            value = 0.0
        elif sdshape < 1e-12:
            value = mean
        elif max_ < 1e-12:
            value = max_
        else:
            while True:
                value = rng.gammavariate(sdshape, mean / sdshape)
                if min_ <= value <= max_:
                    break
    else:
        raise NotImplementedError(
            f"Loi '{prior.law}' non implémentée pour le tirage de valeurs numériques des priors historiques."
        )

    if prior.category in ("N", "T"):
        # DIYABC arrondit à l'entier le plus proche les priors de taille
        # (N) et de temps (T) juste après le tirage -- particuleC.cpp :
        # "if (category<2) value = floor(0.5+value)" (round-half-up),
        # category 0=N, 1=T, 2=A (histparam.category, header.cpp). Le
        # taux d'admixture (A) reste continu, jamais arrondi.
        value = float(math.floor(0.5 + value))

    return value


def _draw_one_group_value(group_prior: GroupPrior, rng: random.Random) -> float:
    """Tire une valeur pour un group prior donné, selon sa loi et ses bornes.
    Quelques gardes-fous pour éviter des erreurs de tirage si le group prior est mal défini.
    Lève ValueError si le group prior n'a pas de loi ou de bornes associées, ou si la loi est GA mais que les priors MEAN et SDSHAPE
    correspondants n'ont pas été tirés avant.
    Lève NotImplementedError si la loi n'est pas encore implémentée via _draw_one_value.
    """
    if group_prior.law is None:
        raise ValueError(
            f"Le group prior {group_prior.name!r} n'a pas de loi associée."
        )
    if group_prior.min is None:
        raise ValueError(
            f"Le group prior {group_prior.name!r} n'a pas de bornes associées."
        )

    if group_prior.law == "GA" and group_prior.mean is None:
        raise ValueError(
            f"Le group prior {group_prior.name!r} a une loi GA mais pas de moyenne associée."
            f"Il faut d'abord tirer la valeau du prior MEAN correspondant du même groupe"
        )

    if group_prior.law == "GA" and group_prior.sdshape is None:
        raise ValueError(
            f"Le group prior {group_prior.name!r} a une loi GA mais pas de moyenne associée."
            f"Il faut d'abord tirer la valeau du prior SDSHAPE correspondant du même groupe"
        )
    bounds = [group_prior.min, group_prior.max, group_prior.mean, group_prior.sdshape]
    prior = Prior(
        name=group_prior.name,
        category="G",  # catégorie fictive pour les group priors
        law=group_prior.law,
        bounds=bounds,
    )
    return _draw_one_value(prior, rng)


def draw_parameter_values(
    priors: list[Prior],
    constraints: list[OrderConstraint],
    seed: int,
    max_attempts: int = 1000,
) -> dict[str, float]:
    """Tire une valeur pour chaque prior, en retirant tant que les
    contraintes d'ordre ne sont pas toutes satisfaites.

    Lève ConstraintsNotSatisfiedError si aucun tirage valide n'est trouvé
    en max_attempts essais -- signe probable d'une configuration de
    contraintes incohérente (bornes de priors incompatibles avec les
    contraintes demandées) plutôt que d'une simple mauvaise chance.
    """
    rng = random.Random(seed)

    for _ in range(max_attempts):
        values = {p.name: _draw_one_value(p, rng) for p in priors}
        if all(c.is_satisfied(values) for c in constraints):
            return values

    raise ConstraintsNotSatisfiedError(
        f"Aucun tirage valide trouvé en {max_attempts} essais. "
        f"Vérifier que les contraintes ({len(constraints)}) sont "
        f"compatibles avec les bornes des priors."
    )


def draw_group_parameter_values(
    group_priors: dict[str, list[GroupPrior]],
    seed: int,
) -> dict[str, dict[str, float]]:
    """
    Tire une valeur pour chaque group prior ou bien des valeurs pour le modèle.

    `seed` est décalé de _GROUP_PRIOR_SEED_OFFSET avant utilisation -- ne
    corrèle jamais ce tirage avec celui de draw_parameter_values, même si
    l'appelant leur passe la même seed de base.

    Point d'attention : ce tirage positionnel entre priors dépendants
    (MEANMU avant GAMMU, etc.) reproduit le comportement de DIYABC, qui lit
    ces lignes dans un ordre fixe sans jamais regarder leur nom (voir
    header.cpp::readHeadersimGroupPrior).
    """
    rng = random.Random(seed + _GROUP_PRIOR_SEED_OFFSET)
    group_priors_values: dict[str, dict[str, float]] = {}
    for group_name in group_priors:
        values = {}
        last_value = None
        for gp in group_priors[group_name]:
            if gp.model:
                continue
            if gp.mean is None:
                gp = dataclasses.replace(
                    gp, mean=last_value
                )  # crée une copie de gp avec la valeur de mean remplacée par last_value
            value = _draw_one_group_value(gp, rng)
            values[gp.name] = value
            last_value = value
        group_priors_values[group_name] = values
    return group_priors_values


def sampling_kappa_per_locus(
    group_prior: GroupPrior,
    k_moy: float,
    n_loci: int,
    check_nloc: bool,
    list_loci: list[LociDescriptionDetailed],
    rng: random.Random,
) -> dict[str, float]:
    """Échantillonne les valeurs de kappa pour chaque groupe de loci en fonction des priorités définies.
    Retourne un dictionnaire avec les noms de groupes comme clés et les valeurs de kappa échantillonnées comme valeurs.
    Pour kappa1, il faut passsé l'argument check_nloc à True pour vérifier le nombre de loci et échantillonner en conséquence.
    Pour kappa2, il faut passsé l'argumentââ check_nloc à False pour
    """

    kappa_values = {}
    if check_nloc:
        if group_prior.sdshape > 0.001 and n_loci > 1:
            group_prior = dataclasses.replace(group_prior, mean=k_moy)
            for locus in list_loci:
                kappa = _draw_one_group_value(group_prior, rng)
                kappa_values[locus.name] = kappa
        else:
            for locus in list_loci:
                kappa_values[locus.name] = k_moy
    else:
        if group_prior.sdshape > 0.001:
            group_prior = dataclasses.replace(group_prior, mean=k_moy)
            for locus in list_loci:
                kappa = _draw_one_group_value(group_prior, rng)
                kappa_values[locus.name] = kappa
        else:
            for locus in list_loci:
                kappa_values[locus.name] = k_moy
    return kappa_values
