"""
Tirage des valeurs numériques des priors, avec retirage si les contraintes
d'ordre (ex: "t4>t3") ne sont pas respectées -- équivalent du mécanisme
"DRAW UNTIL" observé dans header.txt.

Limité pour l'instant à la loi UN (uniforme). Les autres lois (LU, GA --
voir doc DIYABC sur "mean and standard deviation") ne sont pas encore
vérifiées : human n'utilise que UN, donc pas bloquant pour ce POC, mais
à traiter explicitement avant de généraliser à d'autres datasets.

Les priors de catégorie N (taille) et T (temps) sont arrondis à
l'entier le plus proche juste après le tirage, comme DIYABC
(particuleC.cpp, voir _draw_one_value) -- seul le taux d'admixture (A)
reste continu.
"""

import math
import random

from bridge.scenario_types import OrderConstraint, Prior, Scenario


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
    if prior.law == "UN":
        min_, max_, *_ = prior.bounds
        value = rng.uniform(min_, max_)
    else:
        raise NotImplementedError(
            f"Loi '{prior.law}' non implémentée (seule 'UN' est supportée pour "
            f"l'instant). Prior concerné : {prior.name!r}"
        )

    if prior.category in ("N", "T"):
        # DIYABC arrondit à l'entier le plus proche les priors de taille
        # (N) et de temps (T) juste après le tirage -- particuleC.cpp :
        # "if (category<2) value = floor(0.5+value)" (round-half-up),
        # category 0=N, 1=T, 2=A (histparam.category, header.cpp). Le
        # taux d'admixture (A) reste continu, jamais arrondi.
        value = float(math.floor(0.5 + value))

    return value


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
