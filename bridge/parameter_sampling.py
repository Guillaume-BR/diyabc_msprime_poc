"""
Tirage des valeurs numériques des priors, avec retirage si les contraintes
d'ordre (ex: "t4>t3") ne sont pas respectées -- équivalent du mécanisme
"DRAW UNTIL" observé dans header.txt.

Limité pour l'instant à la loi UN (uniforme). Les autres lois (LU, GA --
voir doc DIYABC sur "mean and standard deviation") ne sont pas encore
vérifiées : human n'utilise que UN, donc pas bloquant pour ce POC, mais
à traiter explicitement avant de généraliser à d'autres datasets.
"""

import random

from bridge.scenario_types import Prior, OrderConstraint


class ConstraintsNotSatisfiedError(Exception):
    """Levée si aucun tirage valide n'a été trouvé en max_attempts essais."""


def _draw_one_value(prior: Prior, rng: random.Random) -> float:
    if prior.law == "UN":
        min_, max_, *_ = prior.bounds
        return rng.uniform(min_, max_)

    raise NotImplementedError(
        f"Loi '{prior.law}' non implémentée (seule 'UN' est supportée pour "
        f"l'instant). Prior concerné : {prior.name!r}"
    )


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