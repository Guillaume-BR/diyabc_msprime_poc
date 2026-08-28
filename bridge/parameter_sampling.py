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

from bridge.header_dataclasses import (
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
    """Tire un scénario parmi `scenarios`, pondéré par son `weight`.

    Le poids est le "prior_proba" de
    particuleC.cpp::ParticleC::drawscenario. Reproduit exactement
    l'algorithme C++ : tirage d'un nombre uniforme `ra` dans [0,1),
    puis balayage de la somme cumulée des poids jusqu'à ce qu'elle
    atteigne ou dépasse `ra` -- même logique d'inversion de CDF que
    _draw_single_mutation_edge_child (ancestry_simulation.py),
    vérifiée boundary-compatible avec la boucle C++ ("while ra > sp").

    Ne normalise PAS les poids (comme le C++, qui ne le fait pas non
    plus) : si leur somme est < 1, le DERNIER scénario de la liste sert
    de secours pour tout `ra` au-delà de la somme cumulée -- même
    comportement que la boucle C++, bornée à nscenarios-1.

    Args:
        scenarios: Les scénarios candidats.
        seed: La graine du tirage.

    Returns:
        Le scénario tiré.
    """
    rng = random.Random(seed)
    ra = rng.random()

    cumulative = 0.0
    for scenario in scenarios:
        cumulative += scenario.weight
        if ra <= cumulative:
            return scenario

    return scenarios[-1]  # dernier scénario de secours si somme des poids < ra


def _draw_one_value(prior: Prior, rng: random.Random) -> float:
    """Tire une valeur pour un prior donné, selon sa loi et ses bornes.
    Pour les lois normale, log-normale et gamma, on retire les valeurs hors bornes comme DIYABC.

    Args:
        prior: Le prior à tirer.
        rng: Le générateur aléatoire à utiliser.

    Returns:
        La valeur tirée.

    Raises:
        NotImplementedError: Si la loi n'est pas encore implémentée.
    """
    if prior.min == prior.max:
        # bornes identiques : pas de tirage, valeur fixée
        value = float(prior.min)
    elif prior.law == "UN":  # uniform
        value = rng.uniform(prior.min, prior.max)
    elif prior.law == "LU":  # log-uniform
        value = math.exp(rng.uniform(math.log(prior.min), math.log(prior.max)))
    elif prior.law == "NO":  # normal
        min_, max_, mean, sdshape = prior.min, prior.max, prior.mean, prior.sdshape
        while True:
            value = rng.gauss(mean, sdshape)
            if min_ <= value <= max_:
                break
    elif prior.law == "LN":  # log-gaussian
        min_, max_, mean, sdshape = prior.min, prior.max, prior.mean, prior.sdshape
        while True:
            value = math.exp(rng.gauss(math.log(mean), math.log(sdshape)))
            if min_ <= value <= max_:
                break
    elif prior.law == "GA":  # gamma
        min_, max_, mean, sdshape = prior.min, prior.max, prior.mean, prior.sdshape
        if mean < 1e-12:
            value = 0.0
        elif sdshape < 1e-12:
            value = mean
        elif max_ < 1e-12:
            value = max_
        else:
            while True:
                value = rng.gammavariate(sdshape, mean / sdshape)
                if min_ <= value <= max_:  #
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

    Quelques gardes-fous pour éviter des erreurs de tirage si le group
    prior est mal défini.

    Args:
        group_prior: Le group prior à tirer.
        rng: Le générateur aléatoire à utiliser.

    Returns:
        La valeur tirée.

    Raises:
        ValueError: Si le group prior n'a pas de loi ou de bornes
            associées, ou si la loi est GA mais que les priors MEAN et
            SDSHAPE correspondants n'ont pas été tirés avant.
        NotImplementedError: Si la loi n'est pas encore implémentée
            (propagée depuis _draw_one_value).
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
    prior = Prior(
        name=group_prior.name,
        category="G",  # catégorie fictive pour les group priors
        law=group_prior.law,
        min=group_prior.min,
        max=group_prior.max,
        mean=group_prior.mean,
        sdshape=group_prior.sdshape,
    )
    return _draw_one_value(prior, rng)


def draw_parameter_values(
    priors: list[Prior],
    constraints: list[OrderConstraint],
    seed: int,
    max_attempts: int = 1000,
) -> dict[str, float]:
    """Tire une valeur pour chaque prior, en retirant tant que les
    contraintes d'ordre ne sont pas toutes satisfaites. On reproduit le comportement de DIYABC,
    en gardant le premier tirage satisfaisant toutes les contraintes.

    Args:
        priors: Les priors à tirer.
        constraints: Les contraintes d'ordre à satisfaire (ex: "t4>t3").
        seed: La graine du tirage.
        max_attempts: Le nombre maximal d'essais avant d'abandonner.

    Returns:
        Un dict {nom_prior: valeur}.

    Raises:
        ConstraintsNotSatisfiedError: Si aucun tirage valide n'est
            trouvé en max_attempts essais -- signe probable d'une
            configuration de contraintes incohérente (bornes de
            priors incompatibles avec les contraintes demandées)
            plutôt que d'une simple mauvaise chance.
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
    """Tire une valeur pour chaque group prior ou bien des valeurs pour le modèle.

    `seed` est décalé de _GROUP_PRIOR_SEED_OFFSET avant utilisation -- ne
    corrèle jamais ce tirage avec celui de draw_parameter_values, même si
    l'appelant leur passe la même seed de base.

    Point d'attention : ce tirage positionnel entre priors dépendants
    (MEANMU avant GAMMU, etc.) reproduit le comportement de DIYABC, qui lit
    ces lignes dans un ordre fixe sans jamais regarder leur nom (voir
    header.cpp::readHeadersimGroupPrior).

    Args:
        group_priors: Dict {nom_groupe: [GroupPrior, ...]}.
        seed: La graine de base du tirage (décalée en interne).

    Returns:
        Un dict {nom_groupe: {nom_prior: valeur}}.
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


def sampling_group_local_param(
    group_prior: GroupPrior,
    k_moy: float,
    n_loci: int,
    check_nloc: bool,
    list_loci: list[LociDescriptionDetailed],
    rng: random.Random,
) -> dict[str, float]:
    """Échantillonne le tirage par-locus (second niveau) d'un paramètre de groupe.

    S'applique aussi bien à kappa1/kappa2 (`build_transition_matrix`)
    qu'à mus_rate -- rien de spécifique à kappa dans l'implémentation.

    Args:
        group_prior: Le GroupPrior déjà résolu pour ce groupe (via
            draw_group_parameter_values), dont le `mean` est remplacé
            par `k_moy` avant tirage.
        k_moy: La valeur moyenne du groupe (premier niveau), déjà
            tirée par draw_group_parameter_values.
        n_loci: Le nombre de loci du groupe.
        check_nloc: Si True, un tirage indépendant par locus n'a lieu
            que si `n_loci > 1` en plus de `sdshape > 0.001` (cas de
            kappa1/mus_rate). Si False, seule la condition sur
            `sdshape` est vérifiée (cas de kappa2 -- asymétrie propre à
            DIYABC, reproduite telle quelle).
        list_loci: Les loci du groupe.
        rng: Le générateur aléatoire à utiliser.

    Returns:
        Un dict {nom_locus: valeur} -- soit un tirage indépendant par
        locus, soit `k_moy` répété pour chaque locus.
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


def sample_site_rates(
    p_fixe: float, gams: float, dnalength: int, rng: random.Random
) -> list[float]:
    """Tire mutsit : le taux de mutation relatif par site pour un locus séquence.

    Reproduit header.cpp:707-738 (y compris le "bug" sitefix -- les
    sites fixes sont toujours les premiers de la séquence, pas un
    sous-ensemble aléatoire).

    Args:
        p_fixe: Pourcentage (0-100, pas une fraction 0-1) de sites invariants (`GroupPrior.p_fixe`).
        gams: Forme gamma de l'hétérogénéité de taux par site
            (`GroupPrior.gams`) -- 0.0 est une valeur valide,
            équivalente à un taux uniforme (voir MwcGen::ggamma3).
        dnalength: Longueur du locus.
        rng: Le générateur aléatoire à utiliser.

    Returns:
        La liste `mutsit`, de longueur `dnalength`, normalisée à
        somme 1.
    """
    nb_sites_variables = math.floor(dnalength * (1 - 0.01 * p_fixe) + 0.5)
    nb_sites_fixes = dnalength - nb_sites_variables
    if gams < 1e-12:
        mutsit = [1.0] * dnalength
    else:
        mutsit = [rng.gammavariate(gams, 1.0 / gams) for _ in range(dnalength)]

    # header.cpp:727-738 -- "sitefix" tire des indices aléatoires distincts,
    # MAIS la ligne qui fixe le site utilise l'indice de boucle "i", pas
    # "sitefix[i]" -- donc en pratique, ce sont TOUJOURS les nb_sites_fixes
    # PREMIERS sites (0, 1, 2...) qui sont mis à taux nul, jamais un sous-
    # ensemble aléatoire. On reproduit ce comportement tel quel (fidélité à
    # DIYABC), pas ce que le code semblait vouloir faire.
    # On garde en réserve et on l'active si besoin.
    # sitefix = [0] * (dnalength - nb_sites_variables)
    # for i in range(len(sitefix)):
    #    if i == 0:
    #        sitefix[i]  = rng.randint(1, dnalength)
    #    else:
    #        nouveau = False
    #        while not nouveau:
    #            sitefix[i] = rng.randint(1, dnalength)
    #            nouveau = True
    #            j = 0
    #            while j < i and nouveau:
    #                if sitefix[i] == sitefix[j]:
    #                    nouveau = False
    #                j += 1
    #   mutsit[sitefix[i]] = 0.0
    for i in range(nb_sites_fixes):
        mutsit[i] = 0.0

    total = sum(mutsit)
    mutsit = [x / total for x in mutsit]

    return mutsit
