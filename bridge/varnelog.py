"""
Sortie de bottleneck logistique (Verhulst) dans msprime, alignée sur le patch VARNELOG de DIYABC.

Modèle en générations discrètes (Wright-Fisher) :
    N(k) = NS / (1 + A exp(-r k)),   A = (NS - NF) / NF,   k = 0, 1, 2, ... générations depuis la fondation,
    N(0) = NF exactement (la génération fondatrice).
En temps backward t (0 = échantillonnage, tF = fondation), l'intervalle [t, t+1] a la taille N(tF - t - 1).

mode="harm"  : tailles constantes par morceaux = moyennes harmoniques discrètes exactes,
               MÊME découpage que DIYABC (validation diyabc-MsP vs diyabc-JMC, marqueur par marqueur) ;
mode="exact" : une époque par génération (N(k) exact) jusqu'à l'horizon où N est à 0,1 % de NS.
"""

import math

import msprime
import numpy as np

from bridge.configuration import VARNELOG_K


def verhulst_N(
    s: float | np.ndarray, NF: float, NS: float, r: float
) -> float | np.ndarray:
    """
    Retourne la valeur de la logistique de Verhulst en fonction de différentes paramètres.

    Args:
        s : Temps en générations
        NF : Population des fondateurs
        NS : Population stable
        r : Taux de croissance

    Returns :
        La valeur de la logistique de Verhulst
    """
    return NS / (1.0 + (NS - NF) / NF * np.exp(-r * s))


def _compute_horizon_Dstar(NF: float, NS: float, r: float):
    """
    Calcule l'horizon

    Args:
        NF: population des fondateurs
        NS: population stable
        r: taux de croissance

    Returns:
        Horizon où N est à 0.1% de NS
    """
    A = (NS - NF) / NF
    return max(1.0, math.log(1000.0 * abs(A)) / r) if abs(A) > 1e-12 else 1.0


def _check_varnelog_parameters(NF: float, NS: float, r: float) -> None:
    """
    Vérifie certaines conditions NF >=1, NS >=1 et r>0 sinon lève une erreur.

    Args:
        NF: Taille de la ppulation fondatrice
        NS: Taille de la populaiton stable
        r: Taux de croissande la logistique de Verhulst

    Returns:
        None

    Raises:
        Si une des conditions n'est pas remplie.
    """

    if r <= 0:
        raise ValueError(
            f"Le taux de croissance r={r!r} doit être strictement positif."
        )
    if NF < 1 or NS < 1:
        raise ValueError(
            "La taille des populations fondatrice NF et stable NS doit être supérieure ou égale à 1."
        )


def varnelog_segments(
    t0: float, t1: float, tF: float, NF: int, NS: int, r: float, K: int = VARNELOG_K
) -> list[tuple[float, float, int]]:
    """
    Pour un segment backward [t0, t1] renvoie une liste de (debut, fin, N_harmonique) en temps backward croissant.

    Args:
        t0: Borne inférieure du segment de temps
        t1: Borne supérieure du segment de temps
        tF: Temps de la fondation
        NF: Population des fondateurs
        NS: Population stable
        r: Taux de croissance
        K: Nombre de segments

    Returns:
        Liste de tuples contenant les bornes en temps backward des différentes coupures et la valeur de la logistique sur cet intervalle


    """
    _check_varnelog_parameters(NF, NS, r)
    if t1 > tF:
        raise ValueError(
            f"Le segment [t0={t0!r}, t1={t1!r}] dépasse la fondation tF={tF!r} : il faut t1 <= tF."
        )
    # Constante pour que la logistique de Verhulst aille de NF à NS
    A = (NS - NF) / NF
    # Pour revenir à un temps forward et pas backward
    sa, sb = tF - t1, tF - t0
    # Horizon où N est à 0.1% de NS
    Dstar = _compute_horizon_Dstar(NF, NS, r)
    breaks = {sa, sb}
    if (
        abs(A) > 1e-12
    ):  # coupures à pas égaux en ln N pour avoir le même rapport de N entre chaque niveau
        NDs = NS / (1 + A * math.exp(-r * Dstar))
        # Après K-1, on ne découpe plus car N=NS
        for k in range(1, K):
            Nk = math.exp(
                math.log(NF) + k * (math.log(NDs) - math.log(NF)) / K
            )  # valeur des coupures
            x = (NS / Nk - 1) / A
            if x <= 0:
                continue
            s = int(math.floor(-math.log(x) / r + 0.5))  # arrondi du banquier
            if sa < s < sb:
                breaks.add(
                    s
                )  # on garde les temps se trouvant bien dans le bon intervalle de temps forward
    sD = int(math.ceil(Dstar))
    if sa < sD < sb:
        breaks.add(sD)
    if sa < 1 < sb:  # génération fondatrice isolée (taille NF exactement)
        breaks.add(1)
    sorted_breaks = sorted(breaks)
    out = []
    for k in range(len(breaks) - 1, 0, -1):
        s1, s2 = sorted_breaks[k - 1], sorted_breaks[k]
        if abs(A) <= 1e-12:
            invN = 1.0 / NS
        else:
            # moyenne harmonique DISCRÈTE sur les générations k = s1 .. s2-1
            invN = (
                1
                + A
                * math.exp(-r * s1)
                * (1 - math.exp(-r * (s2 - s1)))
                / ((s2 - s1) * (1 - math.exp(-r)))
            ) / NS
        out.append((tF - s2, tF - s1, max(1, int(math.floor(1 / invN + 0.5)))))
    return out


def add_varnelog(
    demography: msprime.Demography,
    pop: str,
    tL: float,
    tF: float,
    NF: float,
    NS: float,
    r: float,
    K: int = VARNELOG_K,
    mode: str = "harmonic",
) -> None:
    """
    Ajoute la phase logistique de `pop` entre tL et tF (backward). La population doit être créée
    avec initial_size=NS, et sa fondation (add_population_split / mass migration) placée à tF.

    Args:
        demography: demography msprime
        pop: population avec le nom msprime ("pop0")
        tL: Date de début de la phase logistique en temps backward
        tF: Temps des fondateurs
        NF: Population des fondateurs
        NS: Population stable
        r: Taux de croissance de la logistique
        K: Nombre de segments considérés pour la discrétisation
        mode: "harmonic" pour la discrétisation avec moyenne harmonique, "exact" : calcul de Nk directement avec la logistique de Verhulst

    Returns:
        None

    Raises:
        tL doit être strictement inférieur à TF.
        Le mode doit être "harmonic" ou "exact".

    """
    if not (tL < tF):
        raise ValueError("VARNELOG : il faut tL < tF")
    if mode == "harmonic":
        for a, _, Nh in varnelog_segments(tL, tF, tF, NF, NS, r, K):
            # On modifie la population en la laissant stable sur l'intervalle de temps.
            demography.add_population_parameters_change(
                time=a, population=pop, initial_size=Nh, growth_rate=0
            )
    elif mode == "exact":
        _check_varnelog_parameters(NF, NS, r)
        Dh = _compute_horizon_Dstar(NF, NS, r)
        kmax = min(int(tF - tL), int(math.ceil(Dh)))
        for k in range(kmax - 1, -1, -1):  # backward croissant
            Nk = verhulst_N(k, NF, NS, r)
            demography.add_population_parameters_change(
                time=tF - k - 1,
                population=pop,
                initial_size=max(1.0, Nk),
                growth_rate=0,
            )
        if tF - kmax > tL:  # au-delà de l'horizon : N(kmax) ~ NS
            demography.add_population_parameters_change(
                time=tL, population=pop, initial_size=NS, growth_rate=0
            )
    else:
        raise ValueError("mode = 'harmonic' ou 'exact'")
