"""
Construit une msprime.Demography à partir d'un Scenario (header_dataclasses.py)
et d'un dict de valeurs numériques tirées (parameter_sampling.py).

evaluate_expression() est l'équivalent Python de ParticleC::getvalue()
(particuleC.cpp) : transforme une expression texte ("t2-d3", "t1", "0")
en valeur numérique, en utilisant les valeurs déjà tirées des priors.

Événements gérés : SampleEvent, MergeEvent, VarNeEvent, SplitEvent
(admixture) -- couvre les 6 scénarios de human.
"""

import re

import msprime

from bridge.header_dataclasses import (
    MergeEvent,
    SampleEvent,
    Scenario,
    SplitEvent,
    VarNeEvent,
)

# "t2-d3" -> param1="t2", op="-", param2="d3"
# "t1"    -> pas de match -> traité comme un nom de paramètre seul
_EXPR_RE = re.compile(r"^(\w+)([+-])(\w+)$")


def evaluate_expression(expr: str, values: dict[str, float]) -> float:
    """Évalue une expression de temps ou de taille telle qu'elle apparaît
    dans header.txt : un nombre littéral ("0"), un nom de paramètre tiré
    ("t1"), ou une somme/différence de deux noms ("t2-d3", "t2+d3").

    Équivalent de ParticleC::getvalue() en C++.
    """
    match = _EXPR_RE.match(expr)
    if match:
        name1, op, name2 = match.groups()
        v1 = evaluate_expression(name1, values)
        v2 = evaluate_expression(name2, values)
        return v1 + v2 if op == "+" else v1 - v2

    if expr in values:
        return values[expr]

    try:
        return float(expr)
    except ValueError:
        raise ValueError(
            f"Expression '{expr}' non évaluable : ni nombre littéral, "
            f"ni nom de paramètre connu parmi {sorted(values)}"
        ) from None


def build_demography(
    scenario: Scenario, values: dict[str, float]
) -> msprime.Demography:
    """Construit la Demography msprime correspondant au scenario, avec les
    valeurs numériques déjà tirées dans `values`.

    Les populations sont nommées "pop1", "pop2", ... d'après leur indice
    dans header.txt (1-indexed, comme dans le fichier).
    """
    # n_pops = len(scenario.initial_pop_size_exprs)
    demography = msprime.Demography()

    for i, size_expr in enumerate(scenario.initial_pop_size_exprs, start=1):
        demography.add_population(
            name=f"pop{i}",
            initial_size=evaluate_expression(size_expr, values),
            initially_active=True,  # sinon la population ancestrale d'un merge est échantillonnée au temps du merge, pas au présent
        )

    # Les événements sont listés dans header.txt du présent vers le passé
    # (time croissant) -- msprime attend le même ordre pour
    # add_population_split / add_population_parameters_change.
    for event in scenario.events:
        time = evaluate_expression(event.time_expr, values)

        if isinstance(event, SampleEvent):
            # L'échantillonnage est géré séparément, au moment de la
            # simulation (msprime.sim_ancestry(samples=...)), pas ici.
            continue

        if isinstance(event, MergeEvent):
            demography.add_population_split(
                time=time,
                derived=[f"pop{event.derived_pop}"],
                ancestral=f"pop{event.ancestral_pop}",
            )
            continue

        if isinstance(event, VarNeEvent):
            demography.add_population_parameters_change(
                time=time,
                population=f"pop{event.pop}",
                initial_size=evaluate_expression(event.new_size_expr, values),
            )
            continue

        if isinstance(event, SplitEvent):
            # pop_derived disparaît : chaque lignée part vers ancestral_pop1
            # avec probabilité admixture_rate, sinon vers ancestral_pop2 --
            # voir SplitEvent (header_dataclasses.py) pour la sémantique exacte,
            # vérifiée contre particuleC.cpp::ParticleC::split_pop.
            rate = evaluate_expression(event.admixture_rate, values)
            demography.add_admixture(
                time=time,
                derived=f"pop{event.derived_pop}",
                ancestral=[
                    f"pop{event.ancestral_pop1}",
                    f"pop{event.ancestral_pop2}",
                ],
                proportions=[rate, 1 - rate],
            )
            continue

        raise NotImplementedError(
            f"Type d'événement non géré par build_demography : {event!r}"
        )

    demography.sort_events()  ## msprime exige que les événements soient triés par temps croissant
    return demography


def rescale_demography(
    demography: msprime.Demography, factor: float
) -> msprime.Demography:
    """Retourne une COPIE de demography où chaque taille de population est
    multipliée par factor -- ne modifie jamais l'original (une même
    Demography <A> est réutilisée telle quelle pour bâtir la version
    rescalée <X>/<Y>/<M> du même scénario/tirage de valeurs, donc muter en
    place casserait les autres types de locus qui doivent encore utiliser
    la version non rescalée).

    factor : à appeler avec coalescence_coefficient(locus_type,
    sex_ratio) / 2 (voir observed_data.py -- le /2 vient de la conversion
    coeffcoal*N/2 = nombre effectif de copies de gène). Le
    facteur est calculé par l'appelant.

    Deux choses à rescaler :
    - demography.populations : chaque Population a un .initial_size
      (toujours défini, jamais None).
    - demography.events : SEULS les événements PopulationParametersChange
      (générés par add_population_parameters_change, donc par les
      VarNeEvent) ont un .initial_size -- et il peut valoir None (un
      PopulationParametersChange peut ne changer QUE growth_rate, pas la
      taille -- pas notre cas ici, mais à ne pas casser si ça arrive).
      Les autres types d'événements (PopulationSplit, Admixture...) n'ont
      pas d'attribut .initial_size du tout.
    """

    demography_copy = demography.copy()

    for pop in demography_copy.populations:
        pop.initial_size *= factor

    for event in demography_copy.events:
        if hasattr(event, "initial_size") and event.initial_size is not None:
            event.initial_size *= factor

    return demography_copy


def extract_referenced_names(expr: str) -> set[str]:
    """Extrait le ou les noms de paramètres référencés par une expression
    de header.txt, SANS l'évaluer numériquement -- "t2-d3" -> {"t2","d3"},
    "t1" -> {"t1"}, "0" -> set() (un nombre littéral ne référence aucun
    paramètre).

    Utilisé pour déterminer quels paramètres un scénario utilise
    réellement (nécessaire pour filtrer les colonnes du reftable.bin par
    scénario -- voir reftable_loop.write_reftable_bin et notes/
    exploration.md, bug "21 vs 16 paramètres pour le scénario 1").
    """
    match = _EXPR_RE.match(expr)
    if match:
        name1, _, name2 = match.groups()
        return extract_referenced_names(name1) | extract_referenced_names(name2)

    try:
        float(expr)
        return set()  # nombre littéral, pas un nom de paramètre
    except ValueError:
        return {expr}


def get_parameter_names_used_by_scenario(scenario: Scenario) -> set[str]:
    """Collecte l'ensemble des noms de paramètres réellement référencés
    par un scénario : tailles de population initiales, et time_expr /
    new_size_expr de chacun de ses événements.

    C'est ce sous-ensemble (pas la totalité des priors déclarés dans
    header.txt) qui doit constituer les colonnes param[] du reftable.bin
    pour ce scénario -- un scénario donné n'utilise généralement qu'une
    partie des priors globaux (ex: human/header.txt a 21 priors déclarés,
    mais le scénario 1 n'en référence que 16 -- ra/t11/t22/t33/t44
    appartiennent aux scénarios 2-6, pas au scénario 1).
    """
    names: set[str] = set()

    for size_expr in scenario.initial_pop_size_exprs:
        names |= extract_referenced_names(size_expr)

    for event in scenario.events:
        names |= extract_referenced_names(event.time_expr)
        if isinstance(event, VarNeEvent):
            names |= extract_referenced_names(event.new_size_expr)
        if isinstance(event, SplitEvent):
            names |= extract_referenced_names(event.admixture_rate)

    return names
