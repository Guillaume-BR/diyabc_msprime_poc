"""
Construit une msprime.Demography à partir d'un Scenario (scenario_types.py)
et d'un dict de valeurs numériques tirées (parameter_sampling.py).
 
evaluate_expression() est l'équivalent Python de ParticleC::getvalue()
(particuleC.cpp) : transforme une expression texte ("t2-d3", "t1", "0")
en valeur numérique, en utilisant les valeurs déjà tirées des priors.
"""

import re

import msprime

from bridge.scenario_types import Scenario, SampleEvent, MergeEvent, SplitEvent, VarNeEvent

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
        raise ValueError(f"Impossible d'évaluer l'expression {expr!r} : "
                         f"ce n'est ni un nombre, ni un nom de paramètre tiré, "
                         f"ni une somme/différence de deux noms de paramètres.")
    
def build_demography(scenario: Scenario, values: dict[str, float]) -> msprime.Demography:
    """Construit la Demography msprime correspondant au scenario, avec les
    valeurs numériques déjà tirées dans `values`.

    Les populations sont nommées "pop1", "pop2", ... d'après leur indice
    dans header.txt (1-indexed, comme dans le fichier).
    """
    n_pops = len(scenario.initial_pop_size_exprs)
    demography = msprime.Demography()
    for i, size_expr in enumerate(scenario.initial_pop_size_exprs, start=1):
        demography.add_population(name=f"pop{i}", initial_size=evaluate_expression(size_expr, values))
    
    # Les événements sont listés dans header.txt du présent vers le passé
    # (time croissant) -- msprime attend le même ordre pour
    # add_population_split / add_population_parameters_change.
    for event in scenario.events:
        time = evaluate_expression(event.time_expr, values)
        if isinstance(event, SampleEvent):
            # msprime ne gère pas explicitement les échantillonnages
            # (samples) dans la Demography : on les fournit à part
            # dans msprime.sim_ancestry() via l'argument samples.
            continue
        if isinstance(event, MergeEvent): #attention à la dénomination, split en msprime correspond à l'addmixture, merge correspond à l'add_population_split
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
            admixture_rate = evaluate_expression(event.admixture_rate, values)
            demography.add_admixture(
                time=time,
                derived=f"pop{event.derived_pop}",
                ancestral=[f"pop{event.ancestral_pop1}", f"pop{event.ancestral_pop2}"],
                proportions=[admixture_rate, 1 - admixture_rate],
            )
            continue
        
        raise NotImplementedError(f"Type d'événement non géré par build_demography : {event!r}")
    demography.sort_events()
    return demography




