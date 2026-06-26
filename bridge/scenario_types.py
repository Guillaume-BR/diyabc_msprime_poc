from dataclasses import dataclass, field

@dataclass
class MergeEvent:
    time_expr: str
    ancestral_pop: int      # survit (le "a" de merge a b)
    derived_pop: int        # disparaît (le "b" de merge a b)

@dataclass
class VarNeEvent:
    time_expr: str
    pop: int
    new_size_expr: str      # "Nbn3", "N34" — nom de paramètre ou nombre

@dataclass
class SampleEvent:
    time_expr: str          # toujours "0" dans ce qu'on a vu
    pop: int

@dataclass
class Event:
    time_expr: str       # "t2-d3" ou "t1" — texte brut, pas encore évalué
    action: str          # "merge", "varNe", "sample"
    args: list[str]      # ["2", "1"] pour merge, ["3", "Nbn3"] pour varNe

@dataclass
class Scenario:                                    
    index: int
    weight: float
    initial_pop_size_exprs: list[str]
    events: list[Event] = field(default_factory=list)

@dataclass
class Prior:
    """Une ligne de la section 'historical parameters priors'.
    Ex: 'N1 N UN[1000.0,100000.0,0.0,0.0]' ->
        Prior(name="N1", category="N", law="UN", bounds=(1000.0, 100000.0, 0.0, 0.0))
    """

    name: str
    category: str        # "N" (taille), "T" (temps), "A" (taux admixture)
    law: str              # "UN", "LU", "GA"...
    bounds: tuple[float, ...]

@dataclass
class OrderConstraint:
    """Une ligne du type 't4>t3' ou 't431<t32' dans
    'historical parameters priors' : contrainte que doit respecter un
    tirage de valeurs pour être valide.
 
    Observée juste avant 'DRAW UNTIL' dans header.txt -- DIYABC tire les
    valeurs de tous les priors, vérifie ces contraintes, et retire tant
    qu'elles ne sont pas respectées (d'où le nom 'DRAW UNTIL').
 
    Vérifié dans history.cpp (EventC, parsing de 'condition') : 4 opérateurs
    possibles -- '>', '<', '>=', '<='. La contrainte signifie
    "valeur(param1) <operator> valeur(param2)" doit être vraie.
    """
    param1: str      # "t4"
    operator: str    # ">" , "<", ">=", ou "<="
    param2: str      # "t3"
 
    def is_satisfied(self, values: dict[str, float]) -> bool:
        v1, v2 = values[self.param1], values[self.param2]
        if self.operator == ">":
            return v1 > v2
        if self.operator == "<":
            return v1 < v2
        if self.operator == ">=":
            return v1 >= v2
        if self.operator == "<=":
            return v1 <= v2
        raise ValueError(f"Opérateur de contrainte inconnu : {self.operator!r}")

