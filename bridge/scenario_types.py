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
class Scenario:                                    # ← elle est bien là, juste après
    index: int
    weight: float
    initial_pop_size_exprs: list[str]
    events: list[Event] = field(default_factory=list)

@dataclass
class Prior:
    name: str
    category: str        # "N" (taille), "T" (temps), "A" (taux admixture)
    law: str              # "UN", "LU", "GA"
    bounds: tuple[float, ...]