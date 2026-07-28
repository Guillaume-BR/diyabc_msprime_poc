from dataclasses import dataclass, field


@dataclass
class MergeEvent:
    time_expr: str
    ancestral_pop: int  # survit (le "a" de merge a b)
    derived_pop: int  # disparaît (le "b" de merge a b)


@dataclass
class SplitEvent:
    time_expr: str
    derived_pop: int  # disparaît
    ancestral_pop1: int  # survit avec le taux d'admixture
    ancestral_pop2: int  # survit avec le taux 1 - admixture
    admixture_rate: str


@dataclass
class VarNeEvent:
    time_expr: str
    pop: int
    new_size_expr: str  # "Nbn3", "N34" — nom de paramètre ou nombre


@dataclass
class SampleEvent:
    time_expr: str  # toujours "0" dans ce qu'on a vu
    pop: int


@dataclass
class Event:
    time_expr: str  # "t2-d3" ou "t1" — texte brut, pas encore évalué
    action: str  # "merge", "varNe", "sample"
    args: list[str]  # ["2", "1"] pour merge, ["3", "Nbn3"] pour varNe


@dataclass
class Scenario:
    index: int
    weight: float
    initial_pop_size_exprs: list[
        str
    ]  # ["N1", "N2", "N3"] — texte brut, pas encore évalué
    events: list[Event] = field(default_factory=list)


@dataclass
class Prior:
    """Une ligne de la section 'historical parameters priors'.
    Ex: 'N1 N UN[1000.0,100000.0,0.0,0.0]' ->
        Prior(name="N1", category="N", law="UN", bounds=(1000.0, 100000.0, 0.0, 0.0))
    """

    name: str
    category: str  # "N" (taille), "T" (temps), "A" (taux admixture)
    law: str  # "UN", "LU", "GA"...
    bounds: tuple[float, ...]


@dataclass
class GroupPrior:
    """Une ligne de la section 'group priors' du header.txt.
    Ex : group G1 [M]
         MEANMU UN[1e-4,1e-3,5e-4,2]
    -> GroupPrior(
        group="G1",
        ms_or_seq="M",
        name="MEANMU",
        law="UN",
        min=1e-4,
        max=1e-3,
        mean=5e-4,
        sdshape=2,
        model=False,
        name_model = None,
        model_bounds=None)
    """

    group: str  # "G1"
    ms_or_seq: str  # "M" ou "S" ou None si pas précisé
    name: str  # "MEANMU"
    law: str | None  # "UN", "LU", "GA" ou None si pas précisé ou si c'est un model
    min: float | None  # 1e-4 ou None si pas précisé ou si c'est un model
    max: float | None  # 1e-3 ou None si pas précisé
    mean: float | None  # 5e-4 ou None si pas précisé
    sdshape: float | None  # 2 ou None si pas précisé
    model: (
        bool | None
    )  # True si c'est un model, False si c'est une loi, None si pas précisé
    name_model: str | None  # "K2P" ou None si pas précisé ou si c'est une loi
    model_bounds: (
        tuple[float, ...] | None
    )  # (0.0, 0.5, 0.0, 0.0) ou None si c'est une loi


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

    param1: str  # "t4"
    operator: str  # ">" , "<", ">=", ou "<="
    param2: str  # "t3"

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


@dataclass
class LociDescription:
    total_loci: dict[str, int]  # {"A": 10, "X": 5, "Y": 2, "M": 1}
    group: str  # "G1"
    start_index: int  # 0 (déjà converti en 0-based, comme le C++ : prem = N - 1)


@dataclass
class LociDescriptionDetailed:
    """Une ligne de la section 'loci description' du header.txt, au format détaillé
    observé dans sequences-mut ("Locus_M_A_1_ <A> [M] G1 2 40"), pas au format condensé
    observé dans sequences-mut ("5000 <A> G1 from 1"), qui est un format différent
    (header.cpp distingue ces deux cas selon dataobs.filetype).
    """

    locus_name: str  # "Locus_M_A_1_"
    heritage: str  # "A"
    ms_or_seq: str  # "M" ou "S"
    group: str  # "G1"
    motif_size: int | None  # 2
    motif_range: int | None  # 40
    dnalength: int | None  # 100 (pour séquences, sinon None)
