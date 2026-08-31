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
    """Une ligne de la section 'historical parameters priors' de header.txt.

    Ex: 'N1 N UN[1000.0,100000.0,0.0,0.0]' devient
    Prior(name="N1", category="N", law="UN", min=1000.0, max=100000.0,
    mean=0.0, sdshape=0.0).

    Attributes:
        name: Nom du paramètre (ex: "N1").
        category: "N" (taille), "T" (temps) ou "A" (taux admixture).
        law: Loi de tirage ("UN", "LU", "NO", "LN" ou "GA").
        min: Borne basse.
        max: Borne haute.
        mean: Moyenne, utilisée seulement par "NO"/"LN"/"GA" (0.0 et non
            lu sinon).
        sdshape: Écart-type ("NO") ou paramètre de forme de la
            reparamétrisation gamma ("LN"/"GA"), utilisé seulement par
            "NO"/"LN"/"GA" (0.0 et non lu sinon).
    """

    name: str
    category: str
    law: str
    min: float
    max: float
    mean: float
    sdshape: float


@dataclass
class GroupPrior:
    """Une ligne de la section 'group priors' de header.txt (MicroSat/sequences-mut).

    Ex: 'MEANMU UN[1e-4,1e-3,5e-4,2]' sous 'group G1 [M]' devient
    GroupPrior(group="G1", ms_or_seq="M", name="MEANMU", law="UN",
    min=1e-4, max=1e-3, mean=5e-4, sdshape=2, model=False).

    Attributes:
        group: Nom du groupe (ex: "G1").
        ms_or_seq: "M" (MicroSat) ou "S" (séquence), None si pas précisé.
        name: Nom du prior (ex: "MEANMU").
        law: Loi de tirage ("UN", "LU", "GA"...), None si c'est un model.
        min: Borne basse, None si pas précisé ou si c'est un model.
        max: Borne haute, None si pas précisé.
        mean: Moyenne, None si pas précisé.
        sdshape: Écart-type ou paramètre de forme, None si pas précisé.
        model: True si c'est un model, False si c'est une loi, None si pas précisé.
        name_model: Nom du modèle (ex: "K2P"), None si c'est une loi.
        p_fixe: Pourcentage (0-100, pas une fraction 0-1) de sites
            invariants, None si pas précisé.
        gams: Paramètre de forme gamma (hétérogénéité de taux), None si pas précisé.
    """

    group: str
    ms_or_seq: str
    name: str
    law: str | None
    min: float | None
    max: float | None
    mean: float | None
    sdshape: float | None
    model: bool | None
    name_model: str | None
    p_fixe: float | None
    gams: float | None


@dataclass
class OrderConstraint:
    """Contrainte d'ordre entre deux priors, ex: 't4>t3' ou 't431<t32'.

    Trouvée juste avant 'DRAW UNTIL' dans header.txt : DIYABC tire les
    valeurs de tous les priors et retire tant que ces contraintes ne sont
    pas respectées. Vérifié dans history.cpp (EventC, parsing de
    'condition') : 4 opérateurs possibles, '>', '<', '>=', '<='.

    Attributes:
        param1: Nom du premier paramètre (ex: "t4").
        operator: ">", "<", ">=" ou "<=".
        param2: Nom du second paramètre (ex: "t3").
    """

    param1: str
    operator: str
    param2: str

    def is_satisfied(self, values: dict[str, float]) -> bool:
        """Vérifie si un tirage de valeurs respecte la contrainte.

        Args:
            values: Valeurs tirées, indexées par nom de paramètre.

        Returns:
            True si "valeur(param1) <operator> valeur(param2)" est vrai.

        Raises:
            ValueError: Si l'opérateur n'est pas ">", "<", ">=" ou "<=".
        """
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
    loci_counts_by_heritage: dict[str, int]  # {"A": 10, "X": 5, "Y": 2, "M": 1}
    group: str  # "G1"
    start_index: int  # 0 (déjà converti en 0-based, comme le C++ : prem = N - 1)


@dataclass
class LociDescriptionDetailed:
    """Une ligne du format détaillé (un locus par ligne) de 'loci description'.

    Ex: "Locus_M_A_1_ <A> [M] G1 2 40". Utilisé par MicroSat/sequences-mut,
    à distinguer du format condensé ("5000 <A> G1 from 1") -- header.cpp
    distingue ces deux cas selon dataobs.filetype.

    Attributes:
        name: Nom du locus (ex: "Locus_M_A_1_").
        heritage: Type d'hérédité ("A", "H", "X", "Y" ou "M").
        ms_or_seq: "M" (MicroSat) ou "S" (séquence).
        group: Nom du groupe (ex: "G1").
        motif_size: Taille du motif (MicroSat uniquement), None sinon.
        motif_range: Étendue du motif (MicroSat uniquement), None sinon.
        dnalength: Longueur de la séquence (séquences uniquement), None sinon.
    """

    name: str
    heritage: str
    ms_or_seq: str
    group: str
    motif_size: int | None
    motif_range: int | None
    dnalength: int | None
