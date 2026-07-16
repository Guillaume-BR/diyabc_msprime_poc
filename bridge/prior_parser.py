import re

from bridge.scenario_types import OrderConstraint, Prior

_SECTION_START_RE = re.compile(r"^historical parameters priors\s*\(")
_SECTION_END_RE = re.compile(r"^DRAW UNTIL\s*$")

# "N1 N UN[1000.0,100000.0,0.0,0.0]" -> name=N1, category=N, law=UN, bounds_str=...
_PRIOR_LINE_RE = re.compile(r"^(\S+)\s+(\S+)\s+([A-Z]+)\[([^\]]+)\]\s*$")

# "t4>t3", "t431<t32", "t4>=t3", "t4<=t3" -> param1, operator, param2.
# L'ordre des alternatives (>=|<=|>|<) importe : il faut tester les
# opérateurs à deux caractères avant les opérateurs à un caractère, sinon
# "t4>=t3" serait mal découpé en operator=">" + param2="=t3".
_CONSTRAINT_LINE_RE = re.compile(r"^(\S+?)(>=|<=|>|<)(\S+)$")


def _extract_priors_section(header_text: str) -> list[str]:
    """Extrait la section 'historical parameters priors' du texte complet de header.txt,
    et retourne la liste des lignes de cette section (sans les lignes vides)."""
    lines = header_text.splitlines()

    start = next(
        (i for i, line in enumerate(lines) if _SECTION_START_RE.match(line.strip())),
        None,
    )
    if start is None:
        raise ValueError(
            "Impossible de trouver la section 'historical parameters priors' dans le texte fourni"
        )

    end = next(
        (
            i
            for i in range(start + 1, len(lines))
            if _SECTION_END_RE.match(lines[i].strip())
        ),
        None,
    )
    if end is None:
        raise ValueError(
            "Impossible de trouver la fin de la section 'historical parameters priors' dans le texte fourni"
        )

    return [line.strip() for line in lines[start + 1 : end] if line.strip()]


def parse_priors(header_text: str) -> tuple[list[Prior], list[OrderConstraint]]:
    """Extrait les priors et les contraintes d'ordre de header.txt.

    Retourne (priors, constraints). Une ligne qui ne correspond à aucun
    des deux formats connus lève une erreur explicite plutôt que d'être
    silencieusement ignorée.
    """
    priors: list[Prior] = []
    constraints: list[OrderConstraint] = []

    for line in _extract_priors_section(header_text):
        prior_match = _PRIOR_LINE_RE.match(line)
        if prior_match:
            name, category, law, bounds_str = prior_match.groups()
            bounds = tuple(float(b) for b in bounds_str.split(","))
            priors.append(Prior(name=name, category=category, law=law, bounds=bounds))
            continue

        constraint_match = _CONSTRAINT_LINE_RE.match(line)
        if constraint_match:
            param1, operator, param2 = constraint_match.groups()
            constraints.append(
                OrderConstraint(param1=param1, operator=operator, param2=param2)
            )
            continue

        raise ValueError(f"Ligne de la section priors non reconnue : {line!r}")

    return priors, constraints


def is_constant_prior(prior: Prior) -> bool:
    """Détecte si un prior est quasi-dégénéré (min ≈ max), donc en
    pratique une constante déguisée en prior -- DIYABC exclut ces
    paramètres des colonnes du reftable.bin (vérifié indépendamment dans
    readReftable.R et abcranger/readreftable.cpp, voir notes/
    exploration.md et docs/synthese_diyabc_msprime.docx section 5.2).

    Règle exacte (reproduite des deux sources ci-dessus) :
        si maxi != 0.0 : constant si (maxi-mini)/maxi <= 0.000001
        si maxi == 0.0 : jamais considéré comme constant par cette règle
                         (évite une division par zéro -- comportement de
                         readReftable.R, où le test est dans un bloc
                         "if (maxi != 0.0)").

    Ne gère que les priors avec au moins 2 bornes numériques (suffisant
    pour UN/LU/GA, les seules lois rencontrées jusqu'ici -- bounds[0] et
    bounds[1] sont systématiquement min et max dans header.txt).
    """
    mini, maxi = prior.bounds[0], prior.bounds[1]
    if maxi == 0.0:
        return False
    return (maxi - mini) / maxi <= 0.000001
