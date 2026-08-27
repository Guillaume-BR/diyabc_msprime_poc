import re

from bridge.header_dataclasses import GroupPrior, OrderConstraint, Prior

_SECTION_START_RE = re.compile(r"^historical parameters priors\s*\(")
_SECTION_END_RE = re.compile(r"^DRAW UNTIL\s*$")

_SECTION_GROUP_START = re.compile(r"^group priors\s*\(")

# "N1 N UN[1000.0,100000.0,0.0,0.0]" -> name=N1, category=N, law=UN, bounds_str=...
_PRIOR_LINE_RE = re.compile(r"^(\S+)\s+(\S+)\s+([A-Z]+)\[([^\]]+)\]\s*$")

# MEANMU UN[1e-4,1e-3,5e-4,2]
_PRIOR_GROUP_LINE_RE = re.compile(r"^(\S+)\s+([A-Z]+)\[([^\]]+)\]\s*$")

# "t4>t3", "t431<t32", "t4>=t3", "t4<=t3" -> param1, operator, param2.
# L'ordre des alternatives (>=|<=|>|<) importe : il faut tester les
# opérateurs à deux caractères avant les opérateurs à un caractère, sinon
# "t4>=t3" serait mal découpé en operator=">" + param2="=t3".
_CONSTRAINT_LINE_RE = re.compile(r"^(\S+?)(>=|<=|>|<)(\S+)$")


def _extract_historical_priors_section(header_text: str) -> list[str]:
    """Extrait les lignes de la section 'historical parameters priors'.

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        Les lignes de la section, sans les lignes vides.

    Raises:
        ValueError: Si la section ou sa fin est introuvable.
    """
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
            if _SECTION_END_RE.match(lines[i].strip()) or lines[i] == ""
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

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        Un tuple (priors, constraints).

    Raises:
        ValueError: Si une ligne de la section ne correspond à aucun des
            deux formats connus (prior ou contrainte) -- pas d'ignorance
            silencieuse.
    """
    priors: list[Prior] = []
    constraints: list[OrderConstraint] = []

    for line in _extract_historical_priors_section(header_text):
        prior_match = _PRIOR_LINE_RE.match(line)
        if prior_match:
            name, category, law, bounds_str = prior_match.groups()
            min_, max_, mean, sdshape = (float(b) for b in bounds_str.split(","))
            priors.append(
                Prior(
                    name=name,
                    category=category,
                    law=law,
                    min=min_,
                    max=max_,
                    mean=mean,
                    sdshape=sdshape,
                )
            )
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

    Args:
        prior: Le prior à tester.

    Returns:
        True si le prior est quasi-constant selon la règle ci-dessus.
    """
    mini, maxi = prior.min, prior.max
    if maxi == 0.0:
        return False
    return (maxi - mini) / maxi <= 0.000001


def _extract_priors_group_section(header_text: str) -> list[str]:
    """Extrait les lignes de la section 'group priors'.

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        Les lignes de la section, sans les lignes vides.

    Raises:
        ValueError: Si la section ou sa fin est introuvable.
    """
    lines = header_text.splitlines()

    g_start = (
        i for i, line in enumerate(lines) if _SECTION_GROUP_START.match(line.strip())
    )
    start = next(g_start, None)
    if start is None:
        raise ValueError(
            "Impossible de trouver la section 'group priors' dans le texte fourni"
        )

    g_end = (i for i in range(start + 1, len(lines)) if lines[i] == "")
    end = next(g_end, None)
    if end is None:
        raise ValueError(
            "Impossible de trouver la fin de la section 'group priors' dans le texte fourni"
        )
    return [line.strip() for line in lines[start + 1 : end] if line.strip()]


def parse_group_priors(header_text: str) -> dict[str, list[GroupPrior]]:
    """Extrait les priors et models des différents group priors de header.txt.

    Une ligne est soit une loi de prior, soit un model : si elle ne
    correspond pas au format d'une loi, elle est supposée être un model,
    sans validation positive du format 'MODEL ...' -- une ligne réellement
    malformée serait donc mal interprétée silencieusement plutôt que de
    lever une erreur claire (gap connu, pas encore corrigé).

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        {nom_de_groupe: [GroupPrior, ...], ...}.

    Raises:
        ValueError: Si une ligne 'group ...' est malformée, ou si une
            ligne de prior/model apparaît avant toute ligne 'group'.
    """
    group_priors: dict[str, list[GroupPrior]] = {}

    for line in _extract_priors_group_section(header_text):
        # "group G1 [M]" -> group=G1, ms_or_seq=M
        if line.startswith("group"):
            parts = line.split()
            if (
                len(parts) != 3
                or not parts[2].startswith("[")
                or not parts[2].endswith("]")
            ):
                raise ValueError(
                    f"Ligne de la section group priors non reconnue : {line!r}"
                )
            group_name = parts[1]
            ms_or_seq = parts[2][1:-1]  # remove brackets
            group_priors[group_name] = []
            continue

        # "MEANMU UN[1e-4,1e-3,5e-4,2]" -> name=MEANMU, law=UN, bounds_str=...
        prior_match = _PRIOR_GROUP_LINE_RE.match(line)
        if prior_match:
            name, law, bounds_str = prior_match.groups()
            dict_bounds = {"min": None, "max": None, "mean": None, "sdshape": None}
            # On va tester chaque float(bounds_str.split(",")) pour voir si c'est un float ou pas,
            # si c'est pas un float, on renvoie None
            for i in range(len(bounds_str.split(","))):
                try:
                    dict_bounds[list(dict_bounds.keys())[i]] = float(
                        bounds_str.split(",")[i]
                    )
                except ValueError:
                    dict_bounds[list(dict_bounds.keys())[i]] = None
            if not group_priors:
                raise ValueError(
                    f"Ligne de la section group priors avant tout 'group' : {line!r}"
                )
            last_group_name = list(group_priors.keys())[-1]
            group_priors[last_group_name].append(
                GroupPrior(
                    group=last_group_name,
                    ms_or_seq=ms_or_seq,
                    name=name,
                    law=law,
                    min=dict_bounds["min"],
                    max=dict_bounds["max"],
                    mean=dict_bounds["mean"],
                    sdshape=dict_bounds["sdshape"],
                    model=False,
                    name_model=None,
                    p_fixe=None,
                    gams=None,
                )
            )
            continue
        else:
            _, name_model = line.split()[0:2]
            p_fixe, gams = float(line.split()[2]), float(line.split()[3])
            group_priors[last_group_name].append(
                GroupPrior(
                    group=last_group_name,
                    ms_or_seq=ms_or_seq,
                    name=None,
                    law=None,
                    min=None,
                    max=None,
                    mean=None,
                    sdshape=None,
                    model=True,
                    name_model=name_model,
                    p_fixe=p_fixe,
                    gams=gams,
                )
            )
            continue

        raise ValueError(f"Ligne de la section group priors non reconnue : {line!r}")

    return group_priors


def get_parameter_used_by_model(group_prior: GroupPrior) -> tuple[bool, bool]:
    """Détermine les paramètres k1/k2 actifs pour un modèle mutationnel ADN.

    JK -> aucun des deux actif, K2P/HKY -> k1 seul, TN -> les deux.

    Args:
        group_prior: Un GroupPrior de type model (group_prior.model is True).

    Returns:
        (k1_used, k2_used).

    Raises:
        NotImplementedError: Si group_prior n'est pas un model, si
            name_model est absent, ou si le modèle n'est pas géré
            (JK/K2P/HKY/TN).
    """
    if not group_prior.model:
        raise NotImplementedError(
            "Le modèle mutationnel n'est pas défini pour ce locus"
        )
    elif group_prior.name_model is None:
        raise NotImplementedError(
            "Le nom du modèle mutationnel n'est pas défini pour ce locus"
        )
    elif group_prior.name_model == "JK":
        return (False, False)
    elif group_prior.name_model == "K2P" or group_prior.name_model == "HKY":
        return (True, False)
    elif group_prior.name_model == "TN":
        return (True, True)
    else:
        raise NotImplementedError(
            f"Le modèle mutationnel {group_prior.name_model} n'est pas encore implémenté"
        )
