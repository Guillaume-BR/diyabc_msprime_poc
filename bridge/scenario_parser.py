"""
Parseur du texte de header.txt vers les structures définies dans
scenario_types.py.

Aucune valeur numérique n'est calculée ici (voir scenario_types.py pour
la justification). Ce module ne fait que de la transcription texte -> objets.

Vocabulaire géré : sample, merge, varNe, split (admixture) -- couvre les
6 scénarios du dataset human. Voir SplitEvent (scenario_types.py) pour la
sémantique exacte de split, vérifiée contre history.cpp/particuleC.cpp.
"""

import re
import warnings

from bridge.scenario_types import (
    MergeEvent,
    SampleEvent,
    Scenario,
    SplitEvent,
    VarNeEvent,
)

# Capture : "scenario 1 [0.16667] (16)" -> index=1, weight=0.16667, nlines=16
_SCENARIO_HEADER_RE = re.compile(r"^scenario\s+(\d+)\s+\[([\d.]+)\]\s+\((\d+)\)\s*$")


def split_scenario_blocks(header_text: str) -> list[str]:
    """Découpe le texte complet de header.txt en blocs bruts, un par
    scénario, chaque bloc commençant par sa ligne d'en-tête
    'scenario N [poids] (nlignes)' et s'arrêtant juste avant le bloc suivant
    (ou la fin de la section, ex: 'historical parameters priors')."""
    lines = header_text.splitlines()

    # Repère les lignes d'index où démarre chaque scénario
    start_indices = [
        i for i, line in enumerate(lines) if _SCENARIO_HEADER_RE.match(line.strip())
    ]
    if not start_indices:
        raise ValueError(
            "Aucun bloc 'scenario N [...] (...)' trouvé dans le texte fourni"
        )

    # Borne de fin pour le tout dernier scénario : la section des priors,
    # si elle est présente, sinon la fin du fichier.
    priors_section_index = next(
        (
            i
            for i, line in enumerate(lines)
            if line.strip().startswith("historical parameters priors")
        ),
        len(lines),
    )

    blocks = []
    for k, start in enumerate(start_indices):
        end = (
            start_indices[k + 1] if k + 1 < len(start_indices) else priors_section_index
        )
        block = "\n".join(lines[start:end]).strip()
        blocks.append(block)
    return blocks


def parse_scenario_block(block_text: str) -> Scenario:
    """Transforme un bloc brut (en commençant par la ligne 'scenario N [...]')
    en objet Scenario rempli."""
    lines = [line.strip() for line in block_text.splitlines() if line.strip()]

    header_match = _SCENARIO_HEADER_RE.match(lines[0])
    if not header_match:
        raise ValueError(
            f"Première ligne inattendue, pas un en-tête de scénario : {lines[0]!r}"
        )
    index = int(header_match.group(1))
    weight = float(header_match.group(2))

    # Deuxième ligne : tailles de population initiales, ex: "N1 N2 N3 N4"
    initial_pop_size_exprs = lines[1].split()

    events = []
    for line in lines[2:]:
        events.append(_parse_event_line(line))

    return Scenario(
        index=index,
        weight=weight,
        initial_pop_size_exprs=initial_pop_size_exprs,
        events=events,
    )


def _parse_event_line(line: str):
    """Transforme une ligne d'événement, ex: 't1 merge 2 1', en
    SampleEvent / MergeEventn_pops / VarNeEvent selon le mot-clé rencontré.

    Vocabulaire de référence : src-JMC-C++/history.cpp (ScenarioC::read_events).
    """
    tokens = line.split()
    time_expr, action = tokens[0], tokens[1]
    args = tokens[2:]

    if action == "sample":
        # ex: "0 sample 1" -> pop=1
        return SampleEvent(time_expr=time_expr, pop=int(args[0]))

    if action == "merge":
        # ex: "t1 merge 2 1" -> ancestral_pop=2 (survit), derived_pop=1 (disparaît)
        return MergeEvent(
            time_expr=time_expr,
            ancestral_pop=int(args[0]),
            derived_pop=int(args[1]),
        )

    if action == "varNe" or action == "varne":
        # ex: "t2-d3 varNe 3 Nbn3" -> pop=3, new_size_expr="Nbn3"
        return VarNeEvent(
            time_expr=time_expr,
            pop=int(args[0]),
            new_size_expr=args[1],
        )

    if action == "split":
        # ex: "t1 split 3 1 2 r1" -> derived_pop=3 (disparaît), ancestral_pop1=1
        # (reçoit chaque lignée de la pop 3 avec probabilité r1), ancestral_pop2=2
        # (reçoit le complément, 1-r1) -- vérifié dans history.cpp::ScenarioC::
        # read_events (ordre des colonnes pop/pop1/pop2/admixrate) et
        # particuleC.cpp::ParticleC::split_pop (tirage réel : vers pop1 si
        # random() < admixrate, sinon vers pop2).
        return SplitEvent(
            time_expr=time_expr,
            derived_pop=int(args[0]),
            ancestral_pop1=int(args[1]),
            ancestral_pop2=int(args[2]),
            admixture_rate=args[3],
        )

    raise NotImplementedError(
        f"Action '{action}' non gérée par ce parser (vocabulaire connu : "
        f"sample/merge/varNe/split). Ligne : {line!r}"
    )


def parse_header_scenarios(header_text: str) -> list[Scenario]:
    """Point d'entrée principal : header.txt complet -> liste de Scenario.

    Important : seule l'exception NotImplementedError est avalée ici,
    volontairement. Toute autre exception (erreur de parsing réelle, bug)
    doit continuer à se propager normalement.
    """
    blocks = split_scenario_blocks(header_text)
    scenarios = []
    for block in blocks:
        try:
            scenarios.append(parse_scenario_block(block))
        except NotImplementedError as e:
            first_line = block.splitlines()[0]
            warnings.warn(f"Bloc '{first_line}' ignoré : {e}", stacklevel=2)
    return scenarios
