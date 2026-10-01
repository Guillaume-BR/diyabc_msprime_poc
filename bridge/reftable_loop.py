"""Boucle d'itération produisant les nrec "particules" (lignes) d'un futur
reftable.bin : pour chaque particule, un tirage de scénario et de paramètres
distinct, et une simulation msprime complète (calcul des statistiques résumées
via compute_summary_statistics, 100% Python -- plus de subprocess ni de fichier
intermédiaire sur disque).

Parallélisé via ProcessPoolExecutor : chaque particule est indépendante
des autres (son propre tirage, sa propre simulation), donc
embarrassingly parallel.
"""

import os

# Doit s'exécuter AVANT le premier import de numpy (transitif, via
# bridge.pipeline plus bas) : numpy/BLAS lit ces variables une seule fois
# à l'initialisation de son pool de threads, pas à chaque appel. Sans ça,
# chacun des max_workers process de ProcessPoolExecutor essaie d'utiliser
# TOUS les cœurs pour ses propres opérations BLAS -- avec 16 workers sur
# une machine à 16 cœurs, jusqu'à 256 threads se battent pour 16 cœurs
# physiques. Mesuré empiriquement (toy_example3_scenario1, 1000
# particules, filtre MAF actif) : ~30 minutes sans ce fix, <60s avec.
# setdefault() pour ne jamais écraser un réglage déjà choisi explicitement
# par l'appelant.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import re
import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from bridge.ancestry_simulation import prepare_poolseq_observed_reads
from bridge.configuration import _SCENARIO_DRAW_SEED_OFFSET
from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.header_dataclasses import (
    DnaReplayContext,
    MicrosatReplayContext,
    Scenario,
    SnpReplayContext,
)
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    allele_bounds_per_locus,
    base_frequency_by_locus,
    count_individuals_per_sample,
    detect_snp_file_type,
    individual_sexes_per_sample,
    observed_count_sample,
    observed_microsatellites,
    observed_reads,
    observed_sequences,
    parse_maf_ratio,
    parse_mrc_ratio,
    parse_sex_ratio,
)
from bridge.parameter_sampling import draw_scenario
from bridge.pipeline import (
    compute_summary_statistics,
    compute_summary_statistics_dna,
    compute_summary_statistics_dna_from_values,
    compute_summary_statistics_from_values,
    compute_summary_statistics_microsat,
    compute_summary_statistics_microsat_from_values,
    read_header_text,
)
from bridge.prior_parser import (
    get_parameter_used_by_model,
    is_constant_prior,
    parse_group_priors,
    parse_priors,
)
from bridge.scenario_parser import is_serial_scenario, parse_header_scenarios

# _SCENARIO_DRAW_SEED_OFFSET : voir bridge/configuration.py. Décalage
# appliqué à la seed de particule avant de tirer le scénario
# (draw_scenario), pour ne JAMAIS partager la même seed brute avec le
# tirage des paramètres (draw_parameter_values, appelé plus loin dans
# compute_summary_statistics avec seed=seed, sans offset) -- sinon les
# deux tirages, bien qu'indépendants dans l'intention, consomment le
# MÊME premier random.random() sous-jacent (chaque fonction fait son
# propre random.Random(seed) frais), ce qui corrèle artificiellement le
# scénario tiré et la valeur du premier prior déclaré (ex: N1) :
# vérifié empiriquement sur toy_example5 -- scénario 1 (ra petit) donnait
# systématiquement un N1 trop bas, scénario 3 (ra grand) un N1 trop haut.
# Au-delà de _LOCUS_TYPE_SEED_OFFSET (bridge/configuration.py, 0..40_000_000) pour ne
# pas non plus recréer une collision avec ce décalage-là.


@dataclass
class ParticleResult:
    """Le résultat d'une particule : une future ligne du reftable.bin.

    Attributes:
        particle_index: L'index de la particule (0-based).
        scenario_index: L'index 1-based du scénario tiré pour cette
            particule.
        parameter_values: Les valeurs de paramètres historiques tirées,
            {nom: valeur}.
        summary_statistics: Les statistiques résumées calculées,
            {nom_colonne: valeur}.
    """

    particle_index: int
    scenario_index: int
    parameter_values: dict[str, float]
    summary_statistics: dict[str, float]
    group_priors_values: dict[str, float]


# ----------------------------------------------------------------------------
# Pour les fichiers .snp DIYABC : lecture, écriture, rejeux de tirages réels
# ----------------------------------------------------------------------------


def _heritage_types_declared(header_text: str) -> set[str]:
    """Extrait les types d'hérédité déclarés, quel que soit le format des loci.

    Les deux formats de 'loci description' ne portent pas l'information au
    même endroit : le format condensé la donne agrégée
    (`loci_counts_by_heritage`), le format détaillé un locus à la fois
    (`LociDescriptionDetailed.heritage`).

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        L'ensemble des types présents, parmi "A", "H", "X", "Y", "M".
    """
    loci = parse_loci_description(header_text)
    if isinstance(loci, list):
        return {locus.heritage for locus in loci}
    return {
        heritage
        for heritage, count in loci.loci_counts_by_heritage.items()
        if count > 0
    }


# Jetons de la ligne de queue correspondant à un paramètre mutationnel :
# "µseq_2", "k1seq_2", "k2seq_3", "pmic_1", "snimic_1". Les noms de
# statistiques ne peuvent pas collisionner (préfixe en capitales).
_MUTATION_PARAM_TOKEN_RE = re.compile(r"^(µ|mu|k1|k2|p|sni)(seq|mic)_\d+$")

# DIYABC écrit le préfixe mus_rate soit "µ" (U+00B5, cas des `te2`), soit "mu"
# en ASCII (cas de `toy_example1_ms`) -- les deux apparaissent dans de la vraie
# sortie, et le reftable suit toujours son header. On compare donc les jetons
# sous forme normalisée. `group_prior_column_names` produit la variante "µ" ;
# c'est une dette connue, un header écrit en "mu" en reçoit des noms de colonnes
# qui ne correspondent à rien (voir Open work).
# Les DEUX seules orthographes observées dans de la vraie sortie DIYABC :
# "µ" = U+00B5 MICRO SIGN (7 headers de reference/) et "mu" en ASCII
# (toy_example1_ms, 2 headers). Le mu grec U+03BC n'apparaît nulle part et
# n'est DÉLIBÉRÉMENT pas normalisé : un header qui le contiendrait serait
# auto-cohérent avec son propre reftable (bintotxt recopie la ligne de queue),
# donc les deux gardes passeraient, mais `group_prior_column_names` code U+00B5
# en dur et la chaîne de rejeu ne retrouverait jamais la colonne par son nom.
# Mieux vaut que la garde lève sur ce caractère que de l'absorber.
_MU_VARIANTS = ("µ", "mu")


def _normalize_mu(token: str) -> str:
    """Ramène le préfixe mus_rate d'un jeton à la variante "µ"."""
    for variant in _MU_VARIANTS:
        if token.startswith(variant):
            return "µ" + token[len(variant) :]
    return token


def check_header_trailer_line(header_text: str) -> None:
    """Vérifie que la ligne de queue du header est cohérente avec ses priors.

    La dernière ligne de `header.txt` ressemble à de la documentation de
    colonnes mais est **relue comme entrée** par
    `HeaderC::readHeaderAllStat`, qui en dérive
    `nparamhist = jetons - 1 - nstat - nparamut`. Un jeton parasite
    décale ce compte, corrompt l'état interne de DIYABC et produit des
    statistiques fausses de 300 % à 10000 % sans aucun message.

    Cette garde compare les jetons de paramètres mutationnels de la ligne
    de queue à ceux que le header implique réellement
    (`group_prior_column_names`, qui applique
    `get_parameter_used_by_model`). Exemple réel attrapé le 30/09 :
    `reference/toy_example2_ms_dna_50loci_JK` déclare `MODEL JK`, donc
    aucun `k1`, mais sa ligne de queue portait `k1seq_2` et `k1seq_3`
    recopiés de la variante `K2P` -- deux jetons de trop, reftable
    inexploitable.

    Ne vérifie PAS les noms de statistiques : `_historical_columns_order`
    couvre déjà les paramètres historiques lors de l'écriture d'un
    reftable texte, et `stats_group_parser` les colonnes déclarées.

    Args:
        header_text: Texte complet de header.txt.

    Raises:
        ValueError: Si les jetons mutationnels de la ligne de queue
            diffèrent de ceux impliqués par les `group priors`.
    """
    try:
        expected = group_prior_column_names(header_text)
    except ValueError:
        return  # pas de section `group priors` : rien à vérifier (datasets SNP)
    trailer_tokens = header_text.splitlines()[-1].split()
    in_trailer = [
        _normalize_mu(tok)
        for tok in trailer_tokens
        if _MUTATION_PARAM_TOKEN_RE.match(tok)
    ]
    expected = [_normalize_mu(name) for name in expected]
    if in_trailer == expected:
        return
    missing = [t for t in expected if t not in in_trailer]
    unexpected = [t for t in in_trailer if t not in expected]
    raise ValueError(
        "Ligne de queue du header incohérente avec la section `group priors` : "
        f"attendu {expected}, trouvé {in_trailer}"
        + (f" -- en trop : {unexpected}" if unexpected else "")
        + (f" -- manquants : {missing}" if missing else "")
        + ". HeaderC::readHeaderAllStat dérive nparamhist du nombre de jetons de "
        "cette ligne : un écart corrompt silencieusement l'état interne de DIYABC. "
        "Régénérer la ligne de queue pour ce header, puis relancer DIYABC -- "
        "un reftable déjà produit avec ce header ne se rattrape pas."
    )


def check_real_reftable_matches_header(
    header_text: str, real_reftable_path: str | Path
) -> None:
    """Vérifie qu'un vrai reftable DIYABC a été produit avec CE header.

    La première ligne d'un reftable texte porte ses noms de colonnes,
    écrits par `reftable.cpp::bintotxt` depuis la ligne de queue du
    header (`entetehist`, voir "Text-reftable column order" dans
    CLAUDE.md). Comparer les deux détecte donc deux fautes qu'aucun
    horodatage ne distingue : un header **édité après** la génération du
    reftable, et un reftable **produit avec un header corrompu** dont la
    correction ultérieure ne rattrape rien.

    Préférer ce contrôle à une comparaison de dates : DIYABC réécrit
    `headerRF.txt` APRÈS le reftable dans le même run (mesuré : 6 s
    d'écart sur `toy_example2_ms_dna_50loci_TN`), donc un run correct a
    toujours un header plus récent que son reftable.

    Ne compare que les jetons de paramètres mutationnels, pour la même
    raison que `check_header_trailer_line` : ce sont eux que le modèle
    déclaré détermine, et ils ne peuvent pas collisionner avec un nom de
    statistique. Les deux orthographes du préfixe mus_rate (`µ` et `mu`)
    sont normalisées.

    Args:
        header_text: Texte complet de header.txt.
        real_reftable_path: Chemin du reftable réel, au format texte.

    Raises:
        ValueError: Si les deux jeux de jetons diffèrent.
    """
    path = Path(real_reftable_path)
    with path.open() as f:
        first_line = f.readline()
    in_reftable = [
        _normalize_mu(tok)
        for tok in first_line.split()
        if _MUTATION_PARAM_TOKEN_RE.match(tok)
    ]
    in_trailer = [
        _normalize_mu(tok)
        for tok in header_text.splitlines()[-1].split()
        if _MUTATION_PARAM_TOKEN_RE.match(tok)
    ]
    if in_reftable == in_trailer:
        return
    raise ValueError(
        f"{path.name} n'a pas été produit avec ce header : la ligne de queue "
        f"annonce {in_trailer}, le reftable porte {in_reftable}. Soit le header "
        "a été édité après la génération, soit le reftable vient d'un header "
        "corrompu -- dans les deux cas il faut relancer DIYABC, une correction "
        "du header ne rattrape pas un reftable déjà écrit. Ne pas se fier aux "
        "dates : DIYABC écrit headerRF.txt après le reftable."
    )


def raise_if_serial_with_sex_linked_loci(header_text: str) -> None:
    """Refuse la combinaison non implémentée « sériel + loci <X>/<Y> ».

    Les quatre constructeurs d'échantillons sexués
    (`build_sex_stratified_samples_argument` / `build_male_only_samples_
    argument` et leurs jumeaux `_ms_dna`) ne sont PAS sériels-conscients :
    le dispatch par locus écrase les `SampleSet` sériels construits par
    `build_sample_sets_from_scenario`, et ces constructeurs redérivent
    leurs clés msprime depuis la POSITION du bloc observé
    (`f"pop{i}"`, i = i-ème bloc `POP`) au lieu de la population que
    l'événement `sample` désigne. Sur un scénario sériel, `i` dépasse le
    nombre de populations échantillonnées.

    Pourquoi lever plutôt que laisser msprime se plaindre : les deux
    issues possibles sont très inégales, et la plus dangereuse est
    atteignable. Mesuré sur reference/toy_example2_ms_dna, dont le
    scénario 1 déclare 5 populations pour 2 échantillons et dont
    `build_demography` crée `pop1`..`pop5` toutes `initially_active=True`
    -- un `samples={"pop1": 5, "pop2": 5, "pop3": 5}` est accepté SANS
    erreur (15 nœuds au lieu de 10), les échantillons étant rattachés à
    une population ancestrale jamais échantillonnée. Ce n'est que
    lorsque l'indice inventé dépasse l'indice maximal que msprime lève
    un `KeyError` bruyant.

    Appelée au montage du run, une fois, et non par particule : l'échec
    doit précéder la première particule, pas survenir à la 743e. Levée
    dès qu'UN scénario du header est sériel, sans attendre de savoir
    lequel `draw_scenario` tirera.

    Args:
        header_text: Texte complet de header.txt.

    Raises:
        NotImplementedError: Si un scénario est sériel et que le header
            déclare des loci <X> ou <Y>.
    """
    sex_linked = _heritage_types_declared(header_text) & {"X", "Y"}
    if not sex_linked:
        return
    serial_indexes = [
        scenario.index
        for scenario in parse_header_scenarios(header_text)
        if is_serial_scenario(scenario)
    ]
    if not serial_indexes:
        return
    raise NotImplementedError(
        "Échantillonnage sériel non implémenté pour les loci liés au sexe : "
        f"scénario(s) {serial_indexes} échantillonne(nt) deux fois la même "
        f"population, et le header déclare des loci {sorted(sex_linked)}. "
        "Les constructeurs d'échantillons sexués ignorent les SampleSet "
        "sériels et redérivent leurs noms de population depuis la position "
        "du bloc observé ; les échantillons seraient silencieusement "
        "rattachés à une population ancestrale (vérifié sur "
        "toy_example2_ms_dna), sans erreur msprime."
    )


# ── Tirage indépendant de scénario + paramètres, par particule ────────────


def _run_single_particle(
    particle_index: int,
    context: SnpReplayContext,
    scenarios: list[Scenario],
    *,
    num_loci: int | None = None,
    observed_reads_per_locus: list[dict[str, tuple[int, int]]] = None,
    stats_filter: str,
) -> ParticleResult:
    """Calcule une seule particule.

    Fonction top-level (picklable), appelée par chaque worker du
    ProcessPoolExecutor.

    La seed utilisée est dérivée de particle_index, garantissant un
    tirage distinct et reproductible par particule (même particle_index
    -> même résultat, peu importe l'ordre d'exécution des workers).

    IMPORTANT : seed = particle_index + 1, jamais particle_index seul.
    msprime.sim_ancestry rejette explicitement seed=0 (ValueError "seeds
    must be greater than 0 and less than 2^32") -- vérifié empiriquement.
    Donc particle_index=0 (le cas le plus probable, première particule)
    utilise seed=1, pas seed=0.

    Args:
        particle_index: L'index de la particule (0-based).
        context: Le contexte de la simulation.
        scenarios: Les scénarios candidats (chaque particule tire le
            sien).
        num_loci: Voir pipeline.compute_summary_statistics.
        observed_reads_per_locus: PoolSeq uniquement, pré-calculé une
            fois pour toute la boucle (voir run_reftable_simulation).
        stats_filter: "ALL" ou "HEADER", voir
            pipeline.compute_summary_statistics.

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    drawn_scenario = draw_scenario(scenarios, seed + _SCENARIO_DRAW_SEED_OFFSET)

    summary_statistics, parameter_values = compute_summary_statistics(
        context=context,
        scenario_index=drawn_scenario.index,
        num_loci=num_loci,
        seed=seed,
        stats_filter=stats_filter,
        observed_reads_per_locus=observed_reads_per_locus,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=drawn_scenario.index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
        group_priors_values={},
    )


def run_reftable_simulation(
    reference_directory: str | Path,
    scenarios: list[Scenario],
    *,
    num_loci: int | None = None,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Produit nrec particules (lignes de reftable.bin) en parallèle.

    N'écrit rien sur disque par particule (compute_summary_statistics
    est 100% Python, en mémoire).

    Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
    nrec-1), pas dans l'ordre de complétion des workers -- important
    pour la reproductibilité de l'ordre des lignes du reftable final.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .snp observé.
        scenarios: La liste des scénarios candidats (typiquement TOUS
            les scénarios déclarés dans header.txt) : chaque particule
            tire le SIEN au hasard, pondéré par son `weight` (voir
            parameter_sampling.draw_scenario, sémantique vérifiée
            contre particuleC.cpp::ParticleC::drawscenario) -- une même
            particule peut donc finir sur n'importe lequel des
            scénarios de la liste, pas forcément le même pour toutes.
        num_loci: Voir pipeline.compute_summary_statistics.
        nrec: Le nombre de particules à produire.
        stats_filter: "ALL" ou "HEADER", voir
            pipeline.compute_summary_statistics.
        max_workers: Le nombre de process en parallèle (défaut : laissé
            à ProcessPoolExecutor, généralement le nombre de cœurs
            disponibles).

    Returns:
        La liste des ParticleResult, dans l'ordre de particle_index (0
        à nrec-1).
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    snp_filename = header_text.splitlines()[0].strip()
    snp_path = reference_directory / snp_filename
    loci_description = parse_loci_description(header_text)
    snp_file_type = detect_snp_file_type(snp_path)
    counts_per_sample = count_individuals_per_sample(snp_file_path=snp_path)
    sex_ratio = parse_sex_ratio(snp_path)
    maf_ratio = parse_maf_ratio(snp_path)
    mrc_ratio = parse_mrc_ratio(snp_path)
    reads_observed = (
        observed_reads(snp_path, loci_description.loci_counts_by_heritage["A"])
        if snp_file_type == "POOL"
        else None
    )
    sexes_per_sample = (
        individual_sexes_per_sample(snp_file_path=snp_path)
        if snp_file_type == "IND"
        else {}
    )

    context = SnpReplayContext(
        header_text=header_text,
        snp_path=snp_path,
        snp_file_type=snp_file_type,
        loci_description=loci_description,
        counts_per_sample=counts_per_sample,
        sex_ratio=sex_ratio,
        maf_ratio=maf_ratio,
        mrc_ratio=mrc_ratio,
        reads_observed=reads_observed,
        sexes_per_sample=sexes_per_sample,
    )

    results_by_index: dict[int, ParticleResult] = {}

    observed_reads_per_locus = None
    if snp_file_type == "POOL":
        observed_reads_per_locus = prepare_poolseq_observed_reads(context)
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle,
                particle_index,
                context,
                scenarios,
                num_loci=num_loci,
                stats_filter=stats_filter,
                observed_reads_per_locus=observed_reads_per_locus,
            ): particle_index
            for particle_index in range(nrec)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 100 == 0 or done == nrec:
                print(f"  {done}/{nrec} particules complétées")

    return [results_by_index[i] for i in range(nrec)]


# ── Helper partagé (rejeu ET écriture, voir les deux sections suivantes) ──


def _kept_param_names_by_scenario(
    priors: list, scenarios: list[Scenario]
) -> dict[int, list[str]]:
    """Calcule les noms de paramètres à garder, par scénario.

    Pour chaque scénario, la liste (dans l'ordre de déclaration des
    priors) des noms de paramètres à garder : non constants
    (is_constant_prior) ET référencés par CE scénario précis
    (get_parameter_names_used_by_scenario) -- même filtre à deux
    critères que l'ancienne version single-scenario, appliqué
    séparément par scénario.

    Args:
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.

    Returns:
        Un dict {scenario_index: [nom_paramètre, ...]}.
    """
    result = {}
    for scenario in scenarios:
        used_param_names = get_parameter_names_used_by_scenario(scenario)
        result[scenario.index] = [
            p.name
            for p in priors
            if not is_constant_prior(p) and p.name in used_param_names
        ]
    return result


def _historical_columns_order(
    header_line: str, priors: list, scenarios: list[Scenario]
) -> list[str]:
    """Ordre des colonnes de paramètres historiques d'un reftable texte réel.

    Complément de _kept_param_names_by_scenario : elle donne l'ENSEMBLE
    des noms présents pour un scénario, celle-ci donne leur ORDRE. Cet
    ordre est celui de la ligne d'en-tête du fichier lu (`entetehist`,
    c'est-à-dire la ligne trailer du header, reproduite telle quelle par
    reftable.cpp::bintotxt), PAS celui de la section `historical
    parameters priors` -- les deux coïncident sur les headers générés
    par l'interface, mais peuvent diverger sur un header édité à la
    main, ce qui mélabelle alors silencieusement toutes les valeurs
    (trouvé le 2026-09-21 sur toy_example1_ms_modified).

    Args:
        header_line: La première ligne du reftable texte réel
            (`scenario N1 N2 ... <group priors> <stats>`).
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.

    Returns:
        Les noms de colonnes historiques, dans l'ordre du fichier.

    Raises:
        ValueError: Si le nombre de noms lus ne correspond pas à l'union
            des priors non constants utilisés par au moins un scénario
            (typiquement une faute de frappe dans la ligne trailer).
    """
    tokens = header_line.split()[1:]
    prior_names = {p.name for p in priors}
    columns = []
    for token in tokens:
        if token not in prior_names:
            break
        columns.append(token)
    kept = _kept_param_names_by_scenario(priors, scenarios)
    expected = {name for names in kept.values() for name in names}
    if len(columns) != len(expected):
        raise ValueError(
            "Colonnes historiques du reftable incohérentes avec header.txt : "
            f"lues {columns}, attendues {sorted(expected)}, "
            f"différence {sorted(set(columns) ^ expected)}"
        )
    return columns


# ── Écriture du reftable (formats binaire et texte) ────────────────────────


def write_reftable_bin(
    results: list[ParticleResult],
    priors: list,
    scenarios: list[Scenario],
    output_path: str | Path,
) -> None:
    """Écrit un reftable.bin au format binaire DIYABC.

    Vérifié contre reftable.cpp et un vrai reftableRF.bin
    multi-scénario -- voir docs/synthese_diyabc_msprime.docx section 5.

    IMPORTANT -- format à LONGUEUR VARIABLE par ligne, PAS d'union de
    colonnes ni de valeur NA écrite sur disque : chaque ligne écrit
    SEULEMENT nparam[scenario-1] floats de paramètres, ceux de son
    propre scénario -- vérifié empiriquement contre un vrai
    reftableRF.bin multi-scénario (dataset MER modelchoice/PoolSeq) et
    contre reftable.cpp (boucle d'écriture indexée par
    nparam[numscen-1]). La reconstruction en matrice rectangulaire avec
    NA pour les colonnes non concernées est une responsabilité du
    LECTEUR (readReftable.R), jamais de l'écrivain -- ne PAS essayer de
    remplir les colonnes manquantes ici.

    Ne gère PAS les paramètres de mutation (absents de human) -- à
    ajouter (toujours en dernière position, après les paramètres
    démographiques -- voir readReftable.R) si un dataset avec
    microsatellites/séquences est traité plus tard.

    Args:
        results: Les ParticleResult à écrire, une ligne par résultat.
        priors: Les priors déclarés dans header.txt.
        scenarios: La liste de TOUS les scénarios candidats déclarés
            dans header.txt (pas seulement ceux effectivement tirés
            dans `results`) : nscen = len(scenarios), et le numéro de
            scénario écrit par ligne est le numéro 1-indexed du
            header.txt (scenario.index), jamais renuméroté localement
            -- vérifié dans particuleC.cpp::drawscenario et
            reftable.cpp. nparam[i] (nombre de paramètres non
            constants, référencés, du i-ème scénario de `scenarios`)
            pilote directement la taille de chaque enregistrement,
            comme dans reftable.cpp.
        output_path: Chemin où écrire le fichier binaire.

    Raises:
        ValueError: Si results est vide, ou si results contient un
            scenario_index absent de `scenarios`.
    """
    if not results:
        raise ValueError("results est vide : au moins une particule est requise")

    known_indices = {s.index for s in scenarios}
    unknown = {r.scenario_index for r in results} - known_indices
    if unknown:
        raise ValueError(
            f"results contient des scenario_index absents de `scenarios` : {unknown}"
        )

    kept_param_names_by_scenario = _kept_param_names_by_scenario(priors, scenarios)
    stat_names = list(results[0].summary_statistics.keys())

    nrec = len(results)
    nscen = len(scenarios)
    nrecscen = [
        sum(1 for r in results if r.scenario_index == s.index) for s in scenarios
    ]
    group_priors_names = list(results[0].group_priors_values)
    nparam = [
        len(kept_param_names_by_scenario[s.index]) + len(group_priors_names)
        for s in scenarios
    ]
    nstat = len(stat_names)

    with open(output_path, "wb") as f:
        f.write(struct.pack("<i", nrec))
        f.write(struct.pack("<i", nscen))
        for n in nrecscen:
            f.write(struct.pack("<i", n))
        for n in nparam:
            f.write(struct.pack("<i", n))
        f.write(struct.pack("<i", nstat))
        for result in results:
            f.write(struct.pack("<i", result.scenario_index))
            for name in kept_param_names_by_scenario[result.scenario_index]:
                f.write(struct.pack("<f", result.parameter_values[name]))
            for name in group_priors_names:
                f.write(struct.pack("<f", result.group_priors_values[name]))
            for name in stat_names:
                f.write(struct.pack("<f", result.summary_statistics[name]))


def write_reftable_txt(
    results: list[ParticleResult],
    priors: list,
    scenarios: list[Scenario],
    output_path: str | Path,
) -> None:
    """Écrit les résultats au format texte de DIYABC :
    first_records_of_the_reference_table_0.txt.

    Format reproduit depuis particleset.cpp / header.cpp :
      - Ligne 1 : noms de colonnes, chacun centré sur 14 caractères
        (fonction C++ centre(s1, 14)).
      - Lignes suivantes : "%3d  " pour le numéro de scénario, puis
        "  %12.6f" pour chaque paramètre et statistique.

    Note sur le format des paramètres : le C++ distingue categ<2 (%12.0f,
    entiers) vs categ>=2 (%12.3f, flottants), mais cette catégorie n'est
    pas exposée dans les priors Python. On utilise %12.6f uniformément --
    suffisant pour la comparaison statistique des distributions.

    Le texte utilise un jeu de colonnes de paramètres FIXE : l'UNION
    (dans l'ordre de déclaration des priors) des paramètres utilisés par
    au moins un des `scenarios`. Pour une ligne dont le scénario tiré
    n'utilise pas tel paramètre, on écrit :
      - sa valeur RÉELLEMENT TIRÉE si elle est présente dans
        r.parameter_values (cas de `run_reftable_simulation` :
        draw_parameter_values tire TOUS les priors déclarés,
        indépendamment du scénario, voir build_random_demography) ;
      - `nan` sinon (cas de `replay_reftable_simulation` :
        parse_real_reftable_params ne fournit QUE les paramètres du
        scénario propre à chaque ligne, puisque c'est aussi ce que
        DIYABC écrit réellement -- voir sa docstring).
    Jamais une case VIDE dans les deux cas.

    IMPORTANT : ne JAMAIS laisser de case vide ici, même si elle
    représente une valeur non pertinente pour le scénario de la ligne --
    un parseur par espaces (ex: pandas read_csv(sep=r'\\s+'), ou même un
    simple line.split()) ne produit AUCUN token pour une case vide,
    ce qui décale d'une colonne TOUTES les valeurs suivantes sur la
    ligne. Bug découvert empiriquement (comparaison DIYABC/msprime sur
    toy_example5_modif, colonne 'r' non utilisée par le scénario actif
    laissée en blanc -> décalage systématique des statistiques sur
    les 1000 lignes, provoquant des "écarts" massifs et incohérents qui
    n'avaient rien à voir avec un vrai écart de simulation). DIYABC
    lui-même (particleset.cpp, écriture de first_records_of_the_
    reference_table_N.txt) laisse une case vide dans ce cas précis --
    c'est cette case vide, relue naïvement, qui a provoqué le même
    décalage silencieux lors du premier test avec un reftable réel
    multi-scénario (voir parse_real_reftable_params).

    Args:
        results: Les ParticleResult à écrire, une ligne par résultat.
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats (détermine l'union des
            colonnes de paramètres).
        output_path: Chemin où écrire le fichier texte.

    Raises:
        ValueError: Si results est vide.
    """
    if not results:
        raise ValueError("results est vide : au moins une particule est requise")

    kept_param_names_by_scenario = _kept_param_names_by_scenario(priors, scenarios)
    used_by_any = {
        name for names in kept_param_names_by_scenario.values() for name in names
    }
    all_param_names = [p.name for p in priors if p.name in used_by_any]
    all_group_priors_names = list(results[0].group_priors_values)
    stat_names = list(results[0].summary_statistics.keys())

    def _centre(s: str, width: int = 14) -> str:
        """Reproduit centre() de mesutils.cpp : centrage sur width chars."""
        return s.center(width)

    with open(output_path, "w", encoding="utf-8") as f:
        # En-tête : "scenario" + noms des paramètres + noms des stats
        header = (
            _centre("scenario")
            + "".join(_centre(n) for n in all_param_names)
            + "".join(_centre(n) for n in all_group_priors_names)
            + "".join(_centre(n) for n in stat_names)
        )
        f.write(header + "\n")

        # Une ligne par particule
        for r in results:
            line = f"{r.scenario_index:3d}  "
            for name in all_param_names:
                line += f"  {r.parameter_values.get(name, float('nan')):12.6f}"
            for name in all_group_priors_names:
                line += f"  {r.group_priors_values[name]:12.6f}"
            for name in stat_names:
                line += f"  {r.summary_statistics[name]:12.6f}"
            f.write(line + "\n")


def rewrite_real_reftable_txt(
    input_path: str | Path,
    output_path: str | Path,
    priors: list,
    scenarios: list[Scenario],
) -> None:
    """Réécrit un reftable RÉEL de DIYABC en un texte à colonnes de largeur
    FIXE.

    Ex: first_records_of_the_reference_table_0.txt. Remplace les cases
    vides de DIYABC (paramètre non utilisé par le scénario de la ligne,
    voir parse_real_reftable_params) par `nan` -- jamais une case vide.

    Nécessaire pour toute lecture EXTERNE du fichier DIYABC brut avec un
    parseur par espaces générique (ex: `pandas.read_csv(sep=r'\\s+')`,
    utilisé dans les notebooks de comparaison) : sans cette réécriture,
    une ligne dont le scénario n'utilise pas tous les paramètres
    déclarés a MOINS de tokens que la ligne d'en-tête ne le laisse
    penser, ce qui décale silencieusement toutes les colonnes
    suivantes (statistiques comprises) sur cette ligne -- même piège
    que documenté dans parse_real_reftable_params/write_reftable_txt,
    mais ici côté fichier DIYABC lui-même plutôt que côté notre pipeline.

    Le fichier réécrit garde les colonnes de paramètres du fichier
    d'entrée, dans l'ordre de SA ligne d'en-tête (`entetehist`, voir
    _historical_columns_order) -- pas dans l'ordre de déclaration des
    priors, qui peut en différer sur un header édité à la main. Les
    colonnes suivantes (priors de groupe s'il y en a, puis statistiques)
    sont recopiées telles quelles. Le résultat se compare colonne à
    colonne, PAR NOM, avec un reftable_msprime généré par
    run_reftable_simulation/replay_reftable_simulation.

    Args:
        input_path: Chemin du reftable réel brut (format texte).
        output_path: Chemin où écrire le fichier réécrit.
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.
    """
    kept_by_scenario = _kept_param_names_by_scenario(priors, scenarios)

    lines = [line for line in Path(input_path).read_text().splitlines() if line.strip()]
    all_param_names = _historical_columns_order(lines[0], priors, scenarios)
    header_tokens = lines[0].split()
    stat_names = header_tokens[1 + len(all_param_names) :]
    data_lines = lines[1:]

    def _centre(s: str, width: int = 14) -> str:
        return s.center(width)

    with open(output_path, "w", encoding="utf-8") as f:
        header = (
            _centre("scenario")
            + "".join(_centre(n) for n in all_param_names)
            + "".join(_centre(n) for n in stat_names)
        )
        f.write(header + "\n")

        for line in data_lines:
            tokens = line.split()
            scenario_index = int(tokens[0])
            kept = set(kept_by_scenario[scenario_index])
            param_names = [name for name in all_param_names if name in kept]
            n_params = len(param_names)
            param_values = dict(zip(param_names, tokens[1 : 1 + n_params], strict=True))
            stat_values = dict(zip(stat_names, tokens[1 + n_params :], strict=True))

            out_line = f"{scenario_index:3d}  "
            for name in all_param_names:
                value = (
                    float(param_values[name]) if name in param_values else float("nan")
                )
                out_line += f"  {value:12.6f}"
            for name in stat_names:
                out_line += f"  {float(stat_values[name]):12.6f}"
            f.write(out_line + "\n")


# ── Point d'entrée haut niveau (compose tirage indépendant + écriture) ────


def simulate_from_directory(
    test_directory: str | Path,
    context: SnpReplayContext,
    *,
    num_loci: int | None = None,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Point d'entrée pour un sous-dossier de test sous reference/.

    Ex: reference/mon_test/, qui ne contient au départ qu'un header.txt
    (ou headerRF.txt, repli si absent -- voir pipeline.read_header_text)
    et le fichier .snp observé (nommé sur la première ligne du header,
    pas un nom fixe -- voir pipeline.simulate_particle_genotypes).

    Tire les scénarios candidats pondérés par leur `weight` parmi TOUS
    ceux déclarés dans le header (voir parameter_sampling.draw_scenario),
    simule nrec particules, et écrit le résultat dans
    test_directory/reftable_msprime.txt ET .bin -- jamais "reftable.txt",
    pour ne pas être confondu avec le first_records_of_the_reference_
    table_0.txt qu'un vrai run DIYABC produirait dans le même dossier.

    Args:
        test_directory: Le sous-dossier de test.
        num_loci: Voir pipeline.compute_summary_statistics.
        nrec: Le nombre de particules à produire.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle.

    Returns:
        Les ParticleResult (utile pour appeler write_reftable_bin en
        plus, si besoin -- déjà fait ici aussi).
    """
    header_text = context.header_text
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)

    priors, _ = parse_priors(header_text)
    print(f"{len(priors)} priors parsé depuis {test_directory}/header.txt")
    scenarios = parse_header_scenarios(header_text)
    print(f"{len(scenarios)} scénarios parsés depuis {test_directory}/header.txt")

    print(
        f"Simulation de {nrec} particules avec {num_loci if num_loci is not None else 'tous les'} loci)"
    )
    results = run_reftable_simulation(
        reference_directory=test_directory,
        scenarios=scenarios,
        num_loci=num_loci,
        nrec=nrec,
        stats_filter=stats_filter,
        max_workers=max_workers,
    )

    write_reftable_txt(
        results, priors, scenarios, test_directory / "reftable_msprime.txt"
    )
    print(f"Écriture du reftable texte dans {test_directory}/reftable_msprime.txt")

    write_reftable_bin(
        results, priors, scenarios, test_directory / "reftable_msprime.bin"
    )
    print(f"Écriture du reftable binaire dans {test_directory}/reftable_msprime.bin")

    return results


# ── Rejeu des tirages RÉELS de DIYABC (comparaison appariée) ──────────────


def parse_real_reftable_params(
    path: str | Path, priors: list, scenarios: list[Scenario]
) -> list[tuple[int, dict[str, float]]]:
    """Lit un reftable RÉEL produit par DIYABC.

    Ex: first_records_of_the_reference_table_0.txt. Extrait, ligne par
    ligne, le scénario tiré et les valeurs de paramètres RÉELLEMENT
    tirées par DIYABC -- pour les rejouer ensuite côté msprime (voir
    replay_reftable_simulation), afin de comparer les deux simulateurs
    sur EXACTEMENT les mêmes tirages de priors, sans le biais de deux
    tirages indépendants.

    Le nombre de colonnes de paramètres RÉELLEMENT présentes sur une
    ligne dépend du SCÉNARIO de CETTE ligne précise, pas d'un total fixe
    pour tout le fichier : DIYABC (particleset.cpp, écriture de
    first_records_of_the_reference_table_N.txt) parcourt l'UNION de tous
    les noms de paramètres déclarés et, pour un nom NON utilisé par le
    scénario de la ligne courante, écrit une CASE VIDE (aucune valeur,
    juste des espaces) au lieu d'une valeur ou d'un NA -- un parseur par
    espaces ne produit alors AUCUN token pour cette case, ce qui décale
    silencieusement toutes les colonnes suivantes si on suppose un
    nombre de colonnes fixe (union) pour toutes les lignes. Repéré en
    testant un reftable réellement multi-scénario (voir aussi le même
    principe côté écriture dans write_reftable_txt, qui l'évite en
    n'écrivant jamais de case vide).

    On lit donc, pour CHAQUE ligne, les tokens de paramètres
    correspondant SPÉCIFIQUEMENT au scénario de cette ligne, avec deux
    sources distinctes qu'il ne faut pas confondre :
    - l'ENSEMBLE des noms présents vient de _kept_param_names_by_scenario
      (non constant + utilisé par ce scénario, même filtre que
      write_reftable_txt/write_reftable_bin) -- ce qui gère aussi le cas
      single-scénario où certains priors sont devenus constants
      (is_constant_prior) et donc absents des colonnes de sortie ;
    - leur ORDRE vient de la ligne d'en-tête du fichier lu
      (_historical_columns_order), c'est-à-dire de `entetehist` tel que
      reftable.cpp::bintotxt l'imprime -- PAS de l'ordre de déclaration
      des priors dans header.txt. Les deux coïncident sur les headers
      générés par l'interface, mais un header édité à la main peut les
      faire diverger, et l'ordre de déclaration mélabelle alors
      silencieusement toutes les valeurs (2026-09-21,
      toy_example1_ms_modified : la valeur de t32 lue comme t41, etc.).

    Args:
        path: Chemin du reftable réel (format texte).
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.

    Returns:
        Une liste de tuples (scenario_index, {nom_paramètre: valeur}),
        un par ligne du reftable (dans l'ordre du fichier).

    Raises:
        ValueError: Voir _historical_columns_order (ligne d'en-tête
            incohérente avec header.txt).
    """
    priors_kept_by_scenario = _kept_param_names_by_scenario(priors, scenarios)

    lines = [line for line in Path(path).read_text().splitlines() if line.strip()]
    data_lines = lines[1:]
    hist_order = _historical_columns_order(lines[0], priors, scenarios)
    rows = []
    for line in data_lines:
        tokens = line.split()
        scenario_index = int(tokens[0])
        kept = set(priors_kept_by_scenario[scenario_index])
        param_names = [name for name in hist_order if name in kept]
        values = {name: float(tokens[1 + i]) for i, name in enumerate(param_names)}
        rows.append((scenario_index, values))
    return rows


def _run_single_particle_from_values(
    particle_index: int,
    context: SnpReplayContext,
    scenario_index: int,
    values: dict[str, float],
    *,
    num_loci: int | None = None,
    observed_reads_per_locus: list[dict[str, tuple[int, int]]] = None,
    stats_filter: str,
) -> ParticleResult:
    """Variante de _run_single_particle qui NE TIRE AUCUN paramètre.

    Rejoue (scenario_index, values) tels que fournis -- typiquement
    issus de parse_real_reftable_params.

    Args:
        particle_index: L'index de la particule (0-based).
        context: Le contexte de rejouage.
        reference_directory: Le dossier contenant header.txt et le
            fichier .snp observé.
        scenario_index: L'index 1-based du scénario déjà tiré par
            DIYABC pour cette particule.
        values: Les valeurs de paramètres historiques déjà connues,
            {nom: valeur}.
        num_loci: Voir pipeline.compute_summary_statistics_from_values.
        observed_reads_per_locus: PoolSeq uniquement.
        stats_filter: "ALL" ou "HEADER".

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    summary_statistics = compute_summary_statistics_from_values(
        context=context,
        scenario_index=scenario_index,
        values=values,
        num_loci=num_loci,
        seed=seed,
        stats_filter=stats_filter,
        observed_reads_per_locus=observed_reads_per_locus,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=scenario_index,
        parameter_values=values,
        summary_statistics=summary_statistics,
        group_priors_values={},
    )


def replay_reftable_simulation(
    reference_directory: str | Path,
    priors: list,
    scenarios: list[Scenario],
    real_reftable_path: str | Path,
    num_loci: int | None = None,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Rejoue, particule par particule, les tirages de paramètres RÉELLEMENT
    effectués par DIYABC.

    Lit un reftable existant (ex:
    first_records_of_the_reference_table_0.txt) -- au lieu d'en tirer de
    nouveaux indépendamment comme run_reftable_simulation.

    Chaque particule msprime utilise EXACTEMENT le même (N1,N2,N3,ta,
    ts,...) que la particule DIYABC de même particle_index (même ordre
    que les lignes du fichier réel) : tout écart entre les deux résulte
    donc uniquement du moteur de simulation, jamais d'un tirage de prior
    différent -- permet une comparaison appariée ligne à ligne, pas
    seulement une comparaison de distributions agrégées.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .snp observé.
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.
        real_reftable_path: Chemin du reftable réel à rejouer.
        num_loci: Voir pipeline.compute_summary_statistics_from_values.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle.

    Returns:
        Les ParticleResult dans le MÊME ORDRE que les lignes du fichier
        réel -- réutilisable tel quel par write_reftable_txt/
        write_reftable_bin (même type que run_reftable_simulation).
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    snp_filename = header_text.splitlines()[0].strip()
    snp_path = reference_directory / snp_filename
    loci_description = parse_loci_description(header_text)
    snp_file_type = detect_snp_file_type(snp_path)
    counts_per_sample = count_individuals_per_sample(snp_file_path=snp_path)
    sex_ratio = parse_sex_ratio(snp_path)
    maf_ratio = parse_maf_ratio(snp_path)
    mrc_ratio = parse_mrc_ratio(snp_path)
    reads_observed = (
        observed_reads(snp_path, loci_description.loci_counts_by_heritage["A"])
        if snp_file_type == "POOL"
        else None
    )
    sexes_per_sample = (
        individual_sexes_per_sample(snp_file_path=snp_path)
        if snp_file_type == "IND"
        else {}
    )

    context = SnpReplayContext(
        header_text=header_text,
        snp_path=snp_path,
        snp_file_type=snp_file_type,
        loci_description=loci_description,
        counts_per_sample=counts_per_sample,
        sex_ratio=sex_ratio,
        maf_ratio=maf_ratio,
        mrc_ratio=mrc_ratio,
        reads_observed=reads_observed,
        sexes_per_sample=sexes_per_sample,
    )

    # On lit les sorties de diyabc (scénario tiré + valeurs de paramètres RÉELLEMENT tirées) pour
    # les rejouer ensuite côté msprime, afin de comparer les deux simulateurs sur EXACTEMENT
    # les mêmes tirages de priors.
    check_real_reftable_matches_header(header_text, real_reftable_path)
    rows = parse_real_reftable_params(real_reftable_path, priors, scenarios)

    results_by_index: dict[int, ParticleResult] = {}
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_from_values,
                particle_index,
                context,
                scenario_index,
                values,
                num_loci=num_loci,
                stats_filter=stats_filter,
            ): particle_index
            for particle_index, (scenario_index, values) in enumerate(rows)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 100 == 0 or done == len(rows):
                print(
                    f"  {done}/{len(rows)} particules rejouées à partir du reftable réel"
                )
    return [results_by_index[i] for i in range(len(rows))]


# ----------------------------------------------------------------------------
# Pour les séquences ADN : lecture, écriture, rejeux de tirages réels
# ---------------------------------------------------------------------------


def _run_single_particle_dna(
    particle_index: int,
    context: DnaReplayContext,
    scenarios: list[Scenario],
    *,
    stats_filter: str,
) -> ParticleResult:
    """Calcule une seule particule ADN (équivalent DNA de _run_single_particle).

    Fonction top-level (picklable), appelée par chaque worker du
    ProcessPoolExecutor.

    La seed utilisée est dérivée de particle_index, garantissant un
    tirage distinct et reproductible par particule (même particle_index
    -> même résultat, peu importe l'ordre d'exécution des workers).

    IMPORTANT : seed = particle_index + 1, jamais particle_index seul.
    msprime.sim_ancestry rejette explicitement seed=0 (ValueError "seeds
    must be greater than 0 and less than 2^32") -- vérifié empiriquement.
    Donc particle_index=0 (le cas le plus probable, première particule)
    utilise seed=1, pas seed=0.

    Args:
        particle_index: L'index de la particule (0-based).
        context: Le contexte de rejeu pour les séquences ADN.
        scenarios: Les scénarios candidats (chaque particule tire le
            sien).
        stats_filter: "ALL" ou "HEADER".

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    drawn_scenario = draw_scenario(scenarios, seed + _SCENARIO_DRAW_SEED_OFFSET)

    summary_statistics, parameter_values, group_priors_values_nested = (
        compute_summary_statistics_dna(
            context=context,
            scenario_index=drawn_scenario.index,
            seed=seed,
            stats_filter=stats_filter,
        )
    )

    group_priors_values = {
        column: group_priors_values_nested[group][prior]
        for column, group, prior in _group_prior_columns(context.header_text)
    }

    return ParticleResult(
        particle_index=particle_index,
        scenario_index=drawn_scenario.index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
        group_priors_values=group_priors_values,
    )


def run_reftable_simulation_dna(
    reference_directory: str | Path,
    scenarios: list[Scenario],
    *,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Produit nrec particules ADN (lignes de reftable.bin) en parallèle.

    N'écrit rien sur disque par particule (compute_summary_statistics_dna
    est 100% Python, en mémoire).

    Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
    nrec-1), pas dans l'ordre de complétion des workers -- important
    pour la reproductibilité de l'ordre des lignes du reftable final.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        scenarios: La liste des scénarios candidats (typiquement TOUS
            les scénarios déclarés dans header.txt) : chaque particule
            tire le SIEN au hasard, pondéré par son `weight` (voir
            parameter_sampling.draw_scenario, sémantique vérifiée
            contre particuleC.cpp::ParticleC::drawscenario) -- une même
            particule peut donc finir sur n'importe lequel des
            scénarios de la liste, pas forcément le même pour toutes.
        nrec: Le nombre de particules à produire.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle (défaut : laissé
            à ProcessPoolExecutor, généralement le nombre de cœurs
            disponibles).

    Returns:
        La liste des ParticleResult, dans l'ordre de particle_index (0
        à nrec-1).
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename
    list_loci = parse_loci_description(header_text)
    sequences_observed = observed_sequences(mss_path, list_loci)

    context = DnaReplayContext(
        header_text=header_text,
        mss_path=mss_path,
        list_loci=list_loci,
        dna_observed=sequences_observed,
        frequencies_per_locus=base_frequency_by_locus(sequences_observed),
        samples_default=observed_count_sample(mss_path),
        sex_ratio=parse_sex_ratio(mss_path),
    )

    results_by_index: dict[int, ParticleResult] = {}
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_dna,
                particle_index,
                context,
                scenarios,
                stats_filter=stats_filter,
            ): particle_index
            for particle_index in range(nrec)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 100 == 0 or done == nrec:
                print(f"  {done}/{nrec} particules ADN simulées")

    return [results_by_index[i] for i in range(nrec)]


# Rejeu des tirages réels de DIYABC pour les séquences ADN (comparaison appariée)
def _group_prior_columns(header_text: str) -> list[tuple[str, str, str]]:
    """Liste ordonnée de triplets (nom_colonne, nom de groupe, nom du paramètre)
    des noms de colonnes "priors de groupe" d'un vrai
    reftable DIYABC.

    Ex: `µseq_2`, `k1seq_2`, juste après les paramètres historiques et
    avant les colonnes de statistiques sur chaque ligne -- vérifiée
    caractère pour caractère contre la vraie sortie DIYABC de
    `toy_example2_ms_dna` (`µmic_1 pmic_1 snimic_1 µseq_2 k1seq_2
    µseq_3 k1seq_3`, `nparamut=7`).

    Contrairement aux paramètres historiques (`_kept_param_names_by_
    scenario`), ces colonnes NE DÉPENDENT PAS du scénario tiré par la
    particule -- toujours les mêmes colonnes, dans l'ordre de
    déclaration des groupes (`G1`, `G2`, `G3`...) du header, un `mus_rate`
    toujours en premier dans chaque groupe.

    Pour un groupe `[S]` (séquence ADN) : `mus_rate` puis `k1`/`k2` SI ET
    SEULEMENT SI utilisés par le modèle du groupe (`get_parameter_used_
    by_model` -- ex: K2P/HKY n'ont que `k1`, pas de colonne `k2seq_N`).
    Pour un groupe `[M]` (microsat) : `mus_rate` puis `P` puis `SNI`,
    TOUJOURS les trois -- hypothèse basée sur ce qui est observé dans la
    vraie sortie de ce dataset précis, pas sur une règle générale
    vérifiée dans le C++ (MicroSat n'a pas de simulation-side code dans
    ce projet, seulement besoin de savoir COMBIEN de colonnes sauter
    pour atteindre celles des groupes ADN qui suivent).

    Args:
        header_text: Texte complet de header.txt.

    Returns:
            Une liste de tuples (nom_colonne, nom_groupe, nom_paramètre).
    """
    group_priors = parse_group_priors(header_text)
    columns = []
    for group_name, entries in group_priors.items():
        ms_or_seq = entries[0].ms_or_seq
        type_suffix = "seq" if ms_or_seq == "S" else "mic"
        group_number = group_name[1:]  # "G2" -> "2"

        columns.append(
            (f"µ{type_suffix}_{group_number}", group_name, "MEANMU")
        )  # mus_rate

        if ms_or_seq == "S":
            model_entry = next(e for e in entries if e.model)
            k1_used, k2_used = get_parameter_used_by_model(model_entry)
            if k1_used:
                columns.append(
                    (f"k1{type_suffix}_{group_number}", group_name, "MEANK1")
                )
            if k2_used:
                columns.append(
                    (f"k2{type_suffix}_{group_number}", group_name, "MEANK2")
                )
        else:
            columns.append((f"p{type_suffix}_{group_number}", group_name, "MEANP"))
            columns.append((f"sni{type_suffix}_{group_number}", group_name, "MEANSNI"))

    return columns


def group_prior_column_names(header_text: str) -> list[str]:
    """
    Liste des noms sortant dans le reftable DIYABC pour les priors de groupe, dans l'ordre de la ligne d'en-tête du fichier réel.
    Ex: `µseq_2`, `k1seq_2`, juste après les paramètres historiques et avant les colonnes de statistiques sur chaque ligne -- vérifiée caractère pour caractère contre la

    Args:
        header_text: Texte complet de header.txt.

    Returns:
        Une liste de noms de colonnes (str) pour les priors de groupe, dans l'ordre exact de la ligne d'en-tête du reftable DIYABC.
    """
    return [column for column, _, _ in _group_prior_columns(header_text)]


def parse_real_reftable_params_with_group_priors(
    path: str | Path, priors: list, scenarios: list[Scenario], group_priors_names: list
) -> list[tuple[int, dict[str, float], dict[str, float]]]:
    """Variante de parse_real_reftable_params qui lit AUSSI les priors de
    groupe.

    Colonnes `µseq_2`, `k1seq_2`... d'un vrai reftable DIYABC -- une
    fonction séparée plutôt qu'une extension en place, pour ne rien
    risquer sur parse_real_reftable_params et ses appelants SNP déjà
    validés (même choix que _run_single_particle/_run_single_
    particle_from_values : deux fonctions distinctes plutôt qu'une seule
    avec des branches conditionnelles).

    Sur chaque ligne, les colonnes de priors de groupe suivent
    IMMÉDIATEMENT les colonnes de paramètres historiques (elles-mêmes
    variables en nombre selon le scénario tiré par cette ligne -- voir
    priors_kept_by_scenario) et précèdent les colonnes de statistiques.
    L'offset des priors de groupe dépend donc du nombre de paramètres
    historiques de CETTE ligne, et leur ordre suit la ligne d'en-tête du
    fichier, pas la déclaration des priors -- même règle que
    parse_real_reftable_params (voir sa docstring et
    _historical_columns_order).

    Contrairement aux paramètres historiques, les priors de groupe NE
    SONT JAMAIS filtrées par scénario : `group_priors_names` (voir
    group_prior_column_names) est utilisé TEL QUEL pour toutes les
    lignes, quel que soit leur scénario -- vérifié empiriquement sur le
    vrai reftable de toy_example2_ms_dna (`nparamut=7`, un compte
    constant, contrairement à `nparam` qui varie par scénario). Appeler
    `_kept_param_names_by_scenario(group_priors_names, scenarios)` ici
    serait doublement faux : elle attend des objets `Prior` (pas des
    strings -- `is_constant_prior` plante sur `.min`/`.max`), et son filtre
    `get_parameter_names_used_by_scenario` ne reconnaît de toute façon
    que des noms de paramètres historiques, jamais de priors de groupe.

    Args:
        path: Chemin du reftable réel (format texte).
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.
        group_priors_names: Les noms de colonnes de priors de groupe
            (voir group_prior_column_names).

    Returns:
        Une liste de triplets (scenario_index, priors_values,
        group_priors_values) -- les deux dicts de valeurs restent
        SÉPARÉS (pas fusionnés en un seul) car ils alimentent deux
        étapes différentes en aval : priors_values sert à construire la
        démographie, group_priors_values au modèle de mutation ADN
        (k1/k2/mus_rate).
    """
    priors_kept_by_scenario = _kept_param_names_by_scenario(priors, scenarios)

    lines = [line for line in Path(path).read_text().splitlines() if line.strip()]
    data_lines = lines[1:]
    hist_order = _historical_columns_order(lines[0], priors, scenarios)
    rows = []
    for line in data_lines:
        tokens = line.split()
        scenario_index = int(tokens[0])
        kept = set(priors_kept_by_scenario[scenario_index])
        priors_param_names = [name for name in hist_order if name in kept]
        priors_values = {
            name: float(tokens[1 + i]) for i, name in enumerate(priors_param_names)
        }
        group_priors_values = {
            name: float(tokens[1 + len(priors_param_names) + i])
            for i, name in enumerate(group_priors_names)
        }
        rows.append((scenario_index, priors_values, group_priors_values))
    return rows


def _run_single_particle_dna_from_values(
    particle_index: int,
    context: DnaReplayContext,
    scenario_index: int,
    values: dict[str, float],
    group_priors_values: dict[str, float],
    *,
    stats_filter: str,
) -> ParticleResult:
    """Variante de _run_single_particle_dna qui NE TIRE AUCUN paramètre.

    Rejoue (scenario_index, values, group_priors_values) tels que
    fournis -- typiquement issus de
    parse_real_reftable_params_with_group_priors.

    Args:
        particle_index: L'index de la particule (0-based).
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        scenario_index: L'index 1-based du scénario déjà tiré par
            DIYABC pour cette particule.
        values: Les valeurs de paramètres historiques déjà connues,
            {nom: valeur}.
        group_priors_values: Les valeurs de priors de groupe déjà
            connues, {nom_colonne: valeur} (voir
            compute_summary_statistics_dna_from_values).
        stats_filter: "ALL" ou "HEADER".

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    summary_statistics = compute_summary_statistics_dna_from_values(
        context=context,
        scenario_index=scenario_index,
        values=values,
        group_priors_values=group_priors_values,
        seed=seed,
        stats_filter=stats_filter,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=scenario_index,
        parameter_values=values,
        summary_statistics=summary_statistics,
        group_priors_values=group_priors_values,
    )


def replay_reftable_simulation_dna(
    reference_directory: str | Path,
    priors: list,
    group_priors_names: list[str],
    scenarios: list[Scenario],
    real_reftable_path: str | Path,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Rejoue, particule par particule, les tirages RÉELS de DIYABC (équivalent
    ADN de replay_reftable_simulation).

    Lit un reftable réel existant (scénario, paramètres historiques ET
    priors de groupe RÉELLEMENT tirés par DIYABC) et rejoue chaque
    particule côté msprime avec EXACTEMENT les mêmes valeurs -- permet
    une comparaison appariée ligne à ligne, pas seulement une
    comparaison de distributions agrégées.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        priors: Les priors historiques déclarés dans header.txt.
        group_priors_names: Les noms de colonnes de priors de groupe
            (voir group_prior_column_names).
        scenarios: Les scénarios candidats.
        real_reftable_path: Chemin du reftable réel à rejouer.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle.

    Returns:
        Les ParticleResult dans le MÊME ORDRE que les lignes du fichier
        réel.
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename
    list_loci = parse_loci_description(header_text)
    sequences_observed = observed_sequences(mss_path, list_loci)

    context = DnaReplayContext(
        header_text=header_text,
        mss_path=mss_path,
        list_loci=list_loci,
        dna_observed=sequences_observed,
        frequencies_per_locus=base_frequency_by_locus(sequences_observed),
        samples_default=observed_count_sample(mss_path),
        sex_ratio=parse_sex_ratio(mss_path),
    )

    # On lit les sorties de diyabc (scénario tiré + valeurs de paramètres RÉELLEMENT tirées) pour
    # les rejouer ensuite côté msprime, afin de comparer les deux simulateurs sur EXACTEMENT
    # les mêmes tirages de priors.

    check_real_reftable_matches_header(header_text, real_reftable_path)
    rows = parse_real_reftable_params_with_group_priors(
        path=real_reftable_path,
        priors=priors,
        scenarios=scenarios,
        group_priors_names=group_priors_names,
    )

    results_by_index: dict[int, ParticleResult] = {}
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_dna_from_values,
                particle_index,
                context,
                scenario_index,
                values,
                group_priors_values,
                stats_filter=stats_filter,
            ): particle_index
            for particle_index, (
                scenario_index,
                values,
                group_priors_values,
            ) in enumerate(rows)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 100 == 0 or done == len(rows):
                print(
                    f"Rejeu des tirages réels : {done}/{len(rows)} particules terminées"
                )
        return [results_by_index[i] for i in range(len(rows))]


# --------------------------------------------------------------------------
# Pour les microsatellites : lecture, écriture, rejeux de tirages réels
# --------------------------------------------------------------------------


def _run_single_particle_microsat(
    particle_index: int,
    context: MicrosatReplayContext,
    scenarios: list[Scenario],
    *,
    stats_filter: str,
) -> ParticleResult:
    """Calcule une seule particule microsat (équivalent microsat de _run_single_particle).

    Fonction top-level (picklable), appelée par chaque worker du
    ProcessPoolExecutor.

    La seed utilisée est dérivée de particle_index, garantissant un
    tirage distinct et reproductible par particule (même particle_index
    -> même résultat, peu importe l'ordre d'exécution des workers).

    IMPORTANT : seed = particle_index + 1, jamais particle_index seul.
    msprime.sim_ancestry rejette explicitement seed=0 (ValueError "seeds
    must be greater than 0 and less than 2^32") -- vérifié empiriquement.
    Donc particle_index=0 (le cas le plus probable, première particule)
    utilise seed=1, pas seed=0.

    Args:
        particle_index: L'index de la particule (0-based).
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        scenarios: Les scénarios candidats (chaque particule tire le
            sien).
        stats_filter: "ALL" ou "HEADER".

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    drawn_scenario = draw_scenario(scenarios, seed + _SCENARIO_DRAW_SEED_OFFSET)

    summary_statistics, parameter_values, group_priors_values_nested = (
        compute_summary_statistics_microsat(
            context=context,
            scenario_index=drawn_scenario.index,
            seed=seed,
            stats_filter=stats_filter,
        )
    )

    group_priors_values = {
        column: group_priors_values_nested[group][prior]
        for column, group, prior in _group_prior_columns(context.header_text)
    }
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=drawn_scenario.index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
        group_priors_values=group_priors_values,
    )


def run_reftable_simulation_microsat(
    reference_directory: str | Path,
    scenarios: list[Scenario],
    *,
    nrec: int,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Produit nrec particules microsat (lignes de reftable.bin) en parallèle.

    N'écrit rien sur disque par particule (compute_summary_statistics_microsat
    est 100% Python, en mémoire).

    Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
    nrec-1), pas dans l'ordre de complétion des workers -- important
    pour la reproductibilité de l'ordre des lignes du reftable final.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        scenarios: La liste des scénarios candidats (typiquement TOUS
            les scénarios déclarés dans header.txt) : chaque particule
            tire le SIEN au hasard, pondéré par son `weight` (voir
            parameter_sampling.draw_scenario, sémantique vérifiée
            contre particuleC.cpp::ParticleC::drawscenario) -- une même
            particule peut donc finir sur n'importe lequel des
            scénarios de la liste, pas forcément le même pour toutes.
        nrec: Le nombre de particules à produire.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle (défaut : laissé
            à ProcessPoolExecutor, généralement le nombre de cœurs
            disponibles).

    Returns:
        La liste des ParticleResult, dans l'ordre de particle_index (0
        à nrec-1).
    """
    reference_directory = Path(reference_directory)
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename
    list_loci = parse_loci_description(header_text)
    microsat_observed = observed_microsatellites(mss_path, list_loci)

    context = MicrosatReplayContext(
        header_text=header_text,
        mss_path=mss_path,
        list_loci=list_loci,
        microsat_observed=microsat_observed,
        bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
        samples_default=observed_count_sample(mss_path),
        sex_ratio=parse_sex_ratio(mss_path),
    )

    results_by_index: dict[int, ParticleResult] = {}
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_microsat,
                particle_index,
                context,
                scenarios,
                stats_filter=stats_filter,
            ): particle_index
            for particle_index in range(nrec)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 10 == 0 or done == nrec:
                print(f"  {done}/{nrec} particules microsat simulées")

    return [results_by_index[i] for i in range(nrec)]


# Rejeu des tirages réels de DIYABC pour les microsatellites (comparaison appariée)


def _run_single_particle_microsat_from_values(
    particle_index: int,
    context: MicrosatReplayContext,
    scenario_index: int,
    values: dict[str, float],
    group_priors_values: dict[str, float],
    *,
    stats_filter: str,
) -> ParticleResult:
    """Variante de _run_single_particle_microsat qui NE TIRE AUCUN paramètre.

    Rejoue (scenario_index, values, group_priors_values) tels que
    fournis -- typiquement issus de
    parse_real_reftable_params et
    compute_summary_statistics_microsat_from_values.

    Args:
        particle_index: L'index de la particule (0-based).
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        scenario_index: L'index 1-based du scénario déjà tiré par
            DIYABC pour cette particule.
        values: Les valeurs de paramètres historiques déjà connues,
            {nom: valeur}.
        group_priors_values: Les valeurs de priors de groupe déjà
            connues, {nom_colonne: valeur} (voir
            compute_summary_statistics_dna_from_values).
        stats_filter: "ALL" ou "HEADER".

    Returns:
        Le ParticleResult de cette particule.
    """
    seed = particle_index + 1
    summary_statistics = compute_summary_statistics_microsat_from_values(
        context=context,
        scenario_index=scenario_index,
        values=values,
        group_priors_values=group_priors_values,
        seed=seed,
        stats_filter=stats_filter,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=scenario_index,
        parameter_values=values,
        summary_statistics=summary_statistics,
        group_priors_values=group_priors_values,
    )


def replay_reftable_simulation_microsat(
    reference_directory: str | Path,
    priors: list,
    group_priors_names: list[str],
    scenarios: list[Scenario],
    real_reftable_path: str | Path,
    stats_filter: str = "ALL",
    max_workers: int | None = None,
) -> list[ParticleResult]:
    """Rejoue, particule par particule, les tirages RÉELS de DIYABC (équivalent
    microsat de replay_reftable_simulation).

    Lit un reftable réel existant (scénario, paramètres historiques ET
    priors de groupe RÉELLEMENT tirés par DIYABC) et rejoue chaque
    particule côté msprime avec EXACTEMENT les mêmes valeurs -- permet
    une comparaison appariée ligne à ligne, pas seulement une
    comparaison de distributions agrégées.

    Args:
        reference_directory: Le dossier contenant header.txt et le
            fichier .mss observé.
        priors: Les priors historiques déclarés dans header.txt.
        group_priors_names: Les noms de colonnes de priors de groupe
            (voir group_prior_column_names).
        scenarios: Les scénarios candidats.
        real_reftable_path: Chemin du reftable réel à rejouer.
        stats_filter: "ALL" ou "HEADER".
        max_workers: Le nombre de process en parallèle.

    Returns:
        Les ParticleResult dans le MÊME ORDRE que les lignes du fichier
        réel.
    """
    reference_directory = Path(
        reference_directory
    )  # Normalise en path au cas où une str soit passée
    header_text = read_header_text(reference_directory)
    raise_if_serial_with_sex_linked_loci(header_text)
    check_header_trailer_line(header_text)
    mss_filename = header_text.splitlines()[0].strip()
    mss_path = reference_directory / mss_filename
    list_loci = parse_loci_description(header_text)
    microsat_observed = observed_microsatellites(mss_path, list_loci)

    context = MicrosatReplayContext(
        header_text=header_text,
        mss_path=mss_path,
        list_loci=list_loci,
        microsat_observed=microsat_observed,
        bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
        samples_default=observed_count_sample(mss_path),
        sex_ratio=parse_sex_ratio(mss_path),
    )

    print(
        "Lecture et extraction des informations du header.txt et du fichier .mss terminée, lecture des valeurs de diyabc."
    )

    # On lit les sorties de diyabc (scénario tiré + valeurs de paramètres RÉELLEMENT tirées) pour
    # les rejouer ensuite côté msprime, afin de comparer les deux simulateurs sur EXACTEMENT
    # les mêmes tirages de priors.

    check_real_reftable_matches_header(header_text, real_reftable_path)
    rows = parse_real_reftable_params_with_group_priors(
        path=real_reftable_path,
        priors=priors,
        scenarios=scenarios,
        group_priors_names=group_priors_names,
    )

    print("Lecture des tirages réels de DIYABC terminée, lancement du rejeu msprime...")

    results_by_index: dict[int, ParticleResult] = {}
    done = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_microsat_from_values,
                particle_index,
                context,
                scenario_index,
                values,
                group_priors_values,
                stats_filter=stats_filter,
            ): particle_index
            for particle_index, (
                scenario_index,
                values,
                group_priors_values,
            ) in enumerate(rows)
        }

        for future in as_completed(futures):
            particle_index = futures[future]
            results_by_index[particle_index] = future.result()
            done += 1
            if done % 100 == 0 or done == len(rows):
                print(
                    f"Rejeu des tirages réels : {done}/{len(rows)} particules terminées"
                )
    return [results_by_index[i] for i in range(len(rows))]
