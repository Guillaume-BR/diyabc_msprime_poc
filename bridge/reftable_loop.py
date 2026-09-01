"""
Boucle d'itération produisant les nrec "particules" (lignes) d'un futur
reftable.bin : pour chaque particule, un tirage de scénario et de
paramètres distinct, et une simulation msprime complète (calcul des
statistiques résumées via compute_summary_statistics, 100% Python --
plus de subprocess ni de fichier intermédiaire sur disque).

Parallélisé via ProcessPoolExecutor : chaque particule est indépendante
des autres (son propre tirage, sa propre simulation), donc embarrassingly
parallel.
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

import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from bridge.ancestry_simulation import prepare_poolseq_observed_reads
from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.header_dataclasses import Scenario
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import detect_snp_file_type
from bridge.parameter_sampling import draw_scenario
from bridge.pipeline import (
    compute_summary_statistics,
    compute_summary_statistics_dna,
    compute_summary_statistics_dna_from_values,
    compute_summary_statistics_from_values,
    read_header_text,
)
from bridge.prior_parser import (
    get_parameter_used_by_model,
    is_constant_prior,
    parse_group_priors,
    parse_priors,
)
from bridge.scenario_parser import parse_header_scenarios

# Décalage appliqué à la seed de particule avant de tirer le scénario
# (draw_scenario), pour ne JAMAIS partager la même seed brute avec le
# tirage des paramètres (draw_parameter_values, appelé plus loin dans
# compute_summary_statistics avec seed=seed, sans offset) -- sinon les
# deux tirages, bien qu'indépendants dans l'intention, consomment le
# MÊME premier random.random() sous-jacent (chaque fonction fait son
# propre random.Random(seed) frais), ce qui corrèle artificiellement le
# scénario tiré et la valeur du premier prior déclaré (ex: N1) :
# vérifié empiriquement sur toy_example5 -- scénario 1 (ra petit) donnait
# systématiquement un N1 trop bas, scénario 3 (ra grand) un N1 trop haut.
# Au-delà de LOCUS_TYPE_SEED_OFFSET (pipeline.py, 0..40_000_000) pour ne
# pas non plus recréer une collision avec ce décalage-là.
_SCENARIO_DRAW_SEED_OFFSET = 50_000_000


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


# ----------------------------------------------------------------------------
# Pour les fichiers .snp DIYABC : lecture, écriture, rejeux de tirages réels
# ----------------------------------------------------------------------------

# ── Tirage indépendant de scénario + paramètres, par particule ────────────


def _run_single_particle(
    particle_index: int,
    reference_directory: Path,
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
        reference_directory: Le dossier contenant header.txt et le
            fichier .snp observé.
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
        reference_directory=reference_directory,
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

    results_by_index: dict[int, ParticleResult] = {}

    header_text = read_header_text(reference_directory)
    snp_path = reference_directory / header_text.splitlines()[0].strip()
    observed_reads_per_locus = None
    if detect_snp_file_type(snp_path) == "POOL":
        total_loci_poolseq = parse_loci_description(
            header_text
        ).loci_counts_by_heritage["A"]
        observed_reads_per_locus = prepare_poolseq_observed_reads(
            snp_path, total_loci_poolseq
        )
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle,
                particle_index,
                reference_directory,
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
    nparam = [len(kept_param_names_by_scenario[s.index]) for s in scenarios]
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
    stat_names = list(results[0].summary_statistics.keys())

    def _centre(s: str, width: int = 14) -> str:
        """Reproduit centre() de mesutils.cpp : centrage sur width chars."""
        return s.center(width)

    with open(output_path, "w", encoding="utf-8") as f:
        # En-tête : "scenario" + noms des paramètres + noms des stats
        header = (
            _centre("scenario")
            + "".join(_centre(n) for n in all_param_names)
            + "".join(_centre(n) for n in stat_names)
        )
        f.write(header + "\n")

        # Une ligne par particule
        for r in results:
            line = f"{r.scenario_index:3d}  "
            for name in all_param_names:
                line += f"  {r.parameter_values.get(name, float('nan')):12.6f}"
            for name in stat_names:
                line += f"  {r.summary_statistics[name]:12.6f}"
            f.write(line + "\n")


def rewrite_real_reftable_txt(
    input_path: str | Path,
    output_path: str | Path,
    priors: list,
    scenarios: list[Scenario],
) -> None:
    """Réécrit un reftable RÉEL de DIYABC en un texte à colonnes de largeur FIXE.

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

    Le fichier réécrit a EXACTEMENT le même format que celui produit par
    write_reftable_txt (mêmes colonnes de paramètres -- union dans
    l'ordre de déclaration des priors --, dans le même ordre), donc
    directement comparable colonne à colonne avec un reftable_msprime
    généré par run_reftable_simulation/replay_reftable_simulation.

    Args:
        input_path: Chemin du reftable réel brut (format texte).
        output_path: Chemin où écrire le fichier réécrit.
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.
    """
    kept_by_scenario = _kept_param_names_by_scenario(priors, scenarios)
    used_by_any = {name for names in kept_by_scenario.values() for name in names}
    all_param_names = [p.name for p in priors if p.name in used_by_any]

    lines = [line for line in Path(input_path).read_text().splitlines() if line.strip()]
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
            param_names = kept_by_scenario[scenario_index]
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
    pas un nom fixe -- voir pipeline.run_poc_for_directory).

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
    test_directory = Path(test_directory)
    header_text = read_header_text(test_directory)

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

    On lit donc, pour CHAQUE ligne, le nombre de tokens de paramètres
    correspondant SPÉCIFIQUEMENT au scénario de cette ligne
    (kept_by_scenario[scenario_index], même filtre non-constant +
    utilisé-par-ce-scénario que write_reftable_txt/write_reftable_bin),
    pas une union appliquée uniformément -- ce qui gère aussi, en
    particulier, le cas single-scénario où certains priors sont devenus
    constants (is_constant_prior) et donc absents des colonnes de
    sortie.

    Args:
        path: Chemin du reftable réel (format texte).
        priors: Les priors déclarés dans header.txt.
        scenarios: Les scénarios candidats.

    Returns:
        Une liste de tuples (scenario_index, {nom_paramètre: valeur}),
        un par ligne du reftable (dans l'ordre du fichier).
    """
    priors_kept_by_scenario = _kept_param_names_by_scenario(priors, scenarios)

    lines = [line for line in Path(path).read_text().splitlines() if line.strip()]
    data_lines = lines[1:]  # ligne 0 = en-tête

    rows = []
    for line in data_lines:
        tokens = line.split()
        scenario_index = int(tokens[0])
        param_names = priors_kept_by_scenario[scenario_index]
        values = {name: float(tokens[1 + i]) for i, name in enumerate(param_names)}
        rows.append((scenario_index, values))
    return rows


def _run_single_particle_from_values(
    particle_index: int,
    reference_directory: Path,
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
        reference_directory=reference_directory,
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
    """Rejoue, particule par particule, les tirages de paramètres RÉELLEMENT effectués par DIYABC.

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

    # On lit les sorties de diyabc (scénario tiré + valeurs de paramètres RÉELLEMENT tirées) pour
    # les rejouer ensuite côté msprime, afin de comparer les deux simulateurs sur EXACTEMENT
    # les mêmes tirages de priors.
    rows = parse_real_reftable_params(real_reftable_path, priors, scenarios)

    results_by_index: dict[int, ParticleResult] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_from_values,
                particle_index,
                reference_directory,
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
    return [results_by_index[i] for i in range(len(rows))]


# ----------------------------------------------------------------------------
# Pour les séquences ADN : lecture, écriture, rejeux de tirages réels
# ---------------------------------------------------------------------------


def _run_single_particle_dna(
    particle_index: int,
    reference_directory: Path,
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

    summary_statistics, parameter_values = compute_summary_statistics_dna(
        reference_directory=reference_directory,
        scenario_index=drawn_scenario.index,
        seed=seed,
        stats_filter=stats_filter,
    )
    return ParticleResult(
        particle_index=particle_index,
        scenario_index=drawn_scenario.index,
        parameter_values=parameter_values,
        summary_statistics=summary_statistics,
    )


# Rejeu des tirages réels de DIYABC pour les séquences ADN (comparaison appariée)


def group_prior_column_names(header_text: str) -> list[str]:
    """Liste ordonnée des noms de colonnes "priors de groupe" d'un vrai reftable DIYABC.

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
        La liste ordonnée des noms de colonnes de priors de groupe.
    """
    group_priors = parse_group_priors(header_text)

    names = []
    for group_name, entries in group_priors.items():
        ms_or_seq = entries[0].ms_or_seq
        type_suffix = "seq" if ms_or_seq == "S" else "mic"
        group_number = group_name[1:]  # "G2" -> "2"

        names.append(f"µ{type_suffix}_{group_number}")  # mus_rate, toujours présent

        if ms_or_seq == "S":
            model_entry = next(e for e in entries if e.model)
            k1_used, k2_used = get_parameter_used_by_model(model_entry)
            if k1_used:
                names.append(f"k1{type_suffix}_{group_number}")
            if k2_used:
                names.append(f"k2{type_suffix}_{group_number}")
        else:
            names.append(f"p{type_suffix}_{group_number}")
            names.append(f"sni{type_suffix}_{group_number}")

    return names


def parse_real_reftable_params_with_group_priors(
    path: str | Path, priors: list, scenarios: list[Scenario], group_priors_names: list
) -> list[tuple[int, dict[str, float], dict[str, float]]]:
    """Variante de parse_real_reftable_params qui lit AUSSI les priors de groupe.

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
    data_lines = lines[1:]  # ligne 0 = en-tête

    rows = []
    for line in data_lines:
        tokens = line.split()
        scenario_index = int(tokens[0])
        priors_param_names = priors_kept_by_scenario[scenario_index]
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
    reference_directory: Path,
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
        reference_directory=reference_directory,
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
    """Rejoue, particule par particule, les tirages RÉELS de DIYABC (équivalent ADN de replay_reftable_simulation).

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

    # On lit les sorties de diyabc (scénario tiré + valeurs de paramètres RÉELLEMENT tirées) pour
    # les rejouer ensuite côté msprime, afin de comparer les deux simulateurs sur EXACTEMENT
    # les mêmes tirages de priors.

    rows = parse_real_reftable_params_with_group_priors(
        path=real_reftable_path,
        priors=priors,
        scenarios=scenarios,
        group_priors_names=group_priors_names,
    )

    results_by_index: dict[int, ParticleResult] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_single_particle_dna_from_values,
                particle_index,
                reference_directory,
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
    return [results_by_index[i] for i in range(len(rows))]
