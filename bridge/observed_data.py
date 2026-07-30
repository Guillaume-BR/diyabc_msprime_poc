"""
Lecture des données observées DIYABC (fichiers .snp, IndSeq -- format
individu par ligne -- ou PoolSeq -- format pool par population, un locus
par ligne, voir detect_snp_file_type) -- comptage du nombre d'individus
(IndSeq) ou taille haploïde des pools (PoolSeq) par population
(nécessaire pour savoir combien d'échantillons demander à
msprime.sim_ancestry()), lecture du sex-ratio et du sexe par individu
(nécessaires pour les loci <X>/<Y>/<M>, dont la ploïdie et le
coefficient de coalescence dépendent du sexe -- voir ParticleC::calploidy
et DataC::cal_coeffcoal), et lecture des seuils MAF ("<MAF=N%>" ou
"<MAF=hudson>", via parse_maf_ratio) et MRC ("<MRC=N>", via
parse_mrc_ratio, PoolSeq uniquement) consommés par with_maf_filter/
with_maf_filter_shared_ancestry/with_mrc_filter (ancestry_simulation.py)
pour décider s'il faut rejeter-et-resimuler un locus sous le seuil.

Référence : src-JMC-C++/data.cpp (détection du format "IND SEX POP" vs
"POOL", lecture du sex-ratio "<NM=xNF>", du MAF et du MRC en tête de
fichier).

Ce module ne lit PAS les génotypes IndSeq eux-mêmes : on simule des
données artificielles avec msprime, sans jamais réutiliser les
génotypes observés réels dans le pipeline de simulation. Ce n'est PLUS
vrai pour PoolSeq : `observed_reads` lit les comptages de lectures
RÉELLEMENT observés, et leur profondeur totale (`nreads_total`) est
réutilisée telle quelle comme paramètre `n` du tirage binomial des
lectures simulées (voir `simulate_poolseq_reads`,
`ancestry_simulation.py`) -- seule la répartition allèle1/allèle2 est
simulée, jamais la couverture elle-même.
"""

import re
from collections import Counter
from pathlib import Path

from bridge.scenario_types import LociDescriptionDetailed


def _find_header_index(lines: list[str]) -> int:
    """Repère l'index de la ligne d'en-tête 'IND SEX POP' ou 'POOL' parmi les deux
    premières lignes du fichier -- factorisé entre count_samples_per_population
    et individual_sexes_per_population, qui en ont toutes deux besoin.

    L'en-tête peut être précédé ou non d'un commentaire libre en première
    ligne (ex: '<NM=1NF> <MAF=hudson> ...', comportement observé dans
    data.cpp, qui teste les deux cas) : on recherche son index plutôt que
    de supposer sa position, pour ne perdre aucune ligne de données quel
    que soit le cas.

    Lève ValueError si l'en-tête n'est trouvé dans aucune des deux
    premières lignes.
    """

    header_index = next(
        (
            i
            for i in range(min(2, len(lines)))
            if lines[i].split()[:3] == ["IND", "SEX", "POP"]
            or lines[i].split()[:1] == ["POOL"]
        ),
        None,
    )
    if header_index is None:
        raise ValueError(
            f"En-tête 'IND SEX POP' ou 'POOL' non trouvé dans les deux premières "
            f"lignes. Lignes lues : {lines[:2]!r}"
        )
    return header_index


def detect_snp_file_type(snp_file_path: str | Path) -> str:
    """Détecte le type de fichier .snp DIYABC : "INDSEQ" (individus par
    ligne) ou "POOLSEQ" (pools par ligne), à l'aide du header_index
    trouvé par _find_header_index, en lisant la première ligne
    non vide du fichier contenant "IND" ou "POOL".
    """
    lines = Path(snp_file_path).read_text().splitlines()
    header_index = _find_header_index(lines)
    if lines[header_index].split()[0] == "POOL":
        return "POOL"
    elif lines[header_index].split()[0] == "IND":
        return "IND"


def _parse_pool_header_line(lines: list[str], header_index: int) -> dict[str, int]:
    """Parse la ligne d'en-tête POOL du fichier .snp DIYABC, au format
    'POOL POP_NAME:HAPLOID_SAMPLE_SIZE  POP1:200 POP2:200 POP3:200 POP4:200',
    et retourne un dictionnaire {nom_population: taille_haploïde}.

    Ex: pour toy_example4 -> {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}
    """
    header_line = lines[header_index]
    first_pop = header_line.find("POP")
    if first_pop == -1:
        raise ValueError(
            f"En-tête 'POOL' trouvé mais aucune population déclarée dans la ligne : "
            f"{header_line!r}"
        )
    second_pop = header_line.find("POP", first_pop + 1)
    if second_pop == -1:
        raise ValueError(
            f"En-tête 'POOL' trouvé mais aucune population déclarée dans la ligne : "
            f"{header_line!r}"
        )
    counts_by_population = {}
    for part in header_line[second_pop:].split():
        if part.startswith("POP"):
            pop_name, haploid_sample_size = part.split(":")
            counts_by_population[pop_name] = int(haploid_sample_size)
    return counts_by_population


def count_samples_per_population(snp_file_path: str | Path) -> dict[str, int]:
    """Compte le nombre d'individus par population dans un fichier .snp
    DIYABC au format 'IND SEX POP <génotypes...>'n ou 'POOL POP_NAME:HAPLOID_SAMPLE_SIZE'.

    Ex: pour human -> {"ASW": 30, "YRI": 30, ...} ou pour toy_example4 -> {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}.

    IMPORTANT -- garantie d'ordre : le dict retourné préserve l'ordre de
    première apparition des populations dans le fichier (garanti par
    Counter/dict en Python >= 3.7, et vérifié expérimentalement sur
    human). Cet ordre a un sens métier précis : header.txt ne nomme jamais
    les populations (seulement des indices 1,2,3,4) -- le mapping réel,
    vérifié en l'absence de toute référence croisée dans le code C++
    (data.cpp ne relie jamais popname aux indices de scénario), est
    implicite : pop i du scénario = i-ème population dans l'ORDRE
    D'APPARITION de ce fichier. Ne jamais remplacer Counter par un type
    qui ne garantirait pas cet ordre (ex: trier les clés alphabétiquement
    casserait silencieusement ce mapping).
    """
    path = Path(snp_file_path)
    lines = path.read_text().splitlines()
    header_index = _find_header_index(lines)
    type_of_file = detect_snp_file_type(snp_file_path)
    counts = Counter()
    if type_of_file == "IND":
        pop_index = lines[header_index].split().index("POP")

        counts = Counter()
        for line in lines[header_index + 1 :]:
            fields = line.split()
            if not fields:
                continue
            counts[fields[pop_index]] += 1
    elif type_of_file == "POOL":
        counts.update(_parse_pool_header_line(lines, header_index))
    return dict(counts)


def individual_sexes_per_population(
    snp_file_path: str | Path,
) -> dict[str, list[str]]:
    """Lit le sexe de chaque individu, regroupé par population, dans un
    fichier .snp DIYABC au format 'IND SEX POP <génotypes...>'.

    Ex: pour toy_example5 -> {"pop1": ["M", "F", "F", ...], ...}

    Valeurs telles quelles côté DIYABC (data.cpp:702-704) : "M", "F", ou
    "9" (sexe inconnu -- cas de human_snp_all22chr_maf5.snp, où les 120
    individus sont tous "9" puisque le dataset est <A>-only et ne
    renseigne pas le sexe réel). Pas de normalisation en booléen ici :
    c'est à l'appelant (le futur équivalent Python de calploidy) de
    décider quoi faire du cas "9", typiquement lever une erreur si un
    locus <X>/<Y>/<M> est demandé sur des données non sexées.

    Même garantie d'ordre que count_samples_per_population : listes dans
    l'ordre d'apparition des individus dans le fichier, par population
    dans l'ordre de première apparition.
    """
    path = Path(snp_file_path)
    lines = path.read_text().splitlines()
    header_index = _find_header_index(lines)
    fields_header = lines[header_index].split()
    sex_index = fields_header.index("SEX")
    pop_index = fields_header.index("POP")

    sexes_by_population: dict[str, list[str]] = {}
    for line in lines[header_index + 1 :]:
        fields = line.split()
        if not fields:
            continue
        sexes_by_population.setdefault(fields[pop_index], []).append(fields[sex_index])

    return sexes_by_population


def parse_sex_ratio(snp_file_path: str | Path) -> float:
    """Lit le sex-ratio déclaré en tête de fichier .snp DIYABC, au format
    '<NM=xNF> ...' où x = nombre de mâles / nombre de femelles (ex:
    '<NM=0.428571NF>' pour toy_example5).

    Retourne sexratio = x / (1+x) (fraction de mâles), reproduisant
    exactement DataC::readfile (data.cpp:475-486). Comme dans le C++,
    retombe sur 0.5 si le token '<NM=' est absent de la première ligne --
    et, contrairement à _find_header_index, ne cherche JAMAIS que la
    ligne 0 : le C++ ne considère que la toute première ligne du fichier,
    jamais 'IND SEX POP' elle-même.
    """
    path = Path(snp_file_path)
    with path.open() as f:
        first_line = f.readline()

    start = first_line.find("<NM=")
    if start == -1:
        return 0.5

    end = first_line.find("NF>", start + 4)
    ratio_male_to_female = float(first_line[start + 4 : end])
    return ratio_male_to_female / (1.0 + ratio_male_to_female)


def parse_maf_ratio(snp_file_path: str | Path) -> float:
    """Lit le MAF déclaré en tête de fichier .snp DIYABC, au format
    '<MAF=xxx> ...' où xxx = MAF (ex: '<MAF=hudson>' ou '<MAF=0.05>').

    Retourne le seuil MAF, reproduisant exactement DataC::readfile
    (data.cpp:475-497) : 0.0 si le token '<MAF=' est absent de la première
    ligne, ou si xxx n'est pas numérique (ex: 'hudson'), comme le ferait
    atof() en C++ -- 0.0 veut dire "pas de filtre" (algorithme de Hudson
    standard), jamais distingué du cas "MAF=0%" explicite, exactement comme
    côté C++. Contrairement à _find_header_index, ne cherche JAMAIS que la
    ligne 0 : le C++ ne considère que la toute première ligne du fichier.
    """

    path = Path(snp_file_path)
    with path.open() as f:
        first_line = f.readline()

    start = first_line.find("<MAF=")
    if start == -1:
        return 0.0

    end = first_line.find(">", start + 5)
    maf_value = first_line[start + 5 : end]
    try:
        return float(maf_value)
    except ValueError:
        return 0.0


def parse_mrc_ratio(snp_file_path: str | Path) -> float:
    """Lit le MRC déclaré en tête de fichier .snp DIYABC, au format
    '<MRC=xxx> ...' où xxx = MRC ).

    Retourne le seuil MRC, reproduisant exactement DataC::readfile
    (data.cpp:498-508) : 1 si le token '<MRC=' est absent de la première
    ligne, ou si xxx n'est pas numérique (ex: 'hudson'), comme le ferait
    atof() en C++ -- 1 veut dire "pas de filtre"
    """

    path = Path(snp_file_path)
    with path.open() as f:
        first_line = f.readline()

    start = first_line.find("<MRC=")
    if start == -1:
        return 1

    end = first_line.find(">", start + 5)
    mrc_value = first_line[start + 5 : end]
    try:
        return float(mrc_value)
    except ValueError:
        return 1


def population_index_to_name(snp_file_path: str | Path) -> dict[int, str]:
    """Construit le mapping entre l'indice de population utilisé dans
    header.txt (1-indexed : pop1, pop2, ...) et le nom réel de population
    tel qu'il apparaît dans le fichier .snp (ex: "ASW", "YRI"...).

    Ex: {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"} pour human et {1: "POP1", 2: "POP2", 3: "POP3", 4: "POP4"} pour toy_example4.

    Voir la docstring de count_samples_per_population pour la
    justification de ce mapping par ordre d'apparition (header.txt ne
    nomme jamais les populations).
    """
    names_in_order = list(count_samples_per_population(snp_file_path).keys())
    return {i + 1: name for i, name in enumerate(names_in_order)}


def observed_mrc(reads_by_population: dict[str, tuple[int, int]]) -> float:
    """Calcule le MRC observé pour un locus donné, à partir des lectures
    observées par population (dictionnaire {nom_population: (nreads1,
    nreads1+nreads2)}).

    Retourne min(somme reads allèle1, somme reads allèle2) TOUTES
    populations combinées, reproduisant exactement DataC::purgelocMRCPOOLSEQ
    (data.cpp:1087-1093).
    """
    sum_derived = sum(derived for derived, _ in reads_by_population.values())
    sum_total = sum(total for _, total in reads_by_population.values())
    return min(sum_derived, sum_total - sum_derived) if sum_total > 0 else 0.0


def observed_reads(
    snp_file_path: str | Path, num_loci: int | None = None
) -> list[dict[str, tuple[int, int]]]:
    """Lit les lignes de génotypes du fichier .snp DIYABC POOLSEQ, en ignorant
    l'en-tête et les lignes vides. Retourne la liste des lignes de
    comptage de reads (nreads1, nreads1 + nreads2) par population, dans l'ordre d'apparition des populations
    dans le fichier. Chaque tuple contient le nombre de lectures pour l'allèle 1 et
    le nombre total de lectures (allèle 1 + allèle 2) pour cette population.
    On retourne une ligne de la forme POP1: (nreads1, nreads1 + nreads2), POP2: (nreads1, nreads1 + nreads2), ...

    Purge les loci sous le seuil MRC (`<MRC=N>`, via parse_mrc_ratio) --
    reproduit `DataC::purgelocMRCPOOLSEQ` (data.cpp), qui élimine ces
    loci de l'observé AU CHARGEMENT du fichier, avant toute utilisation
    (simulation ou calcul de statobs). Le critère est le même que
    `with_mrc_filter`/`mrcreached` : min(somme reads allèle1, somme
    reads allèle2), TOUTES populations combinées, doit être >= MRC.
    Sans cette purge, les loci quasi-monomorphes (très peu de lectures
    pour l'allèle minoritaire, souvent des erreurs de séquençage) restent
    inclus et faussent silencieusement toutes les statistiques en aval
    -- confirmé empiriquement le 22/07/2026 : 130/130 stats divergaient
    de >1% par rapport au vrai statobs.txt de DIYABC sans cette purge,
    0/130 avec.

    num_loci : si non None, limite le nombre de loci lus à ce nombre
    (utile pour toy_example4, où on ne veut que les 100 premiers loci
    passants le seuil MRC, pour reproduire exactement le statobs.txt de
    DIYABC). Si None, lit tous les loci du fichier.
    """
    if detect_snp_file_type(snp_file_path) != "POOL":
        raise ValueError(
            f"Le fichier {snp_file_path} n'est pas au format POOLSEQ. "
            f"Seul le format POOLSEQ est supporté pour la lecture des reads observés."
        )
    path = Path(snp_file_path)
    lines = path.read_text().splitlines()
    header_index = _find_header_index(lines)
    liste_pop = list(_parse_pool_header_line(lines, header_index).keys())
    mrc = parse_mrc_ratio(snp_file_path)
    rows = []

    # Fonction interne pour vérifier si le locus passe le seuil MRC
    def _passes_mrc(locus_reads: dict[str, tuple[int, int]]) -> bool:
        mrc_observed = observed_mrc(locus_reads)
        return mrc_observed >= mrc

    for line in lines[header_index + 1 :]:
        fields = line.split()
        if not fields:
            continue
        counts_by_population = {}

        for i in range(len(liste_pop)):
            pop_name = liste_pop[i]
            nreads1 = int(fields[2 * i])
            nreads2 = int(fields[2 * i + 1])
            counts_by_population[pop_name] = (nreads1, nreads1 + nreads2)

        if mrc <= 0 or _passes_mrc(counts_by_population):
            rows.append(counts_by_population)
            if num_loci is not None and len(rows) >= num_loci:
                break
    return rows


def coalescence_coefficient(locus_type: str, sex_ratio: float) -> float:
    """Reproduit DataC::cal_coeffcoal (data.cpp:1589-1605) : le
    coefficient qui rescale N (taille de population, prior déclaré dans
    header.txt) dans la formule de temps de coalescence
    (particuleC.cpp:1340 : temps -= coeffcoal * N / n / (n-1) * log(ra)),
    en fonction du type d'héritage du locus et du sex-ratio du dataset.

    locus_type : vocabulaire de LociDescription.total_loci ("A", "H",
    "X", "Y" ou "M" -- voir loci_parser.py), PAS "<A>" avec les chevrons.
    sex_ratio : fraction de mâles, telle que retournée par
    parse_sex_ratio (0.5 = sex-ratio équilibré).

    À l'appel avec n=2 (une paire de lignées), coeffcoal*N/2 donne le
    nombre EFFECTIF de copies de gène dans la population pour ce type de
    locus -- ex: <A> à sex_ratio=0.5 -> 2N (diploïde classique), <X> ->
    1.5N (3/4 de 2N), <Y>/<M> -> 0.5N chacun (1/4 de 2N). Cette fonction
    ne fait PAS cette conversion : elle retourne coeffcoal brut, comme le
    C++, pour rester directement comparable à cal_coeffcoal.
    """

    if locus_type == "A":
        return 16 * sex_ratio * (1 - sex_ratio)
    elif locus_type == "H":
        return 2
    elif locus_type == "X":
        return 18 * sex_ratio * (1 - sex_ratio) / (1 + sex_ratio)
    elif locus_type == "Y":
        return 2 * sex_ratio
    elif locus_type == "M":
        return 2 * (1 - sex_ratio)
    else:
        raise NotImplementedError(
            f"Type de locus inconnu pour le calcul du coefficient de "
            f"coalescence : {locus_type!r}"
        )


# Parsing des séquences ADN observées


def observed_sequences(
    mss_file_path: str | Path, list_loci: list[LociDescriptionDetailed]
) -> dict[str, list[str]]:
    """Lit les séquences ADN observées dans un fichier .mss"""
    lines = Path(mss_file_path).read_text().splitlines()

    g = (i for i, line in enumerate(lines) if line.startswith("POP"))
    first_POP_line_index = next(g, None)

    if first_POP_line_index is None:
        raise ValueError(
            f"Le fichier {mss_file_path} ne contient aucune ligne 'POP'. "
            f"Format de fichier .mss invalide."
        )
    nb_seq = len([locus for locus in list_loci if locus.ms_or_seq == "S"])

    sequences_by_indiv: dict[str, list[str]] = {}
    _MATCH_SEQUENCES = re.compile(r"^\<\[(\S+)\]\>$")
    for line in lines[first_POP_line_index + 1 :]:
        if not line.strip():
            continue
        if line.startswith("POP"):
            continue
        fields = line.split()

        match_counter = 0
        for field, loci in zip(fields[2:], list_loci, strict=True):
            match = _MATCH_SEQUENCES.match(field)
            if not match:
                continue
            match_counter += 1
            sequence = match.group(1).split("][")  # quand il s'agit d'un locus diploïde
            for i in range(len(sequence)):
                sequences_by_indiv.setdefault(loci.name, []).append(sequence[i])

        if match_counter != nb_seq:
            raise ValueError(
                f"Le nombre de séquences observées ({match_counter}) ne correspond pas au nombre de loci séquentiels ({nb_seq})."
            )
    return sequences_by_indiv


def base_frequency_by_locus(
    sequences_by_indiv: dict[str, list[str]],
) -> dict[str, dict[str, float]]:
    """Calcule la fréquence des séquences observées par locus, à partir du
    dictionnaire {nom_locus: [séquence1, séquence2, ...]} retourné par
    observed_sequences. Retourne un dictionnaire {nom_locus: Counter({séquence: fréquence})}.
    """
    base_frequencies_by_locus: dict[str, dict[str, float]] = {}
    for locus_name in sequences_by_indiv:
        sequences = sequences_by_indiv[locus_name]
        n = 0
        n_A, n_C, n_G, n_T = 0, 0, 0, 0
        for sequence in sequences:
            for base in sequence:
                if base not in "ACGTN-":
                    raise ValueError(
                        f"Base inattendue dans la séquence observée pour le locus {locus_name}: {base!r}. "
                        f"Les bases attendues sont A, C, G, T, N ou -."
                    )
                if base in "ACGT":
                    n += 1
                if base == "A":
                    n_A += 1
                elif base == "C":
                    n_C += 1
                elif base == "G":
                    n_G += 1
                elif base == "T":
                    n_T += 1
        if n == 0:
            base_frequencies_by_locus.setdefault(locus_name, {}).update(
                {"pi_A": 0.0, "pi_C": 0.0, "pi_G": 0.0, "pi_T": 0.0}
            )
        else:
            base_frequencies_by_locus.setdefault(locus_name, {}).update(
                {
                    "pi_A": n_A / n,
                    "pi_C": n_C / n,
                    "pi_G": n_G / n,
                    "pi_T": n_T / n,
                }
            )
    return base_frequencies_by_locus
