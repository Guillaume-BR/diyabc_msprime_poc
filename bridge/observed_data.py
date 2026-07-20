"""
Lecture des données observées DIYABC (fichiers .snp, format individu par
ligne) -- comptage du nombre d'individus par population (nécessaire pour
savoir combien d'échantillons demander à msprime.sim_ancestry()), lecture
du sex-ratio et du sexe par individu (nécessaires pour les loci
<X>/<Y>/<M>, dont la ploïdie et le coefficient de coalescence dépendent du
sexe -- voir ParticleC::calploidy et DataC::cal_coeffcoal), et lecture du
seuil MAF ("<MAF=N%>" ou "<MAF=hudson>" en tête de fichier, via
parse_maf_ratio) consommé par with_maf_filter/with_maf_filter_shared_
ancestry (ancestry_simulation.py) pour décider s'il faut rejeter-et-
resimuler un locus sous le seuil.

Référence : src-JMC-C++/data.cpp (détection du format "IND SEX POP",
lecture du sex-ratio "<NM=xNF>" et du MAF en tête de fichier).
Ce module ne lit PAS les génotypes eux-mêmes : on simule des données
artificielles avec msprime, on ne réutilise jamais les données observées
réelles dans le pipeline de simulation.
"""

from collections import Counter
from pathlib import Path


def _find_header_index(lines: list[str]) -> int:
    """Repère l'index de la ligne d'en-tête 'IND SEX POP' parmi les deux
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
        ),
        None,
    )
    if header_index is None:
        raise ValueError(
            f"En-tête 'IND SEX POP' non trouvé dans les deux premières "
            f"lignes. Lignes lues : {lines[:2]!r}"
        )
    return header_index


def count_samples_per_population(snp_file_path: str | Path) -> dict[str, int]:
    """Compte le nombre d'individus par population dans un fichier .snp
    DIYABC au format 'IND SEX POP <génotypes...>'.

    Ex: pour human_snp_all22chr_maf5.snp -> {"ASW": 30, "YRI": 30, ...}

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
    pop_index = lines[header_index].split().index("POP")

    counts = Counter()
    for line in lines[header_index + 1 :]:
        fields = line.split()
        if not fields:
            continue
        counts[fields[pop_index]] += 1

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


def population_index_to_name(snp_file_path: str | Path) -> dict[int, str]:
    """Construit le mapping entre l'indice de population utilisé dans
    header.txt (1-indexed : pop1, pop2, ...) et le nom réel de population
    tel qu'il apparaît dans le fichier .snp (ex: "ASW", "YRI"...).

    Ex: {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"} pour human.

    Voir la docstring de count_samples_per_population pour la
    justification de ce mapping par ordre d'apparition (header.txt ne
    nomme jamais les populations).
    """
    names_in_order = list(count_samples_per_population(snp_file_path).keys())
    return {i + 1: name for i, name in enumerate(names_in_order)}


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
