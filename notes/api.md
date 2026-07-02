
# 📄 ancestry_simulation.py

## `build_samples_argument`

### Signature

```python
build_samples_argument(snp_file_path: str) -> dict[str, int]
```

### Description

Construit l'argument `samples` attendu par msprime.sim_ancestry :
{nom_population_msprime: nombre_d_individus}, où le nom de population
msprime ("pop1", "pop2"...) correspond à l'indice utilisé dans
header.txt, mappé sur le nombre réel d'individus observés pour la
population correspondante (voir observed_data.py pour la
justification du mapping par ordre d'apparition).

---

## `simulate_independent_loci`

### Signature

```python
simulate_independent_loci(demography: msprime.Demography, samples: dict[str, int], num_loci: int, seed: int) -> Iterator[msprime.TreeSequence]
```

### Description

Simule num_loci généalogies indépendantes (un locus SNP = un
réplicat, pas de recombinaison interne ni de liaison entre loci),
sous la démographie donnée.

Retourne un itérateur (pas une liste) : pour 51250 loci, matérialiser
toutes les TreeSequence en mémoire simultanément serait coûteux --
l'appelant doit consommer cet itérateur au fil de l'eau (ex: pour
calculer des statistiques résumées locus par locus).

---

## `_draw_single_mutation_edge_child`

### Signature

```python
_draw_single_mutation_edge_child(ts, rng: random.Random) -> int
```

### Description

Tire le noeud portant la mutation unique, avec probabilité
proportionnelle à la longueur de sa branche -- algorithme de Hudson,
entièrement vectorisé via les tables (pas d'appel branch_length() par
noeud). Valable pour un arbre unique (sequence_length=1).

Chaque edge = une branche (couple parent-enfant) ; edges.child liste
donc tous les noeuds ayant une branche au-dessus d'eux (tous sauf la
racine). Longueur = time[parent] - time[child], calculé en numpy.

Validé empiriquement (proportions observées vs attendues <1% ; valeurs
de statistiques identiques à la version par branch_length() -- voir
notes/exploration.md).

---

## `_draw_single_mutation_node_fast`

### Signature

```python
_draw_single_mutation_node_fast(tree, ts, rng: random.Random) -> int
```

### Description

Version vectorisée : longueur de branche = temps(parent) - temps(noeud),
calculé en numpy sur tous les noeuds d'un coup.

---

## `_draw_single_mutation_node_vectorized`

### Signature

```python
_draw_single_mutation_node_vectorized(ts, rng: random.Random) -> None
```

### Description

Tire le noeud portant la mutation, entièrement en numpy depuis les
tables (pas d'appel branch_length() par noeud). Valable pour un arbre
unique (sequence_length=1, une seule TreeSequence).

---

## `simulate_snp_genotypes`

### Signature

```python
simulate_snp_genotypes(tree_sequences: Iterator[msprime.TreeSequence], seed: int) -> Iterator[dict[str, list[int]]]
```

### Description

Pour chaque TreeSequence (un locus = un arbre indépendant), tire
une mutation UNIQUE selon l'algorithme de Hudson (vectorisé), et
retourne les génotypes (0=ancestral, 1=dérivé) REGROUPÉS PAR
POPULATION.

Voir _draw_single_mutation_edge_child pour l'algorithme de tirage, et
la docstring d'origine pour la justification du modèle (doc DIYABC
section 2.4.3 : exactement une mutation par locus, locus toujours
polymorphe).

---


# 📄 demography_builder.py

## `evaluate_expression`

### Signature

```python
evaluate_expression(expr: str, values: dict[str, float]) -> float
```

### Description

Évalue une expression de temps ou de taille telle qu'elle apparaît
dans header.txt : un nombre littéral ("0"), un nom de paramètre tiré
("t1"), ou une somme/différence de deux noms ("t2-d3", "t2+d3").

Équivalent de ParticleC::getvalue() en C++.

---

## `build_demography`

### Signature

```python
build_demography(scenario: Scenario, values: dict[str, float]) -> msprime.Demography
```

### Description

Construit la Demography msprime correspondant au scenario, avec les
valeurs numériques déjà tirées dans `values`.

Les populations sont nommées "pop1", "pop2", ... d'après leur indice
dans header.txt (1-indexed, comme dans le fichier).

---

## `extract_referenced_names`

### Signature

```python
extract_referenced_names(expr: str) -> set[str]
```

### Description

Extrait le ou les noms de paramètres référencés par une expression
de header.txt, SANS l'évaluer numériquement -- "t2-d3" -> {"t2","d3"},
"t1" -> {"t1"}, "0" -> set() (un nombre littéral ne référence aucun
paramètre).

Utilisé pour déterminer quels paramètres un scénario utilise
réellement (nécessaire pour filtrer les colonnes du reftable.bin par
scénario -- voir reftable_loop.write_reftable_bin et notes/
exploration.md, bug "21 vs 16 paramètres pour le scénario 1").

---

## `get_parameter_names_used_by_scenario`

### Signature

```python
get_parameter_names_used_by_scenario(scenario: Scenario) -> set[str]
```

### Description

Collecte l'ensemble des noms de paramètres réellement référencés
par un scénario : tailles de population initiales, et time_expr /
new_size_expr de chacun de ses événements.

C'est ce sous-ensemble (pas la totalité des priors déclarés dans
header.txt) qui doit constituer les colonnes param[] du reftable.bin
pour ce scénario -- un scénario donné n'utilise généralement qu'une
partie des priors globaux (ex: human/header.txt a 21 priors déclarés,
mais le scénario 1 n'en référence que 16 -- ra/t11/t22/t33/t44
appartiennent aux scénarios 2-6, pas au scénario 1).

---


# 📄 loci_parser.py

## `parse_loci_description`

### Signature

```python
parse_loci_description(header_text: str) -> LociDescription
```

### Description

Extrait la description des loci à partir de header.txt, pour le
format condensé à un seul type d'héritage (cas de human).

Lève NotImplementedError si le format détecté est le format détaillé
(plusieurs lignes, un locus nommé par ligne) ou le format condensé
multi-types -- non nécessaires pour human, à implémenter si on
généralise à un autre dataset.

---

## `rewrite_loci_count`

### Signature

```python
rewrite_loci_count(header_text: str, new_total_loci: int) -> str
```

### Description

Retourne une copie de header_text où le nombre de loci déclaré dans
'loci description' est remplacé par new_total_loci -- nécessaire pour
tester avec un nombre de loci réduit sans avoir à maintenir un
header.txt séparé à la main.

Limité au même format condensé à un seul type que parse_loci_description
(lève NotImplementedError dans les mêmes cas).

---


# 📄 observed_data.py

## `count_samples_per_population`

### Signature

```python
count_samples_per_population(snp_file_path: str | Path) -> dict[str, int]
```

### Description

Compte le nombre d'individus par population dans un fichier .snp
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

L'en-tête 'IND SEX POP' peut être précédé ou non d'un commentaire libre
en première ligne (comportement observé dans data.cpp, qui teste les
deux cas) : on recherche son index plutôt que de supposer sa position,
pour ne perdre aucune ligne de données quel que soit le cas.

Lève ValueError si l'en-tête n'est trouvé dans aucune des deux
premières lignes.

---

## `population_index_to_name`

### Signature

```python
population_index_to_name(snp_file_path: str | Path) -> dict[int, str]
```

### Description

Construit le mapping entre l'indice de population utilisé dans
header.txt (1-indexed : pop1, pop2, ...) et le nom réel de population
tel qu'il apparaît dans le fichier .snp (ex: "ASW", "YRI"...).

Ex: {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"} pour human.

Voir la docstring de count_samples_per_population pour la
justification de ce mapping par ordre d'apparition (header.txt ne
nomme jamais les populations).

---


# 📄 parameter_sampling.py

## `_draw_one_value`

### Signature

```python
_draw_one_value(prior: Prior, rng: random.Random) -> float
```

---

## `draw_parameter_values`

### Signature

```python
draw_parameter_values(priors: list[Prior], constraints: list[OrderConstraint], seed: int, max_attempts: int) -> dict[str, float]
```

### Description

Tire une valeur pour chaque prior, en retirant tant que les
contraintes d'ordre ne sont pas toutes satisfaites.

Lève ConstraintsNotSatisfiedError si aucun tirage valide n'est trouvé
en max_attempts essais -- signe probable d'une configuration de
contraintes incohérente (bornes de priors incompatibles avec les
contraintes demandées) plutôt que d'une simple mauvaise chance.

---


# 📄 pipeline.py

## `build_random_demography`

### Signature

```python
build_random_demography(scenario: Scenario, header_text: str, seed: int) -> tuple[msprime.Demography, dict[str, float]]
```

### Description

Tire des valeurs de paramètres à partir des priors déclarés dans
header_text, puis construit la Demography msprime correspondant à
scenario avec ces valeurs.

Toutes les valeurs de priors du fichier sont tirées (pas seulement
celles utilisées par ce scenario précis) : plus simple, et évite de
casser des contraintes d'ordre qui pourraient porter sur des
paramètres d'autres scénarios.

Retourne (demography, values) -- les valeurs tirées sont renvoyées en
plus de la Demography, car elles seront nécessaires plus tard pour
écrire le reftable.bin (colonnes de paramètres).

---

## `build_random_demography_for_scenario_index`

### Signature

```python
build_random_demography_for_scenario_index(header_text: str, scenario_index: int, seed: int) -> tuple[msprime.Demography, dict[str, float]]
```

### Description

Variante pratique : sélectionne le scénario par son index (1-indexed,
comme dans header.txt) plutôt que de demander un objet Scenario déjà
parsé. Utile pour les tests et l'utilisation interactive.

---

## `run_poc_for_directory`

### Signature

```python
run_poc_for_directory(directory: str | Path, scenario_index: int, num_loci: int, seed: int) -> None
```

### Description

Point d'entrée de haut niveau : équivalent du `-p ./` de DIYABC.

Prend un dossier contenant header.txt et le fichier de données
observées (.snp), et produit num_loci TreeSequence mutées, simulées
sous le scénario demandé.

Le nom du fichier de données est lu sur la PREMIÈRE LIGNE de
header.txt (ex: "human_snp_all22chr_maf5.snp"), pas deviné par
extension -- c'est le contrat du format DIYABC.

Retourne (mutated_tree_sequences, values) : l'itérateur des
TreeSequence mutées, et le dict des valeurs de paramètres tirées
(nécessaires plus tard pour écrire le reftable.bin).

---

## `compute_summary_statistics`

### Signature

```python
compute_summary_statistics(reference_directory: str | Path, scenario_index: int, num_loci: int, seed: int, work_directory: str | Path, general_binary_path: str | Path, stats_filter: str) -> tuple[dict[str, float], dict[str, float]]
```

### Description

Calcule les statistiques résumées sur des données SIMULÉES par
notre pipeline, en utilisant nos formules Python validées
(summary_statistics.py) -- remplace la délégation au binaire C++
(subprocess + fichier .snp intermédiaire).

Retourne (summary_statistics, parameter_values).

---


# 📄 prior_parser.py

## `_extract_priors_section`

### Signature

```python
_extract_priors_section(header_text: str) -> list[str]
```

### Description

Extrait la section 'historical parameters priors' du texte complet de header.txt,
et retourne la liste des lignes de cette section (sans les lignes vides).

---

## `parse_priors`

### Signature

```python
parse_priors(header_text: str) -> tuple[list[Prior], list[OrderConstraint]]
```

### Description

Extrait les priors et les contraintes d'ordre de header.txt.

Retourne (priors, constraints). Une ligne qui ne correspond à aucun
des deux formats connus lève une erreur explicite plutôt que d'être
silencieusement ignorée : contrairement aux événements de scénario, on
n'a pas de raison de s'attendre à du vocabulaire non géré ici pour le
dataset human.

---

## `is_constant_prior`

### Signature

```python
is_constant_prior(prior: Prior) -> bool
```

### Description

Détecte si un prior est quasi-dégénéré (min ≈ max), donc en
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

---


# 📄 reftable_loop.py

## `_run_single_particle`

### Signature

```python
_run_single_particle(particle_index: int, reference_directory: Path, scenario_index: int, num_loci: int, general_binary_path: Path | None, base_work_directory: Path, stats_filter: str) -> ParticleResult
```

---

## `run_reftable_simulation`

### Signature

```python
run_reftable_simulation(reference_directory: str | Path, scenario_index: int, num_loci: int, nrec: int, general_binary_path: str | Path, base_work_directory: str | Path, stats_filter: str, max_workers: int | None) -> list[ParticleResult]
```

### Description

Produit nrec particules (lignes de reftable.bin) en parallèle.

base_work_directory doit déjà exister ; un sous-dossier
"particle_<i>" y est créé pour chacune des nrec particules (donc
nrec sous-dossiers au total -- à nettoyer par l'appelant si besoin,
pas fait automatiquement ici).

Les résultats sont retournés DANS L'ORDRE de particle_index (0 à
nrec-1), pas dans l'ordre de complétion des workers -- important
pour la reproductibilité de l'ordre des lignes du reftable final.

max_workers : nombre de process en parallèle (défaut : laissé à
ProcessPoolExecutor, généralement le nombre de cœurs disponibles).

---

## `write_reftable_bin`

### Signature

```python
write_reftable_bin(results: list[ParticleResult], priors: list, scenario, output_path: str | Path) -> None
```

### Description

Écrit un reftable.bin au format binaire DIYABC (vérifié contre
reftable.cpp, readReftable.R, et abcranger/readreftable.cpp -- voir
docs/synthese_diyabc_msprime.docx section 5).

Limité à un SEUL scénario actif (cohérent avec ce POC, qui ne traite
que le scénario 1 de human) : nscen=1, donc tous les résultats
doivent partager le même scenario_index -- vérifié, lève ValueError
sinon.

Filtre les colonnes de paramètres sur DEUX critères, dans cet ordre :
1. is_constant_prior : exclut les priors quasi-dégénérés (comme
   readReftable.R / abcranger)
2. get_parameter_names_used_by_scenario(scenario) : exclut les
   priors non référencés par CE scénario précis -- correction d'un
   bug découvert empiriquement (readReftable.R levait "indice hors
   limites" : notre code gardait les 21 priors du header.txt entier,
   alors que le scénario 1 n'en référence que 16 -- voir notes/
   exploration.md).

Ne gère PAS les paramètres de mutation (absents de human) -- à
ajouter (toujours en dernière position, après les paramètres
démographiques -- voir readReftable.R) si un dataset avec
microsatellites/séquences est traité plus tard.

---


# 📄 scenario_parser.py

## `split_scenario_blocks`

### Signature

```python
split_scenario_blocks(header_text: str) -> list[str]
```

### Description

Découpe le texte complet de header.txt en blocs bruts, un par
scénario, chaque bloc commençant par sa ligne d'en-tête
'scenario N [poids] (nlignes)' et s'arrêtant juste avant le bloc suivant
(ou la fin de la section, ex: 'historical parameters priors').

---

## `parse_scenario_block`

### Signature

```python
parse_scenario_block(block_text: str) -> Scenario
```

### Description

Transforme un bloc brut (en commençant par la ligne 'scenario N [...]')
en objet Scenario rempli.

---

## `_parse_event_line`

### Signature

```python
_parse_event_line(line: str) -> None
```

### Description

Transforme une ligne d'événement, ex: 't1 merge 2 1', en
SampleEvent / MergeEventn_pops / VarNeEvent selon le mot-clé rencontré.

Vocabulaire de référence : src-JMC-C++/history.cpp (ScenarioC::read_events).

---

## `parse_header_scenarios`

### Signature

```python
parse_header_scenarios(header_text: str) -> list[Scenario]
```

### Description

Point d'entrée principal : header.txt complet -> liste de Scenario.

Les blocs utilisant un vocabulaire pas encore implémenté (ex: 'split',
nécessaire aux scénarios 2/3/5/6 de human) sont sautés avec un
avertissement explicite, plutôt que de faire échouer tout le parsing.

Important : seule l'exception NotImplementedError est avalée ici,
volontairement. Toute autre exception (erreur de parsing réelle, bug)
doit continuer à se propager normalement.

---


# 📄 snp_writer.py

## `_genotypes_to_diploid`

### Signature

```python
_genotypes_to_diploid(haploid_genotypes: list[int]) -> list[int]
```

### Description

Agrège une liste de génotypes haploïdes (une valeur par lignée) en
génotypes diploïdes (0/1/2), en sommant les paires de lignées
consécutives [2i, 2i+1] -- correspondant aux deux copies
chromosomiques d'un même individu (vérifié empiriquement avec
ts.individuals()[i].nodes).

Lève ValueError si le nombre de lignées est impair (incohérent avec
une simulation en ploidy=2).

---

## `write_snp_file`

### Signature

```python
write_snp_file(genotypes_per_locus: list[dict[str, list[int]]], output_path: str | Path) -> None
```

### Description

Écrit un fichier .snp DIYABC à partir de num_loci dicts
{nom_population: [génotypes haploïdes...]}, un par locus (la forme
produite par ancestry_simulation.simulate_snp_genotypes).

Le nom de chaque individu simulé est généré comme "sim_<pop>_<n>"
(ex: "sim_pop1_1", "sim_pop1_2"...). La colonne SEX est fixée à "9"
pour tous les individus -- valeur arbitraire, non confirmée comme
sans impact pour des loci autosomaux <A> (voir notes/exploration.md
pour la justification de cette hypothèse).

genotypes_per_locus doit contenir AU MOINS un locus, et toutes les
populations doivent être présentes et avoir le même nombre de lignées
à chaque locus (cohérence vérifiée par la simulation elle-même, pas
revérifiée ici).

---


# 📄 statobs_parser.py

## `parse_statobs`

### Signature

```python
parse_statobs(statobs_text: str) -> dict[str, float]
```

### Description

Parse le contenu d'un fichier statobsRF.txt/statobs.txt en dict
{nom_colonne: valeur}.

Lève ValueError si le fichier ne contient pas exactement 2 lignes
non vides, ou si le nombre de noms et de valeurs ne correspond pas.

---


# 📄 summary_statistics.py

## `_allele_freq`

### Signature

```python
_allele_freq(haploid_genotypes: list[int]) -> float
```

### Description

Fréquence de l'allèle dérivé (1) dans une population -- équivalent
de locuslist[loc].freq[pop][1] dans le code C++.

---

## `_q1`

### Signature

```python
_q1(haploid_genotypes: list[int]) -> float
```

### Description

Probabilité d'identité par état intra-population, tirage SANS
remise -- formule exacte de sumstat.cpp::q1 (cas SNP, bias=False) :
    q1 = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1))
où y1, y2 = comptes d'allèles 0 et 1 (= freq * n).

---

## `_q2`

### Signature

```python
_q2(haploid_genotypes_a: list[int], haploid_genotypes_b: list[int]) -> float
```

### Description

Probabilité d'identité par état inter-populations -- formule exacte
de sumstat.cpp::q2 (cas SNP) :
    q2 = (y11*y21 + y12*y22) / (n1*n2)

---

## `compute_ML1`

### Signature

```python
compute_ML1(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule ML1p_i pour chaque population i : proportion de loci
monomorphes dans cette population.

Un locus est monomorphe si freq(allele 0) == 0 OU == 1, c'est-à-dire
si sum(genos) == 0 (que des ancestraux) ou sum(genos) == n (que des
dérivés) -- traduit de la condition "freq[pop][0] == 0.0 or == 1.0"
de cal_snfl dans sumstat.cpp.

---

## `compute_ML2`

### Signature

```python
compute_ML2(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

ML2p_i.j : proportion de loci monomorphes ET identiques
(même allèle fixé) dans la paire (pop_i, pop_j).
Référence : cal_snfl(npop=2) -- un locus est fixé sur la paire si
freq0_pop_i == freq0_pop_j ET vaut 0.0 ou 1.0.

---

## `compute_ML3`

### Signature

```python
compute_ML3(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

ML3p_i.j.k : même logique que ML2, sur les triplets de populations.
Référence : cal_snfl(npop=3).

---

## `compute_HW_HB`

### Signature

```python
compute_HW_HB(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule HWm_i (moyenne) et HWv_i (variance) pour chaque population,
et HBm_i.j (moyenne) et HBv_i.j (variance) pour chaque paire.

HW = 1 - q1 (hétérozygotie intra-pop)
HB = 1 - q2 (hétérozygotie inter-pop)

Référence exacte : cal_snhw et cal_snhb dans sumstat.cpp.

---

## `compute_FST1`

### Signature

```python
compute_FST1(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule FST1m_i et FST1v_i pour chaque population i.

Formule exacte lue dans cal_snfsti (sumstat.cpp) :
    FST1m_i = 1 - HWm_i / HBmoy_global
    FST1v_i = HWv_i / (HBmoy_global²)

où HBmoy_global = moyenne de TOUS les HBm (toutes les paires de
populations, pas seulement celles impliquant pop_i) -- confirmé par
le code C++ qui additionne tous les HBm sans filtre de population
(boucle sur grouplist[gr].sumstat, condition curstat == Hbstat).

FST1v n'est PAS une variance empirique locus par locus -- c'est une
formule analytique de propagation d'erreur : HWv_i / HBmoy².
Découverte en relisant précisément cal_snfsti :
    stsnp.mx2 = Hwv / (Hbmoy * Hbmoy)

---

## `_fst_wc`

### Signature

```python
_fst_wc(loci, pops) -> None
```

### Description

Weir & Cockerham vectorisé sur tous les loci.
Retourne (FSTm, FSTv). Formule identique à cal_snfstd, mais toutes
les opérations par-locus sont faites en numpy sur des vecteurs de
longueur n_loci au lieu d'une boucle Python.

---

## `compute_FST2`

### Signature

```python
compute_FST2(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

FST2 via _fst_wc -- même code que FST3/FST4.

---

## `compute_FST3_FST4_FSTG`

### Signature

```python
compute_FST3_FST4_FSTG(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

FST3/FST4/FSTG : Weir & Cockerham généralisé à 3, 4, toutes pops.
Référence : cal_snfstd(npop=3/4/0) dans sumstat.cpp.
Ordre COMB pour FST3/FST4/FSTG.

---

## `compute_NEI`

### Signature

```python
compute_NEI(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire, vectorisé.
NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)
x_prev persiste si n==0.

---

## `compute_AML`

### Signature

```python
compute_AML(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

AMLm et AMLv : coefficient d'admixture ML sur triplets.
Référence : cal_snaml dans sumstat.cpp.
aml = (f3 - f2) / (f1 - f2), clampé à [0,1] si hors bornes.
w=0 si f1==f2 (locus non informatif, exclu de la moyenne pondérée).
Ordre des triplets : HALF (halfsortedbypairs), reproduit empiriquement.
samp[0]=hybride, samp[1]=parent1, samp[2]=parent2.

---

## `compute_F3_F4`

### Signature

```python
compute_F3_F4(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

F3m/F3v : statistique de Patterson f3 sur triplets (HALF).
F4m/F4v : statistique de Patterson f4 sur quadruplets (HALF).
Référence : cal_snf3r et cal_snf4r dans sumstat.cpp (branche SNP).
F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)  où np = nb lignées pop hybride
F4 = (a-b)*(c-d)

---

## `compute_all_statistics`

### Signature

```python
compute_all_statistics(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule toutes les statistiques implémentées et retourne un dict
unifié {nom_stat: valeur} -- même format que parse_statobs().

LIMITES ACTUELLES (à étendre au fil des validations) :
- ML1 seulement (pas ML2/ML3)
- HW, HB (moyenne et variance)
- FST1 seulement (pas FST2/FST3/FST4/FSTG)
- NEI, AML, F3, F4 : NON ENCORE IMPLÉMENTÉS

---
