
# 📄 ancestry_simulation.py

## `build_samples_argument`

### Signature

```python
build_samples_argument(snp_file_path: str) -> dict[str, int]
```

### Description

Construit l'argument `samples` de msprime.sim_ancestry pour un locus <A>.

Le nom de population msprime ("pop1", "pop2"...) correspond à
l'indice utilisé dans header.txt, mappé sur le nombre réel
d'individus observés pour la population correspondante (voir
observed_data.py pour la justification du mapping par ordre
d'apparition).

Un seul appel à count_samples_per_population (pas
population_index_to_name EN PLUS, qui relit et rescanne tout le
fichier .snp pour ne faire que redériver les mêmes clés dans le même
ordre) -- l'indice 1-based se déduit directement de la position dans
ce même dict, garanti dans l'ordre de première apparition (voir sa
docstring).

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Un dict {nom_population_msprime: nombre_d_individus}.

---

## `_sample_sets_from_sexes`

### Signature

```python
_sample_sets_from_sexes(sexes_by_population: dict[str, list[str]]) -> list[msprime.SampleSet]
```

### Description

Construit la liste des SampleSet à partir des sexes par population.

Args:
    sexes_by_population: dict {nom_population: [liste de sexes]}.

Returns:
    Une liste de msprime.SampleSet (2 par population).

---

## `_male_counts_from_sexes`

### Signature

```python
_male_counts_from_sexes(sexes_by_population: dict[str, list[str]]) -> dict[str, int]
```

### Description

Construit le dict {nom_population: nombre_d_individus_mâles} à partir des
sexes par population.

Args:
    sexes_by_population: dict {nom_population: [liste de sexes]}.

Returns:
    Un dict {nom_population: nombre_d_individus_mâles}.

---

## `build_sex_stratified_samples_argument`

### Signature

```python
build_sex_stratified_samples_argument(snp_file_path: str) -> list[msprime.SampleSet]
```

### Description

Construit l'argument `samples` de msprime.sim_ancestry pour un locus <X>.

Contrairement à build_samples_argument (un compte par population,
ploidy uniforme), <X> a besoin d'une ploidy DIFFÉRENTE par individu
selon son sexe (femelles=2 copies, mâles=1 -- voir
ParticleC::calploidy, particuleC.cpp:220-233), donc une liste de
msprime.SampleSet plutôt qu'un simple dict : 2 SampleSet par
population, un pour les femelles (ploidy=2), un pour les mâles
(ploidy=1) -- population= est le nom msprime ("pop1", "pop2"...),
PAS le nom réel du fichier .snp (même traduction que
build_samples_argument, via population_index_to_name).

IMPORTANT -- le ploidy PAR SampleSet ne contrôle QUE le nombre de
lignées regroupées par individu dans le résultat, PAS le taux de
coalescence : cette liste doit être utilisée avec
simulate_independent_loci(..., ploidy=1) et une `demography` déjà
rescalée via rescale_demography(demography,
coalescence_coefficient("X", sex_ratio) / 2) -- vérifié
empiriquement avec le mentor que c'est le ploidy GLOBAL de
sim_ancestry qui interprète initial_size, pas celui des SampleSet.

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    La liste des msprime.SampleSet (2 par population).

Raises:
    ValueError: Si un individu a le sexe "9" (inconnu) -- ex:
        human_snp_all22chr_maf5.snp, où AUCUN individu n'est sexé
        (dataset <A>-only) : on ne peut pas construire un
        échantillonnage <X> dessus, mieux vaut le signaler
        explicitement que de produire un résultat silencieusement
        faux (individual_sexes_per_population laisse ce choix à
        l'appelant, c'est ici qu'il se prend).

---

## `build_male_only_samples_argument`

### Signature

```python
build_male_only_samples_argument(snp_file_path: str) -> dict[str, int]
```

### Description

Construit l'argument `samples` de msprime.sim_ancestry pour un locus <Y>.

Le nom de population msprime ("pop1", "pop2"...) correspond à
l'indice utilisé dans header.txt, mappé sur le nombre réel
d'individus MÂLES observés pour la population correspondante (voir
observed_data.py pour la justification du mapping par ordre
d'apparition).

PAS pour <M> : le mitochondrial est transmis uniquement par les
mères, mais présent et échantillonné chez TOUS les individus
(mâles et femelles), contrairement à <Y> qui n'existe que chez les
mâles -- <M> doit réutiliser build_samples_argument (tout le monde)
avec ploidy=1, pas cette fonction.

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Un dict {nom_population_msprime: nombre_d_individus_mâles}.

Raises:
    ValueError: Si un individu a le sexe "9" (inconnu).

---

## `simulate_independent_loci`

### Signature

```python
simulate_independent_loci(demography: msprime.Demography, samples: dict[str, int] | list[msprime.SampleSet], num_loci: int, seed: int, ploidy: int) -> Iterator[tskit.TreeSequence]
```

### Description

Simule num_loci généalogies indépendantes sous la démographie donnée.

Un locus SNP = un réplicat, pas de recombinaison interne ni de
liaison entre loci. Retourne un itérateur (pas une liste) : pour
51250 loci, matérialiser toutes les TreeSequence en mémoire
simultanément serait coûteux -- l'appelant doit consommer cet
itérateur au fil de l'eau (ex: pour calculer des statistiques
résumées locus par locus).

Args:
    demography: La démographie msprime.
    samples: dict[str, int] (un compte par population, ploidy
        uniforme -- <A>/<M>, voir build_samples_argument) ou
        list[msprime.SampleSet] (ploidy hétérogène par sous-groupe
        au sein d'une population -- <X>, voir
        build_sex_stratified_samples_argument). Les deux formes
        sont transmises telles quelles à msprime.sim_ancestry, qui
        les accepte indifféremment.
    num_loci: Le nombre de généalogies indépendantes à simuler.
    seed: La graine de la simulation.
    ploidy: 2 (défaut) pour <A>, cohérent avec une transmission
        diploïde classique -- chaque "sample individual" de
        `samples` compte pour 2 lignées génomiques. Pour <Y>/<M>,
        passer ploidy=1 avec une `demography` déjà rescalée par
        rescale_demography (voir demography_builder.py) : ces
        loci sont haploïdes (une seule copie de gène transmise),
        et le facteur de rescaling de Ne (coalescence_coefficient,
        observed_data.py) suppose cette combinaison ploidy=1 + Ne
        rescalé, pas ploidy=2 + Ne d'origine. Pour <X>, passer
        aussi ploidy=1 (voir
        build_sex_stratified_samples_argument : c'est le ploidy
        PAR SampleSet, pas ce paramètre global, qui donne 2
        copies aux femelles et 1 aux mâles -- ce paramètre-ci ne
        fixe que le taux de coalescence, via la Demography déjà
        rescalée).

Returns:
    Un itérateur de num_loci TreeSequence indépendantes.

---

## `simulate_shared_ancestry_loci`

### Signature

```python
simulate_shared_ancestry_loci(demography: msprime.Demography, samples: dict[str, int] | list[msprime.SampleSet], num_loci: int, seed: int, ploidy: int) -> Iterator[tskit.TreeSequence]
```

### Description

Simule UNE SEULE généalogie puis la retourne répétée num_loci fois.

Pour <Y>/<M>, dont tous les loci d'un même type partagent la même
généalogie réelle (non-recombinants, transmission uniparentale),
contrairement à <A>/<X> qui tirent un arbre indépendant par locus
(simulate_independent_loci). Reproduit le comportement de
particuleC.cpp:2422-2435 (GeneTreeY/GeneTreeM : premier locus <Y> ou
<M> tire un arbre normalement, tous les suivants COPIENT ce même
arbre -- seule la mutation change d'un locus à l'autre).

IMPORTANT -- ne PAS réimplémenter le tirage de mutation ici :
simulate_snp_genotypes(tree_sequences, seed) lit déjà tree_sequences
au fil de l'eau sans jamais modifier les TreeSequence qu'elle reçoit,
et son rng avance à chaque itération -- lui donner le MÊME objet
TreeSequence répété num_loci fois (au lieu de num_loci objets
différents) suffit à obtenir num_loci mutations indépendantes sur
UNE SEULE généalogie, sans aucune modification de cette fonction
(vérifié empiriquement : 5 répétitions du même arbre -> 5 génotypes
différents).

Args:
    demography: La démographie msprime.
    samples: Même contrat que simulate_independent_loci (voir sa
        docstring) -- cette fonction ne fait que réutiliser
        simulate_independent_loci avec num_loci=1, elle ne
        réinterprète pas ce paramètre.
    num_loci: Le nombre de fois où répéter la généalogie unique.
    seed: La graine de la simulation.
    ploidy: Pour <Y>/<M>, passer ploidy=1 (défaut) avec une `demography` déjà rescalée par
        rescale_demography (voir demography_builder.py) : ces
        loci sont haploïdes (une seule copie de gène transmise),
        et le facteur de rescaling de Ne (coalescence_coefficient,
        observed_data.py) suppose cette combinaison ploidy=1 + Ne
        rescalé, pas ploidy=2 + Ne d'origine.

Returns:
    Un itérateur de num_loci TreeSequence, toutes identiques (le
    même objet Python répété).

---

## `_draw_single_mutation_edge_child`

### Signature

```python
_draw_single_mutation_edge_child(ts: tskit.TreeSequence, rng: random.Random) -> int
```

### Description

Tire le noeud portant la mutation unique de l'algorithme de Hudson.

Probabilité proportionnelle à la longueur de sa branche --
entièrement vectorisé via les tables (pas d'appel branch_length()
par noeud). Valable pour un arbre unique (sequence_length=1).

Chaque edge = une branche (couple parent-enfant) ; edges.child liste
donc tous les noeuds ayant une branche au-dessus d'eux (tous sauf la
racine). Longueur = time[parent] - time[child], calculé en numpy.

Validé empiriquement (proportions observées vs attendues <1% ; valeurs
de statistiques identiques à la version par branch_length() -- voir
notes/exploration.md).

Args:
    ts: La TreeSequence (un seul arbre, sequence_length=1).
    rng: Le générateur aléatoire à utiliser.

Returns:
    L'ID du noeud enfant portant la mutation.

---

## `compute_population_layout`

### Signature

```python
compute_population_layout(ts: tskit.TreeSequence) -> list[tuple[str | None, np.ndarray]]
```

### Description

Calcule le layout (nom de population, IDs des noeuds échantillons) d'une TreeSequence.

Factorisé pour pouvoir être calculé UNE SEULE FOIS et réutilisé sur
plusieurs loci/tentatives qui partagent la même `demography`/
`samples` d'origine -- seule la topologie coalescente varie d'un
réplicat à l'autre, jamais l'assignation des noeuds échantillons aux
populations (vérifié empiriquement, y compris entre réplicats tirés
avec des graines différentes). Voir `simulate_snp_genotypes` (cache
par défaut sur un flux de plusieurs loci) et les boucles de rejet MAF
de `with_maf_filter`/`with_maf_filter_shared_ancestry` (cache
explicite à travers les tentatives, voir notes/exploration.md,
entrée du 20/07/2026).

Args:
    ts: La TreeSequence à inspecter.

Returns:
    La liste des (nom_population, IDs des noeuds échantillons de
    cette population), une entrée par population non vide.

---

## `simulate_snp_genotypes`

### Signature

```python
simulate_snp_genotypes(tree_sequences: Iterator[tskit.TreeSequence], seed: int, population_layout: list[tuple[str | None, np.ndarray]] | None) -> Iterator[dict[str, list[int]]]
```

### Description

Tire une mutation par locus (Hudson) et retourne les génotypes par population.

Pour chaque TreeSequence (un locus = un arbre indépendant), tire
une mutation UNIQUE selon l'algorithme de Hudson (vectorisé), et
retourne les génotypes (0=ancestral, 1=dérivé) REGROUPÉS PAR
POPULATION. Voir _draw_single_mutation_edge_child pour l'algorithme
de tirage, et la docstring d'origine pour la justification du
modèle (doc DIYABC section 2.4.3 : exactement une mutation par
locus, locus toujours polymorphe).

Args:
    tree_sequences: Un itérateur de TreeSequence, un arbre
        indépendant par locus.
    seed: La graine du tirage de mutation.
    population_layout: Voir `compute_population_layout`. Si
        `None` (cas d'un appel unique sur tout un flux de loci,
        ex: chemin `maf=0.0`), calculé UNE SEULE FOIS ici même, au
        premier locus, et réutilisé pour tous les suivants --
        valable car tous les `tree_sequences` d'un même appel
        partagent la même `demography`/`samples` d'origine (mêmes
        réplicats d'un seul appel à simulate_independent_loci/
        simulate_shared_ancestry_loci) : seule la topologie
        coalescente varie d'un locus à l'autre, jamais
        l'assignation des noeuds échantillons aux populations
        (vérifié empiriquement). Si fourni par l'appelant (ex:
        boucles de rejet MAF de `with_maf_filter`/`with_maf_
        filter_shared_ancestry`, qui appellent cette fonction une
        fois PAR TENTATIVE et calculent donc leur propre cache à
        travers les tentatives), utilisé tel quel sans jamais être
        recalculé. Sans ce cache, le redécodage du metadata des
        populations et le refiltrage de ts.samples(population=...)
        à CHAQUE locus représentaient à eux seuls ~20% du temps
        d'une particule sur 5000 loci (voir notes/exploration.md,
        entrée du 20/07/2026) -- le plus gros poste évitable du
        surcoût tskit par locus identifié dans cette investigation.

Returns:
    Un itérateur de dicts {nom_population: [génotype, ...]} (un
    dict par locus).

---

## `observed_maf`

### Signature

```python
observed_maf(locus_genotypes: dict[str, list[int]]) -> float
```

### Description

Calcule la MAF poolée sur toutes les populations pour un locus.

Comme ParticleC::mafreached (min(dérivé, ancestral) / total) -- pas
juste la fréquence dérivée.

Args:
    locus_genotypes: Dict {nom_population: [génotype, ...]} (0 ou
        1) pour un seul locus.

Returns:
    La MAF (fréquence de l'allèle minoritaire).

---

## `with_maf_filter`

### Signature

```python
with_maf_filter(demography: msprime.Demography, samples: dict[str, int] | list[msprime.SampleSet], num_loci: int, maf: float, seed: int, ploidy: int) -> Iterator[dict[str, list[int]]]
```

### Description

Simule des loci SNP indépendants avec filtre MAF.

MAF = minor allele frequency, cf. doc DIYABC section 2.4.3 : si la
fréquence de l'allèle MINORITAIRE (le moins fréquent des deux,
dérivé ou ancestral -- pas forcément le dérivé) est strictement
inférieure à `maf`, on rejette ce locus et on en resimule un
nouveau (nouvelle généalogie + nouvelle mutation, jamais de
recyclage de l'arbre rejeté) jusqu'à obtenir `num_loci` loci
acceptés. Reproduit `ParticleC::mafreached`
(particuleC.cpp:2194-2210).

`maf=0.0` (équivalent DIYABC de `<MAF=hudson>` ou d'un tag absent)
délègue directement à `simulate_independent_loci` +
`simulate_snp_genotypes` avec la même graine pour les deux (comme
le fait déjà chaque branche de `simulate_genotypes_for_locus_type`)
-- comportement et résultats identiques à un appel direct de ces
deux fonctions, pour ne rien changer aux datasets déjà validés qui
n'ont pas de filtre MAF actif (human, toy_example5, ...).

`maf>0.0` : les tentatives sont tirées PAR LOT de
`max(_MAF_BATCH_SIZE, num_loci // 4)` (un seul appel
`simulate_independent_loci(num_replicates=batch_size)` au lieu d'un
appel par tentative individuelle) -- mesuré empiriquement ~5.6x plus
rapide qu'un appel un-par-un sur toy_example3/scenario3/maf=0.05 (voir
_MAF_BATCH_SIZE), le coalescent restant identique (nouvelle
généalogie à chaque rejet, jamais de recyclage d'un arbre rejeté,
comme avant). Le lot scale avec `num_loci` plutôt que d'être fixe :
voir _MAF_BATCH_SIZE pour le détail. La structure population/
échantillons (`population_layout`) ne dépend que de `demography`/
`samples`, jamais de la graine tirée -- elle est donc calculée une
seule fois, à la toute première tentative, et réutilisée pour toutes
les suivantes (voir notes/exploration.md, entrée du 20/07/2026).

Args:
    demography: La démographie msprime.
    samples: Même contrat que simulate_independent_loci.
    num_loci: Le nombre de loci acceptés à produire.
    maf: Le seuil MAF, déjà extrait (ex: via `parse_maf_ratio` sur
        le fichier .snp observé) -- cette fonction ne lit aucun
        fichier, à l'appelant de décider d'où vient le seuil.
    seed: La graine de la simulation.
    ploidy: Transmis tel quel à `simulate_independent_loci` (même
        contrat -- 2 pour <A>, 1 pour <H>/<X> avec une
        `demography` déjà rescalée, voir sa docstring).

Returns:
    Un itérateur de `num_loci` dicts {nom_population:
    [génotype, ...]}, tous au-dessus du seuil MAF.

---

## `with_maf_filter_shared_ancestry`

### Signature

```python
with_maf_filter_shared_ancestry(demography: msprime.Demography, samples: dict[str, int] | list[msprime.SampleSet], num_loci: int, maf: float, seed: int, ploidy: int) -> Iterator[dict[str, list[int]]]
```

### Description

Variante de with_maf_filter pour <Y>/<M> (généalogie partagée).

Contrairement aux loci <A>/<H>/<X> (chaque locus = sa propre
généalogie indépendante), tous les loci <Y> (resp. <M>) d'une même
particule PARTAGENT UNE SEULE généalogie (voir
simulate_shared_ancestry_loci) -- seule la mutation diffère d'un
locus à l'autre.

Reproduit exactement `particuleC.cpp:2424-2495` : le cache
GeneTreeY/GeneTreeM est rempli AVANT le test MAF, donc indépendamment
de son résultat -- la généalogie est tirée UNE SEULE FOIS (au tout
premier appel), et un rejet MAF ne fait jamais redessiner l'arbre,
seulement retirer une nouvelle mutation SUR CE MÊME ARBRE, jusqu'à
obtenir `num_loci` loci acceptés. Voir aussi with_maf_filter (loci
<A>/<H>/<X>), qui redessine au contraire une toute nouvelle
généalogie à chaque rejet -- les deux mécanismes sont réellement
différents côté DIYABC, pas juste une simplification.

Args:
    demography: La démographie msprime (déjà rescalée par
        l'appelant si nécessaire).
    samples: Même contrat que simulate_independent_loci.
    num_loci: Le nombre de loci acceptés à produire.
    maf: Le seuil MAF, déjà extrait (voir with_maf_filter).
        `maf=0.0` délègue directement à
        simulate_shared_ancestry_loci + simulate_snp_genotypes
        avec la même graine pour les deux, comportement identique
        à un appel direct de ces deux fonctions.
    seed: La graine de la simulation.
    ploidy: Transmis tel quel à simulate_independent_loci/
        simulate_shared_ancestry_loci.

Returns:
    Un itérateur de `num_loci` dicts {nom_population:
    [génotype, ...]}, tous au-dessus du seuil MAF.

---

## `simulate_genotypes_for_locus_type`

### Signature

```python
simulate_genotypes_for_locus_type(demography: msprime.Demography, snp_file_path: str, locus_type: str, num_loci: int, seed: int) -> Iterator[dict[str, list[int]]]
```

### Description

Point d'entrée unique de simulation de génotypes SNP, par type de locus.

Choisit la bonne combinaison samples/demography-rescalée-ou-non/
ploidy/fonction de simulation-indépendante-ou-partagée selon
locus_type, puis retourne les génotypes simulés (même contrat de
sortie que simulate_snp_genotypes, qu'on appelle en dernière étape
dans tous les cas -- elle ne dépend jamais de locus_type
elle-même).

sex_ratio n'est PAS un paramètre : il est dérivé automatiquement de
snp_file_path via parse_sex_ratio, comme tout le reste (samples,
sexes par individu) -- l'appelant n'a jamais besoin de le connaître.

Dispatch :
  - "A" : build_samples_argument, demography TELLE QUELLE (pas de
    rescale_demography), ploidy=2, with_maf_filter.
  - "H" : build_samples_argument, demography rescalée par
    coalescence_coefficient("H", sex_ratio) / 2, ploidy=1,
    with_maf_filter.
  - "X" : build_sex_stratified_samples_argument, demography
    rescalée par coalescence_coefficient("X", sex_ratio) / 2,
    ploidy=1, with_maf_filter.
  - "Y" : build_male_only_samples_argument, demography rescalée par
    coalescence_coefficient("Y", sex_ratio) / 2, ploidy=1,
    with_maf_filter_shared_ancestry (arbre unique partagé).
  - "M" : build_samples_argument (TOUT le monde, pas mâles seuls --
    voir la docstring de build_male_only_samples_argument sur ce
    point précis), demography rescalée par
    coalescence_coefficient("M", sex_ratio) / 2, ploidy=1,
    with_maf_filter_shared_ancestry.
  - tout autre locus_type : lever NotImplementedError (même style
    que coalescence_coefficient pour un type inconnu).

MAF : le seuil est lu une fois via parse_maf_ratio(snp_file_path) et
délégué à with_maf_filter ("A"/"H"/"X", généalogie indépendante par
locus) ou with_maf_filter_shared_ancestry ("Y"/"M", généalogie
partagée) -- les deux gèrent elles-mêmes le cas maf=0.0 (pas de
filtre, comportement identique à un appel direct des fonctions
sous-jacentes) et le cas maf>0.0 (boucle de rejet).

"Y"/"M" (MAF quelconque) utilisent with_maf_filter_shared_ancestry,
pas with_maf_filter : ces deux types partagent UNE SEULE généalogie
entre tous leurs loci (simulate_shared_ancestry_loci) -- un rejet MAF
ne redessine jamais l'arbre, seulement la mutation (voir la docstring
de with_maf_filter_shared_ancestry, qui reproduit exactement
particuleC.cpp:2424-2495).

Args:
    demography: La démographie <A> "de base" (construite par
        build_demography, PAS encore rescalée) -- c'est CETTE
        fonction qui décide si/comment la rescaler selon
        locus_type, jamais l'appelant.
    snp_file_path: Chemin du fichier .snp observé.
    locus_type: "A", "H", "X", "Y" ou "M".
    num_loci: Le nombre de loci à simuler.
    seed: La graine de la simulation.

Returns:
    Un itérateur de `num_loci` dicts {nom_population:
    [génotype, ...]}.

Raises:
    NotImplementedError: Si locus_type est inconnu.

---

## `simulate_poolseq_reads`

### Signature

```python
simulate_poolseq_reads(tree_sequences: Iterator[tskit.TreeSequence], observed_reads_per_locus: list[dict[str, tuple[int, int]]], seed: int, population_layout: list[tuple[str | None, np.ndarray]] | None) -> Iterator[dict[str, tuple[int, int]]]
```

### Description

Simule les lectures PoolSeq de chaque locus.

Tire une mutation par locus (même algorithme de Hudson que
simulate_snp_genotypes), puis convertit la proportion de lignées
dérivées de chaque population en un tirage binomial de lectures,
calé sur la profondeur totale RÉELLEMENT observée à ce locus/cette
population (`observed_reads_per_locus`) -- seule la répartition
allèle1/allèle2 est simulée, jamais la couverture elle-même.

Args:
    tree_sequences: Un itérateur de TreeSequence simulées, une par
        locus.
    observed_reads_per_locus: Une liste de dicts {nom_population:
        (nreads_dérivé, nreads_total)} observés, un par locus.
    seed: La graine du tirage de mutation.
    population_layout: Voir `compute_population_layout`. Si
        `None` (cas d'un appel unique sur tout un flux de loci,
        ex: chemin `mrc<=0`), calculé UNE SEULE FOIS ici même, au
        premier locus, et réutilisé pour tous les suivants --
        même principe que `simulate_snp_genotypes`. Si fourni par
        l'appelant (ex: boucle de rejet MRC de `with_mrc_filter`,
        qui appelle cette fonction une fois PAR TENTATIVE et
        calcule donc son propre cache à travers les tentatives),
        utilisé tel quel sans jamais être recalculé.

Returns:
    Un itérateur de dicts {nom_population: (nreads_dérivé,
    nreads_total)} simulés, un par locus.

---

## `_reindex_reads_by_msprime_name`

### Signature

```python
_reindex_reads_by_msprime_name(observed_reads_per_locus: list[dict[str, tuple[int, int]]], snp_file_path: str) -> list[dict[str, tuple[int, int]]]
```

### Description

Reindexe les lectures observées pour utiliser les noms de population msprime.

Args:
    observed_reads_per_locus: Une liste de dicts {nom_population
        réel: (nreads_dérivé, nreads_total)}, un par locus.
    snp_file_path: Chemin du fichier .snp, pour obtenir la
        correspondance des noms de population.

Returns:
    La même liste, avec les clés remplacées par les noms de
    population msprime ("pop1", "pop2"...).

---

## `with_mrc_filter`

### Signature

```python
with_mrc_filter(demography: msprime.Demography, samples: dict[str, int] | list[msprime.SampleSet], num_loci: int, mrc: float, observed_reads_per_locus: list[dict[str, tuple[int, int]]], seed: int, ploidy: int) -> Iterator[dict[str, tuple[int, int]]]
```

### Description

Simule des loci SNP indépendants avec filtre MRC.

MRC = minimum read count. Si le nombre de lectures dérivées est
strictement inférieur à `mrc`, on rejette ce locus et on en
resimule un nouveau. Reproduit le comportement de
`ParticleC::mrc_reached`.

`mrc>0` : les tentatives sont tirées depuis un POOL PARTAGÉ ENTRE TOUS
LES LOCI, pas un pool privé par locus -- contrairement à un batching
naïf "par locus" (un nouveau lot de `_MRC_BATCH_SIZE` généalogies à
chaque `locus_index`, même si ce locus n'a besoin que d'UNE seule
tentative), qui paie un plancher d'un appel `simulate_independent_loci`
PAR LOCUS quel que soit le taux d'acceptation. Ici, un seul flux
continu de généalogies (régénéré par lot de `_MRC_BATCH_SIZE`
uniquement quand épuisé) est consommé par n'importe quel locus qui a
besoin d'une nouvelle tentative -- si la plupart des loci passent dès
le premier tirage (cas courant), un seul lot peut servir des dizaines
de loci au lieu d'un lot par locus. Gain mesuré empiriquement (script
jetable, toy_example4, mrc=5) : ~1.5x supplémentaire par rapport au
batching par-locus, quelle que soit la taille du lot (le partage
compte, pas la taille).

Le compteur `attempt` est GLOBAL et n'est jamais remis à zéro par
locus -- ça élimine par construction le risque de corrélation qui
existait avec l'ancien design par-locus (deux loci ayant besoin du
même nombre de tentatives tiraient alors le même arbre/mutation,
confirmé empiriquement le 22/07/2026) : chaque tentative, tous loci
confondus, consomme une position distincte dans un flux continu,
jamais réutilisée.

Args:
    demography: La démographie msprime.
    samples: Même contrat que simulate_independent_loci.
    num_loci: Le nombre de loci acceptés à produire.
    mrc: Le seuil MRC, déjà extrait via `parse_mrc_ratio` sur le
        fichier .snp observé.
    observed_reads_per_locus: Une liste de dicts {nom_population:
        (nreads_dérivé, nreads_total)} observés, un par locus.
    seed: La graine de la simulation.
    ploidy: Transmis tel quel à `simulate_independent_loci`.

Returns:
    Un itérateur de `num_loci` dicts {nom_population:
    (nreads_dérivé, nreads_total)}, tous au-dessus du seuil MRC.

---

## `prepare_poolseq_observed_reads`

### Signature

```python
prepare_poolseq_observed_reads(snp_file_path: str, num_loci: int) -> list[dict[str, tuple[int, int]]]
```

### Description

Prépare les lectures observées pour la simulation PoolSeq.

Lit le fichier .snp et tronque aux `num_loci` premières entrées.

Args:
    snp_file_path: Chemin du fichier .snp (doit être POOLSEQ).
    num_loci: Le nombre de loci à conserver.

Returns:
    Une liste de dicts {nom_population_msprime: (nreads_dérivé,
    nreads_total)}, un par locus.

---

## `simulate_poolseq_reads_with_mrc_filter`

### Signature

```python
simulate_poolseq_reads_with_mrc_filter(demography: msprime.Demography, snp_file_path: str, seed: int, num_loci: int, observed_reads_per_locus: list[dict[str, tuple[int, int]]]) -> Iterator[dict[str, tuple[int, int]]]
```

### Description

Point d'entrée unique de simulation de lectures pour un fichier PoolSeq.

Pendant de simulate_genotypes_for_locus_type (IndSeq), mais sans
dispatch multi-type : un fichier PoolSeq n'a jamais qu'un seul type
de locus déclaré (`<A>`, cf. `data.cpp:529` -- seule la classe de
locus autosomale diploïde est supportée pour PoolSeq côté DIYABC).
Aucun rescale de Ne n'est nécessaire (même `coeffcoal` que l'IndSeq
`<A>` standard, cf. `data.cpp:1589-1603` -- PoolSeq a `type=15`,
`15 % 5 == 0`, donc tombe dans exactement la même branche que le
cas autosomal diploïde standard).

ploidy=2, PAS ploidy=1 (corrigé le 22/07/2026 -- voir
notes/exploration.md) : DIYABC construit l'arbre de généalogie
PoolSeq en réutilisant tel quel le chemin `<A>` standard --
`HAPLOID_SAMPLE_SIZE` (déclaré dans le fichier .snp, `POOL
pop:N`) est le nombre de COPIES DE GÈNES à échantillonner, mais le
Ne (N1/N2/...) reste interprété en individus DIPLOÏDES, exactement
comme `<A>` IndSeq (`particuleC.cpp:1185`, `data.cpp:970-974`,
confirmé par exploration du code source réel -- aucune branche de
simulation généalogique spécifique à PoolSeq n'existe, seule la
lecture du fichier et le tirage des reads le sont). Passer
`samples=build_samples_argument(...)` (comptes haploïdes) avec
`ploidy=1` -- ce qui a été fait par erreur avant cette correction --
revient à traiter le Ne comme s'il était HAPLOÏDE : la coalescence
devient deux fois plus rapide en unités de générations réelles pour
le MÊME Ne déclaré, ce qui gonfle artificiellement toute
différenciation entre populations (FST/F3/F4 ~60-140% trop élevés,
confirmé empiriquement sur toy_example4 en comparaison appariée
contre un vrai reftable DIYABC) sans affecter les statistiques
mono-population (HW/ML1, qui ne dépendent pas de la vitesse relative
de coalescence entre populations). D'où le //2 ci-dessous : on donne
à msprime un compte d'INDIVIDUS (`ploidy=2` double automatiquement
en lignées), pas un compte de lignées déjà doublé -- le nombre total
de lignées échantillonnées (donc la taille de l'arbre) reste
identique (`HAPLOID_SAMPLE_SIZE`), seule la vitesse de coalescence
relative au Ne change, pour retomber sur la même convention que
l'IndSeq `<A>`.

Compose, dans l'ordre :
  - `parse_mrc_ratio(snp_file_path)` -- seuil MRC (défaut 1 si
    `<MRC=...>` absent, PAS 0 comme pour MAF -- voir
    `parse_mrc_ratio`).
  - `build_samples_argument(snp_file_path)` -- retourne la taille
    HAPLOÏDE du pool par population (cf.
    `count_samples_per_population`/`_parse_pool_header_line`) --
    divisée par 2 ici pour obtenir un compte d'INDIVIDUS diploïdes
    (voir ci-dessus) ; utilisée TELLE QUELLE (non divisée) partout
    ailleurs, notamment comme `pool_sizes` dans
    `summary_statistics.py` (la correction de biais de lecture Q1
    a besoin du vrai `HAPLOID_SAMPLE_SIZE`, pas de sa moitié).
  - `observed_reads(snp_file_path)` -- les lectures RÉELLEMENT
    observées par locus/population, ensuite retraduites vers les
    noms de population msprime (`"pop1"`, `"pop2"`...) via
    `_reindex_reads_by_msprime_name` (les noms réels du fichier .snp
    n'ont aucune raison de coïncider avec cette convention -- voir
    son docstring). Tronquées aux `num_loci` premières entrées :
    c'est CETTE couverture réelle, fixe par emplacement de locus,
    qui sert de paramètre `n` au tirage binomial dans
    `simulate_poolseq_reads`, jamais retirée au hasard (voir
    `with_mrc_filter`/`simulate_poolseq_reads`).
  - `with_mrc_filter(..., ploidy=2)` -- simulation + rejet-et-
    resimule si le critère MRC (min des reads dérivés/ancestraux,
    toutes populations combinées) n'est pas atteint.

Args:
    demography: La démographie de base (PAS encore rescalée --
        comme pour `<A>` en IndSeq, aucun rescale n'est
        nécessaire ici).
    snp_file_path: Chemin du fichier .snp (doit être POOLSEQ).
    seed: La graine de la simulation.
    num_loci: Le nombre de loci à simuler.
    observed_reads_per_locus: Si `None` (défaut), calculé via
        `prepare_poolseq_observed_reads(snp_file_path, num_loci)`.
        Sinon, utilisé tel quel.

Returns:
    Un itérateur de `num_loci` dicts {nom_population:
    (nreads_dérivé, nreads_total)}, tous au-dessus du seuil MRC.

---

## `build_transition_matrix`

### Signature

```python
build_transition_matrix(name_model: str, kappas: tuple[float, float], frequences_by_locus: dict[str, float]) -> np.ndarray
```

### Description

Calcule la matrice de transition (matQ) pour un modèle donné et un locus.

Args:
    name_model: "JK", "K2P", "HKY" ou "TN".
    kappas: (k1, k2) -- k2 ignoré pour "JK"/"K2P"/"HKY".
    frequences_by_locus: Dict {"pi_A": ..., "pi_C": ..., "pi_G":
        ..., "pi_T": ...}, les fréquences de bases observées de ce
        locus (voir observed_data.base_frequency_by_locus).

Returns:
    La matrice 4x4 (ordre A/C/G/T) des probabilités de
    transition, chaque ligne normalisée à somme 1, diagonale nulle.

Raises:
    NotImplementedError: Si name_model est inconnu.

---

## `build_rate_map`

### Signature

```python
build_rate_map(mutsit: list[float], mus_rate: float, dnalength: int) -> msprime.RateMap
```

### Description

Construit le profil de taux de mutation par site d'un locus (msprime.RateMap).

`mus_rate` est un taux moyen PAR SITE (donc le taux total attendu
sur tout le locus est `mus_rate * dnalength`) ; `mutsit` répartit
ce budget total entre les sites (poids relatifs qui somment à 1,
voir parameter_sampling.sample_site_rates -- un site invariant a un
poids de 0). Le taux absolu d'un site donné est donc `mus_rate *
dnalength * mutsit[site]` : on distribue le taux total du locus
site par site, proportionnellement à son poids relatif.

Args:
    mutsit: Le poids de mutation relatif par site, longueur
        dnalength, normalisé à somme 1.
    mus_rate: Le taux de mutation moyen PAR SITE du locus.
    dnalength: La longueur du locus.

Returns:
    La msprime.RateMap correspondante (un taux absolu par site).

Raises:
    ValueError: Si len(mutsit) != dnalength.

---

## `count_loci_per_group`

### Signature

```python
count_loci_per_group(list_loci: list[LociDescriptionDetailed]) -> dict[str, int]
```

### Description

Compte le nombre de loci par groupe.

Args:
    list_loci: La liste des loci détaillés.

Returns:
    Un dict {nom_groupe: nombre_de_loci}.

Raises:
    ValueError: Si un même groupe mélange plusieurs types de loci
        ("M" et "S"), non supporté.

---

## `build_sex_stratified_samples_argument_dna`

### Signature

```python
build_sex_stratified_samples_argument_dna(mss_file_path: str, list_loci: list[LociDescriptionDetailed], locus_name: str) -> list[msprime.SampleSet]
```

### Description

Construit l'argument `samples` pour simulate_independent_loci, stratifié par sexe.

Pour les loci de type "X", on doit échantillonner les individus en fonction de leur sexe.
Cette fonction lit le fichier .snp et construit la liste des SampleSet correspondants.

Args:
    mss_file_path: Chemin du fichier .mss.

Returns:
    Une liste de msprime.SampleSet, stratifiée par sexe.

---

## `build_male_only_samples_argument_dna`

### Signature

```python
build_male_only_samples_argument_dna(mss_file_path: str, list_loci: list[LociDescriptionDetailed], locus_name: str) -> dict[str, int]
```

### Description

Construit l'argument `samples` pour simulate_independent_loci, pour les individus mâles uniquement.

Args:
    mss_file_path: Chemin du fichier .mss.
    list_loci: La liste des loci détaillés.
    locus_name: Le nom du locus pour lequel construire l'argument `samples`.

Returns:
    Un dictionnaire {nom_population: nombre_d_individus_mâles} pour le

---

## `build_group_local_param_per_locus`

### Signature

```python
build_group_local_param_per_locus(header_text: str, seed: int) -> dict[str, tuple[float, float, float]]
```

### Description

Tire k1/k2/mus_rate par locus (hiérarchie à deux niveaux, groupe puis locus).

Pour chaque groupe `[S]` de header.txt : `draw_group_parameter_values`
donne la valeur moyenne par groupe (premier niveau), puis
`sampling_group_local_param` en dérive une valeur par locus (second
niveau, dispersion optionnelle autour de la moyenne du groupe).

Args:
    header_text: Texte complet de header.txt.
    seed: La graine du tirage.

Returns:
    Un dict {nom_locus: (k1, k2, mus_rate)} -- triplet à arité
    fixe pour chaque locus `[S]`, quel que soit le modèle de
    substitution utilisé par son groupe (0.0 pour les kappas non
    utilisés).

---

## `build_matrix_per_locus`

### Signature

```python
build_matrix_per_locus(header_text: str, mss_file_path: str | Path, seed: int) -> dict[str, np.ndarray]
```

### Description

Construit la matrice de transition (matQ) de chaque locus [S].

Pipeline complet header.txt + .mss + seed -> {nom_locus: matQ},
en composant build_group_local_param_per_locus (k1/k2 par locus),
base_frequency_by_locus (pi par locus) et build_transition_matrix.

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    seed: La graine du tirage.

Returns:
    Un dict {nom_locus: matQ} pour chaque locus [S].

---

## `build_rate_map_per_locus`

### Signature

```python
build_rate_map_per_locus(header_text: str, seed: int) -> dict[str, msprime.RateMap]
```

### Description

Construit le profil de taux de mutation (msprime.RateMap) de chaque locus [S].

Args:
    header_text: Texte complet de header.txt.
    seed: La graine du tirage.

Returns:
    Un dict {nom_locus: msprime.RateMap} pour chaque locus [S].

---

## `simulate_dna_mutations`

### Signature

```python
simulate_dna_mutations(tree_sequence: tskit.TreeSequence, transition_matrix: np.ndarray, frequencies: dict[str, float], rate_map: msprime.RateMap, seed: int) -> tskit.TreeSequence
```

### Description

Place les mutations sur une généalogie ADN via msprime.sim_mutations.

Args:
    tree_sequence: La généalogie non mutée (un locus).
    transition_matrix: La matrice de transition (matQ) du locus.
    frequencies: Dict {"pi_A": ..., "pi_C": ..., "pi_G": ...,
        "pi_T": ...}, distribution ancestrale à la racine.
    rate_map: Le profil de taux de mutation par site du locus.
    seed: La graine du tirage de mutation.

Returns:
    La TreeSequence mutée.

---

## `dna_ancestry_parameters_for_heritage`

### Signature

```python
dna_ancestry_parameters_for_heritage(heritage: str, demography: msprime.Demography, sex_ratio: float) -> tuple[msprime.Demography, int]
```

### Description

Détermine la démographie (rescalée ou non) et la ploïdie pour un locus ADN.

Reproduit le même dispatch que `simulate_genotypes_for_locus_type`
côté SNP : "A" utilise la démographie <A> telle quelle en
ploidy=2 ; "H"/"M"/"X"/"Y" la rescalent par
`coalescence_coefficient(heritage, sex_ratio) / 2` en ploidy=1
(mêmes coefficients, mêmes formules, aucune raison structurelle
qu'ils diffèrent entre SNP et séquences ADN -- comp_matQ/
put_mutations opèrent sur le même arbre msprime que le Hudson SNP,
seule la mutation diffère).

Args:
    heritage: "A", "H", "M", "X" ou "Y".
    demography: La démographie <A> de base (PAS encore rescalée).
    sex_ratio: Fraction de mâles (voir observed_data.parse_sex_ratio).

Returns:
    Le tuple (demography, ploidy) à utiliser pour ce locus.

Raises:
    NotImplementedError: Pour tout autre type d'héritage inconnu.

---

## `dna_mutation_simulation_per_locus`

### Signature

```python
dna_mutation_simulation_per_locus(header_text: str, mss_file_path: str | Path, demography: msprime.Demography, seed: int) -> dict[str, tskit.TreeSequence]
```

### Description

Assemble le pipeline complet de simulation ADN, par locus.

Pour chaque locus [S] : généalogie (msprime.sim_ancestry direct,
pas simulate_independent_loci -- sequence_length variable par
locus) + mutation (matQ/RateMap déjà construits).

La démographie et la ploïdie utilisées pour l'arbre de coalescence de
chaque locus dépendent de son type d'héritage (<A>/<H>/<M>/<X>/<Y>, voir
`dna_ancestry_parameters_for_heritage`) -- un groupe peut mélanger des
loci de types différents (ex. toy_example2_ms_dna : G2 <A>, G3 <M>),
donc ce dispatch se fait par locus, jamais une fois pour tout le
dataset.

Pour les loci [S] de type <A>, on tire une graine différente pour chaque
locus (seed + _ANCESTRY_SEED_OFFSET + i), pour que chaque locus <A>
ait sa propre généalogie indépendante tandis que pour les autres loci,
on utilise la graine (seed + _SHARED_M/Y_ANCESTRY_SEED_OFFSET) pour que tous les
loci M ou Y partagent la même généalogie.

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    demography: La démographie <A> de base (PAS encore rescalée).
    seed: La graine de la simulation.

Returns:
    Un dict {nom_locus: TreeSequence mutée} pour chaque locus [S].

---

## `_group_prior_values_from_columns`

### Signature

```python
_group_prior_values_from_columns(group_priors_values: dict[str, float], group_priors: dict) -> dict[str, dict[str, float]]
```

### Description

Reconstruit le dict nested {groupe: {prior: valeur}} depuis les colonnes du reftable réel.

Reshape les colonnes plates du vrai reftable DIYABC (ex:
"µseq_2", "k1seq_2") dans la forme nested que
draw_group_parameter_values produit normalement, pour que
build_group_local_param_per_locus_from_values puisse réutiliser
tel quel le corps de build_group_local_param_per_locus.

Args:
    group_priors_values: Dict {nom_colonne: valeur} tel que lu
        dans le vrai reftable (voir
        reftable_loop.parse_real_reftable_params_with_group_priors).
    group_priors: Dict {nom_groupe: [GroupPrior, ...]} (voir
        prior_parser.parse_group_priors).

Returns:
    Un dict {nom_groupe: {nom_prior: valeur}}, pour les groupes de
    loci ADN ([S]) uniquement -- les groupes MicroSat sont ignorés.

---

## `build_group_local_param_per_locus_from_values`

### Signature

```python
build_group_local_param_per_locus_from_values(header_text: str, group_priors_values: dict[str, float], seed: int) -> dict[str, tuple[float, float, float]]
```

### Description

Variante replay de build_group_local_param_per_locus (tier 1 = valeurs réelles).

Le tirage par-groupe (premier niveau) est remplacé par les valeurs
réellement tirées par DIYABC (`group_priors_values`, via
`_group_prior_values_from_columns`) ; le tirage par-locus (second
niveau, dispersion autour de la moyenne) N'EST PAS remplacé -- il
continue de dépendre de `seed`, car DIYABC n'enregistre pas cette
dispersion dans le reftable, il n'y a donc rien à rejouer pour elle.

Args:
    header_text: Texte complet de header.txt.
    group_priors_values: Dict {nom_colonne: valeur} tel que lu
        dans le vrai reftable.
    seed: La graine du tirage par-locus (second niveau).

Returns:
    Un dict {nom_locus: (k1, k2, mus_rate)}, même contrat que
    build_group_local_param_per_locus.

---

## `build_matrix_per_locus_from_values`

### Signature

```python
build_matrix_per_locus_from_values(header_text: str, mss_file_path: str | Path, group_priors_values: dict[str, float], seed: int) -> dict[str, np.ndarray]
```

### Description

Variante replay de build_matrix_per_locus (premier niveau = valeurs réelles).

Ne tire PAS les k1/k2 moyens par groupe (premier niveau, voir
build_group_local_param_per_locus_from_values) : elle réutilise des
valeurs déjà connues, typiquement les tirages RÉELS d'un reftable
DIYABC existant (voir reftable_loop.
parse_real_reftable_params_with_group_priors) -- pour comparer DIYABC
et msprime sur EXACTEMENT les mêmes valeurs de k1/k2 par groupe, sans
le biais possible de deux tirages indépendants.

Le tirage par-locus (second niveau, la dispersion de k1/k2 autour de
la moyenne du groupe, via sampling_group_local_param à l'intérieur de
build_group_local_param_per_locus_from_values) N'EST PAS remplacé --
il continue d'être tiré depuis `seed`, car cette valeur par-locus
n'est jamais enregistrée dans le vrai reftable DIYABC (seule la
moyenne de groupe l'est) : il n'y a donc rien à rejouer pour lui.

Sinon identique à build_matrix_per_locus (même construction de
matrice de transition via build_transition_matrix, mêmes fréquences
de bases observées).

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    group_priors_values: Dict {nom_colonne: valeur} tel que lu
        dans le vrai reftable.
    seed: La graine du tirage par-locus (second niveau).

Returns:
    Un dict {nom_locus: matQ} pour chaque locus [S], même contrat
    que build_matrix_per_locus.

---

## `build_rate_map_per_locus_from_values`

### Signature

```python
build_rate_map_per_locus_from_values(header_text: str, group_priors_values: dict[str, float], seed: int) -> dict[str, msprime.RateMap]
```

### Description

Variante replay de build_rate_map_per_locus (premier niveau = valeurs réelles).

Ne tire PAS le mus_rate moyen par groupe (premier niveau) :
réutilise group_priors_values, comme
build_matrix_per_locus_from_values -- même principe, voir sa
docstring pour le détail complet.

Le tirage de `mutsit` (sample_site_rates, la dispersion du taux de
mutation par SITE au sein d'un locus) N'EST PAS remplacé -- il
continue de dépendre de `seed` (rng = random.Random(seed +
_SITE_RATE_SEED_OFFSET)), pour la même raison que le tirage
par-locus de k1/k2 : `mutsit` n'est jamais enregistré dans le vrai
reftable DIYABC, il n'y a donc rien à rejouer pour lui non plus.

Args:
    header_text: Texte complet de header.txt.
    group_priors_values: Dict {nom_colonne: valeur} tel que lu
        dans le vrai reftable.
    seed: La graine du tirage par-locus (second niveau) et de
        `mutsit`.

Returns:
    Un dict {nom_locus: msprime.RateMap} pour chaque locus [S],
    même contrat que build_rate_map_per_locus.

---

## `dna_mutation_simulation_per_locus_from_values`

### Signature

```python
dna_mutation_simulation_per_locus_from_values(header_text: str, mss_file_path: str | Path, demography: msprime.Demography, group_priors_values: dict[str, float], seed: int) -> dict[str, tskit.TreeSequence]
```

### Description

Variante replay de dna_mutation_simulation_per_locus (rejeu apparié DIYABC/msprime).

Voir build_matrix_per_locus_from_values pour le principe général :
appelle build_rate_map_per_locus_from_values/
build_matrix_per_locus_from_values au lieu des originales, pour que
les k1/k2/mus_rate moyens par groupe soient ceux RÉELLEMENT tirés
par DIYABC (group_priors_values) plutôt que tirés à nouveau depuis
`seed`.

`demography` reste un paramètre déjà construit par l'appelant, comme
dans la version d'origine -- rien à changer ici pour les paramètres
historiques (N1, ta, ts...), leur propre rejeu "from values" est géré
en amont par pipeline.build_demography_for_scenario_index, réutilisée
telle quelle (générique, ne sait rien de SNP vs ADN).

Pour les loci [S] de type <A>, on tire une graine différente pour chaque
locus (seed + _ANCESTRY_SEED_OFFSET + i), pour que chaque locus <A>/<H>/<X>
ait sa propre généalogie indépendante tandis que pour les loci mitochondriaux,
on utilise la graine (seed + _SHARED_M/Y_ANCESTRY_SEED_OFFSET) pour que tous
les loci M ou Y partagent la même généalogie.

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    demography: La démographie <A> de base (déjà construite par
        l'appelant à partir des valeurs historiques réelles).
    group_priors_values: Dict {nom_colonne: valeur} tel que lu
        dans le vrai reftable.
    seed: La graine du tirage par-locus (second niveau), de la
        généalogie et de la mutation.

Returns:
    Un dict {nom_locus: TreeSequence mutée} pour chaque locus [S],
    même contrat que dna_mutation_simulation_per_locus.

---

## `build_microsat_transition_matrix`

### Signature

```python
build_microsat_transition_matrix(kmin: int, kmax: int, motif_size: int, Pgeom: float, epsilon: float) -> msprime.MatrixMutationModel
```

### Description

Construit la matrice de transition (matQ) d'un locus microsatellite [M].

Args:
    kmin: Nombre minimum d'allèles du locus.
    kmax: Nombre maximum d'allèles du locus.
    motif_size: Taille du motif répétitif du microsatellite.
    Pgeom: Paramètre de mutation du modèle SMM.
    epsilon: Paramètre de précision pour éviter les valeurs nulles.

Returns:
    Modèle de mutation (alleles, root_distribution, transition_matrix) du locus microsatellite.

---

## `build_microsat_local_param_per_locus`

### Signature

```python
build_microsat_local_param_per_locus(header_text, seed) -> dict[str, tuple[float, float]]
```

### Description

Construit les paramètres locaux (mut_rate, Pgeom) de chaque locus [M]. A terme on rajoutera
aussi sni_rate.

Args:
    header_text: Texte complet de header.txt.
    seed: La graine du tirage.

Returns:
    Un dict {nom_locus: (mut_rate,Pgeom)} pour chaque locus microsat.

---

## `build_matrix_microsat_per_locus`

### Signature

```python
build_matrix_microsat_per_locus(header_text: str, mss_file_path: str, seed: int) -> dict[str, msprime.MatrixMutationModel]
```

### Description

Construit la matrice de transition de chaque locus microsatellite [M].

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    seed: La graine du tirage.

Returns:
    Un dict {nom_locus: msprime.MatrixMutationModel} pour chaque locus microsatellite [M].

---

## `microsat_mutation_simulation_per_locus`

### Signature

```python
microsat_mutation_simulation_per_locus(header_text: str, mss_file_path: str | Path, demography: msprime.Demography, seed: int) -> dict[str, tskit.TreeSequence]
```

### Description

Assemble le pipeline complet de simulation microsatellite, par locus.

Pour chaque locus [M] : généalogie (msprime.sim_ancestry direct,
pas simulate_independent_loci -- sequence_length variable par
locus) + mutation (matQ déjà construit).

Args:
    header_text: Texte complet de header.txt.
    mss_file_path: Chemin du fichier .mss.
    demography: La démographie <A> de base (PAS encore rescalée).
    seed: La graine de la simulation.
Returns:
    Un dict {nom_locus: arbre_généalogique} pour chaque locus microsatellite [M].

---


# 📄 demography_builder.py

## `evaluate_expression`

### Signature

```python
evaluate_expression(expr: str, values: dict[str, float]) -> float
```

### Description

Évalue une expression de temps ou de taille de header.txt de manière récursive.

Un nombre littéral ("0"), un nom de paramètre tiré ("t1"), ou une
somme/différence de deux noms ("t2-d3", "t2+d3"). Équivalent de
ParticleC::getvalue() en C++.

Args:
    expr: L'expression à évaluer.
    values: Les valeurs déjà tirées des priors, {nom: valeur}.

Returns:
    La valeur numérique de l'expression.

Raises:
    ValueError: Si expr n'est ni un nombre littéral, ni un nom de
        paramètre connu, ni une somme/différence de tels noms.

---

## `build_demography`

### Signature

```python
build_demography(scenario: Scenario, values: dict[str, float]) -> msprime.Demography
```

### Description

Construit la Demography msprime correspondant au scenario.

Les populations sont nommées "pop1", "pop2", ... d'après leur indice
dans header.txt (1-indexed, comme dans le fichier).

Args:
    scenario: Le scénario parsé (header_dataclasses.Scenario).
    values: Les valeurs numériques déjà tirées des priors, {nom:
        valeur}.

Returns:
    La Demography msprime correspondante, avec ses événements triés
    par temps croissant.

Raises:
    NotImplementedError: Si scenario contient un type d'événement
        non géré.

---

## `rescale_demography`

### Signature

```python
rescale_demography(demography: msprime.Demography, factor: float) -> msprime.Demography
```

### Description

Rescale toutes les tailles de population d'une Demography par factor.

Retourne une COPIE de demography -- ne modifie jamais l'original
(une même Demography <A> est réutilisée telle quelle pour bâtir la
version rescalée <X>/<Y>/<M> du même scénario/tirage de valeurs,
donc muter en place casserait les autres types de locus qui
doivent encore utiliser la version non rescalée).

Deux choses sont rescalées :
- demography.populations : chaque Population a un .initial_size
  (toujours défini, jamais None).
- demography.events : SEULS les événements PopulationParametersChange
  (générés par add_population_parameters_change, donc par les
  VarNeEvent) ont un .initial_size -- et il peut valoir None (un
  PopulationParametersChange peut ne changer QUE growth_rate, pas la
  taille -- pas notre cas ici, mais à ne pas casser si ça arrive).
  Les autres types d'événements (PopulationSplit, Admixture...) n'ont
  pas d'attribut .initial_size du tout.

Args:
    demography: La Demography à rescaler (jamais modifiée).
    factor: Le facteur multiplicatif, à appeler avec
        coalescence_coefficient(locus_type, sex_ratio) / 2 (voir
        observed_data.py -- le /2 vient de la conversion
        coeffcoal*N/2 = nombre effectif de copies de gène). Le
        facteur est calculé par l'appelant.

Returns:
    La copie rescalée de demography.

---

## `extract_referenced_names`

### Signature

```python
extract_referenced_names(expr: str) -> set[str]
```

### Description

Extrait les noms de paramètres référencés par une expression de header.txt.

SANS l'évaluer numériquement -- "t2-d3" -> {"t2","d3"}, "t1" ->
{"t1"}, "0" -> set() (un nombre littéral ne référence aucun
paramètre). Utilisé pour déterminer quels paramètres un scénario
utilise réellement (nécessaire pour filtrer les colonnes du
reftable.bin par scénario -- voir reftable_loop.write_reftable_bin
et notes/exploration.md, bug "21 vs 16 paramètres pour le
scénario 1").

Args:
    expr: L'expression à analyser.

Returns:
    L'ensemble des noms de paramètres référencés (vide pour un
    nombre littéral).s S_A_11

---

## `get_parameter_names_used_by_scenario`

### Signature

```python
get_parameter_names_used_by_scenario(scenario: Scenario) -> set[str]
```

### Description

Collecte les noms de paramètres réellement référencés par un scénario.

    Tailles de population initiales, et time_expr / new_size_expr de
    chacun de ses événements. C'est ce sous-ensemble (pas la totalité
    des priors déclarés dans header.txt) qui doit constituer les
    colonnes param[] du reftable.bin pour ce scénario -- un scénario
    donné n'utilise généralement qu'une partie des priors globaux (ex:
    human/header.txt a 21 priors déclarés, mais le scénario 1 n'en
    référence que 16 -- ra/t11/t22/t33/t44 appartiennent aux scénarios
    2-6, pas au scénario 1).
s S_A_11
    Args:
        scenario: Le scénario parsé (header_dataclasses.Scenario).

    Returns:
        L'ensemble des noms de paramètres utilisés par ce scénario.

---


# 📄 loci_parser.py

## `_find_loci_description_section_index`

### Signature

```python
_find_loci_description_section_index(lines: list[str]) -> int
```

### Description

Repère l'index de la ligne d'en-tête 'loci description (N)'.

Factorisé entre parse_loci_description et rewrite_loci_count, qui en
ont toutes deux besoin.

Args:
    lines: Les lignes de header.txt.

Returns:
    L'index de la ligne d'en-tête.

Raises:
    ValueError: Si la section n'est trouvée nulle part.

---

## `_extract_loci_info_condensed`

### Signature

```python
_extract_loci_info_condensed(text: str) -> tuple[dict[str, int], str, str]
```

### Description

Découpe une ligne de contenu condensée en ses trois composantes.

Une ou plusieurs paires '<n> <type>' suivies de '<groupe> from
<indice>' (single-type ou multi-type indifféremment).

Args:
    text: La ligne de contenu (après la ligne d'en-tête de section).

Returns:
    Un tuple (loci_counts_by_heritage, group, start_1based), où
    loci_counts_by_heritage vaut {type_héritage: compte, ...}.

Raises:
    NotImplementedError: Si le format ne correspond pas.

---

## `parse_loci_description`

### Signature

```python
parse_loci_description(header_text: str) -> LociDescription | list[LociDescriptionDetailed]
```

### Description

Extrait la description des loci à partir de header.txt.

Gère le format condensé (single-type comme human, ou multi-type
comme toy_example5) et le format détaillé un-locus-par-ligne observé
dans sequences-mut.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Un LociDescription pour le format condensé, ou une liste de
    LociDescriptionDetailed (une par locus) pour le format détaillé.

Raises:
    NotImplementedError: Si le format détecté n'est ni l'un ni
        l'autre de ces deux formats.
    ValueError: Propagée si la section 'loci description' est
        introuvable.

---

## `rewrite_loci_count`

### Signature

```python
rewrite_loci_count(header_text: str, new_total_loci: int) -> str
```

### Description

Remplace le nombre de loci déclaré dans 'loci description'.

Nécessaire pour tester avec un nombre de loci réduit sans avoir à
maintenir un header.txt séparé à la main. Contrairement à
parse_loci_description, limité au format condensé à un SEUL type
d'héritage -- pas encore mis à jour pour le multi-type (ex:
toy_example5), non nécessaire pour human.

Args:
    header_text: Texte complet de header.txt.
    new_total_loci: Le nouveau nombre total de loci à déclarer.

Returns:
    Une copie de header_text avec le nombre de loci remplacé.

Raises:
    NotImplementedError: Si la section est au format détaillé
        (multi-ligne), ou si la ligne condensée n'est pas au format
        à un seul type d'héritage.

---


# 📄 observed_data.py

## `_find_header_index`

### Signature

```python
_find_header_index(lines: list[str]) -> int
```

### Description

Repère l'index de la ligne d'en-tête 'IND SEX POP' ou 'POOL'.

Recherchée parmi les deux premières lignes du fichier --
factorisé entre count_samples_per_population et
individual_sexes_per_population, qui en ont toutes deux besoin.
L'en-tête peut être précédé ou non d'un commentaire libre en
première ligne (ex: '<NM=1NF> <MAF=hudson> ...', comportement
observé dans data.cpp, qui teste les deux cas) : on recherche son
index plutôt que de supposer sa position, pour ne perdre aucune
ligne de données quel que soit le cas.

Args:
    lines: Les lignes du fichier .snp.

Returns:
    L'index de la ligne d'en-tête.

Raises:
    ValueError: Si l'en-tête n'est trouvé dans aucune des deux
        premières lignes.

---

## `detect_snp_file_type`

### Signature

```python
detect_snp_file_type(snp_file_path: str | Path) -> str
```

### Description

Détecte le type de fichier .snp DIYABC.

À l'aide du header_index trouvé par _find_header_index, en lisant
la première ligne non vide du fichier contenant "IND" ou "POOL".

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    "IND" (individus par ligne, IndSeq) ou "POOL" (pools par
    ligne, PoolSeq).

---

## `_parse_pool_header_line`

### Signature

```python
_parse_pool_header_line(lines: list[str], header_index: int) -> dict[str, int]
```

### Description

Parse la ligne d'en-tête POOL du fichier .snp DIYABC.

Format : 'POOL POP_NAME:HAPLOID_SAMPLE_SIZE  POP1:200 POP2:200
POP3:200 POP4:200'.

Args:
    lines: Les lignes du fichier .snp.
    header_index: L'index de la ligne d'en-tête POOL.

Returns:
    Un dict {nom_population: taille_haploïde}, ex: pour
    toy_example4 -> {"POP1": 200, "POP2": 200, "POP3": 200,
    "POP4": 200}.

Raises:
    ValueError: Si aucune population n'est déclarée dans la ligne.

---

## `count_samples_per_population`

### Signature

```python
count_samples_per_population(snp_file_path: str | Path) -> dict[str, int]
```

### Description

Compte le nombre d'individus par population dans un fichier .snp DIYABC.

Format 'IND SEX POP <génotypes...>' ou 'POOL
POP_NAME:HAPLOID_SAMPLE_SIZE'.

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

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Un dict {nom_population: effectif}, ex: pour human ->
    {"ASW": 30, "YRI": 30, ...} ou pour toy_example4 ->
    {"POP1": 200, "POP2": 200, "POP3": 200, "POP4": 200}.

---

## `individual_sexes_per_population`

### Signature

```python
individual_sexes_per_population(snp_file_path: str | Path) -> dict[str, list[str]]
```

### Description

Lit le sexe de chaque individu, regroupé par population.

Fichier .snp DIYABC au format 'IND SEX POP <génotypes...>'.
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

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Un dict {nom_population: [sexe, ...]}, ex: pour toy_example5 ->
    {"pop1": ["M", "F", "F", ...], ...}.

---

## `parse_sex_ratio`

### Signature

```python
parse_sex_ratio(snp_file_path: str | Path) -> float
```

### Description

Lit le sex-ratio déclaré en tête de fichier .snp DIYABC.

Format '<NM=xNF> ...' où x = nombre de mâles / nombre de femelles
(ex: '<NM=0.428571NF>' pour toy_example5). Reproduit exactement
DataC::readfile (data.cpp:475-486). Comme dans le C++, retombe sur
0.5 si le token '<NM=' est absent de la première ligne -- et,
contrairement à _find_header_index, ne cherche JAMAIS que la ligne
0 : le C++ ne considère que la toute première ligne du fichier,
jamais 'IND SEX POP' elle-même.

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    sexratio = x / (1+x) (fraction de mâles).

---

## `parse_maf_ratio`

### Signature

```python
parse_maf_ratio(snp_file_path: str | Path) -> float
```

### Description

Lit le Minimum Allele Frequency (MAF) déclaré en tête de fichier .snp DIYABC.

Format '<MAF=xxx> ...' où xxx = MAF (ex: '<MAF=hudson>' ou
'<MAF=0.05>'). Reproduit exactement DataC::readfile
(data.cpp:475-497) : 0.0 si le token '<MAF=' est absent de la
première ligne, ou si xxx n'est pas numérique (ex: 'hudson'),
comme le ferait atof() en C++ -- 0.0 veut dire "pas de filtre"
(algorithme de Hudson standard), jamais distingué du cas "MAF=0%"
explicite, exactement comme côté C++. Contrairement à
_find_header_index, ne cherche JAMAIS que la ligne 0 : le C++ ne
considère que la toute première ligne du fichier.

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Le seuil MAF.

---

## `parse_mrc_ratio`

### Signature

```python
parse_mrc_ratio(snp_file_path: str | Path) -> float
```

### Description

Lit le Minimum Read Count (MRC) déclaré en tête de fichier .snp DIYABC.

Format '<MRC=xxx> ...' où xxx = MRC. Reproduit exactement
DataC::readfile (data.cpp:498-508) : 1 si le token '<MRC=' est
absent de la première ligne, ou si xxx n'est pas numérique (ex:
'hudson'), comme le ferait atof() en C++ -- 1 veut dire "pas de
filtre".

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Le seuil MRC.

---

## `population_index_to_name`

### Signature

```python
population_index_to_name(snp_file_path: str | Path) -> dict[int, str]
```

### Description

Construit le mapping entre l'indice de population de header.txt et son nom réel.

header.txt est 1-indexed (pop1, pop2, ...) ; le nom réel est celui
qui apparaît dans le fichier .snp (ex: "ASW", "YRI"...). Voir la
docstring de count_samples_per_population pour la justification de
ce mapping par ordre d'apparition (header.txt ne nomme jamais les
populations).

Args:
    snp_file_path: Chemin du fichier .snp.

Returns:
    Un dict {indice_1based: nom_population}, ex: {1: "ASW",
    2: "YRI", 3: "CHB", 4: "GBR"} pour human et {1: "POP1",
    2: "POP2", 3: "POP3", 4: "POP4"} pour toy_example4.

---

## `observed_mrc`

### Signature

```python
observed_mrc(reads_by_population: dict[str, tuple[int, int]]) -> float
```

### Description

Calcule le MRC observé pour un locus donné.

Reproduit exactement DataC::purgelocMRCPOOLSEQ (data.cpp:1087-1093).

Args:
    reads_by_population: Dict {nom_population: (nreads1,
        nreads1+nreads2)}.

Returns:
    min(somme reads allèle1, somme reads allèle2) TOUTES
    populations combinées.

---

## `observed_reads`

### Signature

```python
observed_reads(snp_file_path: str | Path, num_loci: int | None) -> list[dict[str, tuple[int, int]]]
```

### Description

Lit les lignes de génotypes du fichier .snp DIYABC POOLSEQ.

Ignore l'en-tête et les lignes vides. Chaque tuple contient le
nombre de lectures pour l'allèle 1 et le nombre total de lectures
(allèle 1 + allèle 2) pour cette population.

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

Args:
    snp_file_path: Chemin du fichier .snp (doit être POOLSEQ).
    num_loci: Si non None, limite le nombre de loci lus à ce
        nombre (utile pour toy_example4, où on ne veut que les 100
        premiers loci passants le seuil MRC, pour reproduire
        exactement le statobs.txt de DIYABC). Si None, lit tous
        les loci du fichier.

Returns:
    La liste des lignes de comptage de reads par population, dans
    l'ordre d'apparition des populations dans le fichier -- chaque
    ligne de la forme {"POP1": (nreads1, nreads1+nreads2),
    "POP2": (nreads1, nreads1+nreads2), ...}.

Raises:
    ValueError: Si le fichier n'est pas au format POOLSEQ.

---

## `coalescence_coefficient`

### Signature

```python
coalescence_coefficient(locus_type: str, sex_ratio: float) -> float
```

### Description

Calcule le coefficient de coalescence en fonction du type d'héritage.

Reproduit DataC::cal_coeffcoal (data.cpp:1589-1605) : le
coefficient qui rescale N (taille de population, prior déclaré dans
header.txt) dans la formule de temps de coalescence
(particuleC.cpp:1340 : temps -= coeffcoal * N / n / (n-1) * log(ra)),
en fonction du type d'héritage du locus et du sex-ratio du dataset.

À l'appel avec n=2 (une paire de lignées), coeffcoal*N/2 donne le
nombre EFFECTIF de copies de gène dans la population pour ce type de
locus -- ex: <A> à sex_ratio=0.5 -> 2N (diploïde classique), <X> ->
1.5N (3/4 de 2N), <Y>/<M> -> 0.5N chacun (1/4 de 2N). Cette fonction
ne fait PAS cette conversion : elle retourne coeffcoal brut, comme le
C++, pour rester directement comparable à cal_coeffcoal.

Args:
    locus_type: Vocabulaire de
        LociDescription.loci_counts_by_heritage ("A", "H", "X",
        "Y" ou "M" -- voir loci_parser.py), PAS "<A>" avec les
        chevrons.
    sex_ratio: Fraction de mâles, telle que retournée par
        parse_sex_ratio (0.5 = sex-ratio équilibré).

Returns:
    Le coefficient coeffcoal brut.

Raises:
    NotImplementedError: Si locus_type est inconnu.

---

## `observed_sequences`

### Signature

```python
observed_sequences(mss_file_path: str | Path, list_loci: list[LociDescriptionDetailed]) -> dict[str, list[str]]
```

### Description

Lit les séquences ADN observées dans un fichier .mss.

Args:
    mss_file_path: Chemin du fichier .mss.
    list_loci: La description détaillée des loci (voir
        loci_parser.parse_loci_description), utilisée pour
        distinguer les loci séquentiels ([S]) des loci
        microsatellites ([M]) et pour nommer les entrées du
        dictionnaire retourné.

Returns:
    Un dict {nom_locus: [séquence, ...]}. Chaque séquence est une
    chaîne de caractères représentant les bases nucléotidiques (A,
    C, G, T, N ou -) pour chaque individu ou pool -- un locus
    diploïde (<A>) produit 2 entrées par individu, un locus
    haploïde (<M>) en produit 1.

Raises:
    ValueError: Si le fichier ne contient aucune ligne 'POP' (indépendamment de la casse), ou
        si le nombre de séquences observées sur une ligne ne
        correspond pas au nombre de loci séquentiels attendu.

---

## `individual_sexes_from_locus_genotype`

### Signature

```python
individual_sexes_from_locus_genotype(mss_file_path: str | Path, list_loci: list[LociDescriptionDetailed], locus_name: str) -> dict[str, list[str]]
```

### Description

Lit le sexe de chaque individu à partir du génotype d'un locus donné.

Args:
    mss_file_path: Chemin du fichier .mss.
    list_loci: La description détaillée des loci (voir
        loci_parser.parse_loci_description), utilisée pour
        distinguer les loci séquentiels ([S]) des loci
        microsatellites ([M]) et pour nommer les entrées du
        dictionnaire retourné.
    locus_name: Le nom du locus à partir duquel déduire le sexe.

Returns:
    Un dict {nom_population: [sexe, ...]}, ex: pour toy_example5 ->
    {"pop1": ["M", "F", "F", ...], ...}.

Raises:
    ValueError: Si le fichier ne contient aucune ligne 'POP' (indépendamment de la casse), ou
        si le locus spécifié n'est pas trouvé dans la liste des loci.

---

## `observed_count_population`

### Signature

```python
observed_count_population(mss_file_path: str | Path) -> dict[str, int]
```

### Description

Compte le nombre d'individus par population dans un fichier .mss.

Args:
    mss_file_path: Chemin du fichier .mss.

Returns:
    Un dict {"pop1": effectif, "pop2": effectif, ...}.

Raises:
    ValueError: Si le fichier ne contient aucune ligne commençant par 'pop' (indépendamment de la casse).

---

## `base_frequency_by_locus`

### Signature

```python
base_frequency_by_locus(sequences_by_locus: dict[str, list[str]]) -> dict[str, dict[str, float]]
```

### Description

Calcule la fréquence des bases observées par locus.

Args:
    sequences_by_locus: Dict {nom_locus: [séquence, ...]} tel que
        retourné par observed_sequences.

Returns:
    Un dict {nom_locus: {"pi_A": ..., "pi_C": ..., "pi_G": ...,
    "pi_T": ...}}.

Raises:
    ValueError: Si une séquence contient une base inattendue
        (autre que A, C, G, T, N ou -).

---

## `observed_microsatellites`

### Signature

```python
observed_microsatellites(mss_file_path: str | Path, list_loci: list[LociDescriptionDetailed]) -> dict[str, list[list[int]]]
```

### Description

Lit les microsatellites observés dans un fichier .mss.

Args:
    mss_file_path: Chemin du fichier .mss.
    list_loci: La description détaillée des loci (voir
        loci_parser.parse_loci_description), utilisée pour
        distinguer les loci séquentiels ([S]) des loci
        microsatellites ([M]) et pour nommer les entrées du
        dictionnaire retourné.

Returns:
    Un dict {nom_locus: [[allèle1, allèle2] / [allèle1] / [allèle2] / [],  ...]} selon le
    contenu du fichier .mss.
    Chaque entrée est une liste de listes représentant les
    allèles observés pour chaque individu.
Raises:
    ValueError: Si le fichier ne contient aucune ligne 'pop', ou
        si le nombre de microsatellites observés sur une ligne ne
        correspond pas au nombre de loci microsatellites attendu.
    ValueError: Si un microsatellite n'est pas au format attendu (3 ou 6 chiffres, ex: '000' ou 'NNNNNN').

---

## `allele_bounds_per_locus`

### Signature

```python
allele_bounds_per_locus(microsatellites_by_locus: dict[str, list[list[int]]], list_loci: list[LociDescriptionDetailed]) -> dict[str, tuple[int, int]]
```

### Description

Calcule les bornes d'allèles observées pour chaque locus microsatellite.

Args:
    microsatellites_by_locus: Dict {nom_locus: [[allèle1, allèle2] / [allèle1] / [allèle2] / [],  ...]} tel que
        retourné par observed_microsatellites.
    list_loci: La description détaillée des loci (voir
                loci_parser.parse_loci_description), utilisée pour
                recupérer motif_size et motif_range pour chaque locus.

Returns:
    Un dict {nom_locus: (min_allèle, max_allèle)}.

Raises:
    ValueError: Si un locus n'est pas trouvé dans la liste des loci ou,
      si le motif_range n'est pas assez large pour contenir tous les allèles observés ou,
      si aucun allèle n'est observé pour un locus.

---


# 📄 parameter_sampling.py

## `draw_scenario`

### Signature

```python
draw_scenario(scenarios: list[Scenario], seed: int) -> Scenario
```

### Description

Tire un scénario parmi `scenarios`, pondéré par son `weight`.

Le poids est le "prior_proba" de
particuleC.cpp::ParticleC::drawscenario. Reproduit exactement
l'algorithme C++ : tirage d'un nombre uniforme `ra` dans [0,1),
puis balayage de la somme cumulée des poids jusqu'à ce qu'elle
atteigne ou dépasse `ra` -- même logique d'inversion de CDF que
_draw_single_mutation_edge_child (ancestry_simulation.py),
vérifiée boundary-compatible avec la boucle C++ ("while ra > sp").

Ne normalise PAS les poids (comme le C++, qui ne le fait pas non
plus) : si leur somme est < 1, le DERNIER scénario de la liste sert
de secours pour tout `ra` au-delà de la somme cumulée -- même
comportement que la boucle C++, bornée à nscenarios-1.

Args:
    scenarios: Les scénarios candidats.
    seed: La graine du tirage.

Returns:
    Le scénario tiré.

---

## `_draw_one_value`

### Signature

```python
_draw_one_value(prior: Prior, rng: random.Random) -> float
```

### Description

Tire une valeur pour un prior donné, selon sa loi et ses bornes.
Pour les lois normale, log-normale et gamma, on retire les valeurs hors bornes comme DIYABC.

Args:
    prior: Le prior à tirer.
    rng: Le générateur aléatoire à utiliser.

Returns:
    La valeur tirée.

Raises:
    NotImplementedError: Si la loi n'est pas encore implémentée.

---

## `_draw_one_group_value`

### Signature

```python
_draw_one_group_value(group_prior: GroupPrior, rng: random.Random) -> float
```

### Description

Tire une valeur pour un group prior donné, selon sa loi et ses bornes.

Quelques gardes-fous pour éviter des erreurs de tirage si le group
prior est mal défini.

Args:
    group_prior: Le group prior à tirer.
    rng: Le générateur aléatoire à utiliser.

Returns:
    La valeur tirée.

Raises:
    ValueError: Si le group prior n'a pas de loi ou de bornes
        associées, ou si la loi est GA mais que les priors MEAN et
        SDSHAPE correspondants n'ont pas été tirés avant.
    NotImplementedError: Si la loi n'est pas encore implémentée
        (propagée depuis _draw_one_value).

---

## `draw_parameter_values`

### Signature

```python
draw_parameter_values(priors: list[Prior], constraints: list[OrderConstraint], seed: int, max_attempts: int) -> dict[str, float]
```

### Description

Tire une valeur pour chaque prior, en retirant tant que les
contraintes d'ordre ne sont pas toutes satisfaites. On reproduit le comportement de DIYABC,
en gardant le premier tirage satisfaisant toutes les contraintes.

Args:
    priors: Les priors à tirer.
    constraints: Les contraintes d'ordre à satisfaire (ex: "t4>t3").
    seed: La graine du tirage.
    max_attempts: Le nombre maximal d'essais avant d'abandonner.

Returns:
    Un dict {nom_prior: valeur}.

Raises:
    ConstraintsNotSatisfiedError: Si aucun tirage valide n'est
        trouvé en max_attempts essais -- signe probable d'une
        configuration de contraintes incohérente (bornes de
        priors incompatibles avec les contraintes demandées)
        plutôt que d'une simple mauvaise chance.

---

## `draw_group_parameter_values`

### Signature

```python
draw_group_parameter_values(group_priors: dict[str, list[GroupPrior]], seed: int) -> dict[str, dict[str, float]]
```

### Description

Tire une valeur pour chaque group prior ou bien des valeurs pour le modèle.

`seed` est décalé de _GROUP_PRIOR_SEED_OFFSET avant utilisation -- ne
corrèle jamais ce tirage avec celui de draw_parameter_values, même si
l'appelant leur passe la même seed de base.

Point d'attention : ce tirage positionnel entre priors dépendants
(MEANMU avant GAMMU, etc.) reproduit le comportement de DIYABC, qui lit
ces lignes dans un ordre fixe sans jamais regarder leur nom (voir
header.cpp::readHeadersimGroupPrior).

Args:
    group_priors: Dict {nom_groupe: [GroupPrior, ...]}.
    seed: La graine de base du tirage (décalée en interne).

Returns:
    Un dict {nom_groupe: {nom_prior: valeur}}.

---

## `sampling_group_local_param`

### Signature

```python
sampling_group_local_param(group_prior: GroupPrior, k_moy: float, n_loci: int, check_nloc: bool, list_loci: list[LociDescriptionDetailed], rng: random.Random) -> dict[str, float]
```

### Description

Échantillonne le tirage par-locus (second niveau) d'un paramètre de groupe.

S'applique aussi bien à kappa1/kappa2 (`build_transition_matrix`)
qu'à mus_rate -- rien de spécifique à kappa dans l'implémentation.

Args:
    group_prior: Le GroupPrior déjà résolu pour ce groupe (via
        draw_group_parameter_values), dont le `mean` est remplacé
        par `k_moy` avant tirage.
    k_moy: La valeur moyenne du groupe (premier niveau), déjà
        tirée par draw_group_parameter_values.
    n_loci: Le nombre de loci du groupe.
    check_nloc: Si True, un tirage indépendant par locus n'a lieu
        que si `n_loci > 1` en plus de `sdshape > 0.001` (cas de
        kappa1/mus_rate). Si False, seule la condition sur
        `sdshape` est vérifiée (cas de kappa2 -- asymétrie propre à
        DIYABC, reproduite telle quelle).
    list_loci: Les loci du groupe.
    rng: Le générateur aléatoire à utiliser.

Returns:
    Un dict {nom_locus: valeur} -- soit un tirage indépendant par
    locus, soit `k_moy` répété pour chaque locus.

---

## `sample_site_rates`

### Signature

```python
sample_site_rates(p_fixe: float, gams: float, dnalength: int, rng: random.Random) -> list[float]
```

### Description

Tire mutsit : le taux de mutation relatif par site pour un locus séquence.

Reproduit header.cpp:707-738 (y compris le "bug" sitefix -- les
sites fixes sont toujours les premiers de la séquence, pas un
sous-ensemble aléatoire).

Args:
    p_fixe: Pourcentage (0-100, pas une fraction 0-1) de sites invariants (`GroupPrior.p_fixe`).
    gams: Forme gamma de l'hétérogénéité de taux par site
        (`GroupPrior.gams`) -- 0.0 est une valeur valide,
        équivalente à un taux uniforme (voir MwcGen::ggamma3).
    dnalength: Longueur du locus.
    rng: Le générateur aléatoire à utiliser.

Returns:
    La liste `mutsit`, de longueur `dnalength`, normalisée à
    somme 1.

---


# 📄 pipeline.py

## `read_header_text`

### Signature

```python
read_header_text(directory: Path) -> str
```

### Description

Lit header.txt si présent, sinon headerRF.txt en repli.

Les deux noms coexistent selon les jeux de données (header.txt =
config initiale fournie par l'utilisateur, headerRF.txt = variante
produite par un run DIYABC réel ; nos sous-dossiers de test n'auront
au départ que l'un des deux).

Args:
    directory: Le dossier contenant header.txt/headerRF.txt.

Returns:
    Le texte complet du fichier trouvé.

---

## `_simulate_genotypes_for_all_locus_types`

### Signature

```python
_simulate_genotypes_for_all_locus_types(demography: msprime.Demography, header_text: str, snp_path: Path) -> Iterator[dict[str, list[int]]]
```

### Description

Simule les génotypes de TOUS les types de locus déclarés dans header_text.

Boucle sur `parse_loci_description(header_text).loci_counts_by_heritage`
(dict[str, int], ex: {"A": 5000} pour human, {"A": 70, "X": 10,
"M": 10, "Y": 10} pour toy_example5), et concatène les génotypes
simulés pour chacun via simulate_genotypes_for_locus_type.

IMPORTANT -- num_loci est un compte PAR TYPE, pas un total : pour un
dataset <A>-only comme human (un seul type déclaré), c'est
rigoureusement identique au comportement actuel (num_loci loci de
type <A>, point). Pour un dataset multi-type comme toy_example5,
si l'on rentre une valeur précise de num_loci,
num_loci loci sont simulés pour CHAQUE type déclaré -- pas les vrais
comptes du header.txt (70/10/10/10) pour lesquels il faut passer
num_loci=None. Le comportement par défaut (num_loci=None) est donc de
simuler le nombre exact de loci déclaré dans header.txt pour chaque
type, ce qui est le plus souvent ce que l'on veut pour un POC ou
un test de validation.

Un seed DISTINCT est dérivé par type de locus via
_LOCUS_TYPE_SEED_OFFSET (voir bridge/configuration.py pour la
justification empirique) -- ne JAMAIS appeler
simulate_genotypes_for_locus_type avec la même seed brute pour
plusieurs types dans cette boucle.

Args:
    demography: La démographie <A> de base.
    header_text: Texte complet de header.txt.
    snp_path: Chemin du fichier .snp observé.
    num_loci: Si None (défaut), simule le nombre exact de loci
        déclaré dans header.txt pour chaque type. Sinon, ce nombre
        de loci pour CHAQUE type déclaré (voir IMPORTANT ci-dessus).
    seed: La graine de base (décalée par type de locus).

Returns:
    Un itérateur des génotypes simulés, tous types de locus
    concaténés -- pas une liste (voir simulate_independent_loci
    pour la justification : ne pas matérialiser 51250 TreeSequence
    en mémoire simultanément).

---

## `_population_names`

### Signature

```python
_population_names(genotypes_list: list[dict[str, list[int]]], snp_path: Path) -> list[str]
```

### Description

Noms de population ("pop1", "pop2"...), dans le même ordre que build_samples_argument.

Dérivés GRATUITEMENT des clés du premier locus déjà simulé
(simulate_snp_genotypes construit ce dict avec exactement les mêmes
noms, voir ancestry_simulation.compute_population_layout) plutôt
que de rescanner le fichier .snp une deuxième fois par particule
(mesuré : ~4% du temps d'une particule sur human, voir
notes/exploration.md, entrée du 20/07/2026).

Args:
    genotypes_list: Les génotypes déjà simulés (voir
        simulate_snp_genotypes), au moins un locus.
    snp_path: Chemin du fichier .snp observé -- utilisé seulement
        en repli si genotypes_list est vide.

Returns:
    La liste des noms de population, dans l'ordre. Repli sur
    build_samples_argument si genotypes_list est vide (num_loci=0,
    cas dégénéré qui n'arrive pas en pratique).

---

## `_filter_statistics`

### Signature

```python
_filter_statistics(summary_stats: dict[str, float], header_text: str, stats_filter: str) -> dict[str, float]
```

### Description

Applique stats_filter ('ALL' ou 'HEADER') à un dict de statistiques déjà calculé.

Factorisé entre compute_summary_statistics et
compute_summary_statistics_from_values (même logique de filtrage,
seule la source des valeurs de paramètres diffère entre les deux).

Args:
    summary_stats: Le dict {nom_colonne: valeur} déjà calculé.
    header_text: Texte complet de header.txt.
    stats_filter: "ALL" (retourne summary_stats tel quel) ou
        "HEADER" (ne garde que les statistiques déclarées dans la
        section 'group summary statistics' de header.txt, dans
        leur ordre de déclaration).

Returns:
    Le dict filtré.

Raises:
    ValueError: Si stats_filter="HEADER" et que header.txt déclare
        une statistique absente de summary_stats (vocabulaire
        obsolète ou non implémenté).
    NotImplementedError: Si stats_filter n'est ni "ALL" ni "HEADER".

---

## `build_random_demography`

### Signature

```python
build_random_demography(scenario: Scenario, header_text: str, seed: int) -> tuple[msprime.Demography, dict[str, float]]
```

### Description

Tire les valeurs de priors puis construit la Demography correspondante.

Toutes les valeurs de priors du fichier sont tirées (pas seulement
celles utilisées par ce scenario précis) : plus simple, et évite de
casser des contraintes d'ordre qui pourraient porter sur des
paramètres d'autres scénarios.

Args:
    scenario: Le scénario parsé (header_dataclasses.Scenario).
    header_text: Texte complet de header.txt.
    seed: La graine du tirage.

Returns:
    Le tuple (demography, values) -- les valeurs tirées sont
    renvoyées en plus de la Demography, car elles seront
    nécessaires plus tard pour écrire le reftable.bin (colonnes de
    paramètres).

---

## `build_random_demography_for_scenario_index`

### Signature

```python
build_random_demography_for_scenario_index(header_text: str, scenario_index: int, seed: int) -> tuple[msprime.Demography, dict[str, float]]
```

### Description

Variante de build_random_demography qui sélectionne le scénario par son index.

1-indexed, comme dans header.txt, plutôt que de demander un objet
Scenario déjà parsé. Utile pour les tests et l'utilisation
interactive.

Args:
    header_text: Texte complet de header.txt.
    scenario_index: L'index 1-based du scénario à utiliser.
    seed: La graine du tirage.

Returns:
    Le tuple (demography, values), même contrat que
    build_random_demography.

Raises:
    ValueError: Si scenario_index ne correspond à aucun scénario
        parsé.

---

## `run_poc_for_directory`

### Signature

```python
run_poc_for_directory(directory: str | Path, scenario_index: int) -> None
```

### Description

Point d'entrée de haut niveau : équivalent du `-p ./` de DIYABC.

Prend un dossier contenant header.txt et le fichier de données
observées (.snp), et produit les génotypes simulés sous le scénario
demandé, pour tous les types de locus déclarés.

Le nom du fichier de données est lu sur la PREMIÈRE LIGNE de
header.txt (ex: "human_snp_all22chr_maf5.snp"), pas deviné par
extension -- c'est le contrat du format DIYABC.

Args:
    directory: Le dossier contenant header.txt et le fichier .snp.
    scenario_index: L'index 1-based du scénario à utiliser.
    num_loci: Voir _simulate_genotypes_for_all_locus_types (compte
        par type, pas un total ; None = comptes réels de header.txt).
    seed: La graine de la simulation.

Returns:
    Le tuple (mutated_tree_sequences, values) : l'itérateur des
    génotypes simulés, et le dict des valeurs de paramètres tirées
    (nécessaires plus tard pour écrire le reftable.bin).

---

## `compute_summary_statistics`

### Signature

```python
compute_summary_statistics(reference_directory: str | Path, scenario_index: int) -> tuple[dict[str, float], dict[str, float]]
```

### Description

Calcule les statistiques résumées SNP/PoolSeq sur des données SIMULÉES.

Utilise nos formules Python validées (summary_statistics.py) --
remplace la délégation au binaire C++ (subprocess + fichier .snp
intermédiaire). Dispatche automatiquement entre le chemin IndSeq
(`run_poc_for_directory` + `compute_all_statistics`) et le chemin
PoolSeq (`simulate_poolseq_reads_with_mrc_filter` +
`compute_all_statistics_poolseq`) selon `detect_snp_file_type`.

Args:
    reference_directory: Le dossier contenant header.txt et le
        fichier .snp observé.
    scenario_index: L'index 1-based du scénario à utiliser.
    num_loci: Voir _simulate_genotypes_for_all_locus_types (IndSeq
        uniquement -- ignoré pour PoolSeq, qui simule toujours tous
        les loci `<A>` déclarés dans header.txt).
    seed: La graine de la simulation.
    work_directory: Gardé pour compatibilité, ignoré.
    general_binary_path: Gardé pour compatibilité, ignoré.
    stats_filter: "ALL" (défaut) : retourne toutes les statistiques
        implémentées (compute_all_statistics), sans filtrage.
        "HEADER" : ne garde, dans l'ordre de déclaration, que les
        statistiques listées dans la section 'group summary
        statistics' de header.txt (voir stats_group_parser.
        parse_requested_statistic_names) -- nécessaire pour que
        reftable_msprime.txt/.bin aient EXACTEMENT les mêmes
        colonnes que le vrai reftable DIYABC (sinon toute
        comparaison colonne-par-nom entre les deux pipelines est
        faussée, comme découvert sur toy_example5_modif :
        'ML3p_1.2.3' calculé par nous mais absent du vrai DIYABC).
    observed_reads_per_locus: PoolSeq uniquement, voir
        simulate_poolseq_reads_with_mrc_filter.

Returns:
    Le tuple (summary_statistics, parameter_values).

Raises:
    ValueError: Si stats_filter="HEADER" et que header.txt déclare
        une statistique qu'on ne sait pas calculer (vocabulaire
        obsolète, ex: human/header.txt -- voir notes/exploration.md).

---

## `build_demography_for_scenario_index`

### Signature

```python
build_demography_for_scenario_index(header_text: str, scenario_index: int, values: dict[str, float]) -> msprime.Demography
```

### Description

Variante de build_random_demography_for_scenario_index qui NE TIRE AUCUNE valeur.

Construit la Demography directement à partir de valeurs de
paramètres déjà connues (ex: reprises telles quelles d'un reftable
DIYABC réel pour servir d'oracle -- voir
reftable_loop.replay_reftable_simulation).

Args:
    header_text: Texte complet de header.txt.
    scenario_index: L'index 1-based du scénario à utiliser.
    values: Les valeurs de paramètres déjà connues, {nom: valeur}.

Returns:
    La Demography correspondante.

Raises:
    ValueError: Si scenario_index ne correspond à aucun scénario
        parsé.

---

## `run_poc_for_directory_with_values`

### Signature

```python
run_poc_for_directory_with_values(directory: str | Path, scenario_index: int, values: dict[str, float]) -> None
```

### Description

Variante de run_poc_for_directory qui prend des valeurs de paramètres déjà connues.

Au lieu d'en tirer de nouvelles -- même contrat par ailleurs
(lecture du nom de fichier .snp sur la première ligne de header.txt,
échantillonnage, simulation).

Args:
    directory: Le dossier contenant header.txt et le fichier .snp.
    scenario_index: L'index 1-based du scénario à utiliser.
    values: Les valeurs de paramètres déjà connues, {nom: valeur}.
    num_loci: Voir _simulate_genotypes_for_all_locus_types.
    seed: La graine de la simulation.

Returns:
    L'itérateur des génotypes simulés (même contrat que
    run_poc_for_directory, sans le dict `values` en plus puisqu'il
    est déjà connu de l'appelant).

---

## `compute_summary_statistics_from_values`

### Signature

```python
compute_summary_statistics_from_values(reference_directory: str | Path, scenario_index: int, values: dict[str, float]) -> dict[str, float]
```

### Description

Variante de compute_summary_statistics qui NE TIRE AUCUNE valeur de prior.

Reprend telles quelles des valeurs de paramètres déjà connues,
typiquement les tirages RÉELS d'un reftable DIYABC existant (voir
reftable_loop.replay_reftable_simulation) -- permet de comparer
DIYABC et msprime sur EXACTEMENT les mêmes tirages de priors, sans le
biais possible de deux tirages indépendants.

Args:
    reference_directory: Le dossier contenant header.txt et le
        fichier .snp observé.
    scenario_index: L'index 1-based du scénario à utiliser.
    values: Les valeurs de paramètres déjà connues, {nom: valeur}.
    num_loci: Voir _simulate_genotypes_for_all_locus_types (IndSeq
        uniquement -- ignoré pour PoolSeq).
    seed: La graine de la simulation.
    stats_filter: "ALL" ou "HEADER", voir compute_summary_statistics.
    observed_reads_per_locus: PoolSeq uniquement, voir
        simulate_poolseq_reads_with_mrc_filter.

Returns:
    Le dict summary_statistics (pas de `values` en retour,
    puisqu'il est déjà connu de l'appelant).

---

## `compute_summary_statistics_dna`

### Signature

```python
compute_summary_statistics_dna(reference_directory: str | Path, scenario_index: int) -> tuple[dict[str, float], dict[str, float]]
```

### Description

Calcule les 13 statistiques résumées ADN (compute_all_statistics_dna)
sur des données SIMULÉES par msprime -- équivalent ADN de
compute_summary_statistics (chemin IND/PoolSeq), pour les datasets
qui déclarent des loci séquence (`[S]`, groupes `G2`/`G3`... de
header.txt) plutôt que des SNP.

Tire les paramètres historiques (N1, ta, ts...) ET les priors de
groupe (k1/k2/mus_rate par groupe ADN, en interne à
dna_mutation_simulation_per_locus) depuis `seed` -- voir
compute_summary_statistics_dna_from_values pour la variante qui
rejoue des valeurs déjà connues plutôt que d'en tirer de nouvelles
(paired comparison avec un vrai reftable DIYABC).

`values` (le second élément du tuple retourné) ne contient QUE les
paramètres historiques, pas les priors de groupe -- dna_mutation_
simulation_per_locus ne renvoie nulle part les valeurs de k1/k2/
mus_rate qu'elle a tirées en interne, donc ce `values` seul ne
suffirait pas à rejouer exactement cette même particule (contrairement
au chemin SNP, où `values` capture tout ce qui a été tiré).

Args:
    reference_directory: dossier contenant header.txt/headerRF.txt
        et le fichier .mss observé (son nom lu sur la première ligne
        du header).
    scenario_index: le scénario à utiliser pour construire la
        démographie (pas de tirage pondéré multi-scénario ici,
        contrairement à reftable_loop.run_reftable_simulation).
    seed: graine de la particule -- dérive toutes les graines
        internes (tirage des paramètres historiques, des priors de
        groupe, des généalogies et mutations par locus).
    stats_filter: "ALL" (toutes les stats implémentées) ou "HEADER"
        (seulement celles déclarées dans header.txt, voir
        compute_summary_statistics pour le détail).

Returns:
    (summary_stats, values) -- summary_stats est le dict {nom_
    colonne_diyabc: valeur} de compute_all_statistics_dna (ex.
    "NSS_2_1"), values est {nom_paramètre_historique: valeur}.

---

## `compute_summary_statistics_dna_from_values`

### Signature

```python
compute_summary_statistics_dna_from_values(reference_directory: str | Path, scenario_index: int, values: dict[str, float], group_priors_values: dict[str, float]) -> dict[str, float]
```

### Description

Variante de compute_summary_statistics_dna qui ne tire AUCUNE valeur de prior.

Reprend telles quelles des valeurs de paramètres déjà connues,
typiquement les tirages RÉELS d'un reftable DIYABC existant (voir
reftable_loop.replay_reftable_simulation) -- permet de comparer
DIYABC et msprime sur EXACTEMENT les mêmes tirages de priors, sans le
biais possible de deux tirages indépendants.

Args:
    reference_directory: Le dossier contenant header.txt et le
        fichier .mss observé.
    scenario_index: L'index 1-based du scénario à utiliser.
    values: Les valeurs de paramètres historiques déjà connues,
        {nom: valeur}.
    group_priors_values: Dict {nom_param: valeur} pour tous les
        groupes ADN déclarés dans header.txt. Ce sont les valeurs
        que dna_mutation_simulation_per_locus aurait tirées en
        interne si on avait appelé la variante "random"
        (compute_summary_statistics_dna) -- elles ne sont pas
        capturées par le dict `values` retourné par cette fonction
        (voir compute_summary_statistics_dna).
    seed: La graine du tirage par-locus (second niveau, généalogie,
        mutation).
    stats_filter: "ALL" ou "HEADER", voir compute_summary_statistics.

Returns:
    Le dict summary_statistics (pas de `values` en retour,
    puisqu'ils sont déjà connus de l'appelant).

---

## `compute_summary_statistics_microsat`

### Signature

```python
compute_summary_statistics_microsat(reference_directory: str | Path, scenario_index: int) -> tuple[dict[str, float], dict[str, float]]
```

### Description

Calcule les statistiques résumées microsat (compute_all_statistics_microsat)
sur des données SIMULÉES par msprime -- équivalent microsat de
compute_summary_statistics (chemin IND/PoolSeq), pour les datasets
qui déclarent des loci microsat (`[M]`, groupes `G4`/`G5`... de
header.txt) plutôt que des SNP.

Tire les paramètres historiques (N1, ta, ts...) ET les priors de
groupe (mus_rate,Pgeom) par groupe microsat, en interne à
microsat_mutation_simulation_per_locus) depuis `seed` -- voir
compute_summary_statistics_microsat_from_values pour la variante qui
rejoue des valeurs déjà connues plutôt que d'en tirer de nouvelles
(paired comparison avec un vrai reftable DIYABC).

`values` (le second élément du tuple retourné) ne contient QUE les
paramètres historiques, pas les priors de groupe -- microsat_mutation_
simulation_per_locus ne renvoie nulle part les valeurs de mut_rate,Pgeom
qu'elle a tirées en interne, donc ce `values` seul ne
suffirait pas à rejouer exactement cette même particule (contrairement
au chemin SNP, où `values` capture tout ce qui a été tiré).

Args:
    reference_directory: dossier contenant header.txt/headerRF.txt
        et le fichier .mss observé (son nom lu sur la première ligne
        du header).
    scenario_index: le scénario à utiliser pour construire la
        démographie (pas de tirage pondéré multi-scénario ici,
        contrairement à reftable_loop.run_reftable_simulation).
    seed: La graine du tirage par-locus (second niveau, généalogie,
        mutation).
    stats_filter: "ALL" (toutes les stats implémentées) ou "HEADER"
        (seulement celles déclarées dans header.txt, voir
        compute_summary_statistics pour le détail).
Returns:
    (summary_stats, values) -- summary_stats est le dict {nom_
    colonne_diyabc: valeur} de compute_all_statistics_microsat (ex.
    "NSS_2_1"), values est {nom_paramètre_historique: valeur}.

---


# 📄 prior_parser.py

## `_extract_historical_priors_section`

### Signature

```python
_extract_historical_priors_section(header_text: str) -> list[str]
```

### Description

Extrait les lignes de la section 'historical parameters priors'.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Les lignes de la section, sans les lignes vides.

Raises:
    ValueError: Si la section ou sa fin est introuvable.

---

## `parse_priors`

### Signature

```python
parse_priors(header_text: str) -> tuple[list[Prior], list[OrderConstraint]]
```

### Description

Extrait les priors et les contraintes d'ordre de header.txt.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Un tuple (priors, constraints).

Raises:
    ValueError: Si une ligne de la section ne correspond à aucun des
        deux formats connus (prior ou contrainte) -- pas d'ignorance
        silencieuse.

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

Args:
    prior: Le prior à tester.

Returns:
    True si le prior est quasi-constant selon la règle ci-dessus.

---

## `_extract_priors_group_section`

### Signature

```python
_extract_priors_group_section(header_text: str) -> list[str]
```

### Description

Extrait les lignes de la section 'group priors'.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Les lignes de la section, sans les lignes vides.

Raises:
    ValueError: Si la section ou sa fin est introuvable.

---

## `parse_group_priors`

### Signature

```python
parse_group_priors(header_text: str) -> dict[str, list[GroupPrior]]
```

### Description

Extrait les priors et models des différents group priors de header.txt.

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

---

## `get_parameter_used_by_model`

### Signature

```python
get_parameter_used_by_model(group_prior: GroupPrior) -> tuple[bool, bool]
```

### Description

Détermine les paramètres k1/k2 actifs pour un modèle mutationnel ADN.

JK -> aucun des deux actif, K2P/HKY -> k1 seul, TN -> les deux.

Args:
    group_prior: Un GroupPrior de type model (group_prior.model is True).

Returns:
    (k1_used, k2_used).

Raises:
    NotImplementedError: Si group_prior n'est pas un model, si
        name_model est absent, ou si le modèle n'est pas géré
        (JK/K2P/HKY/TN).

---


# 📄 reftable_loop.py

## `_run_single_particle`

### Signature

```python
_run_single_particle(particle_index: int, reference_directory: Path, scenarios: list[Scenario]) -> ParticleResult
```

### Description

Calcule une seule particule.

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

---

## `run_reftable_simulation`

### Signature

```python
run_reftable_simulation(reference_directory: str | Path, scenarios: list[Scenario]) -> list[ParticleResult]
```

### Description

Produit nrec particules (lignes de reftable.bin) en parallèle.

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

---

## `_kept_param_names_by_scenario`

### Signature

```python
_kept_param_names_by_scenario(priors: list, scenarios: list[Scenario]) -> dict[int, list[str]]
```

### Description

Calcule les noms de paramètres à garder, par scénario.

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

---

## `write_reftable_bin`

### Signature

```python
write_reftable_bin(results: list[ParticleResult], priors: list, scenarios: list[Scenario], output_path: str | Path) -> None
```

### Description

Écrit un reftable.bin au format binaire DIYABC.

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

---

## `write_reftable_txt`

### Signature

```python
write_reftable_txt(results: list[ParticleResult], priors: list, scenarios: list[Scenario], output_path: str | Path) -> None
```

### Description

Écrit les résultats au format texte de DIYABC :
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
un parseur par espaces (ex: pandas read_csv(sep=r'\s+'), ou même un
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

---

## `rewrite_real_reftable_txt`

### Signature

```python
rewrite_real_reftable_txt(input_path: str | Path, output_path: str | Path, priors: list, scenarios: list[Scenario]) -> None
```

### Description

Réécrit un reftable RÉEL de DIYABC en un texte à colonnes de largeur FIXE.

Ex: first_records_of_the_reference_table_0.txt. Remplace les cases
vides de DIYABC (paramètre non utilisé par le scénario de la ligne,
voir parse_real_reftable_params) par `nan` -- jamais une case vide.

Nécessaire pour toute lecture EXTERNE du fichier DIYABC brut avec un
parseur par espaces générique (ex: `pandas.read_csv(sep=r'\s+')`,
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

---

## `simulate_from_directory`

### Signature

```python
simulate_from_directory(test_directory: str | Path) -> list[ParticleResult]
```

### Description

Point d'entrée pour un sous-dossier de test sous reference/.

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

---

## `parse_real_reftable_params`

### Signature

```python
parse_real_reftable_params(path: str | Path, priors: list, scenarios: list[Scenario]) -> list[tuple[int, dict[str, float]]]
```

### Description

Lit un reftable RÉEL produit par DIYABC.

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

---

## `_run_single_particle_from_values`

### Signature

```python
_run_single_particle_from_values(particle_index: int, reference_directory: Path, scenario_index: int, values: dict[str, float]) -> ParticleResult
```

### Description

Variante de _run_single_particle qui NE TIRE AUCUN paramètre.

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

---

## `replay_reftable_simulation`

### Signature

```python
replay_reftable_simulation(reference_directory: str | Path, priors: list, scenarios: list[Scenario], real_reftable_path: str | Path, num_loci: int | None, stats_filter: str, max_workers: int | None) -> list[ParticleResult]
```

### Description

Rejoue, particule par particule, les tirages de paramètres RÉELLEMENT effectués par DIYABC.

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

---

## `_run_single_particle_dna`

### Signature

```python
_run_single_particle_dna(particle_index: int, reference_directory: Path, scenarios: list[Scenario]) -> ParticleResult
```

### Description

Calcule une seule particule ADN (équivalent DNA de _run_single_particle).

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

---

## `run_reftable_simulation_dna`

### Signature

```python
run_reftable_simulation_dna(reference_directory: str | Path, scenarios: list[Scenario]) -> list[ParticleResult]
```

### Description

Produit nrec particules ADN (lignes de reftable.bin) en parallèle.

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

---

## `group_prior_column_names`

### Signature

```python
group_prior_column_names(header_text: str) -> list[str]
```

### Description

Liste ordonnée des noms de colonnes "priors de groupe" d'un vrai reftable DIYABC.

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

---

## `parse_real_reftable_params_with_group_priors`

### Signature

```python
parse_real_reftable_params_with_group_priors(path: str | Path, priors: list, scenarios: list[Scenario], group_priors_names: list) -> list[tuple[int, dict[str, float], dict[str, float]]]
```

### Description

Variante de parse_real_reftable_params qui lit AUSSI les priors de groupe.

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

---

## `_run_single_particle_dna_from_values`

### Signature

```python
_run_single_particle_dna_from_values(particle_index: int, reference_directory: Path, scenario_index: int, values: dict[str, float], group_priors_values: dict[str, float]) -> ParticleResult
```

### Description

Variante de _run_single_particle_dna qui NE TIRE AUCUN paramètre.

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

---

## `replay_reftable_simulation_dna`

### Signature

```python
replay_reftable_simulation_dna(reference_directory: str | Path, priors: list, group_priors_names: list[str], scenarios: list[Scenario], real_reftable_path: str | Path, stats_filter: str, max_workers: int | None) -> list[ParticleResult]
```

### Description

Rejoue, particule par particule, les tirages RÉELS de DIYABC (équivalent ADN de replay_reftable_simulation).

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

---

## `_run_single_particle_microsat`

### Signature

```python
_run_single_particle_microsat(particle_index: int, reference_directory: Path, scenarios: list[Scenario]) -> ParticleResult
```

### Description

Calcule une seule particule microsat (équivalent microsat de _run_single_particle).

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

---

## `run_reftable_simulation_microsat`

### Signature

```python
run_reftable_simulation_microsat(reference_directory: str | Path, scenarios: list[Scenario]) -> list[ParticleResult]
```

### Description

Produit nrec particules microsat (lignes de reftable.bin) en parallèle.

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

---


# 📄 scenario_parser.py

## `split_scenario_blocks`

### Signature

```python
split_scenario_blocks(header_text: str) -> list[str]
```

### Description

Découpe header.txt en blocs bruts, un par scénario.

Chaque bloc commence par sa ligne d'en-tête
'scenario N [poids] (nlignes)' et s'arrête juste avant le bloc suivant
(ou la fin de la section, ex: 'historical parameters priors').

Args:
    header_text: Texte complet de header.txt.

Returns:
    Un bloc brut par scénario, dans l'ordre du fichier.

Raises:
    ValueError: Si aucun bloc 'scenario N [...] (...)' n'est trouvé.

---

## `parse_scenario_block`

### Signature

```python
parse_scenario_block(block_text: str) -> Scenario
```

### Description

Transforme un bloc brut en objet Scenario rempli.

Args:
    block_text: Bloc brut, en commençant par la ligne 'scenario N [...]'.

Returns:
    Le Scenario correspondant.

Raises:
    ValueError: Si la première ligne n'est pas un en-tête de scénario
        valide ('scenario N [poids] (nlignes)').
    NotImplementedError: Propagée par _parse_event_line si une ligne
        d'événement utilise un mot-clé non géré.

---

## `_parse_event_line`

### Signature

```python
_parse_event_line(line: str) -> None
```

### Description

Transforme une ligne d'événement en objet Event correspondant.

Ex: 't1 merge 2 1' -> MergeEvent(time_expr="t1", ancestral_pop=2,
derived_pop=1). Vocabulaire de référence : src-JMC-C++/history.cpp
(ScenarioC::read_events).

Args:
    line: Une ligne d'événement du bloc scénario (ex: 't1 merge 2 1').

Returns:
    Un SampleEvent, MergeEvent, VarNeEvent ou SplitEvent selon le
    mot-clé rencontré.

Raises:
    NotImplementedError: Si le mot-clé d'action n'est pas dans
        sample/merge/varNe/split.

---

## `parse_header_scenarios`

### Signature

```python
parse_header_scenarios(header_text: str) -> list[Scenario]
```

### Description

Point d'entrée principal : header.txt complet -> liste de Scenario.

Important : seule NotImplementedError est avalée ici, volontairement
(le bloc est ignoré avec un warning) -- toute autre exception (erreur
de parsing réelle, bug) continue de se propager normalement.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Les Scenario parsés avec succès. Un bloc dont le vocabulaire n'est
    pas géré est silencieusement ignoré (warning émis), pas levé.

---


# 📄 snp_writer.py

## `_genotypes_to_diploid`

### Signature

```python
_genotypes_to_diploid(haploid_genotypes: list[int]) -> list[int]
```

### Description

Agrège des génotypes haploïdes en génotypes diploïdes.

Somme les paires de lignées consécutives [2i, 2i+1] --
correspondant aux deux copies chromosomiques d'un même individu
(vérifié empiriquement avec ts.individuals()[i].nodes).

Args:
    haploid_genotypes: Une valeur (0/1) par lignée.

Returns:
    Les génotypes diploïdes (0/1/2), un par individu.

Raises:
    ValueError: Si le nombre de lignées est impair (incohérent
        avec une simulation en ploidy=2).

---

## `write_snp_file`

### Signature

```python
write_snp_file(genotypes_per_locus: list[dict[str, list[int]]], output_path: str | Path) -> None
```

### Description

Écrit un fichier .snp DIYABC à partir de génotypes simulés.

Le nom de chaque individu simulé est généré comme "sim_<pop>_<n>"
(ex: "sim_pop1_1", "sim_pop1_2"...). La colonne SEX est fixée à "9"
pour tous les individus -- valeur arbitraire, non confirmée comme
sans impact pour des loci autosomaux <A> (voir notes/exploration.md
pour la justification de cette hypothèse).

Args:
    genotypes_per_locus: num_loci dicts {nom_population:
        [génotypes haploïdes...]}, un par locus (la forme produite
        par ancestry_simulation.simulate_snp_genotypes). Doit
        contenir AU MOINS un locus, et toutes les populations
        doivent être présentes et avoir le même nombre de lignées
        à chaque locus (cohérence vérifiée par la simulation
        elle-même, pas revérifiée ici).
    output_path: Chemin où écrire le fichier .snp.

Raises:
    ValueError: Si genotypes_per_locus est vide.

---


# 📄 statobs_parser.py

## `parse_statobs`

### Signature

```python
parse_statobs(statobs_text: str) -> dict[str, float]
```

### Description

Parse le contenu d'un fichier statobsRF.txt/statobs.txt.

Args:
    statobs_text: Contenu complet du fichier statobsRF.txt/statobs.txt.

Returns:
    Un dict {nom_colonne: valeur}.

Raises:
    ValueError: Si le fichier ne contient pas exactement 2 lignes
        non vides, ou si le nombre de noms et de valeurs ne
        correspond pas.

---


# 📄 stats_group_parser.py

## `_split_stats_blocks`

### Signature

```python
_split_stats_blocks(header_text: str) -> list[str]
```

### Description

Découpe header.txt en blocs bruts, un par groupe de statistiques.

Chaque bloc commence par sa ligne d'en-tête 'group G1 (N)' et
s'arrête juste avant le bloc suivant (ou la fin de la section, ex:
'scenario').

Args:
    header_text: Texte complet de header.txt.

Returns:
    La liste des blocs bruts, un par groupe, dans l'ordre de
    déclaration.

Raises:
    ValueError: Si la section 'group summary statistics' ou aucun
        bloc 'group Gx (N)' n'est trouvé.

---

## `parse_requested_statistic_names`

### Signature

```python
parse_requested_statistic_names(header_text: str) -> list[str]
```

### Description

Extrait les noms de colonnes de statistiques attendues par header.txt.

Section 'group summary statistics'. Un seul groupe -> colonnes
"STAT_index" (ex: "ML1p_1"). Plusieurs groupes -> colonnes
"STAT_groupe_index" (ex: "ML1p_1_1") -- voir le docstring du
module pour le détail des deux formats.

Args:
    header_text: Texte complet de header.txt.

Returns:
    Les noms de colonnes, dans l'ordre de déclaration.

Raises:
    ValueError: Si une ligne de groupe est inattendue, ou si le
        nombre de statistiques trouvées ne correspond pas au
        compte annoncé par 'group Gx (N)'.

---


# 📄 summary_statistics.py

## `_allele_freq`

### Signature

```python
_allele_freq(haploid_genotypes: list[int]) -> float
```

### Description

Calcule la fréquence de l'allèle dérivé (1) dans une population.

Équivalent de locuslist[loc].freq[pop][1] dans le code C++.

Args:
    haploid_genotypes: Les génotypes (0/1) d'une population, un
        locus.

Returns:
    La fréquence de l'allèle dérivé (nan si population vide).

---

## `_q1`

### Signature

```python
_q1(haploid_genotypes: list[int]) -> float
```

### Description

Calcule la probabilité d'identité par état intra-population.

Tirage SANS remise -- formule exacte de sumstat.cpp::q1 (cas SNP,
bias=False) : `q1 = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1))`, où
y1, y2 = comptes d'allèles 0 et 1 (= freq * n).

Args:
    haploid_genotypes: Les génotypes (0/1) d'une population, un
        locus.

Returns:
    q1 (nan si population de taille <= 1).

---

## `_q2`

### Signature

```python
_q2(haploid_genotypes_a: list[int], haploid_genotypes_b: list[int]) -> float
```

### Description

Calcule la probabilité d'identité par état inter-populations.

Formule exacte de sumstat.cpp::q2 (cas SNP) :
`q2 = (y11*y21 + y12*y22) / (n1*n2)`.

Args:
    haploid_genotypes_a: Les génotypes (0/1) de la population A, un
        locus.
    haploid_genotypes_b: Les génotypes (0/1) de la population B,
        même locus.

Returns:
    q2 (nan si l'une des deux populations est vide).

---

## `_prepare_matrices`

### Signature

```python
_prepare_matrices(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
```

### Description

Construit les matrices (npop, nloci) de comptes et fréquences.

Appelé UNE SEULE FOIS dans compute_all_statistics et transmis via
_mats à toutes les familles de statistiques -- évite de reconstruire
les matrices (npop × nloci) une fois par famille.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population, dans l'ordre voulu
        pour les lignes des matrices.

Returns:
    Le tuple (counts, ns, freq0, freq1) :
    counts[i, l] = nb d'allèles dérivés (1) dans pop i au locus l ;
    ns[i, l] = nb total de lignées dans pop i au locus l ;
    freq1 = counts / ns, freq0 = 1 - freq1.

---

## `_prepare_matrices_poolseq`

### Signature

```python
_prepare_matrices_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
```

### Description

Construit les matrices (npop, nloci) de comptes et fréquences pour POOLSEQ.

Équivalent PoolSeq de _prepare_matrices -- même Returns, mais
`counts`/`ns` viennent des lectures observées (reads), pas de
génotypes individuels.

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population, dans l'ordre voulu
        pour les lignes des matrices.

Returns:
    Le tuple (counts, ns, freq0, freq1) :
    counts[i, l] = nb de lectures dérivées observées dans pop i au
    locus l ; ns[i, l] = nb total de lectures observées dans pop i
    au locus l ; freq1 = counts / ns, freq0 = 1 - freq1.

---

## `_forward_fill`

### Signature

```python
_forward_fill(values: np.ndarray, valid: np.ndarray, fill: float) -> np.ndarray
```

### Description

Propage la dernière valeur valide vue aux positions non valides (forward-fill vectorisé).

Reproduit le comportement de la variable C++ non ré-initialisée
x_prev dans cal_snfstd et cal_snnei. Implémentation : searchsorted
sur les indices valides, O(n log n) tout-numpy, sans boucle Python.

Args:
    values: Les valeurs, certaines invalides.
    valid: Masque booléen, même longueur que values.
    fill: Valeur utilisée pour les positions initiales sans aucune
        valeur valide précédente.

Returns:
    Le tableau avec les positions invalides remplacées par la
    dernière valeur valide vue (ou `fill`).

---

## `_half_sorted_by_pairs`

### Signature

```python
_half_sorted_by_pairs(v: list[int]) -> bool
```

### Description

Filtre HALF de DIYABC (cal_snaml / cal_snf3r / cal_snf4r).

Args:
    v: Une permutation d'indices.

Returns:
    True si v est "half-sorted by pairs".

---

## `_half_arrangements`

### Signature

```python
_half_arrangements(n: int, r: int) -> list[list[int]]
```

### Description

Calcule les arrangements HALF de r éléments parmi n.

Ordre reproduit empiriquement depuis DIYABC (cal_snaml, cal_snf3r,
cal_snf4r).

Args:
    n: Le nombre d'éléments parmi lesquels arranger.
    r: La taille de chaque arrangement.

Returns:
    La liste des arrangements HALF, chacun une liste de r indices.

---

## `compute_ML1`

### Signature

```python
compute_ML1(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

ML1p_i : proportion de loci monomorphes dans la population i.

Un locus est monomorphe si sum==0 (fixé ancestral) ou sum==n (fixé
dérivé).

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Matrices (counts, ns, freq0, freq1) déjà calculées par
        _prepare_matrices (voir compute_all_statistics). Si None,
        calculées ici.

Returns:
    Un dict {"ML1p_i": valeur}.

---

## `compute_ML2`

### Signature

```python
compute_ML2(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

ML2p_i.j : proportion de loci fixés au même allèle dans la paire (i, j).

Référence : cal_snfl(npop=2) -- freq_a == freq_b ∈ {0, 1}.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"ML2p_i.j": valeur}, une entrée par paire de
    populations.

---

## `compute_ML3`

### Signature

```python
compute_ML3(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

ML3p_i.j.k : même logique que ML2, sur les triplets de populations.

Référence : cal_snfl(npop=3).

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"ML3p_i.j.k": valeur}, une entrée par triplet de
    populations.

---

## `compute_HW_HB`

### Signature

```python
compute_HW_HB(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

HWm_i/HWv_i (intra-pop) et HBm_i.j/HBv_i.j (inter-pop).

HW = 1 - q1, HB = 1 - q2 (formules de sumstat.cpp vectorisées) :
q1[i,l] = (y1*(y1-1) + y2*(y2-1)) / (n*(n-1)) [sans remise] ;
q2[i,j,l] = (y1_i*y1_j + y2_i*y2_j) / (n_i*n_j). HWv et HBv
utilisent ddof=1 (validé contre le C++).

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"HWm_i": ..., "HWv_i": ..., "HBm_i.j": ...,
    "HBv_i.j": ...}.

---

## `compute_HW_HB_poolseq`

### Signature

```python
compute_HW_HB_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int], _mats) -> dict[str, float]
```

### Description

Variante PoolSeq de compute_HW_HB.

Correction de biais de lecture Q1 nécessaire côté PoolSeq (la
profondeur de séquençage `pool_sizes` n'est pas le vrai nombre de
copies de gène échantillonnées) -- voir la formule `q1` ci-dessous,
absente du chemin IndSeq.

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool},
        voir observed_data._parse_pool_header_line.
    _mats: Matrices (counts, ns, freq0, freq1) déjà calculées par
        _prepare_matrices_poolseq. Si None, calculées ici.

Returns:
    Un dict {"HWm_i": ..., "HWv_i": ..., "HBm_i.j": ...,
    "HBv_i.j": ...}.

---

## `compute_FST1`

### Signature

```python
compute_FST1(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

FST1m_i = 1 - HWm_i / HBmoy_global, FST1v_i = HWv_i / HBmoy_global².

HBmoy_global = moyenne de TOUS les HBm (toutes paires confondues) --
confirmé dans cal_snfsti (sumstat.cpp), pas seulement les paires de
pop_i. FST1v est une propagation d'erreur analytique, pas une
variance empirique.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"FST1m_i": ..., "FST1v_i": ...}.

---

## `compute_FST1_poolseq`

### Signature

```python
compute_FST1_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int], _mats) -> dict[str, float]
```

### Description

Variante PoolSeq de compute_FST1.

cal_snfsti n'a pas de branche type==3 : elle combine juste des
HW/HB déjà calculés. Duplique donc ici la formule q1/q2 poolseq de
compute_HW_HB_poolseq, exactement comme compute_FST1 duplique déjà
la formule q1/q2 IndSeq plutôt que d'appeler compute_HW_HB -- même
style que l'existant.

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.
    _mats: Voir compute_HW_HB_poolseq.

Returns:
    Un dict {"FST1m_i": ..., "FST1v_i": ...}.

---

## `_fst_wc`

### Signature

```python
_fst_wc(loci, pops, _counts, _ns) -> None
```

### Description

Calcule le FST de Weir & Cockerham, vectorisé sur tous les loci.

Formule identique à cal_snfstd, toutes les opérations par-locus
faites en numpy sur des vecteurs de longueur nloci.

Args:
    loci: Liste de dicts {nom_population: [génotype, ...]}, un
        dict par locus.
    pops: Les noms de population à inclure dans ce calcul.
    _counts: Matrice (len(pops), nloci) de comptes d'allèles
        dérivés, déjà calculée -- passée comme slice de la matrice
        globale depuis compute_FST2/3/4 pour éviter de reconstruire
        les comptes locus par locus pour chaque sous-ensemble. Si
        None, calculée ici.
    _ns: Matrice (len(pops), nloci) de tailles d'échantillon,
        même principe que _counts.

Returns:
    Le tuple (FSTm, FSTv).

---

## `compute_FST2`

### Signature

```python
compute_FST2(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

FST2m_i.j / FST2v_i.j : Weir & Cockerham par paire.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"FST2m_i.j": ..., "FST2v_i.j": ...}.

---

## `_fst_wc_poolseq`

### Signature

```python
_fst_wc_poolseq(pops, pool_sizes, _counts, _ns) -> None
```

### Description

Variante PoolSeq de _fst_wc (cal_snfstd, branche grouplist[gr].type==3).

Calcul en DEUX passes (contrairement à l'IndSeq) : pi1/pi2 (moyennes
pondérées par la profondeur de lecture) doivent être connues avant de
calculer SSP. C_1/C_1_star mélangent la profondeur de lecture (`n`,
variable par locus) et la VRAIE taille du pool (`c`, constante par
population) -- c'est ce mélange qui constitue la correction propre à
PoolSeq (le terme de variance intra-pool supplémentaire).

L'agrégation finale (ratio de sommes num/den + variance via
_forward_fill) est IDENTIQUE à _fst_wc -- confirmé par l'exploration
C++, même code d'agrégation pour les deux types de population.

Args:
    pops: Les noms de population à inclure dans ce calcul.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.
    _counts: Matrice (len(pops), nloci) de lectures dérivées
        (nreads1), slice de la matrice globale -- même contrat que
        _fst_wc.
    _ns: Matrice (len(pops), nloci) de profondeur de lecture
        (nreads_total), même principe que _counts.

Returns:
    Le tuple (FSTm, FSTv).

---

## `compute_FST2_poolseq`

### Signature

```python
compute_FST2_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int], _mats) -> dict[str, float]
```

### Description

Variante PoolSeq de compute_FST2 : FST2m_i.j / FST2v_i.j par paire.

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.
    _mats: Voir compute_HW_HB_poolseq.

Returns:
    Un dict {"FST2m_i.j": ..., "FST2v_i.j": ...}.

---

## `compute_FST3_FST4_poolseq`

### Signature

```python
compute_FST3_FST4_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int], _mats) -> dict[str, float]
```

### Description

Variante PoolSeq de compute_FST3_FST4_FSTG : FST3/FST4 sur triplets/quadruplets (COMB).

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.
    _mats: Voir compute_HW_HB_poolseq.

Returns:
    Un dict {"FST3m_i.j.k": ..., "FST3v_i.j.k": ...} et, s'il y a
    au moins 4 populations, {"FST4m_i.j.k.l": ...,
    "FST4v_i.j.k.l": ...}.

---

## `compute_FST3_FST4_FSTG`

### Signature

```python
compute_FST3_FST4_FSTG(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

FST3/FST4 : Weir & Cockerham sur triplets et quadruplets (COMB).

Ne calcule PAS `FSTG` malgré son nom, seulement FST3/FST4 --
`FSTG` (FST global, toutes populations combinées d'un coup, même
`cal_snfstd` que FST2/3/4 avec npop=0, voir statdefs.cpp:184) n'est
implémenté nulle part dans ce module. Non bloquant en pratique :
aucun header.txt de ce dépôt ne le déclare dans sa section 'group
summary statistics', donc `stats_filter="HEADER"` ne le demande
jamais.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"FST3m_i.j.k": ..., "FST3v_i.j.k": ...} et, s'il y a
    au moins 4 populations, {"FST4m_i.j.k.l": ...,
    "FST4v_i.j.k.l": ...}.

---

## `compute_NEI`

### Signature

```python
compute_NEI(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

NEIm_i.j et NEIv_i.j : distance de Nei (1972) par paire, vectorisée.

`NEI = 1 - (fi*fj + gi*gj) / sqrt(fi²+gi²) / sqrt(fj²+gj²)`. x_prev
persiste si denom==0 (comportement C++ non réinitialisé) --
reproduit via _forward_fill.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"NEIm_i.j": ..., "NEIv_i.j": ...}.

---

## `compute_AML`

### Signature

```python
compute_AML(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

AMLm / AMLv : coefficient d'admixture ML sur triplets HALF.

`aml = (f3-f2)/(f1-f2)` clampé à [0,1]. Les loci non informatifs
(f1==f2, w=0) sont exclus de la moyenne -- équivalent au Welford
pondéré avec w ∈ {0,1}, ce qui réduit à mean/var(ddof=1) sur le
sous-ensemble informatif.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"AMLm_h.p1.p2": ..., "AMLv_h.p1.p2": ...}, une entrée
    par arrangement HALF de 3 populations.

---

## `compute_F3`

### Signature

```python
compute_F3(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

F3m/F3v sur triplets HALF.

`F3 = (f1-f2)*(f1-f3) - f1*(1-f1)/(np-1)` (pop0=hybride,
pop1/2=parents). Tous les loci ont w=1 → mean/var(ddof=1)
directement.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"F3m_i0.i1.i2": ..., "F3v_i0.i1.i2": ...}, une entrée
    par arrangement HALF de 3 populations.

---

## `compute_F3_poolseq`

### Signature

```python
compute_F3_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int], _mats) -> dict[str, float]
```

### Description

Variante PoolSeq de compute_F3 (cal_snf3r, branche grouplist[gr].type==3).

`alpha = ((np*a1p*(a1p-1))/(c1p*(c1p-1)) - a1p/c1p) / (np-1)`
(pop0=hybride), `F3 = alpha + betaBC - betaAB - betaAC`, avec
`beta_XY = (aXp*aYp)/(cXp*cYp)`.

`np` = taille du pool (VRAIE, pas la profondeur de lecture) de la
population hybride -- vient de pool_sizes, pas de `ns`/`_mats`
(contrairement à ns, qui est la profondeur de lecture, variable par
locus). Agrégation identique à compute_F3 (mean/var(ddof=1) simples
sur les loci, pas de ratio de sommes).

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.
    _mats: Voir compute_HW_HB_poolseq.

Returns:
    Un dict {"F3m_i0.i1.i2": ..., "F3v_i0.i1.i2": ...}.

---

## `compute_F4`

### Signature

```python
compute_F4(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str], _mats) -> dict[str, float]
```

### Description

F4m/F4v sur quadruplets HALF.

`F4 = (a-b)*(c-d)`. Tous les loci ont w=1 → mean/var(ddof=1)
directement.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.
    _mats: Voir compute_ML1.

Returns:
    Un dict {"F4m_ia.ib.ic.id": ..., "F4v_ia.ib.ic.id": ...}, une
    entrée par arrangement HALF de 4 populations.

---

## `_genotype_matrix_by_population`

### Signature

```python
_genotype_matrix_by_population(tree_sequence: tskit.TreeSequence) -> dict[str, np.ndarray]
```

### Description

Découpe la matrice de génotypes d'un locus ADN par population.

genotype_matrix() n'est appelé qu'UNE FOIS pour toute la
TreeSequence, puis tranché par population via fancy indexing (pas
de reconstruction par sample).

Args:
    tree_sequence: La TreeSequence mutée d'un locus.

Returns:
    Un dict {nom_pop: matrice (n_sites, n_samples_pop)} --
    convention native de tskit (genotype_matrix() est déjà (sites,
    samples)), pas de transposition.

---

## `_segregating_sites_mask`

### Signature

```python
_segregating_sites_mask(matrix: np.ndarray) -> np.ndarray
```

### Description

Calcule le masque des sites ségrégeants d'une matrice de génotypes.

Factorisé hors de _count_segregating_sites pour être réutilisé par
PSS (_private_segregating_sites_per_locus), qui a besoin du masque
par site, pas seulement du compte agrégé.

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Masque booléen (longueur n_sites) : True si le site n'est pas
    identique chez tous les échantillons.

---

## `_count_segregating_sites`

### Signature

```python
_count_segregating_sites(matrix: np.ndarray) -> int
```

### Description

Compte le nombre de sites polymorphes d'une matrice de génotypes.

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Le nombre de sites ségrégeants.

Raises:
    ValueError: Si matrix.shape[1] == 0 (population sans échantillon).

---

## `compute_NSS`

### Signature

```python
compute_NSS(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule NSS_i (cal_nss1p) : pour chaque population, la moyenne du
nombre de sites ségrégeants sur tous les loci du groupe passé en
argument (un groupe = les TreeSequences des loci séquence d'un même
`group Gx` du header, ex. les 5 loci <A> de G2).

`population_names` fixe explicitement les clés du dict retourné (comme
compute_ML1/compute_HW_HB) -- chaque population attendue a toujours une
valeur (0.0 par défaut, comme le `res = 0.0` du C++), même si
`tree_sequences` est vide, plutôt que d'être silencieusement absente du
résultat.

Suppose que toutes les populations de `population_names` sont présentes
sur tous les loci du groupe (divise par `len(tree_sequences)`, pas par
un décompte par population comme le `nl` du C++ -- lève un KeyError si
ce n'est pas le cas plutôt que d'exclure silencieusement ce locus,
contrairement au C++) -- vérifié vrai sur toy_example2_ms_dna, pas
garanti en général.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `_count_distinct_haplotypes`

### Signature

```python
_count_distinct_haplotypes(matrix: np.ndarray) -> int
```

### Description

Compte le nombre d'haplotypes distincts d'une matrice de génotypes.

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Le nombre d'haplotypes distincts (colonnes uniques).

Raises:
    ValueError: Si matrix.shape[1] == 0 (population sans échantillon).

---

## `compute_NHA`

### Signature

```python
compute_NHA(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule le nombre moyen d'haplotypes distincts par population sur un groupe de loci.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `_pairwise_hamming_distances`

### Signature

```python
_pairwise_hamming_distances(matrix: np.ndarray) -> np.ndarray
```

### Description

Calcule les distances de Hamming par paire d'échantillons.

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Un vecteur 1D de longueur C(n_samples, 2) -- une valeur par
    paire (i, j) avec i < j, pas la matrice carrée
    (n_samples, n_samples) complète.

Raises:
    ValueError: Si matrix.shape[1] == 0 (population sans échantillon).

---

## `compute_MPD`

### Signature

```python
compute_MPD(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule MPD_i (cal_mpd1p) : pour chaque population, la moyenne du
nombre de différences par paire (distance de Hamming) sur tous les
loci du groupe passé en argument (un groupe = les TreeSequences des
loci séquence d'un même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut,
comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

Contrairement à compute_NSS/compute_NHA
(qui divisent par `len(tree_sequences)`), le
dénominateur ici est un compteur PAR POPULATION (comme le `nl` de
cal_mpd1p) : un locus où une population a moins de 2 échantillons
donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
vide), auquel cas ce locus est exclu du calcul pour cette population
-- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
laisser un `nan` (moyenne d'un vecteur vide) empoisonner le résultat.

Suppose quand même que toutes les populations de `population_names`
sont présentes (au moins 1 échantillon) sur tous les loci du groupe
-- lève un KeyError si ce n'est pas le cas, comme les autres
fonctions `compute_*` de ce module.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `compute_VPD`

### Signature

```python
compute_VPD(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule VPD_i (cal_vpd1p) : pour chaque population, la variance du
nombre de différences par paire (distance de Hamming) sur tous les
loci du groupe passé en argument (un groupe = les TreeSequences des
loci séquence d'un même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut,
comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

Contrairement à compute_NSS/compute_NHA
(qui divisent par `len(tree_sequences)`), le
dénominateur ici est un compteur PAR POPULATION (comme le `nl` de
cal_vpd1p) : un locus où une population a moins de 2 échantillons
donne 0 paire (`_pairwise_hamming_distances` retourne un vecteur
vide), auquel cas ce locus est exclu du calcul pour cette population
-- ni ajouté à la somme, ni compté au dénominateur -- plutôt que de
laisser un `nan` (variance d'un vecteur vide) empoisonner le résultat.

Suppose quand même que toutes les populations de `population_names`
sont présentes (au moins 1 échantillon) sur tous les loci du groupe
-- lève un KeyError si ce n'est pas le cas, comme les autres
fonctions `compute_*` de ce module.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_variance}.

---

## `_tajima_constants`

### Signature

```python
_tajima_constants(n_samples: int) -> tuple[float, float, float]
```

### Description

Calcule les constantes a1, e1, e2 pour le D de Tajima (cal_dta1pl,
lignes 1566-1575) à partir du nombre d'échantillons (n_samples) --
ne dépend que de n_samples, jamais des données elles-mêmes.

Args:
    n_samples: Nombre d'échantillons (>= 2).

Returns:
    Le tuple (a1, e1, e2) -- a1 est aussi réutilisé directement dans
    la formule finale du D (S / a1), e1/e2 sont les coefficients de
    la variance sous neutralité. b1, b2, c1, c2, a2 sont des étapes
    intermédiaires purement internes, jamais réutilisées ailleurs.

---

## `_tajima_d_per_locus`

### Signature

```python
_tajima_d_per_locus(matrix: np.ndarray) -> float | None
```

### Description

D de Tajima (cal_dta1pl) pour UNE population sur UN locus.

`D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))` -- pi = MPD (moyenne des
différences par paire), S = NSS (nombre de sites ségrégeants),
a1/e1/e2 = _tajima_constants(n_samples).

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Le D de Tajima, ou `None` si `n_samples < 2` (pas assez
    d'échantillons pour calculer ne serait-ce qu'une paire -- ce
    locus doit être EXCLU de la moyenne du groupe, `OKK = false`
    côté C++). Retourne `0.0` (pas `None`) si `n_samples >= 2` mais
    qu'aucun site n'est ségrégeant (S == 0, le dénominateur de la
    formule est nul) -- ce locus reste INCLUS dans la moyenne du
    groupe avec une valeur de 0.0, contrairement au cas précédent :
    le C++ ne repasse jamais `OKK` à false dans ce cas (lignes
    1579/1594-1598) -- distinction délibérée, pas une
    simplification de notre part.

---

## `compute_DTA`

### Signature

```python
compute_DTA(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule DTA_i (cal_dta1p) : pour chaque population, la moyenne du
D de Tajima (_tajima_d_per_locus) sur tous les loci du groupe passé
en argument (un groupe = les TreeSequences des loci séquence d'un
même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut,
comme le `res = 0.0` du C++), même si `tree_sequences` est vide.

Comme compute_MPD/variance_pairwise_
differences_per_group, le dénominateur est un compteur PAR
POPULATION (le `nl` de cal_dta1p) : un locus où `_tajima_d_per_locus`
retourne `None` (moins de 2 échantillons) est exclu -- ni ajouté à la
somme, ni compté. Un locus où `_tajima_d_per_locus` retourne `0.0`
(0 site ségrégeant) reste, lui, INCLUS dans le compte (voir
_tajima_d_per_locus pour la distinction).

Suppose quand même que toutes les populations de `population_names`
sont présentes (au moins 1 échantillon) sur tous les loci du groupe
-- lève un KeyError si ce n'est pas le cas, comme les autres
fonctions `compute_*` de ce module.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `_private_segregating_sites_per_locus`

### Signature

```python
_private_segregating_sites_per_locus(genotype_matrices: dict[str, np.ndarray], target_pop: str) -> int
```

### Description

Compte les sites ségrégeants "privés" de `target_pop` sur UN locus (cal_pss1p).

Un site ségrégeant privé est ségrégeant dans `target_pop` mais
NULLE PART ailleurs, parmi TOUTES les populations de
`genotype_matrices` (pas seulement celles d'un même groupe -- le
C++ compare à `this->nsample`, le nombre total de populations du
dataset).

Toutes les matrices de `genotype_matrices` viennent du même
`genotype_matrix()` (juste tranchées par colonnes, voir
_genotype_matrix_by_population) -- la ligne `i` désigne donc le MÊME
site physique pour toutes les populations. Contrairement au C++, qui
compare des listes d'indices de sites variables de longueurs
différentes par une recherche imbriquée (`ssa[sample][j] ==
ssa[sa][k]`), on peut donc comparer les masques booléens position par
position directement -- pas de recherche d'égalité nécessaire.

Args:
    genotype_matrices: Dict {nom_population: matrice} pour TOUTES
        les populations du dataset (voir _genotype_matrix_by_population).
    target_pop: La population pour laquelle compter les sites
        privés.

Returns:
    Le nombre de sites ségrégeants privés de target_pop.

---

## `compute_PSS`

### Signature

```python
compute_PSS(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule PSS_i (cal_pss1p) : pour chaque population, la moyenne du
nombre de sites ségrégeants privés (_private_segregating_sites_per_
locus) sur tous les loci du groupe passé en argument.

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. `population_names` DOIT couvrir
TOUTES les populations du dataset (pas seulement celles d'un groupe),
puisque `_private_segregating_sites_per_locus` compare `target_pop` à
toutes les autres populations présentes dans `genotype_matrices`.

Contrairement à compute_NSS et aux autres
fonctions `compute_*` de ce module, le dénominateur ici est
`len(tree_sequences)` SANS AUCUNE exclusion (`nl` s'incrémente sans
condition dans cal_pss1p, ligne 1624 -- pas de garde-fou du tout,
même pas le `samplesize > 0` de NSS/NHA).

Suppose que toutes les populations de `population_names` sont
présentes sur tous les loci du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues (toutes celles du
        dataset).

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `_minor_allele_counts_at_segregating_sites`

### Signature

```python
_minor_allele_counts_at_segregating_sites(matrix: np.ndarray) -> np.ndarray
```

### Description

Compte l'allèle minoritaire à chaque site ségrégeant.

Pour chaque site ségrégeant d'une matrice de génotypes, compte le(s)
allèle(s) le(s) moins fréquent(s) parmi les échantillons (afs,
`nf[jj]` après tri croissant et saut des zéros --
équivalent à `min(comptes des valeurs effectivement présentes)`,
généralisation multi-allélique du "minor allele count" puisqu'un
site ADN peut avoir jusqu'à 4 états A/C/G/T, pas seulement 2 comme
un SNP).

Args:
    matrix: Matrice de génotypes (n_sites, n_samples).

Returns:
    Un vecteur 1D de longueur = nombre de sites SÉGRÉGEANTS (pas
    n_sites) -- un site fixé (une seule valeur parmi les
    échantillons) n'est pas inclus, comme `afs` qui ne pousse rien
    dans `t_afs` quand `jj >= 3` (une seule base présente).

Raises:
    ValueError: Si matrix.shape[1] == 0 (population sans échantillon).

---

## `compute_MNS`

### Signature

```python
compute_MNS(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule MNS_i (cal_mns1p) : pour chaque population, la moyenne,
sur les loci du groupe, de la moyenne (par locus) des comptes
d'allèle minoritaire aux sites ségrégeants
(_minor_allele_counts_at_segregating_sites).

Un locus sans site ségrégeant contribue 0.0 (la boucle C++ sur
`t_afs` ne s'exécute simplement pas -- même effet qu'un vecteur
vide ici). Comme PSS (et contrairement à MPD/VPD/DTA), AUCUNE
exclusion de locus : `nl` = `len(tree_sequences)` sans condition
(cal_mns1p, `nl++` inconditionnel, ligne 1705) -- pas besoin de
compteur par population.

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `compute_VNS`

### Signature

```python
compute_VNS(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule VNS_i (cal_vns1p) : pour chaque population, la moyenne,
sur les loci du groupe, de la variance (par locus) des comptes
d'allèle minoritaire aux sites ségrégeants.

ATTENTION : variance BIAISÉE (ddof=0, division par n, PAS n-1) --
contrairement à compute_VPD (VPD) qui
utilise ddof=1. Vérifié explicitement contre cal_vns1p, ligne 1731 :
`v = (sx2 - sx*sx/a) / a`, pas `/ (a-1)`.

Un locus avec moins de 2 sites ségrégeants contribue 0.0 (`v = 0.0`
explicite, ligne 1734 -- pas assez de points pour une variance).
Comme MNS/PSS, AUCUNE exclusion de locus au niveau du groupe (`nl`
s'incrémente sans condition, ligne 1736) : le `if (v > 0.0) res +=
v` du C++ (ligne 1738) n'exclut rien du dénominateur, il évite
seulement d'ajouter un `v` négatif -- ce qui ne peut mathématiquement
pas arriver pour une vraie variance (`v >= 0` toujours), donc ce
garde-fou n'a aucun effet observable et n'est pas reproduit ici.

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {nom_population: valeur_moyenne}.

---

## `compute_NH2`

### Signature

```python
compute_NH2(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule NH2_ij (cal_nh2p) : pour chaque paire de populations, la
moyenne du nombre d'haplotypes distincts sur tous les loci du groupe
passé en argument (un groupe = les TreeSequences des loci séquence
d'un même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {"i.j": valeur_moyenne}, une entrée par paire de
    populations.

---

## `compute_NS2`

### Signature

```python
compute_NS2(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule NS2_ij (cal_ns2p) : pour chaque paire de populations, la
moyenne du nombre de sites ségrégeants sur tous les loci du groupe
passé en argument (un groupe = les TreeSequences des loci séquence
d'un même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {"i.j": valeur_moyenne}, une entrée par paire de
    populations.

---

## `_mean_pairwise_differences_within_per_locus`

### Signature

```python
_mean_pairwise_differences_within_per_locus(genotype_matrices: dict[str, np.ndarray], pop_a: str, pop_b: str) -> float
```

### Description

Calcule MP2 "within" pour un locus : ratio poolé des sommes, PAS la moyenne des deux MPD.

Additionne les distances de Hamming intra-population de pop_a ET
pop_b (jamais entre les deux), puis divise par le nombre total de
paires des deux côtés -- équivalent à la moyenne des deux MPD
seulement si pop_a et pop_b ont la même taille d'échantillon (voir
compute_MP2).

Args:
    genotype_matrices: Dict {nom_population: matrice}, au moins
        pop_a et pop_b.
    pop_a: Nom de la première population.
    pop_b: Nom de la seconde population.

Returns:
    Le ratio poolé (somme des différences / somme des paires),
    0.0 si aucune des deux populations n'a de paire.

---

## `compute_MP2`

### Signature

```python
compute_MP2(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule MP2_ij (cal_mp2p) : pour chaque paire de populations, la
moyenne du nombre de différences par paire (distance de Hamming)
sur tous les loci du groupe passé en argument (un groupe = les
TreeSequences des loci séquence d'un même `group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {"i.j": valeur_moyenne}, une entrée par paire de
    populations.

---

## `_pairwise_hamming_distances_between`

### Signature

```python
_pairwise_hamming_distances_between(matrix_a: np.ndarray, matrix_b: np.ndarray) -> np.ndarray
```

### Description

Calcule les distances de Hamming par paire entre deux matrices de génotypes.

Contrairement à _pairwise_hamming_distances (matrice unique, une
triangulaire à extraire), ici TOUTE paire (p dans matrix_a, q dans
matrix_b) est valide -- jamais d'auto-comparaison, donc pas de
triangle à extraire.

Args:
    matrix_a: Matrice de génotypes (n_sites, n_samples_a).
    matrix_b: Matrice de génotypes (n_sites, n_samples_b), même
        n_sites que matrix_a.

Returns:
    La matrice 2D (n_samples_a, n_samples_b) des distances de
    Hamming, pas un vecteur aplati.

Raises:
    ValueError: Si matrix_a et matrix_b n'ont pas le même nombre de
        sites, ou si l'une des deux est vide (0 échantillon).

---

## `_mean_pairwise_differences_between_per_locus`

### Signature

```python
_mean_pairwise_differences_between_per_locus(genotype_matrices: dict[str, np.ndarray], pop_a: str, pop_b: str) -> float
```

### Description

Calcule la moyenne des différences par paire (MPB) entre deux populations pour un locus.

Args:
    genotype_matrices: Dict {nom_population: matrice}, au moins
        pop_a et pop_b.
    pop_a: Nom de la première population.
    pop_b: Nom de la seconde population.

Returns:
    La moyenne des distances de Hamming entre chaque échantillon de
    pop_a et chaque échantillon de pop_b.

Raises:
    KeyError: Si pop_a ou pop_b n'est pas dans genotype_matrices.

---

## `compute_MPB`

### Signature

```python
compute_MPB(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule MPB_ij (cal_mpb2p) : pour chaque paire de populations, la
moyenne du nombre de différences par paire (distance de Hamming)
entre les deux populations sur tous les loci du groupe passé en
argument (un groupe = les TreeSequences des loci séquence d'un même
`group Gx` du header).

`population_names` fixe explicitement les clés du dict retourné --
chaque population attendue a toujours une valeur (0.0 par défaut),
même si `tree_sequences` est vide. Suppose que toutes les
populations de `population_names` sont présentes sur tous les loci
du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {"i.j": valeur_moyenne}, une entrée par paire de
    populations.

---

## `compute_HST`

### Signature

```python
compute_HST(tree_sequences: list[tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule HST_ij (cal_fst2p) : pour chaque paire de populations,
une mesure de différenciation type FST à partir des loci du groupe
passé en argument (un groupe = les TreeSequences des loci séquence
d'un même `group Gx` du header).

HST = num/den, où num = somme sur les loci de (Hb_locus - Hw_locus)
et den = somme sur les loci de Hb_locus (Hb = MPB par-locus, Hw =
MPW par-locus) -- un RATIO DE SOMMES accumulées sur tout le groupe,
PAS une moyenne de valeurs par-locus divisée par num_loci (contraire
à compute_MP2/compute_MPB) -- même schéma d'agrégation que _fst_wc
(FST2/FST3/FST4, côté SNP).

`population_names` fixe explicitement les clés du dict retourné --
chaque paire attendue a toujours une valeur (0.0 par défaut, y
compris si `den == 0` pour cette paire), même si `tree_sequences`
est vide. Suppose que toutes les populations de `population_names`
sont présentes sur tous les loci du groupe -- lève un KeyError sinon.

Args:
    tree_sequences: Les TreeSequences mutées du groupe (un locus [S]
        chacune).
    population_names: Les populations attendues.

Returns:
    Un dict {"i.j": valeur}, une entrée par paire de populations.

---

## `compute_all_statistics`

### Signature

```python
compute_all_statistics(genotypes_per_locus: list[dict[str, list[int]]], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule les 130 statistiques résumées SNP (IndSeq).

Les matrices (npop × nloci) de comptes et fréquences sont construites
une seule fois (_prepare_matrices) et transmises à toutes les familles
de statistiques via _mats.

Args:
    genotypes_per_locus: Liste de dicts {nom_population:
        [génotype, ...]}, un dict par locus.
    population_names: Les noms de population.

Returns:
    Un dict {nom_stat: valeur} -- même format que parse_statobs().

---

## `compute_all_statistics_poolseq`

### Signature

```python
compute_all_statistics_poolseq(reads_per_locus: list[dict[str, tuple[int, int]]], population_names: list[str], pool_sizes: dict[str, int]) -> dict[str, float]
```

### Description

Calcule les statistiques résumées SNP pour PoolSeq.

Les matrices (npop × nloci) de comptes et tailles d'échantillon sont
construites une seule fois (_prepare_matrices_poolseq) et transmises
à toutes les familles de statistiques via _mats.

Args:
    reads_per_locus: Liste de dicts {nom_population: (nreads_dérivé,
        nreads_total)}, un dict par locus.
    population_names: Les noms de population.
    pool_sizes: Dict {nom_population: taille_haploïde du pool}.

Returns:
    Un dict {nom_stat: valeur} -- même format que parse_statobs().

---

## `compute_all_statistics_dna`

### Signature

```python
compute_all_statistics_dna(header_text: str, tree_sequences_by_locus: dict[str, tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule les 13 statistiques résumées ADN pour chaque `group Gx`
séquence (`[S]`) du header, et retourne un dict {nom_colonne: valeur}
utilisant les VRAIS noms de colonnes DIYABC (`STAT_<groupe>_<pop-ou-
paire>`, ex. `NSS_2_1`, `NH2_3_1.2`) -- vérifié caractère pour
caractère contre la sortie réelle de `diyabc` sur
`toy_example2_ms_dna` (`STAT_<groupe>_<suffixe>` quand il y a
plusieurs groupes, `STAT_<suffixe>` seul sinon -- même convention
que `stats_group_parser.parse_requested_statistic_names`).

Args:
    header_text: contenu de header.txt/headerRF.txt (pour
        `parse_loci_description`, qui donne le groupe de chaque
        locus).
    tree_sequences_by_locus: {nom_locus: TreeSequence mutée} --
        la sortie de `dna_mutation_simulation_per_locus`. Les loci
        microsat (`ms_or_seq == "M"`) présents dans le header sont
        ignorés ici (pas de code de simulation microsat).
    population_names: toutes les populations du dataset, dans
        l'ordre "pop1".."popN" (leur position dans cette liste,
        pas leur nom, détermine l'indice numérique utilisé dans
        les noms de colonnes).

Returns:
    Un dict {nom_colonne_diyabc: valeur}.

---

## `compute_all_statistics_microsat`

### Signature

```python
compute_all_statistics_microsat(header_text: str, tree_sequences_by_locus: dict[str, tskit.TreeSequence], population_names: list[str]) -> dict[str, float]
```

### Description

Calcule les statistiques résumées microsat pour chaque `group Gx`
microsat (`[M]`) du header, et retourne un dict {nom_colonne: valeur}
utilisant les VRAIS noms de colonnes DIYABC -- vérifié caractère pour
caractère contre la sortie réelle de `diyabc` sur
`toy_example2_ms_dna` (`STAT_<groupe>_<suffixe>` quand il y a
plusieurs groupes, `STAT_<suffixe>` seul sinon -- même convention
que `stats_group_parser.parse_requested_statistic_names`).

Args:
    header_text: contenu de header.txt/headerRF.txt (pour
        `parse_loci_description`, qui donne le groupe de chaque
        locus).
    tree_sequences_by_locus: {nom_locus: TreeSequence mutée} --
        la sortie de `microsat_mutation_simulation_per_locus`. Les loci ADN (`ms_or_seq == "S"`) présents dans le header sont ignorés ici (pas de code de simulation ADN).
    population_names: toutes les populations du dataset, dans
        l'ordre "pop1".."popN" (leur position dans cette liste,
        pas leur nom, détermine l'indice numérique utilisé dans
        les noms de colonnes).
Returns:
    Un dict {nom_colonne_diyabc: valeur}.

---
