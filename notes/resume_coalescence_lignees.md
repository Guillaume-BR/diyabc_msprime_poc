# Coalescence, nombre effectif de copies de gène, et lignées vs loci

Deux notions de base du modèle de coalescence utilisé par msprime, utiles
pour comprendre `coalescence_coefficient`/`rescale_demography`
(`observed_data.py`/`demography_builder.py`) et le rôle du sex-ratio dans
la ploïdie des loci `<X>`/`<Y>`/`<M>`.

## Pourquoi le temps de coalescence dépend du nombre effectif de copies de gène

Le coalescent de Kingman remonte le temps généalogique en fusionnant des
lignées deux par deux. À chaque génération dans le passé, la probabilité
que deux lignées données fusionnent (aient le même parent) est environ
`1/Nₑ`, où `Nₑ` est le nombre de copies du gène qui circulent dans la
population à cette génération.

**Intuition** : plus il y a de copies différentes du gène en circulation,
plus il y a de "parents possibles" parmi lesquels une lignée peut avoir
hérité sa copie — donc moins de chance que deux lignées prises au hasard
partagent un parent récent. À l'inverse, s'il y a peu de copies en
circulation, deux lignées ont statistiquement plus de chances de
"retomber" sur le même ancêtre rapidement.

Avec `k` lignées présentes, il y a `k(k-1)/2` paires possibles, donc le
taux de coalescence total est `~k(k-1)/(2Nₑ)` par génération, et le temps
d'attente moyen avant la prochaine fusion est `~2Nₑ/(k(k-1))` — c'est
exactement la forme de la formule citée dans CLAUDE.md
(`particuleC.cpp:1340`, `temps -= coeffcoal*N/(n*(n-1))*log(ra)`).

**Pourquoi ça dépend du type de locus** : `Nₑ` n'est pas le nombre
d'individus, c'est le nombre de copies du *gène considéré* qui
circulent :

- `<A>` autosomal : chaque individu porte 2 copies → `Nₑ = 2N`.
- `<X>` : les femmes ont 2 copies, les hommes 1 seule → `Nₑ ≈ 1.5N` en
  moyenne.
- `<Y>`/`<M>` : une seule copie par individu, transmise par un seul sexe
  (uniparental, sans recombinaison) → `Nₑ ≈ 0.5N`, beaucoup moins de
  copies en circulation.

Moins de copies (`<Y>`/`<M>`) → coalescence plus rapide, plus de dérive
génétique, moins de diversité — un résultat classique en génétique des
populations (l'ADN mitochondrial coalesce ~4x plus vite que l'autosomal).

C'est ce `Nₑ` (pas `N` brut, le prior de taille de population tiré) qu'il
faut donner à msprime pour que son moteur de coalescence tourne au bon
rythme selon le type de locus — d'où `rescale_demography` : on prend la
`Demography <A>` (bâtie avec `N` brut) et on la multiplie par
`coalescence_coefficient(locus_type, sex_ratio) / 2` pour obtenir `Nₑ`,
sans changer l'algorithme de coalescence lui-même.

## Lignées et nombre de loci : deux axes différents

**Les lignées (`k` dans la formule ci-dessus)** concernent une seule
généalogie, à un seul locus. Elles représentent le nombre de copies de
gène échantillonnées pour ce locus précis, qui diminue au fur et à
mesure qu'on remonte le temps et que des lignées fusionnent (coalescent),
jusqu'à n'en avoir plus qu'une (l'ancêtre commun, MRCA). Par exemple pour
un locus `<A>` avec 30 individus diploïdes échantillonnés dans une
population : `k` démarre à 60 (2 copies par individu) et diminue vers 1.

**Le nombre de loci** (5000 pour `human`, 650 pour `toy_example5`...)
c'est complètement autre chose : c'est le nombre de généalogies
**indépendantes** qu'on simule au total. Le pipeline
(`ancestry_simulation.simulate_independent_loci`) simule un arbre séparé
par locus, sans recombinaison ni lien entre eux — chaque locus a donc sa
propre valeur de `k` qui évolue dans le temps, son propre historique de
coalescence, potentiellement sa propre topologie.

Donc : `k` (nombre de lignées) varie **au sein** d'une généalogie d'un
seul locus (du nombre d'échantillons jusqu'à 1), et le nombre de loci
détermine **combien de fois** on répète toute cette simulation
indépendamment, chacune avec son propre point de départ
`k = n_échantillons`.

**Exemple concret (`human`)** : 4 populations de 30 individus (`ASW`,
`YRI`, `CHB`, `GBR`, vérifié via `count_samples_per_population` sur le
`.snp` réel), soit 120 individus. Pour un locus `<A>` (diploïde,
`ploidy=2`), chaque individu compte pour 2 lignées → `k` démarre à
**240** au temps présent, pour CHACUNE des 5000 généalogies
indépendantes simulées (une par locus `<A>` du dataset) — ce nombre de
départ est le même pour les 5000 loci (même échantillonnage), mais
chaque généalogie évolue ensuite indépendamment (sa propre topologie,
son propre historique de coalescence).

## Généalogie vs topologie : laquelle porte le temps ?

Deux mots qu'on emploie beaucoup dans ce projet et qu'il ne faut pas
confondre — la convention est **l'inverse** de ce qu'on pourrait
intuiter :

- **Topologie** = la structure de branchement de l'arbre (qui est
  parent/enfant de qui), **sans les longueurs de branches, donc sans le
  temps**. C'est la définition explicite de tskit (dont ce projet
  dépend directement) : *"the topology of a tree in a tree sequence
  refers to the relationship among samples ignoring branch lengths"*
  ([tskit topological-analysis
  docs](https://tskit.dev/tskit/docs/stable/topological-analysis.html)).
- **Généalogie** (l'arbre complet, tel que retourné par
  `msprime.sim_ancestry`/`tskit.TreeSequence`) **inclut** le temps :
  chaque noeud porte un `.time`, et la longueur de branche se déduit de
  la différence de temps entre un noeud et son parent
  (`node_times[parents] - node_times[children]`, exactement ce que fait
  `_draw_single_mutation_edge_child` dans `ancestry_simulation.py`).
  tskit ne stocke d'ailleurs jamais les longueurs de branches
  directement, seulement les temps des noeuds — les longueurs de
  branches en sont toujours déduites.

Donc : c'est la **topologie** qui est "sans le temps" (juste la forme
de l'arbre), et c'est la **généalogie complète** qui porte le temps —
pas l'inverse. Cohérent avec des commentaires déjà présents dans le
code (ex: `compute_population_layout`, `"seule la topologie
coalescente varie d'un locus à l'autre"` — on veut dire que la *forme*
de l'arbre change à chaque réplicat, indépendamment du fait que les
temps changent eux aussi).
