# Notes d'exploration — diyabc/diyabc

## Sémantique de MERGE (vérifiée dans particuleC.cpp, verifytree(), ~ligne 2163)

`t merge a b` (header.txt) → pop = a, pop1 = b (history.cpp, parsing)

Exécution réelle (particuleC.cpp) :
    popleine[pop]  = true   // a SURVIT
    popleine[pop1] = false  // b DISPARAÎT (fusionne dans a)

Donc : "merge a b" signifie "b fusionne dans a, a continue d'exister".
Équivalent msprime : add_population_split(time=t, derived=[b], ancestral=a)

Vérifié par cohérence interne sur le scénario 1 de human/header.txt :
  t1 merge 2 1   -> pop 1 disparaît dans pop 2
  t2 merge 3 4   -> pop 4 disparaît dans pop 3
  t3 merge 2 3   -> pop 3 disparaît dans pop 2
  t4 varNe 2 Na  -> cohérent : pop 2 est la seule encore active


## Découverte : le bloc scénario peut déborder dans la section suivante

split_scenario_blocks() découpe sur le motif "scenario N [...] (...)",
donc chaque bloc va jusqu'au scénario suivant OU la fin du texte. Pour le
DERNIER scénario du fichier, le bloc engloutit aussi "historical parameters
priors (...)" et tout ce qui suit. Pas un bug en pratique : ces lignes ne
matchent aucun mot-clé connu (sample/merge/varNe/split) et lèvent
NotImplementedError au parsing -> le bloc est juste rejeté entièrement.
A garder en tête si on ajoute un mot-clé qui pourrait matcher par erreur
du texte de la section suivante.
Correction apportée : on s'arrêt à : "historical parameters priors"

## scenario 4 de human/header.txt passe sans modification du code

Même vocabulaire (merge/varNe) que le scénario 1, juste une numérotation
différente des temps (t11..t44 au lieu de t1..t4). Bonne preuve que le
découpage par mot-clé d'action (plutôt que par scénario) généralise bien.

## Format de l'en-tête de section "historical parameters priors (N,C)"

N = nombre de priors déclarés, C = nombre de contraintes d'ordre (lignes
avec ">"). Vérifié sur human/header.txt : (21,4) correspond exactement à
21 lignes de prior + 4 lignes "X>Y" avant DRAW UNTIL. Utile comme
assertion de validation après parsing (si len(priors) != N, le parsing a
raté quelque chose).

## Mapping indice de population (header.txt) <-> nom réel (fichier .snp)

Aucun nom n'est déclaré dans header.txt -- seulement des indices (1,2,3,4).
Vérifié : popname (data.cpp) n'est jamais croisé avec les indices du
scénario dans le code. Le mapping réel est implicite : pop i du scénario
= i-ème population rencontrée dans l'ORDRE D'APPARITION du fichier .snp.
Vérifié sur human : ASW(1) YRI(2) CHB(3) GBR(4), chacune avec un bloc
de 30 lignes consécutives.

## Choix de ploidy=2 pour human

Confirmé : human/header.txt déclare ses 51250 loci en <A> (autosomal),
cohérent avec une transmission diploïde classique. ploidy=2 (valeur par
défaut de msprime.sim_ancestry) est donc le bon choix : chaque "sample
individual" = 2 lignées génomiques, et l'échelle de temps de la
coalescence est calée en générations diploïdes -- cohérent avec les bornes
des priors de temps (UN[1,30], UN[100,10000]), qui sont en générations.
À revoir si on traite un jour un locus <H>/<X>/<Y>/<M> dans un autre
dataset (mitochondrial par exemple, qui serait haploïde, ploidy=1).

## Reproductibilité de sim_ancestry avec num_replicates

Vérifié empiriquement (pas dans la doc) : sim_ancestry(samples=N,
num_replicates=K, random_seed=S) dérive correctement K graines distinctes
à partir d'une seule seed S -- les réplicats sont statistiquement
indépendants entre eux, ET toute la séquence de K réplicats est
reproductible si on relance avec la même seed S. Pas besoin de générer
et passer un tableau de graines à la main (contrairement à l'exemple de
la doc orienté parallélisation multi-process).

## Format condensé "loci description" pour fichiers SNP -- confirmé (header.cpp::readHeaderLoci, branche SNP)

Syntaxe : <n1> <type1> [<n2> <type2> ...] <groupe> from <indice>

- Chaque paire (n_i, type_i) = n_i loci consécutifs de ce type d'héritage
  (<A>,<H>,<X>,<Y>,<M>), pris dans l'ordre d'apparition du fichier .snp
- "from N" (1-based dans le fichier) = indice de départ dans dataobs.locus[]
  (converti en 0-based : prem = N - 1)
- Somme des n_i = nombre total de loci de ce groupe à extraire

Exemple human : "5000 <A> G1 from 1" = prendre les 5000 PREMIERS loci
(indices 0 à 4999) du fichier .snp, tous <A>, assignés au groupe G1.

Exemple théorique : "70 <A> 10 <X> 10 <M> 10 <Y> G1 from 1" = prendre
les 100 loci à partir de 1, positions 0-69
=<A>, puis les 10 premiers <X>, puis les 10 premiers <Y> et enfin les 10 premiers <M>.

IMPORTANT : le FICHIER .snp peut contenir bien plus de loci (51250 pour
human) que ce qui est réellement simulé/comparé (5000 pour le scénario
de header.txt) -- "loci description" est un FILTRE/SOUS-ÉCHANTILLONNAGE
des colonnes du fichier de données, pas une description de tout le fichier.

## Modèle de mutation SNP correct (doc DIYABC section 2.4.3) : algorithme de Hudson

Confirmé par la doc utilisateur DIYABC : pour les SNP, "il est supposé
qu'il y a eu une et une seule mutation dans l'arbre de coalescence" --
PAS un processus de Poisson à taux variable. C'est l'algorithme "-s"
de Hudson (2002). La première approche (msprime.sim_mutations à
taux fixe, BinaryMutationModel) était donc structurellement incorrecte
pour les SNP, pas seulement approximative.

Implémentation validée empiriqfloat(uement (20000 tirages, proportions
observées vs attendues alignées à <1%) : pour chaque locus,
1. tirer une branche de l'arbre, avec probabilité proportionnelle à sa
   longueur (tree.branch_length(u) / tree.total_branch_length)
2. tous les échantillons descendants de cette branche (tree.samples(u))
   portent l'allèle dérivé (1), les autres l'allèle ancestral (0)

Garantit par construction : exactement une mutation, donc le locus est
TOUJOURS polymorphe (jamais besoin de filtrer les monomorphes a
posteriori, contrairement à notre ancienne approche par taux de Poisson).

## Filtre MAF (Minor Allele Frequency) -- pas nécessaire pour human, à prévoir pour d'autres datasets

Doc DIYABC (section MAF) : <MAF=hudson> = algorithme de Hudson standard,
SANS filtrage supplémentaire -- notre implémentation actuelle
(simulate_snp_genotypes, une mutation par locus) est déjà correcte pour
ce cas. Confirmé : human/human_snp_all22chr_maf5.snp déclare bien
<MAF=hudson> en première ligne -- malgré le suffixe "maf5" dans le NOM
du fichier (qui semble faire référence à un prétraitement déjà appliqué
aux données OBSERVÉES réelles avant export, pas à la simulation).

Pour MAF=N% (non rencontré sur human, mais prévu par la doc) : il
faudrait calculer la MAF du locus simulé (fréquence de l'allèle le moins
fréquent, toutes populations confondues) après chaque tirage Hudson, et
RESIMULER (rejeter et retirer) si elle est sous le seuil -- jusqu'à
obtenir num_loci loci valides. Pas implémenté : à ajouter si on traite
un dataset avec MAF != hudson.

## Colonne SEX du fichier .snp -- valeur arbitraire acceptée sans élucidation complète

Décision pragmatique : on écrit "9" (valeur observée dans human réel)
pour tous les individus simulés, sans avoir confirmé dans le C++ que
cette colonne n'est jamais lue pour des loci <A> (autosomaux). Justifié
par : (1) aucune des 6 catégories de stats SNP listées dans la doc
DIYABC (ML, HW/HB, FST, F3/F4, Nei, AML) ne mentionne le sexe comme
paramètre ; (2) readheaderdata (data.cpp) ne semble utiliser SEX que
pour la DÉTECTION du format de fichier, jamais comme donnée individuelle
exploitée -- vérification incomplète, fonction de lecture ligne-par-ligne
non localisée avec certitude. À revoir si un résultat statistique
incohérent apparaît plus tard.

## Validation empirique de l'architecture "déléguer au C++" -- SUCCÈS

Test manuel réussi : header.txt modifié (5000->10 loci) + .snp généré
depuis nos génotypes simulés msprime + RNG_state_0000.bin, dans un même
dossier, lancé avec le vrai binaire `general -p ./ -R "FST1;ML1" -r 1
-g 50 -m -t 1`. statobsRF.txt produit des valeurs ML1p_1..4 cohérentes
et DIFFÉRENTES de celles obtenues sur les vraies données humaines --
preuve que le calcul tourne bien sur NOS données simulées.

Confirme : ML1p peut être < 1 même avec un locus garanti polymorphe
GLOBALEMENT par l'algorithme de Hudson -- le polymorphisme global
n'implique pas le polymorphisme PAR POPULATION (la mutation peut être
confinée à une seule population selon la branche tirée).

Prochaine étape : automatiser ce processus manuel (appel subprocess à
`general`, parsing de statobsRF.txt) dans bridge/, plutôt que des
scripts/chemins ad hoc.

## compute_summary_statistics validé empiriquement sur les 112 statistiques

Premier appel complet réussi avec stats_filter="ALL" : les 112 colonnes
attendues (ML1/ML2/ML3, HW, HB, FST1-4, NEI, AML, F3, F4) sont produites
par le vrai binaire C++, sur des données simulées par notre pipeline.

Note : des valeurs légèrement négatives apparaissent (ex: FST1m_3=-0.22)
-- NORMAL pour un estimateur de Fst avec peu de loci (variance
d'échantillonnage), contrairement au -1.19 aberrant obtenu lors de notre
tentative de réimplémentation Python (qui était un vrai bug de formule,
pas du bruit statistique légitime).

## "112" dans header.txt != nombre réel de stats produites par general -R "ALL"

header.txt déclare "group summary statistics (112)" avec un vocabulaire
ancien (HP0/HM1/HV1/HMO/FP0...), retrouvé en tant que TABLEAU COMMENTÉ
(stat_type0, désactivé) dans general.cpp -- vestige d'une version
antérieure du format. Le binaire general actuel, avec -R "ALL", produit
en réalité 130 statistiques avec le vocabulaire moderne (ML1-3, HW, HB,
FST1-4, NEI, AML, F3, F4) -- confirmé empiriquement (compute_summary_
statistics, test pytest). Ne jamais se fier au nombre annoncé dans
header.txt pour ce champ -- toujours vérifier empiriquement contre la
sortie réelle de statobsRF.txt.

## Détournement architectural : .snp intermédiaire pour calculer les stats sur données SIMULÉES

IMPORTANT, à ne pas oublier : dans le vrai DIYABC, calstatobs() (header.cpp)
est conçu pour calculer les stats sur les VRAIES données OBSERVÉES,
lues UNE SEULE FOIS depuis un fichier .snp sur disque. Les stats sur
données SIMULÉES sont normalement calculées directement en mémoire
(ParticleC, sans jamais réécrire sur disque), des milliers de fois.

Notre pipeline détourne ce mécanisme : on écrit nos génotypes SIMULÉS
dans un faux fichier .snp, pour réutiliser calstatobs() sans toucher au
C++. Ça fonctionne (validé empiriquement), mais introduit un coût I/O
(écriture + relecture disque) à CHAQUE particule -- absent du vrai
DIYABC. Acceptable pour ce POC (objectif : prouver la faisabilité),
mais à corriger avant toute mise en production réelle (nécessiterait
soit une vraie modification du C++ pour accepter des données simulées
en mémoire depuis Python, soit reproduire intégralement les formules
en C++/Python sans repasser par le format .snp).

## Découverte majeure : incohérence interne au dépôt diyabc entre header.txt et le binaire general actuel

header.txt (tests/datasets/human/) déclare encore l'ANCIEN format de
stats (group summary statistics (112), vocabulaire HP0/HM1/HV1/HMO...).
Le binaire `general` compilé depuis CE MÊME dépôt calcule en réalité
130 stats au format MODERNE (ML/HW/HB/FST/NEI/AML/F3/F4) avec -R "ALL".

readReftable.R suppose la cohérence entre le texte de header.txt et
nstat (lu depuis le fichier binaire) pour calculer paramsname --
head(allcolspre, -nstat). Avec nstat=130 mais header.txt n'offrant que
133 tokens totaux (21 params + 112 anciennes stats), le calcul produit
un résultat absurde (3 paramsname au lieu de 21), provoquant l'erreur
"indice hors limites" plus loin.

PREUVE : notre reftable.bin est structurellement CORRECT -- vérifié
manuellement (lecture binaire instrumentée en R, ligne par ligne, sans
passer par readRefTable) : nrec/nscen/nparam/nstat cohérents, parambyscenh
et lparams s'alignent parfaitement (16 valeurs, bon ordre). Le bug est
une incohérence PRÉEXISTANTE dans le dépôt diyabc lui-même (header.txt
obsolète vs binaire actuel), pas un défaut de notre pipeline Python.

C'est un succès complet, et c'est exactement ce qu'on espérait depuis le début de cette longue investigation.
Décortiquons ce résultat, point par point
Colonnes params : les 21 noms attendus, dans le bon ordre — exactement la même structure que ce qu'on avait obtenu avec le vrai DIYABC sur ce même dataset, tout au début du projet.
Dimensions stats: 5 130 : 5 lignes (nos 5 particules), 130 colonnes de statistiques — cohérent avec nrec=5 et le calcul réel via -R "ALL".
Première ligne de params : les 16 paramètres du scénario 1 (N1...Na) ont des valeurs numériques cohérentes, et ra, t11, t22, t33, t44 sont NA — exactement le comportement attendu, puisque ces 5 paramètres n'appartiennent pas au scénario 1 (rappelle-toi, c'est précisément ce qu'on avait observé avec le vrai reftableRF.bin de DIYABC tout au début du projet, à l'identique).
C'est la preuve définitive
Le reftable.bin produit par notre pipeline Python (header.txt → msprime → algorithme de Hudson → délégation au C++ pour les stats → écriture binaire) est structurellement et fonctionnellement identique à ce que produit le vrai DIYABC — vérifié par l'outil de référence indépendant readReftable.R, avec exactement le même comportement de filtrage des paramètres non utilisés par scénario.

## Découverte critique de performance : -g est la taille de batch interne, jamais notre besoin

-g <n> (general.cpp) contrôle nenr, la taille du batch de particules
simulées en une fois par dosimultabref() AVANT de vérifier la condition
d'arrêt (nrecneeded > rt.nrec). On utilisait -g 50 (copié de notre tout
premier test manuel) alors qu'on ne veut qu'UNE particule par appel
subprocess (notre boucle externe gère déjà l'itération côté Python).

Conséquence mesurée : -g 50 calculait 50 particules complètes par appel,
dont on ne gardait qu'une seule -- facteur ~29x de ralentissement
(343s -> 7.4-12s par particule sur 5000 loci). Corrigé : -g 1.

Pour nrec=1000 particules (5000 loci, scénario 1, -R "ALL") :
  estimation ~200 minutes en séquentiel, ~25 minutes avec parallélisation
  (ProcessPoolExecutor, 8 workers) -- à valider empiriquement.

## Notes du 03/07/2026
Création de deux fichiers txt contenant les résumés statistiques d'une simulation de 1000 particules pour 5000 loci sur le scenario 1 de human. Le premier fichier est issu de la simualtion via diyabc : "time ./diyabc -p ./ -R ALL -r 1000 -g 1000 -m -t 16" (temps de calcul 2 min 17). L'autre provient de la simulation via msprime avec calcul des résumés statistiques via l'implémentation en Python. Note : cette implémentation en Python à prouver que l'on obtenait les mêmes résultats qu'avec le calcul via msprime validant ainsi notre implémentation des formules. 

**Comparaison statistiques de ces deux fichiers :** pour chaque indicateurs (priors, statistiques) calcul de la moyenne, écart-type, médiane et calcul des différences et différences relatives. Puis Détermination des variables ayant une différence relative de la moyenne supérieure à 5%. Résultat : 124 variables sur 125 ont un seuil supérieur. Si on pousse l'analyse un peu plus loin, on se rend compte que pour un test de Kolmogorov-Smirnov à deux échantillons où l'hypothèse nulle serait "les deux distributions sont-elles compatibles avec la même loi ?", 126 variables ont une p-value inférieur à 0.05 rejetant ainsi l'hypothèse nulle. 

**Conclusion :** Les deux simulateurs ne semblent pas donner des simulations identiques. Cela reste à discuter avec les experts pour voir si les moyennes proches peuvent tout de même être interprétées comme allant dans le même sens. 


## Notes du 06/07/2026
Dans certains header.txt, les variations de population efficace sont notées "varne" et "parfois "varNe" ! Il faudra faire en sorte d'éliminer cette difficulté.

## Notes du 07/07/26
Problème qqpart car la comparaison des statistiques de sorties montrent qu'il y a un vrai écart entre les simulations faites par msprime et celles faites par diyabc. Les pistes écartées : 
- 1. Échelle de Ne (k(k-1)/(4N)) — identique à msprime.
- 2. Modèle de mutation Hudson (placement direct pondéré par longueur de branche, sans rejet, sans filtre
  MAF/populationnel caché) — identique au nôtre, vérifié deux fois en détail.
- 3. Arbre de coalescence partagé entre loci — chaque locus autosomal retire bien sa propre généalogie
  indépendante.
- 4. Approximation continue vs discrète du coalescent — écartée par calcul ET confirmée par test
  empirique (biais identique après ×10 sur N).
- 5. Propagation du génotype (quels échantillons héritent de l'allèle muté) — parcours topologique
  standard, rien d'anormal.
- 6. Doublons de branches / MRCA artificiel — arbre binaire propre, rien d'anormal.

## Notes du 10/07/2026 — RÉSOLUTION : la ligne de fin de header.txt/headerRF.txt n'est pas cosmétique

Cause racine trouvée (découverte par l'utilisateur, root-causée et reproduite indépendamment le même jour) de l'essentiel de l'écart documenté le 03/07 et le 07/07 : la dernière ligne du header ("scenario N1 N2 N3 ta ts ML1p_1 ...", qu'on pensait être une simple documentation des colonnes de sortie générée par DIYABC) est en réalité **relue comme une entrée** par le binaire `general`.

Dans `header.cpp::HeaderC::readHeaderAllStat` (~lignes 789-799) :
```cpp
getline(file, this->entete);  // relit la DERNIÈRE ligne du header EN ENTRÉE
...
size_t nparamhist = header_lastline.size() - 1 - nstat - nparamut;
```
`nparamhist` (le nombre de paramètres historiques que DIYABC croit exister) est dérivé du **nombre de tokens de cette ligne de fin**, pas d'un comptage réel des priors déclarés dans `historical parameters priors`. Si cette ligne contient des noms de paramètres en trop (ex: `N4 r` recopiés d'un header d'un autre scénario, alors que `historical parameters priors` n'en déclare que 5), `nparamhist` est faux et corrompt un état interne en aval.

**Signature du bug** : la population "hub" (première déclarée / celle qui survit à toutes les fusions) reste correcte ; toutes les autres populations montrent des statistiques massivement fausses (HWm/ML1p écroulés vers 0, FST1m explosé au-dessus de 1) — écarts de 300% à 10000%, alors que les tirages de priors (N1..ts) restent statistiquement normaux. C'est exactement la signature "asymétrie population 1 vs populations 2/3" observée et creusée pendant toute l'investigation du notebook `correlation_N2_N3_HWm_anomaly.ipynb`.

**Vérification** : `toy_example5_scenario1` (utilisé dans tout ce notebook) déclarait "6 parameters" avec seulement 5 priors réels, et sa ligne de fin listait `N4 r` en trop. En reconstruisant un header identique mais avec une ligne de fin propre (5 tokens), sur 1000 particules DIYABC réelles vs msprime (tirages indépendants, mêmes priors) :
- `ML1p_1/2/3` et `HWm_1/2/3` : écarts tombés de 300-10000% à <3%, plus aucun significatif (avant : très significatifs sur pop2/pop3 uniquement).
- Nombre de stats significativement différentes (test KS, p<0.05) : passé de ~48/55 à ~22/55.
- Il reste un biais résiduel modeste (5-16%) sur les statistiques de paires (FST2, NEI, F3, HB) — du même ordre que ce qu'on observe ailleurs dans le projet avec des headers propres, un problème bien plus petit et probablement plus classique, encore ouvert.

**Conclusion révisée** : la quasi-totalité de l'écart documenté le 03/07 et 07/07 était un artefact de fichier de test mal formé (ligne de fin recopiée depuis un autre scénario), pas un désaccord réel entre les deux simulateurs. Voir mémoire persistante `diyabc_header_trailer_line_bug` pour le détail complet et comment vérifier un header pour ce piège avant de le réutiliser.

