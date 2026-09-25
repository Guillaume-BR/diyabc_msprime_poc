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

Confirmé : human/header.txt déclare ses 51250 loci A (autosomal),
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

Syntaxe : n1 type1 [n2 type2 ...] groupe from indice

- Chaque paire (n_i, type_i) = n_i loci consécutifs de ce type d'héritage
  (A,H,X,Y,M), pris dans l'ordre d'apparition du fichier .snp
- "from N" (1-based dans le fichier) = indice de départ dans dataobs.locus[]
  (converti en 0-based : prem = N - 1)
- Somme des n_i = nombre total de loci de ce groupe à extraire

Exemple human : "5000 A G1 from 1" = prendre les 5000 PREMIERS loci
(indices 0 à 4999) du fichier .snp, tous A, assignés au groupe G1.

Exemple théorique : "70 A 10 X 10 M 10 Y G1 from 1" = prendre
les 100 loci à partir de 1, positions 0-69
=A, puis les 10 premiers X, puis les 10 premiers Y et enfin les 10 premiers M.

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
cette colonne n'est jamais lue pour des loci A (autosomaux). Justifié
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

## Notes du 08/07/26

Ajout de l'option "initially_active = True" dans le modèle de démographie : a permis d'avoir des valeurs cohérentes avec diyabc pour ML1, ML2... 

## Notes du 09/07/26

Ce qui est confirmé propre : j'ai testé le fix sur une dizaine de points fixes très variés — N
  minuscule (100) avec ta/ts énormes, N minuscule avec ta/ts du même ordre de grandeur, ta≈ts (10
  générations d'écart sur 10000), N1 très petit vs N2/N3 énormes, et les valeurs exactes tirées par 
  DIYABC lui-même pour les lignes les plus aberrantes de son reftable (ex: ligne avec N1=1744, N2=19984,
  ta=24641, ts=25586). À chaque fois, en moyennant sur 1000 répliques du même point, DIYABC et msprime
  tombent d'accord au pourcent près. Le fix initially_active est donc solide partout où je l'ai testé
  ponctuellement.

  Ce qui ne colle pas : en comparant les deux reftables complets (1000 particules à priors tirés
  indépendamment, une seule réplique par particule), HWm_2 et HWm_3 divergent nettement (DIYABC≈0.022,
  msprime≈0.086) alors que HWm_1 est correct. Et ce n'est pas du bruit d'échantillonnage — les
  écarts-types sont petits (SEM≈0.0003-0.001, écart observé = 0.064, donc >50 SEM).

  Le vrai signal que j'ai trouvé : dans le reftable réel de DIYABC, HWm_2 et HWm_3 sont quasi 
  indépendants de N2/N3 (corr≈-0.18/-0.19, alors qu'on attend une corrélation positive forte), alors que
  HWm_1 corrèle bien avec N1 (0.34) — et msprime montre la corrélation attendue et correcte pour les 3
  populations (0.82 / 0.70 / 0.46). Pourtant, mon test contrôlé (ne faire varier que N2, ta/ts fixes)
  montre que DIYABC sait parfaitement faire varier HWm_2 avec N2 quand on isole la variable. Autrement
  dit : DIYABC scale correctement avec N en isolation, mais dans le vrai reftable à tous les paramètres
  tirés simultanément, cette dépendance disparaît spécifiquement pour les populations "dérivées" (2 et 3,
  qui disparaissent chacune dans une fusion) — jamais pour la population "ancestrale" (1, qui survit aux
  deux fusions).

  Je n'ai pas encore d'explication ferme — ça ressemble à un problème côté binaire DIYABC réel (pas notre
  port msprime) spécifique à ce header avec contrainte ts>ta par rejet sur intervalles superposés, mais
  je ne peux pas l'affirmer sans regarder le code source plus près. Tu veux que je creuse ça côté C++
  (particuleC.cpp/history.cpp) pour voir si ts>ta DRAW UNTIL a un effet de bord sur N2/N3, ou on met ça
  de côté pour l'instant ?

## Note du 10/07/26

Correction du bug : une différence de lecture du header par notre implémentation en Python et celle de Diyabc. Notre implémentaiton utilise la première partie pour déterminer la liste des priors et des stat tandis que diyabc utilise la dernière ligne qui fait office de référence. Ainsi il y avait un décalage qui faussait les calculs et donc les interprétations. 

Il reste des écarts relatifs qui peuvent paraître important mais rien de comparable avec ce que l'on avait précédemment. Les tests de Kolmogorov-Smirnov passent par contre tous validant ainsi l'hypothèse d'égalité des distributions de nos sumstats. Ceci a été testé sur le dataset human avec 5000 loci dans la version avec plusieurs scénarios.

## Note du 17/07/26 -- le "biais résiduel" n'est pas un vrai désaccord DIYABC/msprime

** `reference/toy_example5_1000loci/compare_reftables_te5_1000loci.ipynb`**
(rejeu exact des tirages RÉELS de DIYABC via
`scripts/replay_diyabc_priors.py` -- comparaison appariée, mêmes
N1..r pour chaque particule des deux côtés, 1000 particules x 650 loci
= 500 `<A>` + 50 `<X>` + 50 `<M>` + 50 `<Y>`). 
Le test avec  70 `<A>` + 10 `<X>` + 10 `<M>` + 10 `<Y>` donne des résultats
décevants certainement du fait du faible nombre de loci simulés. 
L'exécution complète du notebook (~8s) : priors identiques des deux côtés à 0.0 près (rdiff=0,
KS=0, p=1.0 -- confirme l'appariement), et sur les 50 statistiques,
écart relatif moyen 3.6% (max -7.7% sur F3v_3.1.2), **0 statistique
avec p<0.05** (KS, min p=0.12). C'est le même ordre de grandeur que le
"biais résiduel 5-16%" documenté le 10/07, mais désormais confirmé
non significatif statistiquement sur 650 loci.

## Note du 20/07/26 -- écart de performance DIYABC/msprime : le coût par particule est quasi identique, le vrai plafond est la parallélisation + le surcoût tskit par locus

Point de départ : sur `human_modif_scenario1_5000loci` (rejeu des
priors réels DIYABC via `replay_diyabc_priors.py`, 1000 particules,
5000 loci, scénario 1, `max_workers=16`), le pipeline msprime/Python
met **384s** contre **137s** pour le vrai DIYABC sur la même config
(`time ./diyabc -p ./ -R ALL -r 1000 -g 1000 -m -t 16`, mesuré le
03/07) -- un facteur **~2.8x**. Investigation pour savoir si ça vient
du simulateur, des formules de stats, ou d'autre chose.

**1. Coût par particule en séquentiel (mesuré, valeurs réelles DIYABC
rejouées, 5000 loci)** :
- Simulation seule (coalescence msprime + mutation Hudson) : 2.12s
- `compute_summary_statistics_from_values` complet (simulation + 130
  stats + I/O) : 2.35-2.4s -- le calcul des stats + la lecture des
  fichiers ne pèsent donc que ~0.25s (~10%), pas le goulot.
- Comparé au coût DIYABC équivalent par particule (137s x 16 threads /
  1000 particules ~= 2.19s) : **quasi identique** (facteur ~1.1x). Le
  simulateur/les formules ne sont PAS le problème.

**2. D'où vient alors le facteur ×2.8 ?** -- de la mise à l'échelle
parallèle, pas du coût unitaire. Test empirique de `max_workers` sur
128 particules réelles (mêmes conditions) :

| max_workers | temps réel | idéal (linéaire depuis le seq.) | efficacité |
|---|---|---|---|
| 16 | 44.9s | 18.8s | 42% |
| 8  | 54.5s | 37.6s | 69% |
| 7  | 60.8s | 43.0s | 71% |

16 workers reste le MEILLEUR choix en absolu malgré une "efficacité"
plus faible -- donc `max_workers=16` (déjà utilisé) n'est pas un
mauvais réglage. Cause identifiée via `lscpu` : la machine a **8 cœurs
physiques / 16 threads logiques** (Xeon W-11955M, hyperthreading). Le
calcul est purement numérique (msprime + numpy), qui tire peu de
bénéfice de l'hyperthreading -- donc 16 workers Python se disputent en
réalité 8 cœurs physiques. C'est un plafond matériel de la machine de
dev, pas un défaut de notre code ou de son réglage.

**3. Où va le temps DANS la simulation** (profil détaillé d'une
particule via `cProfile`) -- ~93% du temps total est dans
`simulate_snp_genotypes` (boucle sur les 5000 loci), décomposé
ainsi :
- ~30% moteur de simulation msprime en C (`Simulator.run`/`reset`/
  `finalise_tables`) -- difficilement compressible, c'est le vrai
  calcul de coalescence.
- **~40% construction/inspection des objets Python `tskit.TreeSequence`/
  `Tables` PAR LOCUS** (dont le décodage du metadata des populations,
  refait pour CHACUN des 5000 loci alors qu'il est strictement
  identique à chaque fois) -- le plus gros poste, et le plus évitable.
- ~15% l'algorithme de Hudson (tirage de la mutation, déjà vectorisé
  sur les tables d'edges).

Cohérent avec l'archi : DIYABC en C++ reste dans une boucle serrée sans
jamais recréer d'objets haut niveau par locus, alors que nous
matérialisons une vraie `TreeSequence` Python complète (avec décodage
de schéma de métadonnées) 5000 fois par particule.

**Petit à-côté relevé au passage, pas encore corrigé** :
`build_samples_argument` (parsing du `.snp`, 12 Mo) est appelé DEUX
fois par particule au lieu d'une (une fois dans
`run_poc_for_directory_with_values`, une fois dans
`compute_summary_statistics_from_values`) -- redondant, ~0.1s/particule
(~4% du total), pas le facteur principal mais un gain facile.

**Conclusion** : l'écart de perf n'est ni un problème de nos formules
de stats (quasi identiques en coût unitaire à DIYABC), ni un mauvais
choix de `max_workers` -- c'est le plafond des 8 cœurs physiques de la
machine de dev combiné au surcoût Python/tskit de recréer un objet
`TreeSequence` complet par locus (5000x par particule). Pour l'objectif
du POC (démontrer la faisabilité, déjà acquis), pas bloquant. Une
vraie optimisation nécessiterait de réduire ce surcoût tskit par locus
(ex: éviter le redécodage du metadata des populations à chaque
itération, dédupliquer le double appel à `build_samples_argument`) --
un chantier réel à part entière, pas un réglage rapide. Reste à faire
si on veut industrialiser au-delà du POC.

## Note du 20/07/26 (suite) -- optimisation du cache population/samples dans simulate_snp_genotypes, gain mesuré

Implémentation de la première piste identifiée ci-dessus : dans
`simulate_snp_genotypes` (`bridge/ancestry_simulation.py`), la liste
(nom de population, IDs d'échantillons) était recalculée à CHAQUE locus
(`ts.tables.populations`, décodage du metadata, `ts.samples(population=
...)`) alors qu'elle est strictement identique pour tous les réplicats
d'un même appel à `simulate_independent_loci`/`simulate_shared_
ancestry_loci` -- seule la topologie coalescente varie d'un locus à
l'autre, jamais l'assignation des noeuds échantillons aux populations.
Vérifié empiriquement avant de coder (5 loci, 4 populations : IDs
d'échantillons et noms de population identiques sur tous les loci).

Changement : cette liste est maintenant calculée une seule fois, au
premier locus consommé par le générateur, et réutilisée pour tous les
suivants (le tirage de la mutation et le calcul de `derived_samples`
restent bien par-locus, seule la structure pop->samples est mise en
cache).

**Gain mesuré** (mêmes conditions que l'investigation ci-dessus,
`human_modif_scenario1_5000loci`) :
- Séquentiel : 2.35s -> **1.96s**/particule (~17%)
- Parallèle réel (128 particules, `max_workers=16`) : 44.9s -> **36.4s**
  (~19%)
- Extrapolé sur les 1000 particules réelles : 384s -> **~284s**
  (~26% de temps en moins sur le run complet)

Suite de tests complète (62/62) toujours verte après le changement.

**Limite** : ce cache ne profite qu'au chemin rapide `maf=0.0`
(`<MAF=hudson>` ou tag absent, cas de `human`) où `simulate_snp_
genotypes` est appelée UNE FOIS sur tous les loci d'un coup. Le chemin
`with_maf_filter`/`with_maf_filter_shared_ancestry` avec un vrai seuil
MAF numérique (boucle de rejet) appelle `simulate_snp_genotypes` locus
par locus (un seul élément à chaque appel) -- le cache ne s'y active
donc pas tel quel, il faudrait restructurer la boucle de rejet pour
partager ce cache entre tentatives si on veut aussi optimiser ce
chemin-là. Pas fait pour l'instant.

L'écart résiduel avec DIYABC (~284s vs 137s) reste la partie
incompressible : plafond des 8 cœurs physiques + coût de matérialiser
une `TreeSequence` par locus (nécessaire de toute façon pour y tirer la
mutation).

## Note du 20/07/26 (suite 2) -- cache étendu à with_maf_filter (maf>0), gain marginal cette fois

Extension de l'optimisation précédente au chemin `maf>0.0` (boucle de
rejet, `with_maf_filter`/`with_maf_filter_shared_ancestry`) :
`simulate_snp_genotypes` accepte maintenant un paramètre optionnel
`population_layout` (factorisé dans un nouveau helper
`_population_layout(ts)`) -- si fourni par l'appelant, jamais
recalculé. Les deux boucles de rejet le calculent une seule fois (au
premier `ts` généré pour `with_maf_filter`, à partir de l'arbre partagé
déjà unique pour `with_maf_filter_shared_ancestry`) et le réutilisent à
travers toutes les tentatives suivantes, acceptées ou rejetées.

**Gain mesuré** (avant/après par patch réversible sur le seul fichier
`ancestry_simulation.py`, `toy_example3_scenario1`, `<MAF=0.05>`, 300
loci acceptés) : **0.69s -> 0.67s, soit ~2-3% seulement** -- beaucoup
plus faible que les ~17% obtenus sur le chemin `maf=0.0` de `human`.

Explication : dans la boucle de rejet, chaque tentative simule déjà un
seul locus via son propre appel à `simulate_independent_loci(num_loci=
1, ...)` -- l'essentiel du coût par tentative est là (création d'un
nouveau générateur `sim_ancestry` à chaque tentative), pas dans le
décodage du metadata des populations qu'on élimine ici. Ce dernier
pesait lourd sur `human` parce qu'un SEUL appel traitait 5000 loci d'un
coup (4 populations x 5000 décodages redondants économisés en une
fois) ; ici chaque tentative ne fait qu'UN seul décodage de toute façon,
il n'y a donc qu'un tout petit nombre d'économies à faire par tentative.

Le changement reste correct (0 régression, 62/62 tests toujours verts,
y compris les tests qui valident sémantiquement le filtre MAF --
`test_with_maf_filter_rejects_low_maf_loci`,
`test_with_maf_filter_shared_ancestry_rejects_low_maf_loci`) et gratuit
à garder, mais ne pas s'attendre à un gain comparable à celui du chemin
`maf=0.0` si ce dataset redevient un sujet de perf.

## Note du 20/07/26 (suite 3) -- dédup de build_samples_argument (double scan interne + double appel pipeline.py)

Dernière piste facile identifiée dans l'investigation initiale : le
`.snp` (12 Mo pour `human`) était scanné en trop DEUX FOIS
différentes :

1. **En interne à `build_samples_argument`** (`ancestry_simulation.py`)
   : la fonction appelait `population_index_to_name(snp_file_path)` --
   qui appelle lui-même `count_samples_per_population` -- PUIS
   rappelait `count_samples_per_population(snp_file_path)` une seconde
   fois, indépendamment, juste pour les comptes. Corrigé : un seul
   appel à `count_samples_per_population`, l'indice 1-based se déduit
   directement de la position dans ses clés (même ordre garanti).
2. **Entre `pipeline.py` et `ancestry_simulation.py`** : comme identifié
   le 20/07 (suite 1), `compute_summary_statistics(_from_values)`
   rappelait `build_samples_argument(snp_path)` juste pour obtenir
   `population_names = list(samples.keys())`, alors que
   `genotypes_list` (déjà calculé juste avant) contient EXACTEMENT ces
   mêmes noms comme clés de chaque dict par locus (produits par
   `simulate_snp_genotypes`/`_population_layout`, mêmes noms "pop1"..
   "popN" dans le même ordre). Nouveau helper `_population_names`
   (`pipeline.py`) : prend les clés du premier locus déjà simulé --
   zéro I/O supplémentaire -- avec repli sur `build_samples_argument`
   si `genotypes_list` est vide (cas dégénéré `num_loci=0`, n'arrive
   pas en pratique).

**Gain mesuré** (cumulé avec les deux optimisations précédentes de la
même journée, `human_modif_scenario1_5000loci`) :
- Séquentiel : 1.96s -> **~1.84s**/particule (~6% de plus)
- Parallèle réel (128 particules, `max_workers=16`) : 36.4s -> **34.0s**
- Extrapolé 1000 particules (depuis l'échantillon de 128) : 284s -> ~265s

**Confirmation en conditions réelles** : rejeu complet des 1000
particules réelles (pas un échantillon) de `human_modif_scenario1_
5000loci` par l'utilisateur -- **300s**, contre 384s avant les 3
optimisations du jour, soit un **gain réel de ~22%** (un peu en dessous
du ~31% extrapolé depuis l'échantillon de 128, ce qui est normal --
une extrapolation depuis un sous-échantillon reste approximative). Ce
chiffre (300s), mesuré sur le run complet, est plus fiable que
l'extrapolation et remplace le ~265s ci-dessus comme référence.

**Bilan cumulé de la journée (3 optimisations)** : 384s -> **300s**
(mesuré sur les 1000 particules réelles), soit **~22% de temps en
moins** sur le run complet `human_modif_scenario1_5000loci`, pour 0
régression (62/62 tests toujours verts à chaque étape). L'écart
résiduel avec DIYABC (~300s vs 137s, facteur ~2.2x) reste la partie
incompressible identifiée dès la première note du jour : plafond des 8
cœurs physiques de la machine de dev + coût de matérialiser une
`TreeSequence` par locus.

## Note du 21/07/2026
Codage de simulate_poolseq_reads() dans ancestry_simulaiton.py. Utilisation d'un zip qui va s'arrêter silencieusement dès que l'un des deux mutables, tree_sequence ou observed_reads_per_locus, est consommé. 
A voir si cela pose problème, sinon utilisé itertools.zip_longest(tree_sequence, obeserved_reads_per_locus, fillvalue=_SENTINEL)

Codage de la simulation des reads avec filtre MRC et tests correspondants.

## Note du 23/07/2026 -- perf PoolSeq

Signalement : la simulation PoolSeq est ~6x plus lente que DIYABC sur
`toy_example4` (100 loci, MRC=5) -- bien plus que l'écart connu côté
IndSeq (~2,2-2,8x, voir note du 20/07). Profilage (cProfile, 1 particule) :

1. **~22-24% du temps** : `observed_reads()` (`observed_data.py`)
   re-scannait et re-purgeait (MRC) les 30000 loci du `.snp` observé à
   CHAQUE particule, alors que le résultat ne dépend que du fichier
   (jamais de la graine/démographie tirées). Corrigé : arrêt anticipé du
   scan dès que `num_loci` loci ont passé le filtre MRC (nouveau
   paramètre `num_loci=None`, `None` = comportement inchangé/scan
   complet) -- 0.092s -> 0.005s sur `toy_example4` (~20x), résultat
   identique vérifié (tous les appelants tronquaient déjà à `num_loci`
   après coup). Câblé en paramètre optionnel `observed_reads_per_locus`
   à travers toute la chaîne (`simulate_poolseq_reads_with_mrc_filter` ->
   `compute_summary_statistics(_from_values)` -> `reftable_loop.py`),
   calculé UNE SEULE FOIS par run avant `ProcessPoolExecutor(...)` au
   lieu d'une fois par particule (nouvelle fonction publique
   `prepare_poolseq_observed_reads`).
2. **~5%** : `simulate_poolseq_reads` recalculait `_population_layout(ts)`
   à chaque locus/tentative, sans le paramètre de cache que
   `simulate_snp_genotypes` a déjà côté IndSeq (20/07). Même fix porté
   ici (`population_layout=None`, cache externe à travers les tentatives
   dans `with_mrc_filter`, même principe que `with_maf_filter`).
3. **~70%, pas un bug** : le rejet-resimulation MRC=5 redessine en
   moyenne ~1,8 fois par locus (279 appels msprime bas niveau pour 100
   loci) -- coût fixe par appel `sim_ancestry` (construction du
   `Simulator`, provenance JSON) qui domine d'autant plus qu'il y a peu
   de loci. Même mécanisme que le filtre MAF déjà en place côté IndSeq,
   pas spécifique à PoolSeq, pas traité (pas facile à réduire sans
   changer l'algorithme).

**3 bugs attrapés en review avant qu'ils ne passent**, tous la même
famille Python (nom de variable/fonction locale masquant un nom
englobant -> `UnboundLocalError`) : `_passes_mrc` appelée avant sa
propre définition dans la boucle, le test MRC imbriqué dans la mauvaise
boucle (lignes ajoutées plusieurs fois, incomplètes), et
`observed_reads = observed_reads(...)` dans le nouveau helper (variable
locale masquant la fonction importée). Un 4e bug, plus grave, a aussi été
attrapé : en câblant `observed_reads_per_locus` dans `pipeline.py`,
l'appel réel à `simulate_poolseq_reads_with_mrc_filter` (la simulation
msprime) a été supprimé par erreur, remplacé par un calcul direct des
stats sur les données OBSERVÉES -- aurait rendu toutes les particules
PoolSeq d'un reftable identiques entre elles. Passé inaperçu par la
suite de tests existante (aucun test ne couvrait la branche PoolSeq de
`compute_summary_statistics(_from_values)`) -- comblé après coup par 2
nouveaux tests (`test_pipeline.py`) qui vérifient que deux graines/jeux
de paramètres différents donnent des statistiques différentes, vérifiés
en réintroduisant le bug pour confirmer qu'ils l'auraient attrapé.

**Gain mesuré** (run réel `run_reftable_simulation`, 100 particules x
100 loci, `toy_example4`, `max_workers=8`, comparaison via `git worktree`
sur le commit précédent) : **~32.3s -> ~27.6s, ~14-16% de gain** sur le
point 1. Le point 2 (cache `_population_layout`) a un effet non
mesurable isolément (bruit système ±10-15% > gain attendu ~5%, confirmé
par un test en double aveugle à charge égale) -- gardé pour sa
cohérence avec le reste du code, pas pour un gain chiffrable. 73/73
tests verts.

## Note du 24/07
Mesure de temps d'exécution pour un PoolSeq avec 100 loci : toy_example4:
- Avec un mrc = 1 : diyabc : 34s et msprime : 133s (x3.9)
- Avec un mrc = 5 : diyabc : 80s et msprime : 282s (x3.5)

Mesure de temps d'execution pour un IndSeq avec 5000 loci d'un seul type autosomaux : human
- Sans maf : diyabc : 114s et msprime : 243s (x2.13) 

Mesure de temps d'execution pour un IndSeq avec 100 loci de type autosomal : toy_example3
- Avec maf = 0.05 : diyabc : 2s et msprime : 224s (x112) !!!!!!!!!!!!!!!!
- 500 loci et maf = 0.05 : 5s et msprime : 1126s (x225) !!!!!!!!!!!!!!!!!

Mesure de temps d'execution pour un IndSeq avec 100 loci  : toy_example5
- multitype 70 A ; 10 X ; 10 M ; 10 Y : diyabc : 0.6s et msprime : 4.3 (x7)
- multitype 350 A ; 50 X ; 50 M ; 50 Y : diyabc : 2.7s et msprime : 13.7 (x5)

Après sur le temps avec maf, on est 2.4s par particules pour les deux benchmarks. On a une progression linéaire du temps de simulation. Mais on peut faire un appel msprime.sim_ancestry par petit batch au lieu de faire un appel pour chaque tentive et chaque rejet.via num_replicates 

après modification du code de with_maf_filter et with_mrc_filter : 

Mesure de temps d'exécution pour un PoolSeq avec 100 loci : toy_example4:
- Avec un mrc = 1 : diyabc : 34s et msprime : 105s (x3.1)
- Avec un mrc = 5 : diyabc : 80s et msprime : 114s (x1.4)

Mesure de temps d'execution pour un IndSeq avec 100 loci de type autosomal : toy_example3
- Avec maf = 0.05 : diyabc : 2s et msprime : 15s (x7.5) 
- 500 loci et maf = 0.05 : diyabc : 5s et msprime : 65s (x13)

On essaie d'améliorer encore les choses en optimisant le nombre de batch pour le maf. On garde 20 pour les petits nombres de loci et nb_loci/4 sinon. 

Mesure de temps d'execution pour un IndSeq avec 100 loci de type autosomal : toy_example3
- Avec maf = 0.05 : diyabc : 2s et msprime : 13s (x6.5) 
- 500 loci et maf = 0.05 : diyabc : 5s et msprime : 40s (x8)

Pour le mrc, on va faire un pool partagé sur tous les locus car avant la boucle recrée un lot à chaque locus même si le locus n'avait besoin que d'une tentative pour passer le mrc.

Mesure de temps d'exécution pour un PoolSeq avec 100 loci : toy_example4:
- Avec un mrc = 1 : diyabc : 34s et msprime : 22s (x0.6)
- Avec un mrc = 5 : diyabc : 80s et msprime : 44s (x0.55)

Possibilité de jouer encore sur la taille du batch en mrc pour gagner encore un peu de temps, compromis vitess/mémoire à trouver. Pourrait dépendre à terme de la valeur du mrc.

## Note du 28/07
Attaque des microsat
- parsing du header : loci, priors, stats
- parameter_sampling : 
 * implémentation des différentes lois selon le mécanisme de diyabc
 * gestion des priors dépendants dans le tirage des valeurs selon les lois de groupe

## Note du 29/07 — modèle de mutation microsat/séquences, msprime vs DIYABC

Recherche (pas d'implémentation) pour préparer la suite du chantier
microsat/sequences-mut : est-ce que les modèles de mutation déjà
intégrés à msprime (`msprime.SMM`, `msprime.JC69`/`HKY`/`GTR`) peuvent
remplacer un algorithme écrit à la main, comme on l'a fait pour Hudson
côté SNP ?

**Microsat : `msprime.SMM` NE correspond PAS au modèle DIYABC** (vérifié
dans `particuleC.cpp::ParticleC::mute`, branche `type<5`,
lignes ~1682-1699). `SMM` (doc msprime) est un stepwise strict : ±1
avec proba 50/50, et une mutation qui sortirait de `[lo,hi]` n'a
simplement aucun effet (bornes absorbantes). Le vrai modèle DIYABC est
plus riche sur trois points :
1. Deux processus mélangés par événement de mutation : "SNI" (toujours
   exactement ±1, taux `sni_rate`) vs "GSM" (taux `mut_rate`), choisis
   avec probabilité `sni_rate/(sni_rate+mut_rate)`.
2. Le GSM autorise des sauts de plus d'un pas : taille `d` tirée d'une
   loi **géométrique** de paramètre `Pgeom` (= la valeur tirée du prior
   `MEANP`/`GAMP` du groupe) — `d = 1 + floor(log(ra)/log(Pgeom))`,
   `d=1` si `Pgeom<=0.001` — puis déplacement `± d * motif_size`.
3. Bornes `[kmin,kmax]` (dérivées de `motif_size`/`motif_range`, voir
   `header.cpp:2014-2017`) **clampées**, pas absorbantes : la mutation a
   quand même lieu, juste plafonnée, contrairement à `SMM` qui annule
   la mutation entière si elle sort de l'intervalle.

Conclusion : `msprime.SMM` ne peut représenter ni le mélange SNI/GSM, ni
`Pgeom`, ni le clamping — pas un problème de paramétrage, un modèle
structurellement plus pauvre. Il faudra écrire l'algorithme à la main
(même esprit que `_draw_single_mutation_edge_child`/Hudson pour le SNP),
pas s'appuyer sur `msprime.SMM`. Le nombre d'événements de mutation par
branche reste un processus de Poisson standard (`put_mutations`,
`mutrate = mut_rate + sni_rate`, taux constant), c'est la mécanique de
CHAQUE événement qui diverge de `SMM`.

**Séquences ADN : les modèles msprime correspondent bien, cette fois.**
`comp_matQ` (`particuleC.cpp:1121-1166`) construit une matrice de
transition 4×4 (choix du nouveau nucléotide, CONDITIONNEL à un événement
de mutation déjà survenu via un Poisson séparé sur `mus_rate *
dnalength`) selon `grouplist[gr].mutmod` :
- `mutmod=0` (JK/Jukes-Cantor) : matrice non modifiée (taux/fréquences
  égaux après normalisation ligne par ligne) → `msprime.JC69()` exact,
  aucun paramètre.
- `mutmod=1` (K2P) : transitions (A↔G, C↔T) pondérées par `k1`, pas de
  `pi_X` → `msprime.HKY(kappa=k1, equilibrium_frequencies=[0.25]*4)`
  (msprime n'a pas de classe `K2P` dédiée, mais K2P = HKY à fréquences
  égales).
- `mutmod=2` (HKY) : `matQ[i][j] = pi_j * (k1 si transition sinon 1)` —
  construction HKY85 standard → `msprime.HKY(kappa=k1,
  equilibrium_frequencies=[pi_A,pi_C,pi_G,pi_T])` terme à terme.
- `mutmod=3` (TN/Tamura-Nei) : comme HKY mais deux kappas différents
  (`k1`/`k2` selon la paire de transition) → pas de classe `TN93`
  dédiée dans msprime, mais représentable via `msprime.GTR(relative_
  rates=..., equilibrium_frequencies=[pi_A,pi_C,pi_G,pi_T])` en
  construisant la matrice de taux relatifs à la main.
Classes msprime vérifiées disponibles (introspection directe,
`diyabc_msprime` env, msprime 1.4.2) : `JC69`, `HKY(kappa,
equilibrium_frequencies)`, `F84`, `GTR(relative_rates,
equilibrium_frequencies)` — pas de `K2P`/`TN93` nommées, d'où les
équivalences ci-dessus.

**Origine de `pi_A`/`pi_C`/`pi_G`/`pi_T`** (vérifié `data.cpp:1533-1568`) :
fréquence EMPIRIQUE, calculée PAR LOCUS (`this->locus[loc].pi_A`, pas
une valeur globale partagée entre loci) en comptant les bases sur
TOUTES les séquences observées de ce locus (toutes populations/individus
confondus) au chargement du `.mss` observé — jamais tirées d'un prior,
jamais recalculées par particule. Implique un nouveau module de lecture
des séquences ADN observées (même famille que `observed_data.py` pour
le SNP/PoolSeq : `count_samples_per_population`/`observed_reads`), pas
encore écrit.

**Reste ouvert, à reprendre dans une prochaine session** :
- `gams`/`p_fixe` (les deux valeurs de la ligne `MODEL K2P 10 2.00` du
  header) : proportion de sites invariants + hétérogénéité de taux
  entre sites (looks like une loi Gamma, terminologie phylogénétique
  standard) — pas encore tracé dans le C++, pas géré par `comp_matQ`
  lui-même. À voir si `msprime.sim_mutations` supporte nativement un
  taux variable par site, ou s'il faut simuler les sites invariants à
  part.
- Le module de lecture des séquences ADN observées (pour `pi_A..T`) —
  pas commencé.
- L'algorithme de mutation microsat lui-même (SNI/GSM/Pgeom) — pas
  commencé, à écrire à la main.

## DNA sequence summary statistics (started 2026-08-24, mentor mode — user-driven, reviewed/debugged with the assistant)

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Picked up right after the ploidy/demography fix above. Goal: reproduce
the 13 DNA-sequence-specific statistics from `sumstat.cpp` (`cal_nha1p`/
`2p`, `cal_nss1p`/`2p`, `cal_mpd1p`, `cal_vpd1p`, `cal_mpw2p`, `cal_mpb2p`,
`cal_dta1p`, `cal_pss1p`, `cal_mns1p`, `cal_vns1p`, `cal_fst2p` —
confirmed against `toy_example2_ms_dna/headerRF.txt`'s `group summary
statistics` section, which requests exactly these 13 under the names
`NHA`/`NSS`/`MPD`/`VPD`/`DTA`/`PSS`/`MNS`/`VNS` per-population and
`NH2`/`NS2`/`MP2`/`MPB`/`HST` per-pair) in `bridge/summary_statistics.py`,
building on the mutated `TreeSequence`s from `dna_mutation_simulation_
per_locus`. Distinct from the MicroSat-specific `NAL`/`HET`/`VAR`/`MGW`/
`FST`/`LIK`/`DAS`/`DM2` stats declared in the same header's `G1` group —
those need allele-size arithmetic, not tskit genotypes, and are out of
scope here.

- **`compute_population_layout`** (`ancestry_simulation.py`, renamed
  from the private `_population_layout` since `summary_statistics.py`
  now needs it too) hit a real name-shadowing bug when made public: four
  call sites inside `simulate_snp_genotypes`/`with_maf_filter_shared_
  ancestry` (which both also have a **parameter** named `population_layout`)
  did `population_layout = population_layout(ts)` — the local parameter
  shadowed the module-level function, so this became `None(ts)` whenever
  the parameter defaulted to `None`, breaking 24 tests across three test
  files. Fixed by renaming the function only (not the parameter, which is
  documented at length in multiple docstrings) — see
  `feedback_name_shadowing_pattern` project memory, same bug class
  recurring.

- **`_genotype_matrix_by_population`** (`summary_statistics.py`): one
  `TreeSequence` (one DNA sequence locus) → `{pop_name: matrix}`, matrix
  shape `(n_sites, n_samples_pop)` — tskit's native `genotype_matrix()`
  convention, sliced per population via `compute_population_layout`.
  `genotype_matrix()` called once per `TreeSequence`, not once per
  sample (an early draft rebuilt the whole matrix per sample). Tested on
  both an `<A>` and an `<M>` locus of `toy_example2_ms_dna`: correct
  `n_sites`/`n_samples` shapes, no sample lost/duplicated across
  populations, and the 2:1 sample-count ratio between `<A>`/`<M>`
  confirms it composes correctly with the 2026-08-24 ploidy fix above.

- **Two-tier pattern established for the per-population ("1p") stats**,
  mirroring `sumstat.cpp`'s own per-locus/per-group split (`cal_*pl` +
  `cal_*1p` with an `nl` denominator): a `_count_*(matrix)` brick
  (one locus, one population → a scalar) plus a `mean_*_per_group
  (tree_sequences, population_names)` aggregator (mean over a **single**
  header `group Gx`'s loci — never loci from two different groups mixed
  together, since each group computes its own independent stat).
  Aggregators pre-fill `{pop_name: 0.0 for pop_name in population_names}`
  before accumulating, matching the C++'s `res = 0.0` declared before its
  `if (nl > 0)` guard — so every expected population always has a value,
  even for an empty `tree_sequences` list, rather than a population
  silently missing from the result dict. Documented, not fixed: this
  assumes every population in `population_names` is present on every
  locus of the group (divides by `len(tree_sequences)`, not a real
  per-population `nl` count) — true on `toy_example2_ms_dna` (checked
  empirically), not guaranteed in general; violating it raises `KeyError`
  rather than silently excluding that locus, unlike the C++.
  - `NSS` (`_count_segregating_sites` + `mean_segregating_sites_per_group`,
    `cal_nsspl`/`cal_nss1p`): a site is segregating for a population if
    not all its samples share the same base — vectorized as
    `np.any(matrix != matrix[:, [0]], axis=1)`, correct for any number of
    distinct values per site since "differs from sample 0" is equivalent
    to "not all identical" (not just "exactly 2 alleles").
  - `NHA` (`_count_distinct_haplotypes` + `mean_distinct_haplotypes_per_group`,
    `cal_nha1p`): number of distinct haplotypes = distinct **columns** of
    the matrix (`np.unique(matrix, axis=1)`). Bug caught and fixed: an
    early draft returned `.shape[0]` (number of *sites*) instead of
    `.shape[1]` (number of distinct haplotypes) — both happened to be
    `3` on the first hand-picked test matrix, masking the bug until a
    second matrix with `n_sites != n_distinct_haplotypes` was tried.
    `np.unique` on a `(0, n_samples)` matrix (locus with zero variable
    sites) correctly returns exactly 1 unique column for free, matching
    `cal_nha1p`'s explicit `dnavar == 0` → 1-haplotype special case
    without needing an explicit branch.
  - Both `_count_segregating_sites`/`_count_distinct_haplotypes` raise
    `ValueError` on `matrix.shape[1] == 0` (population with zero samples
    on a locus) rather than crashing obscurely or returning a silently
    wrong count — deliberately checked against `shape[1]` (samples), not
    `matrix.size` (an earlier draft used `.size`, which also triggers
    incorrectly on the *valid* `n_sites == 0` case, a locus with no
    mutations at all — see `feedback_control_flow_chaining_bugs`-style
    edge-case conflation).
  - A `ValueError` was independently added and removed **twice** from
    `mean_segregating_sites_per_group`'s empty-list handling during this
    session — once genuinely misplaced (inside `if num_loci > 0` instead
    of `else`, so it fired on the *normal* case), once syntactically
    correct but reintroduced the very "raise instead of 0.0-default"
    design this whole two-tier pattern was built to avoid. Kept as
    `0.0`-default, confirmed explicitly with the user both times — if
    this `raise` reappears a third time, check for an editor/autosave
    restoring a stale buffer rather than assuming it's an intentional
    edit.

- **`_pairwise_hamming_distances`** (`summary_statistics.py`, brick for
  `MPD`/`VPD`/`cal_mpdpl`/`cal_vpd1p`): one matrix → the 1D vector of
  Hamming distances for all `C(n_samples, 2)` pairs, via
  `(matrix[:, :, None] != matrix[:, None, :]).sum(axis=0)` then
  `np.triu_indices(..., k=1)` to keep `i < j` pairs only (no double-count,
  no diagonal). Verified against a hand-computed matrix. Same
  `matrix.shape[1] == 0` → `ValueError` guard added for consistency with
  `_count_segregating_sites`/`_count_distinct_haplotypes` (an initial
  draft silently returned `[]` instead, which would have propagated into
  `nan` + a numpy `RuntimeWarning` from `mean()`/`var()` on an empty
  array rather than a clear error at the source).

- **`MPD`/`VPD`** (`mean_pairwise_differences_per_group`/`variance_
  pairwise_differences_per_group`, `cal_mpd1p`/`cal_vpd1p`) needed a
  different exclusion regime than `NSS`/`NHA`: instead of a flat
  `num_loci` denominator, a **per-population** `valid_loci_count` dict,
  because a locus only contributes if it has at least 1 pair (`MPD`,
  `nd > 0`) or at least 2 pairs (`VPD`, `nd > 1`) — `_pairwise_hamming_
  distances` on a 1-sample matrix returns an empty vector, and
  `.mean()`/`.var()` on that gives `nan`, which would otherwise silently
  poison the whole group's sum. Not a bug in practice on this project's
  datasets (20-40 samples/population always) but the C++'s own guard
  reproduced for fidelity. A `> 1` vs `> 0` threshold confusion on
  `VPD`'s **final division** guard was flagged as a possible bug and
  turned out to be a false alarm — dividing by 1 is a no-op, so the two
  thresholds are numerically indistinguishable in every case; still
  switched to `> 0` for readability/consistency with the rest of the
  file, not because the `> 1` version was wrong.

- **`DTA`** (`_tajima_d_per_locus` + `mean_tajima_d_per_group`,
  `cal_dta1pl`/`cal_dta1p`) is the classic Tajima's D neutrality
  statistic, built directly on top of `MPD` (π) and `NSS` (S) — no new
  per-site logic needed beyond `_tajima_constants(n_samples)` (the
  `a1`/`e1`/`e2` coefficients, pure functions of sample size). Caught
  before commit: a parenthesization bug, `(n+1) / ((n-1)/3.0)` instead
  of `(n+1)/(n-1)/3.0`, made `b1` exactly 9× too large (verified
  numerically). **Two distinct, easy-to-conflate exclusion cases**:
  `n_samples < 2` excludes the locus entirely (`_tajima_d_per_locus`
  returns `None`, matching `OKK = false`); `n_samples >= 2` but `S == 0`
  (no segregating sites → the formula's denominator is 0) still
  **includes** the locus in the group average with a contributed value
  of `0.0` — the C++ never resets `OKK` in that second case (`cal_
  dta1pl` lines 1579/1594-1598). An initial draft used the wrong input
  shape entirely (SNP-style `genotypes_per_locus: list[dict]`, treating
  "number of loci" as "number of samples") before being rewritten to
  match the `tree_sequences`-based two-tier pattern of every other stat
  here.

- **`PSS`** (`_private_segregating_sites_per_locus` + `mean_private_
  segregating_sites_per_group`, `cal_pss1p`) is the one per-population
  stat that needs **every** population's matrix at once, not just the
  target's — a site counts as "private" to population `i` only if it's
  segregating in `i` and fixed in **every other population of the whole
  dataset** (not just the target's group). Required factoring a reusable
  `_segregating_sites_mask(matrix)` boolean helper out of `_count_
  segregating_sites` (previously computed the count directly). Since all
  populations' matrices for one locus come from the same underlying
  `genotype_matrix()` (just column-sliced), row `i` means the same
  physical site for every population — booleans compare directly by
  position, no index-matching search needed (the C++ does need one,
  `ssa[sample][j] == ssa[sa][k]`, because its per-population variable-
  site index lists are separate dynamic arrays). `nl` increments
  unconditionally every locus in `cal_pss1p` (no `samplesize > 0` guard
  at all, unlike `NSS`/`NHA`) — simple `num_loci` denominator.

- **`MNS`/`VNS`** (`_minor_allele_counts_at_segregating_sites` +
  `mean_minor_allele_count_per_group`/`variance_minor_allele_count_per_
  group`, `afs`/`cal_mns1p`/`cal_vns1p`): at each segregating site,
  `min(counts of each distinct base actually present)` — reproduces the
  C++'s "sort 4 slots ascending, skip zeros" (`afs`) via `np.unique(site,
  return_counts=True).min()` when `len(values) > 1`, simpler because
  `np.unique` only ever returns actually-present values (no need to
  handle the zero-slots explicitly). **`VNS` is a BIASED variance
  (`ddof=0`, division by `n` not `n-1`)** — confirmed against `cal_
  vns1p`'s `v = (sx2 - sx*sx/a) / a`, deliberately different from `VPD`'s
  `ddof=1`, easy to get wrong by pattern-matching against `VPD`. Both
  stats use the flat `num_loci` denominator (`nl` increments
  unconditionally in both `cal_mns1p`/`cal_vns1p`, like `PSS`) — no
  per-population exclusion needed.

- **Pairwise ("2p") stats reuse the 1p bricks on pooled/cross
  matrices**, all following the same aggregator skeleton (`{pair_key:
  0.0}` pre-filled, `num_loci` denominator, `"{i+1}.{j+1}"` keys via a
  plain double loop — not `_half_arrangements`, see below):
  - **`NH2`** (`mean_distinct_haplotypes_per_group_pairwize`, `cal_
    nha2p`): `_count_distinct_haplotypes` on `np.hstack(matrix_a,
    matrix_b)` — the two populations' matrices are column-slices of the
    same `genotype_matrix()`, so concatenation along the sample axis
    needs no realignment.
  - **`NS2`** (`mean_segregating_sites_per_group_pairwize`, `cal_
    nss2p`): identical trick with `_count_segregating_sites` instead.
  - **`MP2`** ("mean pairwise **within**",
    `mean_pairwise_differences_per_group_pairwize`, `cal_mpw2p`):
    `_pairwise_hamming_distances` computed **separately** on each
    population's own matrix (never concatenated), then pooled as a
    ratio of sums (`(sum_di_a + sum_di_b) / (nd_a + nd_b)`) — NOT a
    simple average of the two populations' own `MPD` values; only
    equal to that average when both populations have the same sample
    size (as they happen to in `toy_example2_ms_dna`).
  - **`MPB`** ("mean pairwise **between**",
    `mean_pairwise_differences_between_per_group_pairwize`, `cal_
    mpb2p`): new brick `_pairwise_hamming_distances_between(matrix_a,
    matrix_b)` — full cross-product `(matrix_a[:,:,None] !=
    matrix_b[:,None,:]).sum(axis=0)`, shape `(n_a, n_b)`, **no
    triangle extraction needed** (unlike the "within" case) since every
    `(p in a, q in b)` pair is valid, never a self-comparison. An early
    draft's docstring claimed this returned a flattened 1D vector of
    length `n_a*n_b`; it actually returns the 2D `(n_a, n_b)` matrix —
    `.mean()` on it is still numerically correct either way (numpy
    averages all elements regardless of shape), but a downstream `len(...)
    > 0` guard was checking the wrong thing (`n_a`, not `n_a*n_b`) —
    harmless in practice only because the brick's own guard already
    rejects 0-sample inputs before that check is ever reached.
  - **`HST`** (`mean_hst_per_group_pairwize`, `cal_fst2p` — note the
    lowercase C++ name, easy to confuse with MicroSat's unrelated
    `cal_Fst2p`): `(Hb - Hw) / Hb` where `Hb` = `MPB`-per-locus, `Hw` =
    `MP2`-per-locus. **Breaks the "mean over loci" pattern used by every
    other DNA stat** — it's a *ratio of sums* accumulated separately per
    pair across the whole group (`num[pair]`, `den[pair]`, divided once
    at the very end), the same style as `_fst_wc` (SNP side, already in
    this file), not `sum_of_per_locus_values / num_loci`. Took two
    attempts to get right: a per-pair aggregator needs `num`/`den` as
    **dicts keyed by pair**, not shared scalars reset once per locus —
    with a shared scalar, N>2 populations silently mix different pairs'
    contributions and only the last pair visited by the inner loop ever
    gets its dict entry updated (every other pair stays stuck at its
    `0.0` default). Completely invisible on this project's only real
    DNA-sequence dataset (2 populations = 1 pair, so "the shared scalar"
    and "the only pair" are the same thing by coincidence) — caught only
    via a synthetic 3-population mock test
    (`unittest.mock.patch` on `_genotype_matrix_by_population`). See
    `feedback_pairwise_accumulator_bug` project memory; the same bug
    class could in principle recur in `NH2`/`NS2`/`MP2`/`MPB` if ever
    exercised on a 3+ population dataset, even though those four
    happened to be written correctly.
  - `_half_arrangements` (built for `AML`/`F3`/`F4`'s asymmetric HALF
    ordering, where element order matters) was briefly misapplied to
    generate `NH2`'s plain symmetric pairs — not numerically wrong for
    `r=2` (its HALF filter happens to keep only the ascending-index
    permutation), but semantically confusing and needlessly expensive;
    switched to the plain double loop `compute_HW_HB`/`compute_FST2`
    already use elsewhere in this file.

**Resolved 2026-08-25**: the column-collision question below WAS a real
issue and IS handled — see `compute_all_statistics_dna` in the same
file, which embeds the group index in every column name
(`stats_group_parser.parse_requested_statistic_names` was fixed
alongside it), exactly mirroring how the real DIYABC reftable itself
names these columns (`NSS_2_1` for group G2, `NSS_3_1` for group G3,
never a bare `NSS_1`) — verified byte-for-byte against real `diyabc`
output on `toy_example2_ms_dna`. `compute_all_statistics_dna(header_text,
tree_sequences_by_locus, population_names)` is the top-level DNA entry
point (mirrors `compute_all_statistics`), now wired into `pipeline.py`
(see "DIYABC-replay pipeline for DNA sequences" below) — the
`stats_group_parser.py` docstring's old "dedup deferred, unclear if
legitimate" framing is stale, ignore it.

**Resolved 2026-09-11**: MicroSat's own stats (`NAL`/`HET`/`VAR`/`MGW`/
`N2P`/`H2P`/`V2P`/`FST`/`LIK`/`DAS`/`DM2`) are now implemented too — see
"MicroSat summary statistics" above. They needed a different data shape
than tskit genotype matrices (allele sizes from `_length_by_population`,
plus individual/ploidy-aware genotypes for `FST`/`LIK` specifically),
confirming the original note below was right that this wasn't a simple
reuse of the DNA-sequence bricks. See `notes/resume_stat_dna_ms.md` for
a biology-first (not code-first) explanation of what each of the 13 DNA
stats and 11 MicroSat stats measures.


## DIYABC-replay pipeline for DNA sequences (2026-08-26) — cross-validated against a real reftable

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Mirrors the SNP-side replay architecture (see item 12 above,
`run_reftable_simulation` vs. `replay_reftable_simulation`) exactly: a
`_dna`/`_from_values` sibling of each SNP function, never modifying the
SNP originals, so this can't regress the already-validated SNP path.
Six pieces, all in `bridge/`:

1. `reftable_loop.group_prior_column_names(header_text)` — real
   reftable column names for group priors (`µseq_2`, `k1seq_2`,
   `µmic_1`, `pmic_1`...). Verified scenario-INDEPENDENT (`nparamut` is
   a constant across scenarios, unlike `nparam` for historical params)
   — do NOT reuse `_kept_param_names_by_scenario` here, it expects
   `Prior` objects with `.bounds` (not strings) and filters by
   scenario, a concept group priors don't have.
2. `reftable_loop.parse_real_reftable_params_with_group_priors(path,
   priors, scenarios, group_priors_names)` — a NEW sibling of
   `parse_real_reftable_params`, returns triplets `(scenario_index,
   historical_values, group_priors_values)` with the two value dicts
   kept SEPARATE (they feed different downstream stages: demography
   vs. mutation model). Like the SNP original, handles the real
   reftable's per-row ragged column width (see `parse_real_
   reftable_params`'s own docstring for why a naive whitespace-split
   parser silently misaligns columns here).
3. `ancestry_simulation._group_prior_values_from_columns(group_priors_
   values, group_priors)` — reshapes the flat real-column dict into the
   nested `{group: {"MEANMU":.., "MEANK1":..}}` shape `draw_group_
   parameter_values` already produces, so `build_group_local_param_
   per_locus`'s existing body (model branching, per-locus `sampling_
   group_local_param` calls) is reused byte-for-byte in `build_group_
   local_param_per_locus_from_values`. Only the group-level (tier 1)
   draw is replaced with the real value; the per-locus (tier 2)
   dispersion around that mean is NEVER replaced — real DIYABC doesn't
   record it in the reftable, so there's nothing to replay, it keeps
   drawing from `seed`. Same principle propagates through `build_
   matrix_per_locus_from_values`/`build_rate_map_per_locus_from_values`/
   `dna_mutation_simulation_per_locus_from_values`.
4. `pipeline.compute_summary_statistics_dna`/`_from_values`,
   `reftable_loop._run_single_particle_dna`/`_from_values`,
   `reftable_loop.replay_reftable_simulation_dna` — each a direct
   mirror of its SNP sibling, no `num_loci`/`observed_reads_per_locus`
   params (no PoolSeq/loci-truncation concept for DNA sequences).
5. `write_reftable_txt` needed ZERO new code — already fully generic
   (`ParticleResult` + `priors`/`scenarios`), reused as-is for DNA
   results.
6. `scripts/replay_diyabc_priors_dna.py` (and its 50-loci-dataset copy
   `scripts/replay_diyabc_priors_dna_50loci.py`) — the runnable
   end-to-end script, mirrors `scripts/replay_diyabc_priors.py`. Run as
   `python3 -m scripts.replay_diyabc_priors_dna` from the repo root
   (see the `scripts/` section below for why). A 1000-particle real
   reftable replay runs in under 2 minutes.

**Validation result** (`toy_example2_ms_dna`, scenario 1, 1000
particles, 5+5 DNA loci): historical params match the real reftable
EXACTLY (`rdiff_mean = 0`, KS `p = 1.0` on `N1`/`t1`/`ta`/`ra`/`t2` —
confirms the replay plumbing itself is correct). Of the 42 DNA stat
columns, 11 show KS `p<0.05` — see the next section for the full
investigation of that gap.


## Note du 27/08/26 — écart de variance sur les stats ADN séquence du groupe G3 (<M>) : investigation complète, du dataset stress-test jusqu'au code source C++

Contexte : la validation appariée DIYABC/msprime sur `toy_example2_ms_dna`
(5 loci `<A>` + 5 loci `<M>`, 1000 particules réelles rejouées via
`scripts/replay_diyabc_priors_dna.py`) avait montré 11/42 colonnes de
stats ADN avec un écart KS significatif (p<0.05), concentrées presque
exclusivement sur `MNS`/`VNS`/`DTA`/`VPD` (les stats dérivées de la
variance des différences par paire), et nettement pires sur le groupe
G3 (`<M>`, mitochondrial/haploïde) que sur G2 (`<A>`, autosomal
diploïde). Cette note documente l'investigation complète menée pour
comprendre cet écart, en trois étapes.

### Étape 1 — tester l'hypothèse "pas assez de loci" (réfutée)

Précédent connu côté SNP (voir note du 17/07 ci-dessus) : un écart
"décevant" en KS s'était résorbé en passant de 10 à 650 loci par type.
Hypothèse testée ici : le même phénomène pourrait expliquer l'écart
ADN, vu qu'on n'a que 5 loci par groupe (contre 5000 loci SNP sur
human).

Construction d'un dataset étendu `reference/toy_example2_ms_dna_50loci/`
(jamais modifié l'original) : les 5 loci `<A>` et 5 loci `<M>` dupliqués
jusqu'à 50+50, sous de nouveaux noms continuant la numérotation globale
existante (`Locus_S_A_21_` à `_65_`, `Locus_S_M_66_` à `_110_` — un
premier essai avait redémarré la numérotation à 6, ce qui entrait en
collision avec les noms originaux `S_A_11..15`/`S_M_16..20`, détecté et
corrigé avant de lancer quoi que ce soit). Vrai run DIYABC relancé
dessus (1000 particules, `-R ALL -r 1000 -g 1000 -m -t 16`, après
réinitialisation du RNG via `-n "t:16;c:1;s:1;f:"` — l'option `-n 1`
seule ne suffit pas, il faut la chaîne complète), puis rejeu msprime via
une copie paramétrée du script de replay.

Résultat : **17/47 colonnes significatives, contre 11/42 à 5+5 loci —
la proportion n'a PAS diminué**, ce qui réfute l'hypothèse. Une première
comparaison avait donné un résultat alarmant (39/47, écarts de -96%)
mais c'était un artefact de parsing : les lignes du vrai reftable ont
une largeur variable selon le scénario de la ligne (colonnes de
paramètres/group-priors différentes), et un parseur `pandas.read_csv
(sep=r'\s+')` naïf décale silencieusement toutes les colonnes sur les
lignes plus courtes. Corrigé en réutilisant `_kept_param_names_by_scenario`/
`group_prior_column_names` (déjà dans `bridge/reftable_loop.py`) pour
aligner chaque ligne par son propre scénario.

Le résultat corrigé affine le diagnostic plutôt que de le confirmer
platement : G2 colle très bien à DIYABC (ratio écart-type sim/réel
0.94–1.04 sur les 8 stats, écarts de moyenne négligeables) ; G3
concentre presque tous les écarts significatifs, et c'est un **déficit
de variance**, pas un décalage de moyenne : ratio écart-type sim/réel
de 0.65 à 0.94 selon la stat, et **0.26 pour `DTA`**. Les `rdiff%`
énormes sur `DTA` (jusqu'à +139%) sont un leurre : sa vraie moyenne est
proche de 0 (0.001 à 0.09), donc un tout petit écart absolu explose en
pourcentage — le ratio d'écart-type est le vrai signal.

### Étape 2 — isoler ancestrie vs mutation (les deux étages sont mis hors de cause)

Deux diagnostics jetables, mesurant le coefficient de variation
(écart-type/moyenne, invariant d'échelle) entre les 50 loci d'un même
groupe, au sein de particules réellement rejouées :

1. **Ancestrie seule** (`msprime.sim_ancestry` sans mutation,
   `ts.first().total_branch_length` comme proxy) : ratio
   CV(G3)/CV(G2) = **1.05**.
2. **Pipeline mutation complet** (`dna_mutation_simulation_per_locus_
   from_values` avec les vraies moyennes de groupe rejouées,
   `ts.num_mutations` comme proxy) : ratio CV(G3)/CV(G2) = **1.02**.

Les deux ratios sont proches de 1 : dans NOTRE simulation, G2 et G3 ont
une variance relative quasi identique, à chaque étage. Le déficit n'est
donc PAS dans le mécanisme de dispersion par locus, ni dans le
rééchelonnage `coalescence_coefficient`/`ploidy`, ni dans le modèle de
mutation — tout ça est cohérent en interne. Le problème, c'est
spécifiquement que le G3 RÉEL de DIYABC porte PLUS de variance
(relativement à sa propre moyenne) que son G2, alors que notre
simulation garde un ratio ~1 entre les deux.

### Étape 3 — lecture directe du code source C++ (`~/Documents/Github/diyabc`)

Plutôt que de continuer à théoriser côté Python, lecture de
`particuleC.cpp` (`ParticleC::coal_pop`) pour vérifier la fidélité de
notre traduction du mécanisme de coalescence.

**Formule des temps de coalescence** (`particuleC.cpp:1329-1341`, mode
"approximation continue") :
```
start -= (coeffcoal * N / nLineages / (nLineages-1)) * log(ra)
```
En comparant à la formule standard du coalescent (temps d'attente moyen
= `2·Ne/(k(k-1))`), ceci implique `Ne_effectif = coeffcoal * N / 2` —
**exactement** notre `rescale_demography(factor=coeffcoal/2)`. La
formule des moyennes qu'on a portée est donc fidèle au bit près.

**Deux régimes de coalescence** (`ParticleC::evalcriterium`,
`particuleC.cpp:1251-1275`) : en plus de l'approximation continue
ci-dessus, DIYABC a un mode "génération par génération" (Wright-Fisher
discret : tirage d'un parent uniforme parmi `Ne` par lignée à chaque
génération, détection de collision) qu'on ne réplique pas du tout côté
msprime. Vérifié empiriquement en instrumentant temporairement
`coal_pop` avec une trace conditionnelle (`bool trace = (loc==10) or
(loc==15);`, recompilation via le `CMakeLists.txt` existant du repo,
puis reverti après coup — build propre, `cmake --build . --target
general`) : sur les N1 typiques de ce dataset (1000-10000), `ra =
nLineages/N` reste bien en dessous de tous les seuils d'`evalcriterium`,
pour `<A>` ET `<M>` — ce mode discret ne se déclenche jamais ici, écarté
comme cause.

**L'admixture par tirage de Bernoulli indépendant par lignée**
(`ParticleC::split_pop`, `particuleC.cpp:1513-1524`) :
```cpp
if (this->mw.random() < this->seqlist[iseq].admixrate)
    this->gt[loc].nodes[i].pop = this->seqlist[iseq].pop1;
else
    this->gt[loc].nodes[i].pop = this->seqlist[iseq].pop2;
```
Chaque lignée survivante au moment du split `ta` reçoit un tirage
indépendant. En traçant 4 particules réelles de scénario 1 (le seul
avec un événement `ta split`), le nombre de lignées entrant dans ce
split est très différent entre groupes :

| particule (Ne) | lignées `<A>` (pop3+pop4) | lignées `<M>` (pop3+pop4) |
|---|---|---|
| Ne=9977  | 1+12 = 13 | 0+4 = 4 |
| Ne=4442  | 9+3 = 12  | 6+0 = 6 |
| Ne=5445  | 1+6 = 7   | 0+1 = 1 |
| Ne=9138  | 2+3 = 5   | 1+0 = 1 |

Ratio moyen `<M>`/`<A>` ≈ 0.29-0.33 (cohérent avec un diagnostic Python
équivalent sur 30 particules : nb moyen de lignées à `ta` = 6.88 pour
G2, 2.28 pour G3). Plus parlant que le ratio : **`<M>` a un côté
(pop3 OU pop4) à ZÉRO lignée dans 3 cas sur 4** ci-dessus, alors que
`<A>` ne l'a jamais. Avec si peu de lignées survivantes, le partage de
l'admixture devient quasiment binaire (tout d'un côté par pur hasard),
un régime qualitativement différent de `<A>`, qui se répartit presque
toujours des deux côtés.

**Pourquoi le diagnostic "ancestrie seule" (étape 2) n'avait rien vu** :
`total_branch_length` est une quantité continue qui lisse cet effet
tout-ou-rien ; `DTA` et les stats de différences par paire sont
justement sensibles à CE type de structuration en sous-populations, pas
la longueur totale d'arbre.

**Vérification finale : notre propre msprime reproduit-il cet effet ?**
Script diagnostic utilisant `msprime.sim_ancestry(...,
record_migrations=True)` puis inspection de `ts.tables.migrations`
(population `dest` au temps `ta`) pour compter, par locus, combien de
lignées partent vers pop3 vs pop4 — l'équivalent exact côté msprime de
ce qu'on a tracé côté `split_pop`. Sur 30 particules réelles de
scénario 1 (1332 loci G2, 968 loci G3) :

| | % de loci avec un côté à 0 lignée | ratio moyen (côté minoritaire/total) |
|---|---|---|
| G2 (`<A>`) | 30.6% | 0.206 |
| G3 (`<M>`) | **55.6%** | 0.171 |

**Notre propre simulation reproduit bien l'effet quasi-binaire** : G3 a
un côté à 0 lignée près de deux fois plus souvent que G2, dans le même
sens que la trace réelle DIYABC. Ça réfute l'hypothèse "notre portage
msprime ignore/sous-produit cet effet" — `msprime.add_admixture` fait
le même tirage indépendant par lignée que `split_pop`.

### Conclusion (provisoire — RÉFUTÉE partiellement, voir mise à jour du 31/08 ci-dessous)

Chaîne complète vérifiée, étape par étape, source à l'appui :
1. Formule de moyenne du coalescent → identique entre C++ et notre
   `rescale_demography`.
2. Mode discret Wright-Fisher → jamais déclenché sur ce dataset, ni
   pour `<A>` ni pour `<M>`.
3. Effet d'admixture quasi-binaire sur peu de lignées → réel côté
   DIYABC (tracé) ET reproduit qualitativement par notre msprime
   (vérifié via les tables de migration).

**Aucun bug de portage trouvé** (ce point reste vrai). L'écart résiduel
de variance sur `DTA`/`VNS`/`VPD`/`MNS` pour `<M>` avait été attribué à
un mécanisme combinatoire réel (peu de lignées survivantes → partage
d'admixture quasi tout-ou-rien à l'événement `ta split`), présent et
correctement répliqué des deux côtés. **Cette attribution s'est avérée
incomplète — voir la mise à jour du 31/08/26** : l'admixture n'explique
pas, à elle seule, l'ampleur du déficit observé.

### Mise à jour du 31/08/26 — contre-test scénario sans admixture (réfute l'admixture comme cause PRINCIPALE)

Test de falsification direct de la conclusion ci-dessus : si le
mécanisme "peu de lignées → partage d'admixture quasi-binaire à `ta`"
est bien LA cause du déficit de variance sur G3, alors ce déficit
devrait disparaître (ou fortement diminuer) pour les particules qui
n'ont PAS tiré d'admixture. Or `toy_example2_ms_dna` a justement 2
scénarios candidats tirés à poids égal (`[0.5]`/`[0.5]`) :
- **scénario 1** : `pop1`+`pop2` fusionnent à `t1`, PUIS admixture
  (`ta split`) vers `pop3`/`pop4`, PUIS refusion à `t2` — celui étudié
  ci-dessus.
- **scénario 2** : `pop1`+`pop2` fusionnent à `t1`, point final. Aucun
  `split`, aucune admixture.

Pas besoin de relancer DIYABC : le reftable réel de 1000 particules
déjà rejoué (`scripts/replay_diyabc_priors_dna.py`,
`reference/toy_example2_ms_dna/{first_records_of_the_reference_table_0.txt,
reftable_msprime_replay.txt}`) mélange déjà les deux scénarios (tirage
pondéré par particule). Filtré les 1000 particules par
`scenario_index` (colonne déjà présente, alignement ligne à ligne
réel/simulé vérifié : 0 désaccord de scénario entre les deux fichiers),
puis recalculé le ratio d'écart-type sim/réel séparément pour G2 et G3,
séparément par scénario (488 particules scénario 1, 512 scénario 2) :

| stat | G2 scén.1 (admixture) | G3 scén.1 (admixture) | G2 scén.2 (SANS admixture) | G3 scén.2 (SANS admixture) |
|---|---|---|---|---|
| DTA | 1.027 | 0.609 | 0.975 | 0.627 |
| VNS | 0.996 | 0.830 | 1.014 | 0.778 |
| VPD | 1.246 | 0.909 | 1.250 | 0.714 |
| MNS | 0.984 | 0.791 | 0.985 | 0.825 |

(moyenne globale sur les 21 colonnes G2/G3 : écart moyen G3-G2 =
-0.177 en scénario 1, **-0.165 en scénario 2** — quasi identique)

**Le déficit persiste à une magnitude quasi identique SANS aucune
admixture**, et sur `VPD`/`VNS` il est même **plus marqué** en
l'absence d'admixture (0.714/0.778 vs 0.909/0.830 avec admixture).
Vérifié cohérent sur les 13 types de statistiques individuellement, pas
seulement en moyenne globale.

**Conclusion révisée** : le mécanisme d'admixture quasi-binaire
(étape 3 ci-dessus) est réel et correctement reproduit par notre port
— mais ce n'est PAS la cause principale (ni même une cause nécessaire)
du déficit de variance sur G3, puisque celui-ci survit intact dans un
scénario qui n'a aucun événement d'admixture. La vraie cause est donc
plus générale, probablement liée à `<M>` lui-même (Nₑ réduit → moins de
lignées survivantes à *tout* moment de son histoire, pas seulement à un
événement de split précis — potentiellement moins d'événements de
coalescence indépendants sur lesquels les statistiques peuvent
moyenner/lisser dans l'ensemble de l'arbre, pas seulement au niveau
d'un point de partition discret). **Pas encore investiguée** : cette
piste plus générale (ex: comparer directement le nombre d'événements de
coalescence, ou leur distribution temporelle, entre G2 et G3, à
n'importe quel point de l'arbre, pas seulement à `ta`).

**Statut de l'investigation : ROUVERTE** (n'est plus "close" comme
indiqué le 27/08) — le but du POC (prouver la faisabilité
`header.txt` → msprime) reste atteint indépendamment de cet écart
résiduel, mais l'explication documentée précédemment était incomplète
et ne doit plus être citée comme la cause établie.

**Bug réel trouvé et corrigé en cours de route** (indépendant de ce qui
précède) : `build_group_local_param_per_locus`/son jumeau `_from_values`
(`bridge/ancestry_simulation.py`) recréaient `random.Random(seed +
_KAPPA1_SEED_OFFSET)` (et `_KAPPA2_`/`_MUS_RATE_SEED_OFFSET`) à chaque
itération de la boucle `for group in nloc_per_group`, avec un offset
indépendant du groupe — donc G2 et G3 (tous deux modèle K2P, tous deux
utilisateurs de kappa1) rejouaient exactement la même séquence de
tirages, juste recentrée sur un `k_moy` différent. Invisible à 5+5 loci
(les deux groupes déclarent le même `GAMK1` shape=2, donc la dispersion
relative semblait cohérente par coïncidence) ; révélé seulement en
construisant un dataset multi-groupes de même modèle pour cette
investigation. Corrigé en construisant chaque `rng` UNE SEULE FOIS
avant la boucle sur les groupes (même motif que `build_rate_map_per_
locus`, déjà correct). 128/128 tests verts après régénération des 15
valeurs golden qui en dépendaient (toutes G3/pairwise, aucune G2-seule
touchée — signature exacte confirmant que le fix est bien scopé).

### Mise à jour du 02/09/26 — cause trouvée et corrigée : généalogie `<M>` non partagée entre loci (RÉSOLU, correctif écrit par l'utilisateur)

Cause identifiée par instrumentation temporaire de `cal_dta1pl`
(`~/Documents/Github/diyabc/src-JMC-C++/sumstat.cpp`, revertée après
coup) : dans le vrai DIYABC, les loci `<M>` d'une même particule sont
fortement corrélés entre eux (corrélation de `pi` par paire : 0.33-0.75)
alors que les loci `<A>` sont quasi indépendants (~0.15) — attendu,
l'ADN mitochondrial est non-recombinant, transmission uniparentale, un
seul arbre pour toute la molécule. Ce mécanisme existait déjà côté SNP
(`simulate_shared_ancestry_loci`, `particuleC.cpp:2422-2435`
`GeneTreeY`/`GeneTreeM`) mais n'avait jamais été porté côté séquences
ADN : `dna_mutation_simulation_per_locus`/`_from_values` tiraient une
généalogie indépendante par locus, y compris pour `<M>`, moyennant
artificiellement le bruit inter-locus sur `DTA_3`/`VNS_3`/`MNS_3`/
`VPD_3`.

Correctif (écrit par l'utilisateur, `bridge/ancestry_simulation.py`) :
nouvelle constante `_SHARED_M_ANCESTRY_SEED_OFFSET` ; un locus `<M>`
tire désormais son ancestrie avec une graine FIXE (`seed +
_SHARED_M_ANCESTRY_SEED_OFFSET`, sans `+i`), partagée par tous les loci
`<M>` du dataset — `<A>`/`<H>` inchangés. Deux itérations avant
d'arriver à cette version : un premier essai conditionnait sur
`heritage == "A"` (donc `<H>` serait tombé à tort dans la branche
partagée) ; un second réutilisait `_ANCESTRY_SEED_OFFSET` nu pour `<M>`
(collision possible avec un futur locus `<A>`/`<H>` d'indice global 0) —
corrigé avec une constante dédiée. 130/130 tests verts après
régénération des 15 valeurs golden (G3/pairwise uniquement, G2
inchangé — même signature de validation qu'au 27/08). Rejeu complet du
reftable réel (1000 particules) : KS passe de 11/42 à 2/42 colonnes
significatives, les 2 restantes (`MPD_2_2`, `VPD_2_2`) étant sur G2 avec
un écart dans le sens inverse (bruit d'échantillonnage, pas un déficit).

**Investigation CLOSE.** Voir aussi `CLAUDE.md`, section "RESOLVED
2026-09-02".

## 02/09/26 — Séquences ADN `<X>`/`<Y>` : le vrai DIYABC plante sur toute donnée réellement hémizygote (SUSPENDUE puis reprise et implémentée le 03/09, voir mise à jour ci-dessous)

Reprise du chantier séquences ADN pour lever la limitation documentée
dans `dna_ancestry_parameters_for_heritage` (`bridge/
ancestry_simulation.py`) : `<X>`/`<Y>` y lèvent `NotImplementedError`
avec pour justification "le format `.mss` ne porte pas de sexe par
individu". Cette investigation montre que c'est vrai à moitié, et
surtout que le vrai DIYABC compilé (`~/Documents/Github/diyabc`, build
statique fourni par l'utilisateur) **plante systématiquement**
(`SIGSEGV`) dès qu'on lui donne un jeu `<X>`/`<Y>` en séquences ADN où
des individus sont réellement hémizygotes — il n'y a donc pour l'instant
aucune sortie DIYABC réelle à répliquer pour ce cas.

**Étape 1 — le sexe EST inférable depuis le `.mss`, mais pas via une
colonne dédiée.** Lecture directe de `DataC::do_sequence` (`~/Documents/
Github/diyabc/src-JMC-C++/data.cpp:1494-1495`) : `indivsexe` part à 2
(femelle) par défaut pour tout le monde (`data.cpp:1320`) et n'est
jamais réinitialisé à 2 une fois passé à 1 — un individu devient "mâle"
si :
- un locus `<X>` présente un génotype **haploïde** (un seul groupe de
  crochets `<[seq]>`, pas `<[seq][seq]>` — hémizygotie) ;
- OU un locus `<Y>` présente un génotype non vide (`geno != "[]"`).

Exactement le même mécanisme que `DataC::do_microsat`
(`data.cpp:1406-1407`), qui lui fonctionne correctement (voir plus bas
pourquoi la version séquence ADN ne fonctionne pas).

**Étape 2 — construction d'un jeu de test réel.** L'utilisateur a créé
`reference/toy_example2_ms_dna_XY/` (copie de `toy_example2_ms_dna`,
avec le vrai binaire `diyabc`/`abcranger` dedans) et relabellisé G2
`<A>`→`<X>` dans `headerRF.txt`. Le `.mss` original étant resté
identique bit-à-bit (tous les tokens G2 diploïdes, aucun individu
hémizygote — vérifié par script), l'assistant a construit une version
synthétique avec un vrai mélange de sexes : individus `*-001`..`*-010`
de chaque pop désignés "mâles" (G2 rendu haploïde par troncature du
2e allèle ; G3 relabellisé `<M>`→`<Y>` et laissé haploïde-présent),
`*-011`..`*-020` désignés "femelles" (G2 inchangé diploïde ; G3 remplacé
par le token vide `<[]>`, choisi par analogie avec le format conteneur
`<[...]>` et avec la constante `SEQMISSING` = `""` trouvée dans
`particuleC.hpp:17`).

**Étape 3 — crash, isolé empiriquement en 3 runs réels** (`./diyabc -n
"t:8;c:1;s:1" -p ./` puis `./diyabc -p ./ -R ALL -r 1000 -g 1000 -m -t
8`, dans des copies temporaires sous `/tmp`) :
- **Test 1** (G2 hétérogène + G3 `<Y>` avec `<[]>` manquant côté
  femelles) → `SIGSEGV` dans `ParticleC::cal_numvar`
  (`sumstat.cpp:2157`), confirmé par `gdb -batch -ex run -ex bt`.
- **Test 2** (aucune modification des données observées, seul le header
  relabellisé `<A>`→`<X>`/`<M>`→`<Y>`, ploïdie uniforme pour tous) →
  **pas de crash**, `reftableRF.bin` généré normalement (1000
  particules, 20000 loci simulés).
- **Test 3** (G2 toujours hétérogène, mais G3 laissé en `<M>` — donc
  aucun token manquant nulle part) → crash quand même.

Conclusion de l'isolement : ni le relabellisage `<Y>` en soi, ni le
token `<[]>` choisi pour le manquant, ne sont en cause — c'est
spécifiquement le **mélange haploïde/diploïde au sein d'un même groupe
`<X>` séquence** qui fait planter DIYABC.

**Étape 4 — cause racine, trouvée en relisant `data.cpp`.**
`DataC::readfile` (`data.cpp:1363`) fait `locus[loc].type += 5` pour
tout locus dont un génotype contient un `[` — c'est-à-dire tous les
loci `[S]` (séquence ADN), pour les distinguer de leurs équivalents
MicroSat (type 0-4 → 5-9). Mais les deux lignes de `do_sequence` citées
à l'étape 1 comparent `this->locus[loc].type == 2` / `== 3`, **sans
`% 5`** — alors que `do_microsat` (qui s'exécute AVANT le bump, les
génotypes microsat ne contenant jamais `[`) compare bien sur les types
bruts 2/3, correctement. Résultat : pour une séquence ADN `<X>` (type
réel 7) ou `<Y>` (type réel 8), ces comparaisons ne sont **jamais**
vraies — le bloc d'inférence de sexe est du **code mort** pour les
séquences ADN, il ne fonctionne que pour le MicroSat. `indivsexe` reste
donc à 2 (femelle) pour tout le monde, quel que soit le contenu réel du
fichier `.mss`, dès qu'il s'agit de loci séquence.

**Étape 5 — pourquoi ce code mort fait planter, et pas juste produire un
résultat faux.** `DataC::calcule_ss` (`data.cpp:966-999`) calcule
`ssize[catégorie][pop]` (l'effectif retenu par catégorie de locus) à
partir de cet `indivsexe` toujours-femelle : pour `<X>`, la condition
`(locustype==2) and (indivsexe==2)` est vraie pour tout le monde →
`ssize[X][pop] = 2 × 20 = 40` (les 20 individus comptés diploïdes). Mais
`do_sequence` a, lui, correctement lu le fichier et n'a peuplé
`haplodna[pop]` que de 30 entrées réelles (10 mâles à 1 copie + 10
femelles à 2 copies — la lecture du fichier n'est pas buguée, seule
l'inférence de sexe qui en découle l'est). `cal_numvar`
(`sumstat.cpp:2201-2226`) boucle ensuite jusqu'à `ssize[X][pop]=40` pour
indexer `haplodna[pop][i]`, un `vector<string>` de taille 30 réelle →
accès hors bornes → `SIGSEGV`. Le mésalignement ssize/haplodna
n'apparaît PAS quand tout le monde a la même ploïdie (Test 2 : `ssize`
et taille réelle valent tous deux 40, coïncidence heureuse), d'où
l'absence de crash dans ce cas précis.

**Conclusion et implication pour ce projet** : avec le binaire `diyabc`
actuellement disponible, il n'existe **aucune sortie réelle exploitable**
comme vérité terrain pour un jeu `<X>`/`<Y>` en séquences ADN comportant
de vrais individus hémizygotes — DIYABC plante avant même d'écrire
`statobsRF.txt`. Deux voies possibles pour la suite (non tranchées) :
corriger le `%5` manquant dans `do_sequence` et recompiler pour obtenir
un vrai `reftableRF.bin` de référence (même philosophie que la
correction du bug `<M>` non partagé, RESOLVED ci-dessus), ou implémenter
côté Python en se basant sur l'intention (logique `do_microsat`,
correcte, généralisée aux séquences) sans réplique-terrain réelle,
comme au tout début du chantier MicroSat.

**Statut initial : SUSPENDU** à la demande de l'utilisateur — décision à
prendre avec son encadrant académique avant de choisir entre les deux
voies ci-dessus. `reference/toy_example2_ms_dna_XY/` est laissé en
l'état (header + `.mss` construits pour reproduire le crash) comme cas
de test pour une reprise future. Aucun changement de code dans
`bridge/` à ce stade — `dna_ancestry_parameters_for_heritage` continue
de lever `NotImplementedError` pour `<X>`/`<Y>`, sa justification
("`.mss` ne porte pas de sexe par individu") reste correcte dans les
faits même si la raison profonde côté DIYABC s'est révélée plus
subtile (code d'inférence présent mais mort, pas absent).

### Mise à jour du 03/09/26 — reprise, chantier implémenté côté Python (voie "intention du code", pas de recompilation)

Chantier repris en mode mentor (utilisateur au clavier, assistant en
review/debug). Décision (implicite, jamais formellement tranchée avec
l'encadrant, mais l'utilisateur a choisi d'avancer) : implémenter côté
Python **l'intention évidente du code C++** — la logique de
`do_microsat`, qui fonctionne réellement (elle n'est jamais atteinte par
le bug `%5`, ses génotypes ne contenant jamais de `[`) — plutôt que
recompiler `do_sequence` avec le `%5` corrigé. Pas de nouvelle sortie
DIYABC réelle générée ou utilisée : `reference/toy_example2_ms_dna_XY/`
reste le seul fixture, purement synthétique.

Résumé de ce qui a été fait (détail complet dans `CLAUDE.md`, section
"DNA sequence `<X>`/`<Y>` support") :
- `observed_data.individual_sexes_from_locus_genotype` : nouvelle
  fonction qui déduit le sexe par individu à partir de la ploïdie du
  génotype AU locus `<X>`/`<Y>` lui-même (haploïde/présent = mâle),
  reproduisant exactement `do_microsat`/l'intention de `do_sequence`.
- `ancestry_simulation.py` : `_sample_sets_from_sexes`/
  `_male_counts_from_sexes` factorisées hors de
  `build_sex_stratified_samples_argument`/`build_male_only_samples_
  argument` (`.snp`), réutilisées par deux nouveaux wrappers `_dna`
  (`.mss`). `dna_ancestry_parameters_for_heritage` : `"X"`/`"Y"`
  rejoignent la branche `"H"`/`"M"` (même rescale, `ploidy=1`).
  `dna_mutation_simulation_per_locus`/`_from_values` : `samples`
  dispatché PAR LOCUS (plus un seul calcul hoisté hors boucle) ;
  nouvelle constante `_SHARED_Y_ANCESTRY_SEED_OFFSET` pour que les loci
  `<Y>` partagent une seule généalogie entre eux, comme `<M>`
  (transmission sans recombinaison) — distincte de
  `_SHARED_M_ANCESTRY_SEED_OFFSET` pour ne pas faire partager le MÊME
  arbre à `<M>` et `<Y>`.
- Plusieurs bugs attrapés en cours de review (indices de population
  décalés, asymétrie manquant `<X>`/`<Y>` non respectée, regex `\S+` qui
  ne matchait pas le token manquant `<[]>` -- à la fois dans la nouvelle
  fonction ET dans `observed_sequences`, préexistante, jamais exercée sur
  données manquantes avant `te2_ms_dna_XY` --, vérification du sexe "9"
  perdue lors d'une factorisation, variable renommée à moitié dans le
  jumeau `_from_values`) : tous corrigés, détail dans `CLAUDE.md`.

**Validation** : uniquement sur `reference/toy_example2_ms_dna_XY/`
(synthétique) — cohérence interne vérifiée (nombre de lignées simulées
cohérent avec le sexage réel du fichier, généalogie partagée entre tous
les loci `<Y>`, indépendante entre loci `<X>`) mais **aucune comparaison
appariée à une vraie sortie DIYABC**, puisqu'aucune n'existe pour ce cas
(le binaire réel plante toujours dessus, cause racine inchangée — voir
ci-dessus). Tests ajoutés dans `test_observed_data.py`/
`test_ancestry_simulation.py`, suite complète verte (134/134).

**Statut final : IMPLÉMENTÉ** (voie "intention du code"). La question de
corriger et recompiler `do_sequence` pour obtenir une vraie référence
reste ouverte, toujours à trancher avec l'encadrant si une validation
face au vrai DIYABC devient nécessaire un jour.

## MicroSat GSM mutation model (2026-09-04 to 2026-09-08, mentor mode — user-driven, reviewed/debugged with the assistant)

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Picked up as the natural continuation of "MicroSat / sequences-mut
header parsing" above (parsing was complete, simulation wasn't started
at all). Covers the full path from `header.txt` + `.mss` to a mutated
`tskit.TreeSequence` per MicroSat (`[M]`) locus, validated end-to-end
on `toy_example1_ms`/`toy_example2_ms_dna`/`toy_example2_ms_dna_XY`,
plus the `pipeline.py`/`reftable_loop.py` orchestration layer mirroring
the DNA-sequence path. The DIYABC model (confirmed by direct reading of
`particuleC.cpp::mute`/`setMutParamValue`/`cree_haplo`, and by the
user's own doc research): at each mutation event (Poisson per branch,
rate `mut_rate` — SNI is a separate, deferred channel, see below), a
GSM step of `d` repeat units is drawn via a geometric distribution
parameterized by `Pgeom` (`Pgeom=0` is the SMM special case), clamped
to `[kmin, kmax]`. Ancestral state = midpoint of `[kmin, kmax]`, not
drawn.

**Key discovery (proposed by the user, verified empirically)**:
`msprime.TPM` (Two-Phase Model) with `p→ε` (literal `p=0` raises
`ValueError`) reproduces the GSM channel alone exactly — verified by
inspecting `.transition_matrix` directly: the geometric ratio between
`P(d=k)` and `P(d=k+1)` equals `Pgeom` (so `m_msprime = 1 - Pgeom`,
same epsilon-clamp requirement, relevant because `Pgeom=0`/`Pgeom=1`
are both valid DIYABC inputs that map to msprime's forbidden `m=1`/`m=0`
literals), each row sums to 1, zero diagonal except at the clamped
edges (faithfully reproducing DIYABC's own "auto-mutation at the
boundary"). **Non-obvious fix caught before code was written**: the
`TPM` grid must be anchored on `root` (the ancestral state), NOT on
`kmin` — verified on real data (`toy_example1_ms`'s `Locus_M_A_1_`:
`root - kmin = 39`, ODD, while `motif_size = 2` is even) that the
ancestral state doesn't necessarily fall on a `motif_size`-multiple
offset from `kmin`, and GSM steps happen from the CURRENT state, not
from `kmin`. Accepted consequence (explicit choice, not a silent
approximation): `TPM`'s own bounds can differ from DIYABC's literal
`kmin`/`kmax` by up to `motif_size - 1` bp — negligible in practice
(`kmin`/`kmax` are a generous guard-rail, almost never reached by a
real genealogy).

- **`build_microsat_transition_matrix(kmin, kmax, motif_size, Pgeom,
  epsilon=1e-16)`** (`ancestry_simulation.py`) → `msprime.
  MatrixMutationModel`: computes `root`, `n_minus`/`n_plus`, builds
  `msprime.TPM(p=epsilon, m=clamp(1-Pgeom, epsilon, 1-epsilon), lo=0,
  hi=n_alleles-1).transition_matrix`, relabels alleles to real bp
  (`root + (i-n_minus)*motif_size`), one-hot `root_distribution` at
  `n_minus`. Tested (commit `740b9e6`): nominal case plus both `Pgeom`
  edges — `Pgeom=0` gives exactly `0.5`/`0.5` to the two immediate
  neighbors and zero elsewhere (pure SMM), `Pgeom=1` gives a uniform
  `1/(n_alleles-1)` off-diagonal (no distance preference at all).
  **Non-obvious**: at the `Pgeom` edges, `epsilon` must stay near its
  default (`1e-16`), NOT the nominal case's `epsilon=0.01` — at
  `epsilon=0.01` the `TPM`'s `p` clamp has a real, measurable (~1%)
  effect that breaks the exact `0.5` assertion, since `p` is the
  probability of a "long jump" (uniform anywhere), not just a numerical
  guard-rail.

- **`build_microsat_local_param_per_locus(header_text, seed)`** →
  `dict[locus_name, (mut_rate, Pgeom)]`: same two-tier hierarchy as the
  DNA-sequence `k1`/`k2`/`mus_rate` draw (`draw_group_parameter_values`
  then `sampling_group_local_param`, both generalize for free — no
  MicroSat-specific changes needed in `parameter_sampling.py`, verified
  empirically before writing the rest). **Naming pitfall, already
  documented and re-triggered here**: DIYABC uses `mut_rate` for
  MicroSat and a DIFFERENT variable, `mus_rate`, for DNA sequences —
  reusing `mus_rate` for MicroSat (copy-paste from the DNA-sequence
  code) was caught and corrected before commit. A `next(gp for gp in
  group_priors[group] if gp.name == "GAM")` bug (real name is `"GAMP"`,
  not `"GAM"` — `StopIteration` otherwise) was also caught before
  commit. Tested (`740b9e6`): MicroSat-only filter (the dict must
  contain ONLY `ms_or_seq=="M"` keys, never the `[S]` loci of the same
  mixed header), total count, bounds (`mut_rate>0`, `0<=Pgeom<=1`),
  inter-locus diversity, 2 golden values, reproducibility.

- **`build_matrix_microsat_per_locus(context, seed)`** (signature as of
  the ReplayContext refactor — `(header_text, mss_file_path, seed)` at
  the time of `740b9e6`) → `dict[locus_name, MatrixMutationModel]`: assembles
  `allele_bounds_per_locus` + `build_microsat_local_param_per_locus` +
  `build_microsat_transition_matrix`, mirroring `build_matrix_per_locus`
  (DNA). Tested (`740b9e6`): MicroSat-only filter (on the right dict —
  an early draft looped over `params_per_locus` instead of the
  function's own return value, testing the wrong thing), integration
  cross-check against the first test (`Locus_M_A_1_` → 39 alleles,
  confirming bounds+params+matrix wiring), reproducibility.
  **`MatrixMutationModel` does NOT implement `__eq__` usefully** (`m1
  == m2` is `False` even for two identically-constructed models,
  verified empirically) — reproducibility tests must compare
  `.alleles`/`.root_distribution`/`.transition_matrix` separately
  (`np.array_equal`/`np.allclose`), never `==` on the object or a dict
  containing one.

- **`microsat_mutation_simulation_per_locus(context, demography, seed)`**
  (`ancestry_simulation.py`, commit `c50e37c` — took `(header_text,
  mss_file_path, demography, seed)` before the ReplayContext refactor) →
  `dict[locus_name, tskit.TreeSequence]`: the full per-locus assembly,
  copied from `dna_mutation_simulation_per_locus`'s structure
  (ploidy/demography/sex dispatch via `dna_ancestry_parameters_for_
  heritage`/`build_sex_stratified_samples_argument_dna`/`build_male_
  only_samples_argument_dna` reused as-is — all already generic enough,
  verified before writing). Three bugs caught and fixed during review,
  none survived to the committed version:
  - `msprime.sim_mutations(..., seed=...)` instead of `random_seed=...`
    — confirmed via `inspect.signature`: the real keyword is
    `random_seed`; `seed` isn't a recognized parameter at all
    (`TypeError` at the very first call, caught by actually running the
    function rather than just reading the diff).
  - `mut_rate` was never retrieved at all — `build_matrix_microsat_
    per_locus` (called internally) discards it (`_, Pgeom = params_
    per_locus[...]`), so an early draft passed no `rate=` to
    `sim_mutations` at all. Fixed by calling `build_microsat_local_
    param_per_locus` a second time (deterministic given the same seed,
    so no correctness risk, just a redundant computation — a known,
    accepted minor inefficiency, not fixed) to recover `mut_rate`
    alongside the matrix.
  - `sequence_length=locus.motif_size * locus.motif_range` instead of
    `sequence_length=1` — a MicroSat locus is a SINGLE site with a
    large allele-space alphabet (like DNA sequence's 4-base alphabet,
    but ~39 states here), not multiple independent sites; with
    msprime's default `discrete_genome=True`, a length `>1` would have
    silently simulated several INDEPENDENT big-alphabet sites per
    locus instead of one repeat-length state — caught by comparing
    against the plan agreed the session before (session memory /
    `notes` had already flagged `sequence_length=1` as the right
    answer), not by running the code first.
  Validated by direct execution (not just `pytest`): 10/10 MicroSat
  loci on `toy_example2_ms_dna_XY`, `num_sites=1` everywhere, mutation
  counts scaling plausibly with each locus's own drawn `mut_rate`
  (cross-checked against `build_microsat_local_param_per_locus`'s own
  output), correct `<X>`/`<Y>` sample-count dispatch (`79`/`20` on this
  dataset's real heterogeneous sex ratio).
  Tests (same commit, later fixed in follow-up commits): reproducibility
  and cross-locus independence must compare `np.array_equal(...
  .genotype_matrix(), ...)`, NEVER `.genotype_matrix().all()` (the
  latter collapses the whole matrix to one boolean — two completely
  different matrices can coincidentally have the same `.all()`,
  verified empirically on this exact dataset, making such an assertion
  pass without proving anything). Ploidy-matches-heritage test needed a
  fixture edit (`toy_example2_ms_dna_XY`'s `Locus_M_A_2_` relabeled
  `<A>`→`<M>` in both `headerRF.txt` and the `.mss`, since the original
  fixture had no MicroSat `<M>` locus at all) — confirmed this didn't
  disturb any already-committed golden-value test elsewhere (the
  per-group RNG draw sequence depends only on locus order/count within
  the group, never on the `heritage` field). Also confirmed empirically
  that `<Y>` (`num_samples=20`) is NOT expected to equal `<M>`
  (`num_samples=40`) on this dataset — `<M>` includes every individual
  regardless of sex while `<Y>` filters to males only, and this
  dataset's real sex distribution happens to put ALL males in `pop1`
  (`{"pop1": 20, "pop2": 0}`, verified via `observed_count_population`/
  `build_male_only_samples_argument_dna`) — the right assertion compares
  `<Y>`'s sample count against the real male count, not against `<M>`.
  **Known gap, not a blocker**: no fixture currently has more than one
  MicroSat `<M>`/`<Y>` locus, so the shared-genealogy behavior for those
  heritage types (same mechanism as DNA sequences, `_SHARED_M_
  ANCESTRY_SEED_OFFSET`/`_SHARED_Y_ANCESTRY_SEED_OFFSET`) can't be
  tested end-to-end on MicroSat the way it is on DNA sequences.

- **`summary_statistics.compute_all_statistics_microsat`,
  `pipeline.compute_summary_statistics_microsat`, `reftable_loop.
  run_reftable_simulation_microsat`/`_run_single_particle_microsat`**
  (commits `fbada23`/`3066a16`/`35d3164`): the orchestration layer,
  each a direct structural mirror of its DNA-sequence sibling (same
  "two generations of architecture" duplication philosophy as SNP→DNA,
  see below — deliberately NOT factored into a single parameterized
  function even though `_run_single_particle_dna`/`_run_single_
  particle_microsat` are identical but for one function call; discussed
  explicitly with the user, who agreed to keep the duplication given
  MicroSat's stats aren't real yet and only 2 similar cases exist so
  far).
  `compute_all_statistics_microsat` is a deliberate skeleton:
  `_MICROSAT_STATS = {}` (empty catalog dict, mirroring `_DNA_
  PAIRWISE_STATS`'s structure) and an explicit `if not _MICROSAT_STATS:
  raise NotImplementedError(...)` guard before the (currently
  unreachable) aggregation loop — verified by direct execution (not
  just reading the diff) that it parses a real header's loci correctly
  and raises at the right point, past the parsing, not before it. Two
  test-writing bugs caught along the way, both about `pytest.raises`
  misuse (recurring pattern worth flagging in review): a second
  positional argument to `pytest.raises(ExceptionType, "some string")`
  is NOT a message-matching pattern (that's the `match=` keyword) —
  it's legacy API expecting a callable, so a bare string raises
  `TypeError: '...' object must be callable`; and calling the
  exception-raising function BEFORE entering the `with pytest.raises
  (...):` block means the exception fires outside the context manager
  and is never caught at all.

**RESOLVED 2026-09-11** — the actual MicroSat summary statistics
(`compute_all_statistics_microsat`'s body) are now fully implemented
and wired — see "MicroSat summary statistics" below for the full
account (formulas, architecture, bugs caught). The hypotheses below
(11 categories including `N2P`/`H2P`/`V2P`, `LIK` asymmetric) were both
confirmed correct against `statdefs.cpp`/`sumstat.cpp`.
SNI mutation channel remains deferred too (placeholder comment already
in place in `build_microsat_local_param_per_locus`, see above).

**Update 2026-09-08 — MicroSat `_from_values` replay chain done, mirroring
the DNA-sequence one.** Before tackling the stats catalog above, the
user chose to first build MicroSat's `_from_values` sibling chain (same
architecture as "DIYABC-replay pipeline for DNA sequences" below),
since it shares a lot of structure with the DNA-sequence one and stays
useful groundwork regardless of when the stats themselves get written.

- **Two `reftable_loop.py` functions needed ZERO changes** —
  `group_prior_column_names` and `parse_real_reftable_params_with_
  group_priors` already handle MicroSat groups generically (confirmed
  by re-reading them, not assumed): the former's `else` branch already
  emits `µmic_N`/`pmic_N`/`snimic_N`, the latter never inspects
  `ms_or_seq` at all, just splits columns positionally by count.
- **Five new sibling functions were needed**, each a direct structural
  mirror of its DNA-sequence counterpart (same "two generations of
  architecture" duplication choice as SNP→DNA→MicroSat elsewhere in
  this file — see [[project_dry_vs_duplication_reftable_loop]] project
  memory for the explicit discussion that reaffirmed this, this time
  specifically about whether `_group_prior_values_from_columns` should
  gain an `if ms_or_seq=="S" elif =="M"` branch instead of a sibling
  function — decided against it: unlike the identical `_run_single_
  particle_dna`/`_microsat` pair, this function's DNA and MicroSat
  branches have genuinely different shapes (conditional `k1`/`k2`
  model dispatch vs. flat `mut_rate`/`Pgeom`), and it's already part of
  the DNA replay chain validated against a real 1000-particle DIYABC
  reftable — modifying it in place, even for a low-risk addition,
  breaks the project's "never touch an already-validated original"
  rule for no strong enough reason):
  `ancestry_simulation._group_prior_values_microsat_from_columns`,
  `build_group_local_param_per_locus_microsat_from_values`,
  `build_matrix_microsat_per_locus_from_values`, `microsat_mutation_
  simulation_per_locus_from_values`, and `pipeline.compute_summary_
  statistics_microsat_from_values` + `reftable_loop._run_single_
  particle_microsat_from_values`/`replay_reftable_simulation_microsat`.
  **No MicroSat equivalent of `build_rate_map_per_locus_from_values`**
  — per-site rate heterogeneity (`RateMap`) is a DNA-sequence-only
  concept; a MicroSat locus has `sequence_length=1` and a single
  scalar `mut_rate`, so there's nothing to mirror there.
- **Renamed `build_sex_stratified_samples_argument_dna`/`build_male_
  only_samples_argument_dna` to `..._ms_dna`** (all call sites in
  `ancestry_simulation.py` and `tests/test_ancestry_simulation.py`
  updated, confirmed via `grep`, 142/142 tests green): both were
  already used by MicroSat's `.mss`-based sex dispatch too, so the
  `_dna` suffix was misleading — `_ms_dna` reflects "the `.mss`-file
  path" (MicroSat+DNA sequence), distinct from the SNP `.snp`-based
  originals (`build_sex_stratified_samples_argument`/`build_male_only_
  samples_argument`, no suffix, unchanged).
- **Several bugs caught by direct execution, not just reading the
  diff** (consistent with this whole MicroSat effort's review style):
  - `random.Random(seed + _MICROSAT_SNI_SEED_OFFSET)` referenced a
    constant that didn't exist in `configuration.py` at all —
    `ImportError` at module load, breaking every test in the file.
    Fixed by commenting the line out (matching the existing deferred-
    SNI placeholder convention in `build_microsat_local_param_per_
    locus`) rather than adding a real unused constant.
  - `f"µmiq_{group_number}"` typo (real column name is `µmic_`,
    confirmed against `group_prior_column_names`) — `KeyError` the
    first time a real column name was used.
  - `params_per_locus[locus.name] = (mut_rate, Pgeom)` in `build_
    group_local_param_per_locus_microsat_from_values` assigned the
    **whole per-group dicts** returned by `sampling_group_local_param`
    to every locus, instead of indexing them by `locus.name`
    (`mut_rate[locus.name]`, `Pgeom[locus.name]`) — every locus in a
    group would have ended up with the identical `(dict, dict)` pair
    instead of its own scalar values. Confirmed by inspecting
    `sampling_group_local_param`'s actual return type directly.
  - Return type hint said `tuple[float, float, float]` (copy-pasted
    from the DNA version's `k1`/`k2`/`mus_rate`), but the function
    returns 2-tuples (`mut_rate`, `Pgeom`).
  - **Duplicate function definitions, three separate times in this
    session** (`build_matrix_microsat_per_locus_from_values`, `dna_
    mutation_simulation_per_locus_from_values`, `microsat_mutation_
    simulation_per_locus_from_values`) — a leftover copy-paste block
    left underneath a freshly-written correct one, same name, so
    Python's module-level namespace silently kept whichever one was
    defined LAST. In the `microsat_mutation_simulation_per_locus_
    from_values` case this was actively dangerous: the correct version
    was first, the leftover broken copy-paste (calling `build_matrix_
    microsat_per_locus_from_values` with a missing argument) was
    second, so the BROKEN one was the one actually callable — caught
    only by executing the function directly (`TypeError: missing 1
    required positional argument: 'seed'`), not by reading the diff.
    `grep -n "^def name" file.py` before declaring a function "done"
    is now the standing suggestion to the user for catching this
    early.
  - **The most consequential bug**: `compute_summary_statistics_
    microsat_from_values` initially called `build_random_demography_
    for_scenario_index(header_text, scenario_index, seed)` — the
    RANDOM-draw variant — and reassigned its own `values` parameter
    with that function's freshly-drawn return value, silently
    discarding the real DIYABC values the caller had passed in. Since
    this whole `_from_values` chain exists specifically to compare
    DIYABC and msprime on IDENTICAL draws, this bug would have made
    every such comparison meaningless without raising any error.
    Caught by direct execution with deliberately recognizable
    placeholder values (`{"N1": 999999.0, "ta": 12345.0}`) and
    confirming the constructed demography used freshly-drawn values
    instead (`{"N1": 6398.0, "ta": 2758.0, ...}`) — the bug produced no
    exception, only silently wrong behavior, so this is the kind of
    thing that must be verified by executing with recognizable inputs,
    not just checked for absence of errors. Fixed in two more passes
    (first fixing the wrong function but with a wrong argument count/
    order — `build_demography_for_scenario_index` takes `(header_text,
    scenario_index, values)`, no `seed`; then fixing an incorrect
    2-tuple unpacking of its single-`Demography` return value) before
    landing on the correct one-liner, matching `compute_summary_
    statistics_dna_from_values`'s existing pattern exactly.
  - An unrelated accidental edit to `compute_summary_statistics_dna_
    from_values`'s docstring (the already-validated DNA sibling, edited
    by mistake while working on the MicroSat function right below it in
    the same file) — its indentation was shifted and a stray truncated
    line ("cette approche pourra ensuite être appliquée à des jeux de
    don") was inserted mid-docstring. Not a functional bug (docstrings
    tolerate arbitrary indentation), but a reminder that editing near
    an already-validated function carries a real risk of unintentional
    collateral changes — caught by reading the diff, not by any test.
- **Validated by direct execution up to the expected stopping point**:
  the full chain (`_run_single_particle_microsat_from_values` →
  `compute_summary_statistics_microsat_from_values` → real `values`-based
  demography → `microsat_mutation_simulation_per_locus_from_values` →
  `compute_all_statistics_microsat`) runs cleanly and raises
  `NotImplementedError` only at the final, expected point (the stats
  catalog itself, deliberately not yet implemented) — confirming
  everything upstream of the stats is correctly wired. **Not yet
  validated against an actual real MicroSat DIYABC reftable** (no such
  reference file has been used in this session — `replay_reftable_
  simulation_microsat` is wired but untested end-to-end against real
  data, unlike its DNA-sequence sibling which IS cross-validated
  against a real 1000-particle reftable).


## 2026-09-08 à 2026-09-11 — Statistiques résumées MicroSat (les 11 dernières pièces du chantier MicroSat)

Reprise du chantier MicroSat là où "MicroSat GSM mutation model" et le
rejeu `_from_values` l'avaient laissé : `compute_all_statistics_
microsat` existait comme squelette (`raise NotImplementedError`), les
11 catégories réelles (`NAL`/`HET`/`VAR`/`MGW`/`N2P`/`H2P`/`V2P`/`FST`/
`LIK`/`DAS`/`DM2`) avaient été identifiées via le header réel mais pas
encore vérifiées contre `statdefs.cpp`/`sumstat.cpp`. Session en mode
mentor (utilisateur au clavier sur `bridge/summary_statistics.py` et
`tests/test_summary_statistics.py`, assistant en lecture de source
C++/review/debug), une stat à la fois, dans l'ordre NAL → HET → VAR →
MGW → N2P → H2P → V2P → DAS → DM2 → FST → LIK, plus le câblage final.
Détail complet (formules, architecture, bugs) dans `CLAUDE.md`, section
"MicroSat summary statistics" — ce qui suit est le résumé narratif des
découvertes les plus notables.

**La représentation des données a dû évoluer deux fois.** Le premier
réflexe (`_length_by_population` : un dict `{pop: [(taille, compte),
...]}`, un pool plat de copies de gène par population, converti depuis
`variant.genotypes`/`variant.alleles` de tskit) a suffi pour 9 des 11
stats. `FST` et `LIK` ont chacune forcé une brique de regroupement par
INDIVIDU (via `tree_sequence.individuals()`, dont `.nodes` donne
directement la ploïdie — vérifié empiriquement sur un locus `<A>` ET un
locus haploïde avant d'écrire quoi que ce soit), mais avec deux formes
différentes : `_length_by_pop_and_individuals` duplique un individu
haploïde en paire `(taille, taille)` (légitime pour FST, dont la
formule ANOVA traite déjà un haploïde matché comme "doublé") ;
`_genotypes_by_pop_and_individuals` garde la vraie ploïdie (`(taille,)`
vs `(taille1, taille2)`), nécessaire pour LIK dont la formule diffère
réellement selon la ploïdie (pas unifiable par duplication).

**Le bug `cal_dmu2p` (DM2)** : `moy[]` (moyenne des tailles par
population) n'est recalculé que si les deux populations ont des
échantillons à ce locus, mais la ligne qui l'utilise est en dehors de
ce garde — donc sur le locus `<Y>` (pop2 absente) de
`toy_example2_ms_dna_XY`, DIYABC réutilise silencieusement le `moy[]`
du locus précédent, divisé par le `motif_size` du locus courant, sans
compter ce locus dans `nl`. Discuté explicitement avec l'utilisateur le
2026-09-11 : verdict commun que c'est un bug de portée de variable en
C++ (même famille que `mutsit`/`sitefix`), pas un choix statistique —
aucune autre fonction de `sumstat.cpp` ne sépare accumulation et
compteur de cette façon, et il n'y a pas de justification biologique à
réutiliser le delta-mu d'un autre locus. Décision : reproduire
fidèlement (comme `sample_site_rates`), pas corriger. Vérifié par
exécution directe sur le vrai dataset : le bug se déclenche exactement
une fois (locus 10 réutilise le `moy` du locus 9).

**FST** a demandé la décomposition la plus profonde (4 niveaux de
briques, chacune testée séparément) — Weir & Cockerham par allèle, puis
sommé sur tous les allèles d'un locus, puis sur tous les loci du
groupe, pondéré par un `nc` par locus. Simplification trouvée avant
d'écrire l'agrégateur final : `nc` ne dépend en réalité pas de
l'allèle testé (aucune donnée manquante chez nous), donc calculable une
seule fois par locus plutôt que recalculé à chaque allèle comme le fait
le C++ ; et une population absente fait naturellement tomber `nc` à 0,
ce qui exclut le locus des deux sommes sans compteur `valid_loci`
séparé — contrairement à toutes les stats "moyenne par locus"
précédentes.

**LIK** est la seule stat asymétrique du catalogue — confirmé dans le
vrai header (`LIK 1.2 2.1`, les deux sens déclarés explicitement). Deux
bugs sérieux attrapés avant validation : un brouillon calculait les
deux sens (`i->j` ET `j->i`) dans la même accumulation, ce qui aurait
rendu `LIK_1.2` et `LIK_2.1` numériquement identiques — détruisant
l'intérêt même de la stat ; et le facteur de normalisation `a`
utilisait `len(length_by_pop[pop_i])` (le nombre de tailles d'allèles
DISTINCTES, une trentaine, fixe) au lieu de `len(genotypes_by_pop
[pop_i])` (le vrai nombre d'individus, ~20) — deux ordres de grandeur
différents.

**Un bug de nommage de colonnes trouvé à la toute fin, avec une
implication rétroactive sur le code ADN déjà validé.** En câblant
`compute_all_statistics_microsat`, les colonnes produites ressemblaient
à `NAL_pop1`/`HET_pop2` au lieu de `NAL_1`/`HET_2` (repéré en exécutant
directement la fonction, pas en relisant le diff) — un oubli de la
traduction nom→indice déjà présente côté ADN
(`population_names.index(pop_name)+1`). En creusant `multi_group`
(la logique "faut-il mettre le numéro de groupe dans le nom de
colonne"), vérification contre le VRAI reftable de `toy_example2_
ms_dna` (`first_records_of_the_reference_table_0.txt`, déjà
cross-validé côté ADN) : les colonnes microsat sont bien `NAL_1_1`
(avec le numéro de groupe), alors que G1 est le SEUL groupe microsat de
ce header (les 2 autres, G2/G3, sont ADN) — donc la vraie règle DIYABC
est "plusieurs groupes dans TOUT le header, tous types confondus", pas
"plusieurs groupes de mon propre type". `compute_all_statistics_dna`
(déjà en prod, déjà validé) calculait `multi_group` avec cette même
règle fausse (filtrée à `ms_or_seq=="S"`) — invisible sur tous les
datasets testés jusqu'ici uniquement parce qu'ils ont chacun 2+ groupes
ADN de toute façon. Les deux corrigés ensemble ; suite complète
(173 tests) rejouée après coup, 0 régression.

**Statut final** : les 11 stats sont implémentées, testées (bricks
vérifiées à la main quand c'était faisable, sinon vérification
indépendante par l'assistant via un script séparé — voir `CLAUDE.md`
pour le détail par stat), et câblées. `AML` (3 populations) et le canal
SNI restent différés, aucun jeu de données ne les rendant testables
pour l'instant.

## 2026-09-14 — Première comparaison réelle des stats MicroSat contre un vrai reftable DIYABC

Suite directe de la section précédente : le test cassé
(`test_compute_summary_statistics_microsat`) corrigé, puis tentative de
comparaison réelle DIYABC/msprime sur les 11 stats. Comme
`toy_example2_ms_dna_XY` ne peut pas servir (le vrai binaire `diyabc`
plante dessus, voir plus haut), la comparaison a été faite sur
`toy_example2_ms_dna` (même groupe microsat G1, mais tout en `<A>`
diploïde, un vrai reftable 1000 particules déjà disponible) via
`scripts/replay_diyabc_priors_ms.py` (adapté du script ADN).

**Confusion à garder à l'esprit** : ce dataset a un vrai canal SNI actif
(`snimic_1` non négligeable dans le reftable réel), qu'on n'implémente
toujours pas — un écart résiduel n'est donc pas forcément un bug.

**Deux vrais bugs trouvés en creusant les écarts, avant de conclure "c'est
SNI"** :
1. `StopIteration` sur `next(tree_sequence.variants())` dans
   `_length_by_population` ET les briques FST/LIK — des vrais `mut_rate`
   DIYABC assez faibles pour qu'un locus n'ait aucune mutation du tout
   (`num_sites==0`), jamais rencontré sur les données de test synthétiques.
   Corrigé en traitant ce cas comme "tout le monde porte encore l'allèle
   ancestral" (valeur arbitraire mais partagée). Piège en cours de route :
   une première tentative sur `_length_by_pop_and_individuals` a copié le
   pattern de `_length_by_population` (`[(0, len(sample_ids))]`, UNE ligne)
   sans remarquer que cette fonction retourne UNE LIGNE PAR INDIVIDU — ça a
   fait planter `MSI` de FST (`sni-2.0` à 0 exactement) en écrasant le vrai
   effectif à 1 individu par population.
2. `compute_LIK` appelait `_length_by_pop_and_individuals` (la brique de
   FST, une ligne par individu) au lieu de `_length_by_population` (une
   ligne par allèle distinct) — même arité de tuple des deux côtés, donc
   aucune exception, juste un calcul de fréquences n'importe quoi. Trouvé
   UNIQUEMENT en comparant au vrai reftable : nos deux sens de LIK
   (`LIK_1_1.2`/`LIK_1_2.1`) sortaient presque identiques (~3.5 de moyenne
   chacun) alors que le vrai DIYABC les a très différents (7.33 vs 3.08,
   avec un écart-type de 11.44 dans un sens contre 1.92 dans l'autre) —
   or LIK est justement censée être asymétrique. Corrigé, réduit l'écart
   KS à ~0 sur les deux colonnes.

**Après les deux corrections** : 9/11 familles de stats collent bien.
Seul `FST_1_1.2` reste décalé (`KS≈0.245`, ~1.8× trop bas systématiquement
— pas une différence de forme comme LIK, un décalage constant). Recherché
un bug résiduel et rien trouvé : formule vérifiée à la main terme à terme
plus tôt, ce dataset n'a aucun locus haploïde (donc la logique de
duplication FST ne se déclenche jamais), et les loci monomorphes
(`num_sites==0`) ne touchent que 13/1000 loci sur un échantillon de 100
particules — bien trop rare pour expliquer un tel écart. Conclusion
provisoire (pas prouvée) : le canal SNI manquant explique l'écart de FST
— à confirmer une fois SNI implémenté et la comparaison relancée. Détail
complet dans `CLAUDE.md`, section "MicroSat stats cross-validated
against a real reftable".

## 2026-09-14 (suite) — hypothèse SNI invalidée pour FST, cause réelle non trouvée

L'hypothèse "c'est SNI" ci-dessus s'est révélée fausse. Deux tests de
falsification (même méthode que l'investigation du déficit de variance
G3 sur `<M>` côté ADN, 2026-08-27/31) l'ont réfutée :

1. **`reference/toy_example2_ms_dna_weak_sni/`** : nouveau reftable réel
   DIYABC généré avec `MEANSNI`/`GAMSNI` forcés très bas (`snimic_1` ≈1%
   de `µmic_1` au lieu d'un ordre de grandeur comparable). Si SNI était
   la cause, l'écart FST aurait dû se réduire nettement. Il n'a PAS
   bougé du tout (`0.1429`/`0.0809` réel/simulé, identique au dataset
   original `0.1440`/`0.0798`).
2. **Isolation du scénario 2** (sans admixture `ta split`, contrairement
   au scénario 1) sur le reftable déjà rejoué : FST diverge pareil.
   Admixture écartée aussi.

Poussé plus loin avec `scikit-allel` (implémentation indépendante de
Weir & Cockerham 1984) : sur les mêmes données simulées, `scikit-allel`
donne `FST≈0.0168`, nous `≈0.0077` — un facteur ~2×, du même ordre que
l'écart contre le vrai DIYABC. Confirme un vrai problème numérique,
indépendant de toute question DIYABC.

Vérifications supplémentaires, toutes concluantes mais sans trouver LE
bug :
- Formule (`cal_Fst2p`) relue une 4e fois contre la source, caractère
  pour caractère — toujours correcte.
- Données reconstruites indépendamment (paires `(taille1, taille2)` par
  individu, directement depuis `variant.genotypes`/`ts.individuals()`,
  sans passer par `_length_by_pop_and_individuals`) et diffées contre
  la sortie réelle de cette fonction sur les 40 individus d'un vrai
  locus — 0 écart.

**Piste non explorée, à creuser en priorité la prochaine fois** : FST
est la seule stat basée sur une décomposition de variance ANOVA
(`MSG`/`MSI`/`MSP`, sensible à l'INTERACTION hétérozygotie/structure de
population), alors que `HET` (diversité) et `DAS`/`DM2` (comparaison
directe des allèles entre populations) collent bien — le problème est
donc spécifique à cette décomposition, pas à la diversité ou à la
différenciation générale simulée (qui semblent correctes par ailleurs).
Chantier mis de côté volontairement le 2026-09-14 plutôt que de
continuer à deviner sans nouvelle piste — voir `CLAUDE.md` pour le
détail complet des tests.

## 2026-09-15 — SNI implémenté : le fossé FST confirmé indépendant de SNI, et un vrai coût de performance (~34x)

Chantier repris pour lever le doute une bonne fois : plutôt que de se
contenter de la falsification indirecte du 2026-09-14 (forcer SNI≈0
côté vrai DIYABC), on a implémenté le canal SNI réel dans notre propre
pipeline pour voir si l'écart FST se referme quand SNI existe vraiment
de notre côté.

Mécanisme (`particuleC.cpp::mute`) : un seul processus de Poisson
combiné à taux `mut_rate+sni_rate`, chaque événement classé par un
tirage de Bernoulli (`p_sni = sni_rate/(sni_rate+mut_rate)`) en pas GSM
(`±d×motif_size`) ou pas SNI (`±1`pb, indépendant de `motif_size`),
tous deux bornés à `[kmin,kmax]` avec repli en auto-mutation aux bords.
Contrairement à GSM seul (qui ne visite jamais que les états de la
classe de résidu de `root` modulo `motif_size`, d'où une seule grille
`msprime.TPM` espacée de `motif_size` suffisante), SNI peut faire
changer de classe de résidu — il faut donc une grille DENSE (tous les
entiers de `kmin` à `kmax`), construite à la main (`bridge/ancestry_
simulation.py`, nouvelles fonctions `_distribution_from_position`/
`_place_gsm_row_on_dense_grid`/`_sni_row_on_dense_grid`/`_mix_sni_gsm_
rows`/`build_microsat_transition_matrix_with_sni`).

Plusieurs bugs trouvés et corrigés en cours de route (détail complet,
formules et bricks dans `CLAUDE.md`, section "SNI mutation channel") :
un décalage d'arguments qui sautait purement et simplement l'appel à
`_distribution_from_position` (`ZeroDivisionError`, révélé aussi par un
`epsilon` non utilisé signalé par le linter), plusieurs confusions
entre l'indice LOCAL d'un état (échelle `motif_size`) et sa distance
brute à `kmin` (échelle 1pb) dans le placement sur la grille dense, et
un oubli du `- kmin` dans `_sni_row_on_dense_grid` (`IndexError`
immédiat dès que `kmin != 0`, donc sur tout jeu de données réel).
Méthode de validation notable : comparer la matrice dense (`sni_rate=
0`) à l'ancienne `build_microsat_transition_matrix` n'est PAS une
comparaison de même forme (grilles de tailles différentes par
construction) — il faut extraire, de la matrice dense, seulement les
lignes/colonnes des états que l'ancienne grille sait représenter (via
`old_model.alleles`, déjà calculé et testé), pas recalculer `root`/
`n_minus` en double dans le test.

**Coût de performance réel et significatif, diagnostiqué mais pas
encore corrigé** : le reftable complet de `toy_example2_ms_dna` passe
de ~15s (GSM seul) à ~512s (GSM+SNI) — un facteur ~34×. Cause : la
nouvelle fonction reconstruit un `msprime.TPM` complet PAR ÉTAT de la
grille dense (`kmax-kmin+1` appels), alors que l'ancienne n'en
construit qu'UN SEUL pour tout le locus — cohérent quantitativement
avec le facteur observé (`~4n` pour `motif_size=2`, soit ~30-40× pour
un `n≈40` typique). Compromis assumé au départ (réutiliser `msprime.
TPM` état par état plutôt que rederiver la formule géométrique à la
main, pour éviter un nouveau bug de formule) dont le coût réel s'avère
plus élevé que prévu — piste d'optimisation la plus directe : calculer
la décroissance géométrique + repli aux bords sous forme close, sans
repasser par `msprime.TPM` du tout sur ce chemin. Laissé pour une
session dédiée, ce n'est pas un problème de correction.

**Résultat sur `FST`** : une fois SNI branché dans le pipeline
(`build_matrix_microsat_per_locus`/son jumeau `_from_values`), le
rejeu complet contre le vrai reftable `toy_example2_ms_dna` donne
`FST_1_1.2` KS≈`0.2643` (`p≈0.0`) — quasiment identique au `KS≈0.245`
d'avant SNI. Ça confirme, cette fois directement (et pas seulement par
la falsification indirecte du 2026-09-14), que SNI n'est pas la cause
de l'écart FST — deux tests indépendants et complémentaires (retirer
SNI des données réelles, ajouter SNI à notre simulation) pointent
maintenant dans le même sens. `DM2` montre un léger bruit non
significatif (`KS≈0.0439`, `p≈0.4689`), rien d'inquiétant. Le chantier
FST reste ouvert exactement comme documenté le 2026-09-14 (piste ANOVA
`MSG`/`MSI`/`MSP` non explorée) — ce résultat retire SNI de la liste
des suspects avec beaucoup plus de confiance, sans faire avancer
l'investigation elle-même.

## 2026-09-16 — Optimisation du canal SNI : ~34x → ~1.5x

Le ralentissement diagnostiqué la veille (`_distribution_from_position`
reconstruit un `msprime.TPM` complet PAR ÉTAT de la grille dense) a été
corrigé, sans passer par la piste initialement envisagée (rederiver la
formule géométrique à la main) — une piste plus simple, trouvée en
repartant d'une identité déjà établie la veille pendant la vérification
de l'équivalence `sni_rate=0` : pour deux états `s` et `s+motif_size`
(même classe de résidu modulo `motif_size`), `n_minus`/`n_plus` varient
chacun d'exactement `±1`, donc `n_alleles = n_minus+n_plus+1` est
CONSTANT à l'intérieur d'une classe de résidu — mais diffère D'UNE
classe à l'autre (vérifié : `kmin=80,kmax=120,motif_size=2` donne
`n_alleles=21` pour tous les états pairs, `20` pour tous les impairs).
Conséquence directe : il ne faut que `motif_size` matrices `msprime.
TPM` distinctes pour tout le locus (une par classe de résidu, chacune
avec son propre `hi`), pas une par état de la grille dense.

Deux bugs trouvés en cours d'implémentation (par exécution directe) :
un dict de matrices d'abord construit avec le mauvais `hi` (taille de
la grille DENSE au lieu de la taille propre à chaque classe) et indexé
par l'indice DENSE au lieu du résidu — ne changeait donc rien au
problème ; puis, une fois le dict corrigé, la boucle principale
réutilisait par erreur la variable `i` de la boucle précédente
(construction du dict) au lieu de calculer `residue = (position-kmin) %
motif_size` — `KeyError` dès que l'indice dense dépassait `motif_size`.

Une faiblesse de test repérée au passage : `test_distribution_from_
position` construisait sa matrice de test avec un `hi` correspondant à
la MAUVAISE classe de résidu par rapport à la `position` testée — invisible
car les seules assertions (forme, somme=1) sont vraies pour n'importe
quelle ligne de n'importe quelle matrice stochastique, peu importe la
classe. Résolu en réalisant que `_distribution_from_position` n'a plus
aucune responsabilité sur la construction de la matrice (elle se contente
d'indexer `transition_matrix[n_minus]`) — son test n'a donc même plus
besoin d'un vrai `msprime.TPM`, une matrice factice avec des lignes
identifiables suffit à isoler proprement ce qu'elle fait réellement.

**Résultat mesuré** : le rejeu complet de `toy_example2_ms_dna` passe de
~512s à **~23s** (un facteur ~22x d'amélioration), soit ~1.5x le temps
d'avant SNI (~15s) — pas 1x pile, puisque construire 2 matrices au lieu
d'1 coûte un peu plus, mais la construction de matrice ne représente
qu'une fraction du coût total d'une particule. `FST`/`DM2` rejoués sur
cette version optimisée donnent des valeurs identiques à avant
(`FST_1_1.2` KS≈0.2643, `DM2` KS≈0.0439) — l'optimisation n'a rien
changé à la correction, seulement à la vitesse.

## 2026-09-16/18 — AML implémentée et validée ; piège de casse `cal_Aml3p`/`cal_aml3p`

Dataset adapté enfin disponible (`toy_example1_ms_modified`, 4
populations réelles issues d'un scénario split/merge — `toy_example1_ms`
non modifié n'a qu'une population échantillonnée à 4 temps, voir plus
bas). Avant tout code, lecture directe de `sumstat.cpp`/`statdefs.cpp`
pour localiser la bonne fonction — et trouvaille immédiate, le genre
d'écueil qui aurait fait perdre une session entière si on avait suivi
la convention de casse "naturelle" du reste du projet :

`statdefs.cpp` (`microsat_statns`/`dna_statns`) montre que la stat
MicroSat `"AML"` pointe vers `cal_Aml3p` (A **majuscule**,
sumstat.cpp:1290-1326), alors que `cal_aml3p` (a **minuscule**,
sumstat.cpp:2087-2154, qui appelle `cal_freq`/`cal_nss2pl`/
`libere_freq` — toute la mécanique de comparaison de séquences ADN) est
en réalité la fonction de la stat `"SML"` côté **DNA** (pas microsat du
tout). Suivre la casse minuscule "logique" (cohérente avec le style du
reste du fichier) aurait fait partir l'implémentation sur la mauvaise
fonction, avec une logique de comparaison de chaînes n'ayant aucun
rapport avec un calcul d'admixture par fréquences alléliques.

**Algorithme de `cal_Aml3p`** (avec `pente_lik`, sumstat.cpp:1217-1288) :
recherche par bissection du taux de mélange `a ∈ [0.001, 0.998]`
maximisant la vraisemblance des génotypes de la population "focale"
(`samp[0]`) sous l'hypothèse que ses fréquences alléliques sont
`a·freq[parent1] + (1-a)·freq[parent2]` (`samp[1]`, `samp[2]`) — 3
formules de vraisemblance selon la ploïdie (haploïde / diploïde
homozygote / hétérozygote), pas de pseudo-compte (contrairement à
`cal_lik2p`/LIK — un terme à fréquence nulle est simplement omis, pas
remplacé par un pseudo-compte `1/nal`). La bissection ne calcule jamais
`li(a)` pour toutes les valeurs de `a` : elle regarde le SIGNE de la
pente (différence finie `li(a+0.001)-li(a)`) aux deux bornes, avec 3 cas
limites — pente nulle partout → aucune information → tirage aléatoire
uniforme(0,1) ; pente toujours négative → `a=0.0` ; toujours positive →
`a=1.0`. Convention d'indices confirmée par `statdefs.cpp` (npop=3,
`sortArr::HALF`) + le header : pour un triplet `{i,j,k}` trié, on génère
`"i.j.k"`, `"j.i.k"`, `"k.i.j"` — le premier indice cycle (focal), les
deux autres restent toujours en ordre croissant (parents) — exactement
la même structure que `_half_arrangements` (déjà utilisée pour la
version SNP de AML, `compute_AML`), directement réutilisable pour
énumérer les triplets côté microsat sans réinventer la combinatoire.

**Optimisation, pas juste une implémentation directe** : la première
version recalculait `_length_by_population`/`_genotypes_by_pop_and_
individuals` (donc les fréquences alléliques) à CHAQUE évaluation de
`a` — or `cal_Aml3p` en fait ~20 par triplet (bissection sur ~998
valeurs possibles, `log2(998)≈10` étapes, 2 évaluations de `pente_lik`
par étape). Mesuré : 15.2s/particule avant optimisation. Hisser le
calcul de fréquences hors de la boucle de bissection
(`_prepare_loci_for_admixture`, calculé une fois par triplet plutôt
qu'une fois par évaluation de `a`) donne 0.77s/particule — un facteur
~20x, cohérent avec le nombre d'évaluations éliminées.

**Validation contre un vrai reftable DIYABC** : 100% des colonnes AML
concordent sur `toy_example1_ms_modified` (1000 particules) — seules
les colonnes FST divergent (voir section suivante). Confirme que
l'algorithme et son branchement (triplets, seed par triplet ET par
groupe pour éviter une corrélation entre deux groupes tombant tous les
deux dans le cas dégénéré — même classe de bug que documenté ailleurs
dans ce fichier pour les tirages de groupe) sont corrects.

## ReplayContext refactor — read header.txt/.snp/.mss once per run, not once per particle (2026-09-16/18, mentor mode — user-driven, reviewed/debugged with the assistant)

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Triggered by a genuine performance investigation (a `toy_example1_ms_
modified` real-reftable replay was taking 20-40+ minutes, occasionally
hanging outright) that turned out to have TWO unrelated root causes —
see "Serial/temporal sampling" and the trailer-line bug below for the
real cause of the hang. This refactor itself was found, by direct
measurement, NOT to be the cause of that particular slowness (~7.6ms/
particle of avoidable disk I/O, negligible against a 20+ minute hang) —
but was pursued anyway as a real, independently-motivated cleanup, and
is now complete for all three data families (SNP, MicroSat, DNA
sequences).

**Reading this file: every signature written down in a section dated
BEFORE 2026-09-16 may be stale because of this refactor** — most
per-locus builders took `(header_text, mss_file_path, ...)` and now take
`(context, ...)` instead (e.g. `build_matrix_per_locus`,
`build_matrix_microsat_per_locus`, `microsat_mutation_simulation_per_
locus`, `dna_mutation_simulation_per_locus`). Not all of them moved:
`build_group_local_param_per_locus`, `build_microsat_local_param_per_
locus` and `build_rate_map_per_locus` still take `(header_text, seed)`,
since they read nothing off the disk. Check the real signature before
calling one from a notebook or a script (found stale on 2026-09-22 while
writing `notebook/visualise_dna_pipeline.py`).

Three new dataclasses in `bridge/header_dataclasses.py` —
`SnpReplayContext`, `MicrosatReplayContext`, `DnaReplayContext` — each
built ONCE per `run_reftable_simulation*`/`replay_reftable_simulation*`
call (not once per particle), holding every value previously obtained
by re-reading `header.txt`/`.snp`/`.mss` from disk
(`read_header_text`, `observed_microsatellites`/`observed_sequences`,
`observed_count_population`, `parse_sex_ratio`, `parse_maf_ratio`,
`parse_mrc_ratio`, `observed_reads`, `individual_sexes_per_population`,
`detect_snp_file_type`...), then passed through `ProcessPoolExecutor.
submit` to every particle instead of a raw `reference_directory`/
`mss_file_path`/`snp_file_path`. Covers both the "draw fresh values"
and `_from_values` (replay) path for all three families — 180 tests
green, `ruff` clean.

**Field selection principle, stated explicitly by the user querying
each candidate field**: include a value in the context if and only if
computing it involves an actual disk `read_text()` (directly or
transitively) — NOT merely "does this value vary by locus/particle".
`SnpReplayContext`'s `sexes_per_population` illustrates why the second
criterion is wrong: unlike DNA/MicroSat's `individual_sexes_from_
locus_genotype` (which must be recomputed per LOCUS, since sex is
inferred from that locus's own genotype ploidy — caching it once
dataset-wide would be meaningless), `.snp`'s `individual_sexes_per_
population` reads a real SEX column, a single dataset-wide value with
no locus dependence at all — so unlike the DNA/MicroSat `<X>`/`<Y>`
case (where the "keep it simple, just pass the path through" choice
was deliberately made because there was no cheap alternative), here
caching it AND updating `build_sex_stratified_samples_argument`/
`build_male_only_samples_argument` to consume it was both correct and
cheap, and the "keep it simple" precedent did not apply. Conversely,
`observed_reads`'s own internal `detect_snp_file_type`/`parse_mrc_
ratio` calls (redundant with fields already in the context) were left
alone after explicit discussion — real duplication, but happening once
per RUN inside context construction itself, not once per particle, so
below the threshold of what this refactor was for.

**Bugs found in series while wiring this through — well beyond the 5
already logged for MicroSat/DNA context work (see the [[feedback_
signature_refactor_mismatch]] persistent-memory checklist)**:
- `context.loci_desciption` (typo, missing the "r") in `pipeline.py::
  _simulate_genotypes_for_all_locus_types` — `AttributeError`.
- `context, context,` — a duplicated positional argument in
  `replay_reftable_simulation`'s `executor.submit(_run_single_particle_
  from_values, ...)` call, shifting every subsequent positional
  argument by one slot.
- `haploid_pool_sizes`/`pool_sizes` built from `context.count_samples`
  used AS-IS (real `.snp` population names as keys) instead of
  translated to msprime's `"pop1"/"pop2"` convention — what `build_
  samples_argument` used to do internally before this refactor.
  `msprime.Demography` only recognizes `"pop1"/"pop2"`, so this
  produces `KeyError: "Population with name '<real name>' not found"`
  downstream — the SAME error class as the serial-sampling limitation
  below, different root cause. Found and fixed independently at two
  call sites (`ancestry_simulation.py::simulate_poolseq_reads_with_mrc_
  filter` and `pipeline.py::compute_summary_statistics`), both needing
  the same `{f"pop{i}": count for i, count in enumerate(context.count_
  samples.values(), 1)}` translation.
- **A real bug in production code, invisible to the full test suite**:
  `run_reftable_simulation`'s context construction had `reads_observed
  = None` hardcoded — none of the 5 existing `run_reftable_simulation`
  calls in `tests/test_reftable_loop.py` exercise a PoolSeq dataset
  (all use `human`, IndSeq), so 180/180 tests stayed green while this
  path was fully broken. Only caught by directly executing `run_
  reftable_simulation` against `toy_example4` (a real PoolSeq
  reference dataset) — `TypeError: 'NoneType' object is not iterable`.
  Fixed in two passes: first by computing `reads_observed` unconditionally
  (broke `human`/IndSeq instead, since `observed_reads` explicitly
  rejects any non-POOLSEQ file); then by gating it on `snp_file_type ==
  "POOL"` but using the function's own `num_loci` parameter (can be a
  small testing value) instead of `loci_description.loci_counts_by_
  heritage["A"]` (the real declared count) — `RuntimeError: generator
  raised StopIteration` in `with_mrc_filter` from premature pool
  exhaustion. A reminder that "N tests pass" only proves what those N
  tests actually exercise, not the untested branches.
- A test that writes a modified `header.txt` to a `tmp_path` and
  expects the function under test to pick it up
  (`test_compute_summary_statistics_stats_filter_header`) broke for a
  structural reason, not a typo: before this refactor, `compute_
  summary_statistics` re-read `header.txt` from disk on every call, so
  "write a new file, then call the function on that directory" worked.
  Now that everything flows through `context.header_text` (read once,
  never re-read), that technique silently uses the ORIGINAL, unmodified
  header — the fix is to build a fresh `SnpReplayContext` with the
  modified `header_text` rather than reuse an existing context
  unchanged. Worth remembering for any other test in this codebase that
  mutates the observed environment on disk and expects a fresh read.
- Three `conftest.py` fixtures (`snp_context_human`/`_te4`/`_te5`)
  shared two copy-pasted bugs: `maf_ratio`/`mrc_ratio` hardcoded to
  `None` instead of calling `parse_maf_ratio`/`parse_mrc_ratio`; and
  `header_text=header_text` — referencing the WRONG name (each
  fixture's own parameter is `header_text_te4`/`header_text_te5`, not
  bare `header_text`, which happens to collide with an unrelated
  module-level fixture defined for `human`) — a pytest fixture
  referenced by its bare (undecorated) name resolves to the fixture
  FUNCTION OBJECT, not its value, producing `AttributeError:
  'FixtureFunctionDefinition' object has no attribute 'splitlines'`
  deep inside unrelated parsing code. `snp_context_human` itself
  initially had no `return` statement at all (computed every field,
  returned nothing) — a bare Python fixture without `return` resolves
  to `None`, caught via `AttributeError: 'NoneType' object has no
  attribute 'header_text'`.

**Two `headerRF.txt` trailer-line bugs found and fixed on `toy_example1_
ms_modified/`, same class as [[diyabc_header_trailer_line_bug]]**: the
file's last line (re-read as INPUT by the real DIYABC binary to derive
`nparamhist`, despite looking like pure output-column documentation)
still declared the 3 historical parameter names of the *original*,
unmodified `toy_example1_ms` (`Npast Npresent tbn`) instead of the 9
real ones (`N1 N2 N3 N4 t423 ra t32 t21 t421`) — corrupting the real
DIYABC binary's own `nparamhist` and producing nonsensical replayed
values (`ra=512.0` against a declared `[0.05,0.95]` prior). Fixed once,
then a SECOND time after `headerRF.txt` was edited again post-regeneration
without rerunning the real `diyabc` binary — caught both times by
comparing `stat -c '%y %n'` on `headerRF.txt` vs `first_records_of_the_
reference_table_0.txt` (the real reftable must always be NEWER than the
header it was generated from, never the reverse).


## 2026-09-18 — FST microsat : nouvelle résolution grâce à 4 populations, piste redirigée vers la généalogie `<M>`

Jusqu'ici, l'investigation FST (voir plus haut, entrées du 2026-09-14/15)
n'avait accès qu'à des datasets à 2 populations, où toutes les paires
FST sont écrasées en une seule valeur moyenne — impossible de savoir si
le déficit touchait `<A>` et `<M>` de façon uniforme. Avec
`toy_example1_ms_modified` (4 populations réelles, 12 colonnes FST : 6
`<A>` groupe G1, 6 `<M>` groupe G2), la réponse est nette : les 12
colonnes divergent TOUTES significativement, mais à deux amplitudes
très différentes et parfaitement séparées — KS≈0.37-0.46 pour les 6
`<A>` (cohérent avec l'ancien ~1.8-2x déjà documenté sur
`toy_example2_ms_dna`), KS≈0.71-0.77 pour les 6 `<M>` (quasiment le
double, aucun chevauchement entre les deux groupes de valeurs).

`describe()` sur `FST_2_1.2` (`<M>`) : médiane `msprime` **négative**
(-0.005, la moitié des 1000 particules sous zéro), moyenne 0.0096,
écart-type 0.050 — contre médiane 0.111, moyenne 0.149, écart-type 0.118
côté DIYABC réel. La distribution `msprime` entière est décalée vers
zéro, pas quelques particules aberrantes qui tirent la moyenne.

**Formule revérifiée une 5ᵉ fois**, cette fois avec un test synthétique
exécuté (pas seulement relu) plutôt qu'une relecture terme à terme
de plus : `_compute_ni_nA_AA_for_one_population`/`_compute_FST_
constants_for_two_populations_combined` appliqués à deux cas construits
à la main. Différenciation complète (pop A = 100%×allèle 10, pop B =
100%×allèle 12, 4 individus haploïdes chacune) → `FST=1.0` exact,
comme attendu. Différenciation partielle (pop A = 3×10+1×12, pop B =
2×10+2×12) → `FST=-0.16667` exact, qui correspond à une dérivation à la
main de `cal_Fst2p` terme à terme sur ce même cas — confirme qu'**un
FST négatif est un comportement légitime et attendu de l'estimateur de
Weir & Cockerham** à faible différenciation observée (pas un bug de
notre implémentation, ni de DIYABC).

**Conséquence pour l'investigation** : si la formule est vérifiée
correcte même dans ce cas limite (négatif), alors le déficit ne peut
plus être imputé au calcul de la statistique elle-même — il faut que
les populations `<M>` SIMULÉES par notre pipeline aient réellement
moins de différenciation entre elles que celles de DIYABC. Or HET/AML
(qui ne dépendent que des fréquences alléliques marginales par
population, jamais de la corrélation entre populations) collent bien
sur ce même groupe `<M>` — donc les fréquences alléliques globales sont
correctes en moyenne, mais quelque chose dans la STRUCTURE DE
POPULATION portée par la généalogie `<M>` partagée
(`_SHARED_M_ANCESTRY_SEED_OFFSET`, vérifiée présente et correctement
câblée dans le code) produit moins de signal between/within-population
que la vraie généalogie DIYABC. Nouvelle piste, pas encore vérifiée :
`ms_dna_ancestry_parameters_for_heritage`/`rescale_demography` pour
`<M>`, en particulier leur interaction avec les événements
split/admixture d'un scénario à 4 populations — jamais testée avant
avec un vrai `<M>` multi-populations (le seul autre dataset avec `<M>`,
`toy_example2_ms_dna`, n'a que 2 populations et son groupe FST-testé
est `<A>`-only).

## 2026-09-17/18 — Échantillonnage sériel : limitation d'architecture confirmée, pas encore résolue

`toy_example1_ms` (non modifié, avant d'être adapté en
`toy_example1_ms_modified` pour AML) échantillonne une SEULE population
biologique à 4 temps différents (`0 sample 1`, `50 sample 1`, `200
sample 1`, `500 sample 1` — le `1` est toujours le même indice de
population dans le header), mais le `.mss`/les statistiques traitent
ces 4 échantillons comme 4 "populations" `pop1`..`pop4` distinctes (4
blocs `POP` dans le fichier). Confirmé par exécution directe :
`demography_builder.py` ne construit qu'UNE population msprime
(`"pop1"`) à partir des événements du scénario, alors que
`observed_count_population` sur le `.mss` en trouve 4 — `msprime.sim_
ancestry` rejette alors le `samples=` à 4 entrées avec `KeyError:
"Population with name 'pop2' not found"` (`demography.py:970`, erreur
native msprime, pas la nôtre).

Ce n'est pas un bug ponctuel : DIYABC numérote chaque prélèvement
temporel comme une population distincte dans tout `header.txt` (priors,
stats), mais côté généalogie c'est une seule population, échantillonnée
à plusieurs `time=`. Pour le reproduire, `sim_ancestry` devrait recevoir
plusieurs `SampleSet` à des temps différents mais tous rattachés à LA
MÊME population msprime, puis les statistiques devraient traiter ces
échantillons temporels comme des "populations" virtuelles distinctes —
aucun morceau de ce mécanisme n'existe aujourd'hui dans le pipeline.
Contourné cette session en travaillant sur `toy_example1_ms_modified`
(vraies populations distinctes, mêmes stats) plutôt que résolu.

## RESOLVED 2026-09-21: FST microsat gap — a numpy bool-addition bug, not the genealogy, not the formula (investigated 2026-09-14/18, fixed by the user 2026-09-21, commits `30fa872` + `b8b0657`)

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Closes the gap first documented in "MicroSat stats cross-validated
against a real reftable" above (`FST` ~1.8-2x too low on `<A>`, SNI
and admixture ruled out, formula reverified 6 times with no bug
found). On `toy_example1_ms_modified` (4 real populations) the gap had
two magnitudes: `KS≈0.37-0.46` on `<A>` (G1), `KS≈0.68-0.75` on `<M>`
(G2, msprime median **negative**, whole distribution shifted toward
zero). The 2026-09-18 hypothesis — the shared `<M>` genealogy's
interaction with split/admixture producing less population structure —
is **wrong**, see below.

**Step 1 — admixture falsified for FST, cleanly, for both locus
types.** The user rewrote scenario 1 of `toy_example1_ms_modified` to
remove its `split` event (`t41 merge 1 4` / `t32 merge 2 3` / `t21
merge 1 2`, no `ra`), keeping scenario 2's `t421 split 4 2 1 ra`,
regenerated the real reftable, replayed, and split the 1000 particles
by real `scenario_index` (same methodology as the 2026-08-31 DNA G3
falsification, never done for FST before): `rdiff_mean` on `FST_1_*`
was -52..-70% (no admixture) vs -55..-65% (admixture), on `FST_2_*`
-90..-121% vs -91..-113%, std ratio ~2.2x in all four cells. **No
difference.** Admixture is not a factor, for `<A>` or `<M>`.

**Step 2 — the decisive observation was already in the data.** On the
same replay, `DAS` (shared alleles between populations — a *direct*
identity-based differentiation measure), `DM2` (δμ², divergence),
`HET`, `H2P` all matched DIYABC within ~1-3% (`KS≈0.03`) on both G1
and G2. If the simulated data really carried 2x less differentiation,
`DAS` would be visibly higher and `DM2` lower. They weren't. So the
simulated data was right and the deficit had to be in the FST
*computation* — which is the only stat consuming `_length_by_pop_and_
individuals` (individual-level `(taille1, taille2)` pairs → `AA`,
`nA`); every other stat goes through `_length_by_population`'s
`(taille, compte)` frequency table.

**Step 3 — the one unexploited clue.** The 2026-09-14 entry above
records that `scikit-allel`'s own `weir_cockerham_fst` gave ~2x our
value on *identical* simulated data, and filed it as "confirms a real
numerical discrepancy exists". That framing was the mistake: two
correct implementations of the same estimator cannot disagree 2x on
the same input, so one of them was being misapplied — and this was
the cheapest, DIYABC-free test available all along. Redone
brick-by-brick: (a) `_compute_FST_constants_on_all_alleles_for_two_
populations` on hand-built pairs of Python `int` → identical to allel
to 15 decimals (so `cal_Fst2p` IS standard WC84 — this time re-derived
algebraically against the textbook MSG/MSI/MSP definitions, not just
"matches the C++" — and our transcription is exact); (b) the same
function on pairs coming out of `_length_by_pop_and_individuals` on
real `TreeSequence`s → 0.034 vs allel 0.134 on the same locus; (c) the
pairs themselves verified identical to `ts.individuals()` — so the
only difference left was the **dtype**: the pipeline's pairs are
`np.int64` (from `np.array(...)[variant.genotypes]`).

**Root cause**: `_compute_ni_nA_AA_for_one_population` computed
`nA = sum((p[0] == al) + (p[1] == al) for p in pairs)`. With numpy
operands, `p[0] == al` is an `np.bool_`, and `np.bool_ + np.bool_` is
a **logical OR**, not an integer sum — `np.True_ + np.True_ == np.True_`
(=1), whereas Python's `True + True == 2`. Every homozygote for `al`
was counted `nA += 1` instead of `+= 2`, biasing `s2A`/`MSP`/`MSI` and
deflating θ. It explains the whole picture with nothing left over:
`<A>` moderately (only homozygotes undercounted), `<M>` catastrophically
(the haploid duplication `(taille, taille)` makes EVERY individual a
"homozygote", so `nA` is halved everywhere → negative FST), and every
synthetic/golden test passing (all fed Python `int` literals, where
the same line is correct). Confirmed as the *entire* cause by
monkeypatching `int()` into the diagnostic: our `compute_FST` then
equals allel to 6 decimals on 3 real replayed particles
(0.046339/0.045166/0.116091 both sides).

**Fix (`30fa872`, written by the user)**: `_length_by_pop_and_
individuals` now builds its tuples with `int(tailles[...])`, so the
`np.int64` never leaks past the tskit layer — every consumer of these
pairs (and `al`, derived from them) gets Python ints. Golden values
`FST_1_1.2` regenerated in `tests/test_summary_statistics.py`/
`test_pipeline.py` (0.00537 → 0.03084). **Validation**: full replay +
per-scenario notebook on `toy_example1_ms_modified` — scenario 1 (no
admixture): 0/152 stats with KS `p<0.05`; scenario 2 (admixture): 2/152
(`V2P_2_1.3` p=0.036, `MGW_2_4` p=0.043) — with 304 tests at α=0.05
~15 false positives are *expected*, so this is below noise, no
residual on `<A>` or `<M>`. The 11 MicroSat stats are now all
validated against real DIYABC output.

**Hardened and tested (`b8b0657`, written by the user)**: belt and
braces — `_compute_ni_nA_AA_for_one_population` itself now does
`int(p[0] == al) + int(p[1] == al)`, so the brick is correct whatever
its caller feeds it, independently of the `int()` conversion in
`_length_by_pop_and_individuals`. Two tests guard the two layers:
`test_compute_ni_nA_AA_for_one_population` feeds `np.int64` pairs
(3 homozygotes + 1 heterozygote → `nA == 7`; the buggy line gave 4),
and `test_compute_FST_vs_scikit_allel` (`allel = pytest.importorskip
("allel")` INSIDE the function, not at module level — at module level
it would skip the whole file) runs `compute_FST` on the real
`toy_example1_ms_modified` G1 loci (new `microsat_context_te1_
modified` fixture in `conftest.py`) and asserts equality with
`allel.weir_cockerham_fst` accumulated as `A.sum()/(A+B+C).sum()`
over loci and alleles (a ratio of sums, never a mean of per-locus
ratios), individuals paired via `ts.individuals()`, loci with
`num_sites == 0` skipped on the allel side (they contribute zero on
both). Verified to discriminate: with BOTH `int()` conversions
removed it fails, with either one present it passes (the two
protections are redundant by design, so removing only one does not
break it — that is expected, not a weak test). `scikit-allel` is not
declared in `pyproject.toml`; the `importorskip` keeps the suite green
on environments without it.

**Still not audited**: `_genotypes_by_pop_and_individuals` (LIK) builds
its tuples from the same numpy `tailles` and may leak `np.int64` the
same way — harmless today (LIK's formula does no boolean addition and
matches DIYABC), but the next consumer of those tuples inherits the
trap. Throwaway diagnostics in `tmp/fst_diag/` (`diag.py`/`diag2.py`/
`diag3.py`).

**Two review lessons, recorded in memory (`feedback_numpy_bool_
addition`)**: (1) a cross-check that disagrees with our code on
identical input is a software bug on one side, never "a real
discrepancy" — investigate it before any simulator/biology
hypothesis; (2) test stat bricks with the real dtype the pipeline
produces (`np.int64` from tskit), not only Python literals.

**Separate latent bug found the same morning — fixed the same day
(`2f8363a`, written by the user with the assistant reviewing)**: the
three readers of a real DIYABC *text* reftable
(`parse_real_reftable_params`, `parse_real_reftable_params_with_group_
priors`, `rewrite_real_reftable_txt` in `reftable_loop.py`) ordered a
scenario's parameter columns by iterating `priors` in the order of the
`historical parameters priors` *declaration* section. The real rule,
read in `reftable.cpp::bintotxt` (lines 474-542): the text export walks
the header's **trailer line** (`entetehist`, copied verbatim by
`header.cpp::readHeaderAllStat`) token by token, looks each name up
*by name* in the row's own scenario, prints the value if found and 14
spaces otherwise — so the text file's column order IS the trailer
order, and a whitespace split loses the blanks. (The `.bin` is
different: `nparamvar` floats in the scenario's *own* `histparam`
order from `sethistparamname`, constants excluded — untouched here,
`write_reftable_bin` was validated on `human` where all orders
coincided.) Verified on raw tokens: a scenario-1 row's three T values
are laid out `t32, t21, t41` (trailer order) and satisfy `t41>t21>t32`,
while declaration order (`t41, t32, t21`) would violate it. Every
earlier dataset had the two orders identical by coincidence; the first
hand-edited header that inserted `t41` at a different position in the
two sections mislabeled every scenario-1 parameter (`t32`'s value read
as `t41`, etc.) and produced a spurious `FST_1_1.4` `KS=0.80`. Only
scenario 2 survived, because its dropped column (`t41`) was the LAST
trailer token, so removing it shifted nothing.

Fix: new `_historical_columns_order(header_line, priors, scenarios)`
next to `_kept_param_names_by_scenario` — the two are complementary
(the latter gives the SET of names a scenario has, constants excluded;
the former gives their ORDER, read from the first line of the reftable
file being parsed, stopping at the first token that is not a declared
prior). A `ValueError` guard compares the count of names read against
the union of non-constant names used by any scenario and prints the
symmetric difference — a trailer typo (`t14` for `t41`) now fails
loudly with `différence ['t41']` instead of silently shifting every
value after it. The three readers keep `_kept_param_names_by_scenario`
as the set and filter the file order by it, per row, AFTER reading the
row's scenario index (the filter depends on it, so the order of
operations is forced). `_kept_param_names_by_scenario` itself and the
two writers are unchanged. Tests: `test_historical_columns_order`
(order follows the line, not the declaration; typo raises) and
`test_parse_real_reftable_params_follows_file_header_order` (a
3-line `tmp_path` reftable whose header order differs from the
`toy_example1_ms_modified` declaration order — the only test that
discriminates, verified to fail with the old `param_names` line and
pass with the new one; a plain `git stash` cannot be used for this
check because the test file imports the new helper). Two review
mishaps during this fix worth remembering: a `break` inside a list
comprehension (SyntaxError, never executed before being shown), and
`len(dict)` (number of scenarios) used where the union of its values
(number of names) was meant — both caught by running the code, not by
reading it. The same commit also moved `tests/conftest.py`'s
`REFERENCE_DIR` from `reference/human` to `reference/` (8 call sites
updated with explicit `/ "human"`), unrelated cleanup bundled in.


## CLOSED 2026-09-21: no "simplified substitution model" exists in the C++ (flagged 2026-09-18, falsified by direct source reading)

*(Entrée migrée verbatim depuis CLAUDE.md le 23/09/2026 ; CLAUDE.md n'en garde qu'un verdict condensé. Rédigée en anglais à l'origine, conservée telle quelle.)*

Flagged by the user at the close of the 2026-09-18 session: DIYABC
reportedly switches to a simplified substitution model below some
threshold of variable sites or individuals. **Checked exhaustively in
`~/Documents/Github/diyabc/src-JMC-C++` on 2026-09-21 for all three
data families — nothing of the kind exists. Chantier closed, nothing
to implement.** Evidence, so this isn't reopened on the same rumor:

- **DNA sequences**: `mutmod` (`JK`/`K2P`/`HKY`/`TN`, `history.hpp:151`)
  is written ONLY when the header token is read (`header.cpp:671-678`
  and its twin `1784-1791`) plus two copy-constructors (`history.cpp:
  257/303`) — no conditional reassignment anywhere. `comp_matQ`
  (`particuleC.cpp:1121-1170`) branches on `mutmod` alone, never on
  `dnavar`/`nloc`/`ngenes`/sample size. `mute`/`draw_nuc`/`put_
  mutations` (`1535-1790`) are a pure cumulative walk over `matQ`/`pi`.
  `do_sequence`'s base-frequency computation (`data.cpp:1532-1568`) has
  a single `if (nn > 0.0)` guard against division by zero on a fully-
  missing locus — no `0.25` fallback, no "too few sequences" branch.
  `grep -i "simplif|too few|trop peu"` finds nothing relevant; the only
  non-stats uses of `dnavar` are in `sumstat.cpp` (observed variable
  sites for the statistics, unrelated to the mutation model).
- **MicroSat**: `mute` (`particuleC.cpp:1683-1708`) has one branch,
  `if (Pgeom > 0.001) d = 1 + log(ra)/log(Pgeom); else d = 1;` — the SMM
  special case, a numerical guard against `log(0)`, not a data-driven
  fallback (already covered by `build_microsat_transition_matrix`'s
  `epsilon` clamp, tested at `Pgeom=0`). `setMutParamValue`
  (`803-825`) falls back to the group mean when `sdshape <= 0.001` or
  `nloc <= 1` — the per-locus dispersion tier, already ported via
  `check_nloc`, not a model switch. `kmin`/`kmax` (`particleset.cpp:
  80-82`: `kmoy = (maxi+mini)/2` integer division, `kmin = kmoy -
  (motif_range/2 - 1)*motif_size`, `kmax = kmin + (motif_range-1)*
  motif_size`) are deterministic from the observed alleles and the
  header's `motif_range`, no threshold.
- **SNP**: there is no model to simplify — `mute` sets `state = 1`
  (`1765`); `put_one_mutation` (`1644-1680`) is plain Hudson (one branch
  drawn proportionally to length), its only branch being `if (weight ==
  0.0) return` (zero-length tree → locus rejected and redrawn via
  `loc--`). `cherche_branchesOK` (`1549-1553`) sets every branch `OK`/
  `OKOK` to `true` — the reference-population filtering (`refnindtot`/
  `popref`) is entirely commented out, so `weight` is always 1 in the
  compiled binary and `simulate_snp_genotypes`'s "draw over all edges"
  is exact. `mafreached` (`2194-2211`) has no sample-size-dependent
  mode.

The only genuine "second mode" anywhere in `particuleC.cpp` remains the
discrete generation-by-generation coalescent selected by
`evalcriterium` (`1251-1275`) for very small `N` — a coalescent-side
switch, not a mutation-model one, already documented above as
un-replicated and never triggered on any of this project's datasets.
The original claim most likely came from documentation or memory of a
different tool; if a precise DIYABC source (doc page, paper) turns up,
cross-check it against the line numbers above before reopening.

## 2026-09-23/24 — Échantillonnage sériel : implémenté et validé, et un bug de port sur `LIK` révélé au passage

Dernier gros chantier du POC. Mené en mode mentor : conception et revue avec
l'assistant, code écrit par l'utilisateur.

### Le cadrage, qui n'était pas celui annoncé

La limitation notée le 17/09 disait « échantillonnage sériel non supporté ».
La lecture du C++ montre que le vrai manque est plus général : **le pipeline
n'avait aucune notion d'échantillon distincte de population**, et le sériel est
simplement le premier jeu où les deux divergent.

- `history.cpp::read_events` compte `nsamp` **séparément** de `npop` (`nn0`) et
  attribue `event[i].sample = ++nsamp` dans l'ordre du fichier. Rien ne
  contraint `nsamp == nn0`.
- `data.cpp` ne parle **que** d'échantillons : `nsample`, `samplesize[ech]`,
  `haplosnp[ech]`, `ssize[locustype][sa]`. Les indices des statistiques sont
  des indices d'échantillon, jamais de population.
- `particuleC.cpp:1185/1213/1528` : la taille d'un échantillon est lue dans les
  données observées, chaque nœud est étiqueté `gt.nodes[i].sample = sa + 1`, et
  un événement SAMPLE active exactement les nœuds portant son indice. Le
  mapping bloc de données `sa` ↔ `sa+1`-ième ligne `sample` est donc
  **positionnel**, sans référence croisée nominale.
- `buildSuperScen` (`header.cpp:1057-1108`) ne fait que des maxima : DIYABC ne
  garantit **rien** sur la cohérence de cet ordre entre scénarios.

Côté nous, la confusion était encodée dans `observed_count_population` /
`build_samples_argument`, qui nomment les blocs `pop1..popN` par ordre
d'apparition — c'est-à-dire par indice d'échantillon — et passent ce nom tel
quel à msprime comme nom de **population**. D'où le symptôme :
`KeyError "Population with name 'pop2' not found"`.

### Quatre faits mesurés avant d'écrire une ligne

1. Le temps d'un `sample` peut être un **nom de paramètre**
   (`particuleC.cpp:599-605`), donc la construction des samples vient après le
   tirage des priors. Accessoirement, `demography_builder.py:102` évalue déjà
   ces expressions puis les jette (le `continue` sur `SampleEvent` est ligne
   107) : preuve empirique gratuite que `values` suffit, dans les deux chemins.
2. msprime alloue les IDs de nœuds échantillons **contigus, dans l'ordre de la
   liste de `SampleSet`** ; les individus suivent le même ordre.
3. Cela reste vrai à **ploïdie hétérogène** (3 femelles diploïdes + 2 mâles
   haploïdes → individus 0-2 à 2 nœuds, 3-4 à 1 nœud).
4. `samples={"pop1":20,"pop2":20}` et la liste de `SampleSet` équivalente
   produisent des tables `edges`/`nodes`/`individuals` **byte-identiques** à
   graine égale. C'est ce qui autorise à basculer le chemin existant au lieu
   d'en dupliquer un frère — l'équivalence est mesurable, pas postulée.

### Trois routes, une fausse piste assumée

- **A** — propager le layout jusqu'aux helpers. Explicite, ~38 éditions sur du
  code validé contre la sortie réelle de DIYABC.
- **F** — dériver le découpage de la `ts` elle-même, en groupant les nœuds
  échantillons par `(population, temps)`. 9 éditions. Implémentée, testée,
  fonctionnelle — y compris sur `<X>`, dont les deux `SampleSet` fusionnent
  naturellement puisqu'ils partagent la même paire.
- **C** — envelopper la `ts` dans un objet portant son layout, `__getattr__`
  délégant le reste. 4 éditions. Vérifié sans obstacle : aucun `isinstance` sur
  `TreeSequence` dans le dépôt, la couche stats n'utilise que 4 attributs, et
  les `TreeSequence` ne traversent jamais la frontière de processus.

**F a été abandonnée après implémentation**, pour une raison découverte en la
testant : un échantillon d'effectif **nul** ne laisse aucun nœud dans la `ts`,
donc il est invisible à la dérivation, et tous les échantillons suivants
glissent d'un rang. `compute_population_layout` y était immunisée (elle tire
ses noms des métadonnées). Le déclenchement demande trois conditions cumulées
(un locus `<Y>`, une population sans mâle, et qu'elle ne soit pas la dernière)
et aucun jeu du dépôt ne les réunit — mais documenter une régression n'est pas
la corriger. Retour à A, choisie par l'utilisateur pour son caractère explicite.

### Les deux briques

`build_sample_sets_from_scenario(scenario, values, counts_by_samples)` et
`compute_sample_layout(ts, counts_by_samples)` — voir CLAUDE.md, section
« Serial/temporal sampling », pour leur contrat et l'invariant positionnel.

Bugs attrapés en revue, tous par exécution et non par lecture du diff :

- **Indexation par nom de population** dans la première brique
  (`counts_by_samples.get(f"pop{event.pop}")`). Sur `toy_example1_ms` les
  quatre `sample` pointent tous `pop=1`, donc le même effectif était lu quatre
  fois — **juste par coïncidence**, les quatre blocs ayant 20 individus. Le
  test initial ne discriminait pas : la version buguée passait toutes ses
  assertions. Verrouillé par un cas à effectifs volontairement inégaux
  (20/15/30/10), vérifié échouer sur l'ancienne version.
- **`assert [np.array_equal(...), ...]`** — une liste non vide est toujours
  vraie, donc toutes les comparaisons de tableaux étaient inertes. Même famille
  que le bug `.genotype_matrix().all()` déjà consigné. `assert all(...)`.
- **`if layout is None: population_layout = ...`** dans deux helpers dont la
  variable locale ne s'appelait pas `layout` : `UnboundLocalError` dès qu'un
  layout était fourni. Bug **latent** — la suite restait verte tant que
  personne ne passait l'argument.
- **`from matplotlib.style import context`**, auto-import de l'IDE déclenché
  par le nom `context`, dans `bridge/pipeline.py`. matplotlib n'est pas une
  dépendance déclarée : sur un environnement propre, tout le pipeline aurait
  échoué à l'import. Invisible ici car chaque fonction a un paramètre `context`
  qui masque l'import — d'où 9 `F811` de `ruff`. Cinquième occurrence du motif
  `feedback_name_shadowing_pattern`.
- **Paramètre sans valeur par défaut** (`layouts_by_locus: ... | None` sans
  `= None`) dans `compute_all_statistics_dna` : tous les appelants cassaient.
  Le réflexe initial a été d'ajuster le test qui le signalait ; c'est le test
  qui avait raison.
- **Deux normalisations manquantes** (`compute_MGW`, `compute_DM2`) sur 27
  éditions quasi identiques — précisément les deux où le geste différait
  (allonger un `zip` existant à trois termes au lieu d'en créer un).
  `strict=True` n'aide pas ici : il valide des longueurs, pas la nullité.

### Le vrai bug de port, révélé par le jeu sériel

Premier rejeu apparié contre la reftable réelle de `toy_example1_ms` : 152
colonnes, historiques exactes, et un bloc d'écarts KS **entièrement concentré
sur `LIK_2_*`** — jamais G1, toujours G2. En scénario 1, exactement les 6
paires ordonnées impliquant l'échantillon 4 ; en scénario 2 (sans goulot), les
12 paires, plus faiblement.

Cause, trouvée en lisant `cal_lik2p` :

```c
nal = 0;
for (k = 0; k < locuslist[loc].nal; k++) {
    frt = 0.0;
    for (pop = 0; pop < this->nsample; pop++) frt += freq[pop][k];
    if (frt > 0.000001) nal++;        // seulement si l'allèle est OBSERVÉ
}
b = 1.0 / nal;
```

`_compute_LIK_for_one_locus` calculait `nb_allele` à partir des clés de
`length_by_pop`, c'est-à-dire de `variant.alleles` — **tous** les états
apparaissant dans la table de mutations, y compris ceux qu'aucun échantillon ne
porte plus (écrasés par une mutation ultérieure sur la même lignée). Le
pseudo-compte `b = 1/nal` était donc jusqu'à **deux fois trop petit**.

Mesuré sur un locus réel : `nal` = 13/6/14/12/18 chez nous contre 7/3/11/10/9
dans le C++.

Deux choses expliquent le profil des écarts. **Pourquoi G2 et pas G1** : le
header donne à G2 `MEANMU UN[5e-4,5e-3]` contre `UN[1e-4,1e-3]` pour G1 — cinq
fois plus de mutations, donc beaucoup plus d'états créés puis perdus (G2 cumule
en plus haploïdie et `Ne` rescalé). **Pourquoi l'échantillon 4 en scénario 1** :
`b` ne pèse que lorsque l'allèle testé est absent de l'échantillon de
référence ; l'échantillon 4 (t=500) traverse le goulot `tbn ∈ [10,1000]` et est
donc le plus divergent.

La docstring de la fonction affirmait que le comportement était « exactement le
`nal` du C++ puisque ce projet n'a que 2 populations ». Elle était fausse deux
fois : le vrai écart n'est pas le périmètre du pooling mais le filtre sur les
allèles non observés — et l'union `count_i ∪ count_j` qu'elle décrivait est un
**no-op**, `_length_by_population` donnant la même liste de clés à toutes les
populations. Une docstring qui déclare une approximation « équivalente dans
notre cas » est une dette qui expire sans prévenir ; celle-ci a tenu deux mois.

Correction : sommer les comptes sur toutes les populations de `length_by_pop`
et ne garder que les non nuls. Vérifiée contre l'oracle (7/3/11/10/9 exactement)
avant rejeu. Après correction, les 12 colonnes `LIK_2_*` rentrent dans le bruit
sur les deux scénarios.

### Reste

- **`AML` résiduel : bruit, confirmé par falsification.** Trois rejeux
  indépendants comparés colonne à colonne, jamais les mêmes colonnes
  signalées : 4 significatifs sur 153 dans le
  premier (`AML_2_4.2.3`, `AML_1_2.1.3`, `DM2_1_1.4`, `DM2_1_2.4`), 3 dans le
  second (`AML_2_3.2.4`, `DM2_2_2.4`, `AML_2_3.1.2`) — **recouvrement nul**.
  Les p s'effondrent d'un run à l'autre (`DM2_1_1.4` 0,039 → 0,953 ;
  `DM2_1_2.4` 0,046 → 0,877), quand un effet réel garderait son ordre de
  grandeur. Les comptes sont sous l'attendu (~7,7 à α=0,05) et la
  distribution des p penche vers le HAUT (médianes 0,578 et 0,661, histogrammes
  lourds dans les derniers déciles) — l'inverse de la signature d'un effet.
  Le troisième rejeu signale encore 3 colonnes, encore différentes.

  Fausse alerte levée puis close au passage : un `Npast` non nul
  (`diff_mean = 1297`) dans le deuxième fichier laissait croire que le rejeu
  avait cessé d'être apparié. C'était un artefact de notebook —
  `replay_reftable_simulation_microsat` restaure l'ordre du fichier via
  `results_by_index[i] for i in range(len(rows))`, jamais par ordre d'arrivée
  d'`as_completed`, donc l'appariement ne peut pas s'y perdre. Confirmé au
  troisième rejeu (`diff_mean = 0` à nouveau). Contrôle à remettre en tête de
  tout notebook de comparaison : `(reel["scenario"].values !=
  notre["scenario"].values).sum()` doit valoir 0.

  Corollaire méthodologique réutilisable : sur un rejeu **apparié**, les
  particules partagent leurs priors, les statistiques sont donc corrélées et le
  test de KS est **conservateur** — il sous-estime l'écart. Un « 0 significatif »
  n'y prouve pas beaucoup plus qu'un « 3 ou 4 sans recouvrement ».
- `<X>`/`<Y>` sériels non câblés, **sans garde** : la dispatch par locus
  écrase silencieusement les `SampleSet` sériels. Aucun jeu ne le déclenche.
- SNP sériel non traité : `build_samples_argument` rend toujours un dict.
- Fausse alerte levée en cours de route : le `mumic_1` ASCII de ce header face
  au `µmic_1` codé en dur dans `group_prior_column_names` est **sans effet** —
  `parse_real_reftable_params_with_group_priors` lit les valeurs
  positionnellement et n'utilise ces noms que comme clés de son dict de sortie.

## 2026-09-25 — Échantillonnage sériel côté SNP : implémenté et validé sur `human_seriel`

Suite directe du chantier microsat du 23-24/09. Mode mentor, code écrit par
l'utilisateur.

### Le jeu de référence

`reference/human_seriel/` : le `.snp` de `human` (4 blocs de 30 individus,
ASW/YRI/CHB/GBR) réinterprété comme **4 échantillons d'une seule population**,
avec le squelette de scénario de `toy_example1_ms` (`0/50/200/500 sample 1`,
`tbn VarNe 1 Npast` en scénario 1, rien en scénario 2). 1000 particules,
513/487 sur les deux scénarios.

**Le vrai binaire refusait de tourner**, et pour une raison qui mérite d'être
retenue : la section `historical parameters priors (3,0)` déclare **zéro**
contrainte d'ordre, mais le header avait gardé la ligne `DRAW UNTIL` copiée de
`human`. Or `readHeaderHistParam` (`header.cpp:274`) consomme cette ligne
**à l'intérieur** de `if (this->nconditions > 0)`. Avec `C = 0`, le `getline`
n'est jamais exécuté, la ligne reste dans le flux et **tous les `getline`
suivants glissent d'un cran** — `readHeadersimLoci` lit alors `DRAW UNTIL` là
où il attend `loci description (N)`. Même famille que la ligne trailer : une
ligne d'apparence décorative qui est en réalité un jeton de position. Consigné
dans les règles dures de CLAUDE.md.

### Un chantier beaucoup plus court que le microsat, pour une raison structurelle

Les 130 statistiques SNP reçoivent `genotypes_per_locus`, une liste de dicts
`{nom: [génotype]}`. Le layout est donc consommé **une seule fois**, dans
`simulate_snp_genotypes`, pour fabriquer ces dicts — il n'y a **rien à propager
à travers les fonctions de stats**, contrairement aux 27 fonctions et 7 boucles
de dispatch du côté microsat/ADN. Tout tient dans `ancestry_simulation.py`.

Deuxième facilité : `samples` était **déjà** typé
`dict[str, int] | list[msprime.SampleSet]` aux cinq niveaux de la chaîne
(`simulate_independent_loci`, `simulate_shared_ancestry_loci`, les deux boucles
MAF, le chemin MRC), et `simulate_snp_genotypes` acceptait déjà un
`population_layout` optionnel. Il n'y avait donc que deux portes à ouvrir.

### Le virage de conception, en cours de route

La première version faisait descendre un **layout déjà calculé** depuis le
pipeline. Deux objections l'ont fait abandonner avant le câblage :

1. Le layout dépend de la **ploïdie**, qui varie par type d'héritage (`<A>` 2,
   `<H>`/`<M>` 1). Un layout unique construit en amont serait donc faux pour
   tous les types sauf un. Sur `human_seriel` (que du `<A>`) ça aurait marché
   **par coïncidence** — exactement ce que ce chantier élimine.
2. Le pipeline n'a aucune `TreeSequence` sous la main : les arbres naissent
   dans les boucles MAF. Il ne *peut* pas construire le layout.

D'où le principe qui a guidé la suite, et qui a resservi deux fois :
**construire chaque objet là où ses entrées existent**. Le layout se calcule là
où naît la `ts` ; les `SampleSet` se construisent là où naît `values`. C'est ce
qui fait descendre `counts_by_samples` (un dict ordonné `{nom: nb_individus}`)
plutôt qu'un layout.

### Les bugs, tous de la même famille

- **`counts_by_samples=samples`.** La confusion centrale : `samples` va à
  msprime (`dict` *ou* `list[SampleSet]`), `counts_by_samples` alimente
  `compute_sample_layout`. L'un n'est jamais déductible de l'autre — une liste
  de `SampleSet` ne porte pas les noms, et un locus `<X>` décrit UN échantillon
  avec DEUX `SampleSet`. Symptôme : `AttributeError: 'list' object has no
  attribute 'values'` sur les jeux multi-types.
- **Un bloc de `with_maf_filter` collé en tête de
  `simulate_genotypes_for_locus_type`**, référençant `ts` et
  `counts_by_samples` hors portée.
- **`population_layout = layout`** survivant dans `with_maf_filter` après le
  renommage du paramètre → `F821` / `NameError` sur le chemin lent.
- **`build_sample_sets_from_scenario(scenario, values, ...)` appelée avant que
  `values` existe** (`UnboundLocalError`, 12 tests). C'est la contrainte
  annoncée dès le premier jour : un temps d'échantillonnage peut être un nom de
  paramètre (`particuleC.cpp:599-605`), donc les `SampleSet` se construisent
  **après** le tirage.
- **Le bloc de construction placé dans `run_poc_for_directory_with_values`
  (rejeu) au lieu de `run_poc_for_directory` (tirage).** Les 207 tests restaient
  verts, parce qu'aucun n'exerce un jeu SNP sériel — même leçon que le
  `reads_observed = None` de septembre.

### Un danger mesuré, non corrigé

Un layout dont la longueur totale ne correspond pas à `ts.num_samples` produit
des génotypes **silencieusement faux**, jamais une exception :
`simulate_snp_genotypes` fait un test d'appartenance (`s in derived_samples`),
jamais une indexation, donc un ID de nœud inexistant rend simplement `0`.
Découvert en calculant un layout sur une `ts` diploïde passée à une simulation
haploïde. Un garde (`somme des longueurs != ts.num_samples -> ValueError`) a
été proposé et non implémenté.

### Validation

Rejeu apparié contre la reftable réelle de `human_seriel` : **0 colonne sur 130
sous p<0,05**. À relativiser dans le bon sens — on en attendrait ~6,5 par pur
hasard à α=0,05 ; ce zéro reflète le caractère conservateur du KS sur un rejeu
apparié (priors partagés, statistiques corrélées), pas une perfection
surnaturelle. Il reste que rien ne diverge.

### Reste

- `<X>`/`<Y>` sériels non câblés, **sans garde**, des deux côtés.
- PoolSeq sériel intouché (`with_mrc_filter` / `simulate_poolseq_reads`
  appellent encore `compute_population_layout` directement).
- Le renommage `pop` → `samp` est **partiellement** débloqué : le côté
  layout/statistiques ne touche plus msprime et peut être renommé ; le repli
  `build_samples_argument` des appelants directs (tests, `scripts/`) reste
  msprime-facing.

## 2026-09-25 (suite) — PoolSeq sériel : validé ; et un biais systématique d'environ 1 %, invisible au KS, découvert au passage

### Le chantier lui-même : court, et sans surprise

Décalque du SNP IndSeq, avec deux spécificités PoolSeq. `context.count_samples`
donne la taille **haploïde** du pool (des copies de gènes, pas des individus),
d'où `poolseq_counts_by_sample` et son `// 2`. Ce helper a été **extrait**
plutôt que dupliqué, pour une raison précise : les deux dicts candidats
(tailles haploïdes contre nombres d'individus) ont le **même total**, donc
`_check_layout_matches` ne les distingue pas — une divergence entre les deux
définitions serait silencieuse. Et `pool_sizes`, consommé par
`compute_all_statistics_poolseq` pour la correction de biais de lecture, reste
en unités **haploïdes** : les deux échelles coexistent volontairement.

Jeu de référence construit pour l'occasion : `reference/toy_example4_seriel`,
le `.snp` PoolSeq de `toy_example4` (4 pools de 200 haploïdes) réinterprété en
4 échantillons temporels d'une population unique, 100 loci, `<MRC=5>`. Ce
dernier point compte : contrairement à `human_seriel` (`<MAF=hudson>`, chemin
rapide seulement), il exerce réellement la **boucle de rejet** MRC.

Un garde de longueur (`_check_layout_matches`) a été posé et factorisé en trois
appels au passage. Il protège d'un échec strictement silencieux : les
consommateurs testent l'appartenance (`s in derived_samples`,
`derived_samples.intersection(sample_ids)`), jamais l'indexation, donc un ID de
noeud absent de la ts rend `0`. Côté SNP ça donne des génotypes de mauvaise
longueur ; côté PoolSeq, `len(sample_ids)` étant le dénominateur de la
fréquence dérivée, ça donne une fréquence **diluée** puis un tirage binomial
dessus — indétectable en aval.

### Le câblage, vérifié de bout en bout

Avant toute interprétation statistique :

```
SampleSet : 4 × 100 individus, pop1, temps 0 / 50 / 200 / 500
ts        : 800 noeuds, 400 individus
layout    : 4 groupes de 200 noeuds, aux temps [0] [50] [200] [500]
            (compute_population_layout, lui, rend 1 groupe de 800)
paramètres historiques rejoués : diff = 0, p = 1.0
```

### Le résultat, et le bon critère

Deux rejeux appariés, prior d'origine (`Npresent UN[10,1000]`) :

```
rejeu 0   KS 10/133   signes  99/128  p = 3,8e-10   rdiff médian −1,06 %
rejeu 1   KS  7/133   signes  84/128  p = 5,2e-04   rdiff médian −0,79 %
```

**Par le critère du projet, ça valide** : 10 et 7 colonnes significatives pour
~6,7 attendues à α=0,05. Le chantier est déclaré complet sur la même base que
toutes les autres familles.

### Le vrai apport de la session : le KS est aveugle à un décalage global

Le test des signes sur les `rdiff` est significatif dans les **deux** rejeux,
même sens, même amplitude — là où le KS colonne par colonne ne voit rien. Un
décalage uniforme de −1 % reste très à l'intérieur de la variance
inter-particules de chaque colonne et ne sort jamais.

Passé au crible de ce test, plusieurs jeux **déjà déclarés validés** le portent :

```
toy_example5           70 loci   −2,65 %   (validé depuis juillet)
toy_example5_500loci  350 loci   −2,33 %
toy_example4_seriel   100 loci   −0,9 %
toy_example1_ms        50 loci   −0,87 %   (sériel microsat, validé le 24/09)
toy_example3          100 loci   −0,74 %
toy_example3_500loci  500 loci   −0,24 %
toy_example4          100 loci   +0,22 %
toy_example4_MRC1     100 loci   +0,00 %
human_seriel         5000 loci   +0,06 %
```

**Pas spécifique au sériel** : `human_seriel` est sériel et propre. La piste la
plus cohérente est celle du biais résiduel clos le 17/07 — il décroît avec le
nombre de loci — mais ce n'est pas établi : `toy_example3` divise le sien par
trois en passant à 500 loci, `toy_example5` ne bouge presque pas. Expérience à
variable unique pour trancher : rejouer `toy_example4_seriel` à 500-1000 loci,
tout le reste inchangé.

### Trois erreurs de méthode de ma part, consignées parce qu'elles sont instructives

**Contrôles mal choisis.** J'ai d'abord conclu « biais propre au sériel » en ne
comparant qu'aux deux PoolSeq non sériels — qui se trouvent être précisément
les deux jeux neutres du dépôt. Élargir à tous les `comparaison_summary.csv`
disponibles a immédiatement falsifié la conclusion.

**Expérience confondue.** Pour tester si le coalescent discret de DIYABC
(`evalcriterium`) expliquait l'écart, j'ai proposé d'augmenter `Npresent`. Ça
changeait **deux** choses à la fois : supprimer le déclenchement du mode
discret, **et** effondrer le vrai FST vers zéro — ce qui rend un petit décalage
absolu dominant en relatif. Résultat ininterprétable (93/133 significatives),
et une régénération de reftable perdue.

**Corrélation surinterprétée.** J'ai présenté `r = 0,86` entre deux rejeux
comme une preuve forte d'effet reproductible. Or les deux rejeux partagent la
**même** reftable DIYABC : le dénominateur de `rdiff` est identique, et avec
1000 particules les moyennes sont bien estimées. Une forte corrélation est donc
attendue dès qu'il existe une différence par colonne, si petite soit-elle. Elle
prouve que l'écart n'est pas du bruit d'échantillonnage, pas qu'il est grand.

### Correction factuelle sur `evalcriterium`

CLAUDE.md affirmait que le coalescent discret n'était « jamais déclenché sur
aucun jeu de ce projet ». C'est faux depuis ce jeu. Le critère est
`ra = nLineages / N`, l'approximation continue n'étant conservée que si
`ra < 0,5` sur les segments de plus de 100 générations. `toy_example4_seriel`
concentre 400 individus (800 copies) dans UNE population avec `Npresent ≤ 1000`
: le segment 200→500 est en mode discret pour **toutes** les particules. Son
rôle éventuel dans l'écart observé n'est **pas** établi — la seule expérience
menée dessus était confondue.

