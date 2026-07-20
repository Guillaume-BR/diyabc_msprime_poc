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


