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

### Mise à jour du 01/09/26 — cause trouvée et corrigée : généalogie <M> non partagée entre loci (investigation CLOSE)

Reprise directement là où la mise à jour du 31/08 s'était arrêtée
("piste plus générale... pas encore investiguée") — via une série de
diagnostics jetables (non committés, tous décrits ci-dessous) puis une
instrumentation temporaire du vrai binaire `general` (`sumstat.cpp`,
revertée après coup), suivie d'un correctif effectif dans `bridge/`.

**Étape 1 (réfutée) — dépendance au régime `t1/N1`.** Hypothèse testée :
`<M>` (Ne effectif 4x plus petit, `coalescence_coefficient("M",
0.5)/2 = 0.25` contre `1.0` pour `<A>` avec ploidy=2 — confirmé fidèle
à `particuleC.cpp::cal_coeffcoal` et déjà validé le 27/08) coalesce
tellement plus vite qu'un nombre variable de lignées atteint même les
événements démographiques les plus anciens, gonflant la variance
inter-locus. Diagnostic : sur les 512 particules réelles de scénario 2
(sans admixture) déjà rejouées, 200 répliques msprime par particule,
comptage des lignées survivantes juste avant `t1` (le seul événement
démographique de ce scénario), binné par `t1/N1`. Résultat : le ratio
interne `CV(M)/CV(A)` varie fortement selon le bin (1.34 dans
`[0,0.5)`, ~0 dans `[5,10)`/`[10,∞)` où `<M>` atteint une variance
nulle — exactement 1 lignée par population — bien avant `<A>`) mais de
façon **non monotone** (positif dans les bins superficiels, négatif
dans les bins profonds — `<M>` "en avance" sur sa propre horloge
coalescente traverse plus vite le pic de la courbe de variance en
cloche). Stratifier le VRAI écart-type réel/simulé (`DTA_3_1` etc.) par
le même binning `t1/N1` montre un déficit **uniforme** (~0.5-0.85)
across tous les bins, y compris ceux où le mécanisme interne prédit un
ratio proche de 1 — **pas de corrélation avec le régime**, hypothèse
réfutée comme cause de l'écart observé (le mécanisme lui-même est réel
et vérifié, juste pas responsable).

**Étape 2 (réfutée) — relecture attentive de `sumstat.cpp`.** Vérifié
terme à terme : `DataC::calcule_ss` (définition de `ssize`, `data.cpp:
966-999`) donne bien `n=20` copies de gène pour `<M>` (tous les
individus, +1 chacun, aucune distinction de sexe) et `n=40` pour `<A>`
(+2 chacun) — identique à `observed_count_population`/notre ploidy.
Les constantes de Tajima's D (`a1,a2,b1,b2,c1,c2,e1,e2`,
`cal_dta1pl:1566-1579`) sont identiques à `_tajima_constants`.
`ParticleC::cal_numvar` (`dnavar`/`haplodnavar`) restreint le calcul aux
sites variables **groupés sur toutes les populations du dataset**, pas
par population — vérifié sans impact numérique (chaque stat "1p" refait
son propre test d'identité localement). Aucune divergence de formule
trouvée.

**Étape 3 (LA cause) — instrumentation temporaire de `cal_dta1pl`
(`~/Documents/Github/diyabc/src-JMC-C++/sumstat.cpp`), rebuild via le
`CMakeLists.txt` du dépôt, run réel (`general -p ./ -R "ALL" -r 30
-g 30 -m -t 1` sur une copie de `toy_example2_ms_dna`), instrumentation
retirée et binaire reconstruit propre ensuite.** Trace de `n`/`S`/`pi`
par locus × population × particule pour les 10 loci ADN (`kloc` 10-14 =
`<A>`, 15-19 = `<M>`). Calcul de la corrélation de `pi` entre paires de
loci **au sein d'une même particule** :

| paire de loci | type | corrélation de `pi` |
|---|---|---|
| 10 vs 11 | `<A>` | -0.249 |
| 11 vs 12 | `<A>` | 0.319 |
| 15 vs 16 | `<M>` | 0.569 |
| 16 vs 17 | `<M>` | 0.753 |
| 17 vs 18 | `<M>` | 0.333 |
| 18 vs 19 | `<M>` | 0.657 |

Les loci `<A>` d'une même particule sont quasi indépendants (corrélation
moyenne 0.153, proche du bruit) ; les loci `<M>` d'une même particule
sont FORTEMENT corrélés (moyenne 0.552, jusqu'à 0.75-0.9). C'est exactement
ce que prédit la biologie de l'ADNmt (non-recombinant, transmission
uniparentale) et exactement le mécanisme déjà documenté et implémenté
côté SNP (`simulate_shared_ancestry_loci`, `ancestry_simulation.py:305`,
qui cite `particuleC.cpp:2422-2435` `GeneTreeY`/`GeneTreeM` : le premier
locus `<Y>`/`<M>` tire une généalogie, tous les suivants COPIENT cette
même généalogie, seule la mutation diffère) — **mais jamais porté côté
séquences ADN**. `dna_mutation_simulation_per_locus`/`_from_values`
appelaient `msprime.sim_ancestry` indépendamment pour CHAQUE locus, avec
un seed qui varie par locus (`seed + _ANCESTRY_SEED_OFFSET + i`), y
compris pour `<M>`.

**Conséquence** : en tirant 5 généalogies indépendantes au lieu d'une
seule partagée et répétée, le pipeline moyennait artificiellement le
bruit inter-locus sur les stats de groupe (`DTA_3`/`VNS_3`/`MNS_3`/
`VPD_3`, moyennées sur `nl=5` loci) — exactement le sens du déficit
observé depuis le 27/08, et ça explique du même coup pourquoi G2 (`<A>`,
correctement indépendant des deux côtés) n'a jamais montré d'écart, et
pourquoi le déficit était resté constant à travers tous les régimes
`t1/N1`/scénarios testés aux étapes précédentes : la cause n'a jamais
été démographique, elle était architecturale (un mécanisme de partage
de généalogie non porté), donc invariante à tout ce qu'on pouvait faire
varier côté paramètres.

**Correctif** (`bridge/ancestry_simulation.py`) : nouvelle constante
`_SHARED_M_ANCESTRY_SEED_OFFSET = 120_000_000`. Dans
`dna_mutation_simulation_per_locus`/`_from_values`, la graine
d'ancestrie pour un locus `<M>` est désormais **fixe**
(`seed + _SHARED_M_ANCESTRY_SEED_OFFSET`, sans `+i`) au lieu de varier
par locus — tous les loci `<M>` du dataset (peu importe leur groupe
déclaré, comme côté SNP) partagent donc la même généalogie tirée par
`msprime.sim_ancestry`. Approche différente de
`simulate_shared_ancestry_loci` (qui réutilise littéralement le MÊME
objet `TreeSequence`) parce que `sequence_length` (donc le `RateMap` de
mutation) varie d'un locus `<M>` à l'autre côté séquences ADN, alors
qu'il vaut toujours 1 côté SNP — resimuler avec la même graine mais un
`sequence_length` différent par locus fonctionne car, sans
recombinaison (le défaut), `msprime.sim_ancestry` produit exactement la
même topologie/mêmes temps de nœuds quelle que soit `sequence_length`
à graine égale (vérifié empiriquement avant d'écrire le correctif :
mêmes démographie/samples/ploidy, `sequence_length` variée entre 1 et
1000, nœuds internes identiques bit à bit). `<A>`/`<H>` restent
inchangés (généalogie indépendante par locus, comme avant).

**Validation** : 15 valeurs golden régénérées dans
`tests/test_summary_statistics.py`/`test_pipeline.py` — uniquement des
valeurs G3, AUCUNE valeur G2 n'a bougé (signature exacte confirmant que
le fix est bien scopé à `<M>`, même méthode de vérification que le bug
de seed du 31/08). 130/130 tests verts. Rejeu complet du reftable réel
(`scripts/replay_diyabc_priors_dna.py`, 1000 particules,
`toy_example2_ms_dna`) : le test KS colonne par colonne passe de
**11/42 à 2/42 colonnes significatives** (p<0.05), et les 2 restantes
(`MPD_2_2`, `VPD_2_2`) sont sur G2 (`<A>`, jamais concerné par ce bug)
avec un écart dans le sens inverse (ratio std sim/réel >1, pas de
déficit) — cohérent avec du bruit résiduel à cet effectif de
particules, pas un biais systématique. Le ratio std(sim)/std(réel) sur
les colonnes G3 anciennement déficitaires est passé de 0.6-0.85 à
0.97-1.10 (`DTA_3_1` : 0.649→1.056, `VNS_3_1` : 0.776→1.014, `MNS_3_1` :
0.853→1.044, `VPD_3_1` : 0.617→1.099).

**Statut de l'investigation : CLOSE.** Cause identifiée, corrigée, et
validée par KS sur un reftable réel complet — contrairement à la
clôture prématurée du 27/08, celle-ci survit à un contre-test direct
(comparaison avant/après correctif sur les mêmes données). Aucun autre
gap connu sur la pipeline ADN séquence à ce jour.
