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
les 100 premiers loci (indices 0-99), répartis comme : positions 0-69
=<A>, 70-79=<X>, 80-89=<M>, 90-99=<Y>.

IMPORTANT : le FICHIER .snp peut contenir bien plus de loci (51250 pour
human) que ce qui est réellement simulé/comparé (5000 pour le scénario
de header.txt) -- "loci description" est un FILTRE/SOUS-ÉCHANTILLONNAGE
des colonnes du fichier de données, pas une description de tout le fichier.

## Modèle de mutation SNP correct (doc DIYABC section 2.4.3) : algorithme de Hudson

Confirmé par la doc utilisateur DIYABC : pour les SNP, "il est supposé
qu'il y a eu une et une seule mutation dans l'arbre de coalescence" --
PAS un processus de Poisson à taux variable. C'est l'algorithme "-s"
de Hudson (2002). Notre précédente approche (msprime.sim_mutations à
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