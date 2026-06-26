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

Confirmé : human/header.txt déclare ses 5000 loci en <A> (autosomal),
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

## Modèle de mutation SNP de human : non élucidé précisément, simplifié pour le POC

Découverte importante : human/header.txt n'a AUCUNE section "group priors"
(pas de MEANMU/MEANSNI/MEANP), contrairement à sequences-mut. Le format
"5000 <A> G1 from 1" n'active aucune des branches [M]/[S]/[P] de
header.cpp (qui testent ss[2] contre ces marqueurs -- ici ss[2]="G1",
donc aucune branche ne matche). "from 1" est un texte libre/commentaire,
absent de tout parseur C++ (vérifié par grep sur tout le dépôt).

Conséquence : put_mutations() (particuleC.cpp) utiliserait mutrate =
mut_rate + sni_rate pour ce type de locus (<5), mais ces valeurs ne sont
jamais déclarées dans header.txt pour human -- donc soit des valeurs par
défaut codées en dur existent ailleurs dans le C++ (non trouvées), soit
human suit un algorithme de simulation SNP distinct du modèle de mutation
à taux fixe (un article tiers, arxiv 2501.17107, mentionne un "algorithme
de simulation SNP" spécifique à DIYABC-RF v1.0, sans en détailler le
mécanisme).

DÉCISION POC : reporté. On utilise un modèle de mutation msprime simplifié
et raisonnable (taux fixe choisi à la main), sans prétendre reproduire
l'algorithme exact de DIYABC pour les SNP. À creuser sérieusement avant
toute comparaison statistique fine avec le reftableRF.bin de référence.

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