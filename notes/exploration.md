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