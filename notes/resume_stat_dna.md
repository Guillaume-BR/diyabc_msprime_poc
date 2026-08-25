# Ce que mesurent les statistiques ADN

13 statistiques résumées calculées sur des séquences ADN (groupes `[S]`
du header, ex. `toy_example2_ms_dna`), dont 8 se calculent pour chaque
population séparément, et 5 pour chaque paire de populations. Elles
servent toutes à décrire, à partir d'un échantillon de séquences, à quel
point une population est génétiquement diverse en elle-même, et à quel
point deux populations sont différenciées l'une de l'autre.

## Statistiques d'UNE population (diversité interne)

**NSS — nombre de sites ségrégeants**
Sur toute la longueur de la séquence, combien de positions varient
d'un individu à l'autre dans l'échantillon ? Un site "ségrégeant" est
un site où on n'observe pas la même base partout. C'est la mesure la
plus brute du polymorphisme : plus il y a de sites qui varient, plus la
population est diverse à ce locus.

**NHA — nombre d'haplotypes distincts**
Combien de versions différentes de la séquence entière observe-t-on
dans l'échantillon ? Si tout le monde a exactement la même séquence :
1 haplotype. Si chaque individu a sa propre combinaison unique de
variants : autant d'haplotypes que d'individus. Contrairement à NSS
(qui compte site par site), NHA regarde la séquence comme un tout —
deux individus peuvent partager tous leurs sites variables sauf un et
compter comme deux haplotypes différents.

**MPD — diversité nucléotidique (π)**
La mesure la plus connue en génétique des populations. On prend toutes
les paires possibles d'individus de l'échantillon, on compte pour
chaque paire le nombre de sites où ils diffèrent, et on fait la
moyenne. Répond à la question : "si je prends deux individus au hasard
dans cette population, à quel point sont-ils différents en moyenne ?"

**VPD — variance de la diversité nucléotidique**
Est-ce que toutes les paires d'individus sont à peu près également
distantes les unes des autres (VPD faible), ou est-ce qu'il y a une
grande hétérogénéité — certaines paires très proches, d'autres très
éloignées (VPD élevée) ? Une variance élevée peut par exemple trahir
la présence de deux sous-groupes distincts cachés dans l'échantillon.

**DTA — D de Tajima**
Un test de neutralité classique. Il existe deux façons indépendantes
d'estimer la diversité génétique théorique d'une population : à partir
de MPD (la diversité moyenne par paire) ou à partir de NSS (le nombre
de sites variables). Sous un modèle simple — population de taille
stable, pas de sélection — les deux devraient donner à peu près le
même résultat, et DTA ≈ 0. Un DTA négatif suggère un excès de mutations
rares (typique d'une population en expansion récente, ou sous sélection
purificatrice). Un DTA positif suggère un déficit de variants rares
(typique d'un goulot d'étranglement démographique récent, ou d'une
sélection balancée).

**PSS — sites ségrégeants privés**
Parmi les sites variables de cette population, combien sont
**spécifiques** à elle — c'est-à-dire ségrégeants ici, mais parfaitement
fixés (pas de variation du tout) dans TOUTES les autres populations du
jeu de données ? Une population avec beaucoup de sites privés porte des
mutations qui lui sont propres, pas partagées — signe d'isolement
génétique ou de divergence ancienne.

**MNS — compte moyen de l'allèle minoritaire**
Sur les sites variables, on regarde à chaque fois quelle est la base la
moins fréquente et combien d'individus la portent (par exemple, si 38
individus ont "A" et 2 ont "G" à un site, l'allèle minoritaire est "G",
compté 2 fois). MNS est la moyenne de ces comptes sur tous les sites
variables. Un MNS bas veut dire que les variants rares sont vraiment
rares (souvent un ou deux individus) — signe de mutations récentes.
Un MNS plus élevé veut dire que même le variant "minoritaire" est
partagé par une bonne partie de l'échantillon — variation plus ancienne
ou plus équilibrée.

**VNS — variance du compte de l'allèle minoritaire**
Est-ce que la rareté des variants est homogène d'un site à l'autre
(VNS faible), ou très inégale — certains sites avec un variant
extrêmement rare, d'autres avec un variant presque à 50/50 (VNS
élevée) ?

## Statistiques PAR PAIRE de populations (différenciation)

Les 5 suivantes reprennent les mêmes idées que ci-dessus, mais
appliquées à deux populations à la fois, pour mesurer leur
**différenciation** plutôt que la diversité de chacune séparément.

**NH2** — comme NHA, mais en regroupant les échantillons des deux
populations ensemble dans un seul pool avant de compter les haplotypes
distincts. Si les deux populations partagent largement les mêmes
haplotypes, NH2 restera proche du nombre d'haplotypes d'une population
seule. Si elles ont des répertoires d'haplotypes complètement
différents, NH2 se rapprochera de la somme des deux — signe de
divergence.

**NS2** — comme NH2, mais pour NSS : on regroupe les échantillons des
deux populations et on compte les sites ségrégeants dans ce pool
combiné. Un site peut très bien être fixé dans pop `i` seule ET dans
pop `j` seule, tout en étant ségrégeant une fois les deux réunies (si
elles n'ont pas la même base fixée chacune de leur côté) — NS2 est donc
toujours ≥ au NSS de chaque population prise séparément.

**MP2 — diversité "au sein", pool des deux populations**
Comme MPD, mais en ne comparant jamais un individu de pop `i` à un
individu de pop `j` — on prend les paires à l'intérieur de pop `i`
d'un côté, à l'intérieur de pop `j` de l'autre, et on combine les deux
résultats (pondérés par le nombre de paires de chacune). Une façon de
résumer "à quel point les individus sont différents entre eux, quand on
reste à l'intérieur d'une même population" — sur les deux populations à
la fois.

**MPB — diversité "entre" les deux populations**
Le miroir de MP2 : cette fois on ne compare QUE des individus de
populations différentes (un de pop `i` avec un de pop `j`, jamais deux
individus de la même population), et on fait la moyenne sur toutes ces
paires croisées. Répond à : "si je prends un individu au hasard dans
chaque population, à quel point sont-ils différents en moyenne ?"

**HST — différenciation génétique (type FST)**
Compare MPB (diversité ENTRE les populations) à MP2 (diversité AU SEIN
des populations) : `HST = (MPB - MP2) / MPB`. Si les deux populations
sont indifférenciées, un individu pris au hasard dans l'une n'est pas
plus différent d'un individu de l'autre que de quelqu'un de sa propre
population — MPB ≈ MP2, et HST ≈ 0. Plus les populations ont divergé,
plus MPB dépasse MP2, et HST se rapproche de 1. C'est la statistique la
plus directement interprétable comme "à quel point ces deux populations
sont-elles génétiquement distinctes ?".

## État d'avancement (2026-08-25)

Les 13 statistiques (NSS, NHA, MPD, VPD, DTA, PSS, MNS, VNS par
population ; NH2, NS2, MP2, MPB, HST par paire) sont toutes
implémentées et testées dans `bridge/summary_statistics.py`. Reste à
faire : validation contre un vrai reftable DIYABC (pas encore généré
pour `toy_example2_ms_dna`), câblage dans `pipeline.py`/
`reftable_loop.py`, et la question ouverte des noms de colonnes
identiques entre groupes G2/G3 du header.
