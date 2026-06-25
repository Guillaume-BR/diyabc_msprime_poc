# POC — Remplacement du simulateur DIYABC par msprime

## Objectif
Démontrer la faisabilité de remplacer le simulateur génétique de DIYABC
(C++, particuleC.cpp::dosimulpart) par msprime, sur le cas le plus simple
possible : le scénario 1 du dataset `human` (4 populations, fusions
uniquement, pas de split/admixture).

## Critère de succès du POC
Produire, à partir du `header.txt` du dataset `human`, restreint au
scénario 1 :
1. Une `msprime.Demography` correcte (vérifiée par construction manuelle)
2. Une simulation de coalescence + mutation SNP
3. Des statistiques résumées comparables (même ordre de grandeur,
   même structure) à celles du `reftableRF.bin` de référence dans
   `reference/human/`

## Ce que ce POC NE couvre PAS (volontairement, pour l'instant)
- Les scénarios 2 à 6 (split/admixture, numérotation t11..t44)
- Les microsatellites (GSM) — risque connu, traité séparément
- L'intégration avec le binaire C++ (sous-processus) — étape suivante
- Les marqueurs liés au sexe / mitochondriaux

## Structure
- `reference/` — fichiers générés par le DIYABC historique (NE PAS MODIFIER)
- `bridge/` — code Python de traduction header → msprime → stats
- `tests/` — tests automatisés
- `notes/` — journal des découvertes faites en explorant le C++ source
- `docs/` — documents de synthèse

## Comment a été générée la référence
Depuis `diyabc/tests/datasets/human/` (dépôt diyabc/diyabc) :
    general -p ./ -R "FST1;ML1" -r 100 -g 50 -m -t 8

## État au 25/06/2026

Fait et testé (pytest tests/ -v, tout vert) :
1. Demography msprime correcte pour le scénario 1 — ✅
2. Simulation de coalescence + mutation SNP — ✅ (modèle de mutation
   simplifié, PAS une reproduction fidèle de l'algorithme DIYABC réel
   pour les SNP, qui reste non élucidé — voir notes/exploration.md)

Reste à faire :
3. Statistiques résumées comparables au reftableRF.bin de référence —
   prochaine étape : lire sumstat.cpp pour comprendre le calcul exact
   de FST1/ML1, puis les reproduire à partir des TreeSequence msprime
   (probablement via tskit.TreeSequence.Fst() ou un calcul manuel).

Modules du pipeline (bridge/) :
- scenario_types.py, scenario_parser.py — parsing scénarios ✅
- prior_parser.py, parameter_sampling.py — priors + tirage sous contraintes ✅
- demography_builder.py — construction Demography msprime ✅
- observed_data.py — comptage échantillons + mapping pop ✅
- ancestry_simulation.py — coalescence + mutation ✅
- pipeline.py — orchestration de tout ce qui précède ✅