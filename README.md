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