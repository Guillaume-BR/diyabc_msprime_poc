# POC — Remplacement du simulateur DIYABC par msprime

## Objectif
Démontrer la faisabilité de remplacer le simulateur génétique de DIYABC
(C++, `particuleC.cpp::dosimulpart`) par `msprime` : un pipeline
`header.txt` → `msprime.Demography` → coalescence+mutation → statistiques
résumées, entièrement en Python, produisant un `reftable.bin`
structurellement et statistiquement équivalent à celui du vrai DIYABC.

## État actuel (2026-09-16)

**Côté SNP** (format condensé de `loci description`) : COMPLET et validé
contre le vrai DIYABC sur `human`, `toy_example3`/`4`/`5` :
- Tous les types d'héritage `<A>/<H>/<X>/<Y>/<M>`
- Événements split/admixture et scénarios multiples (tirage pondéré)
- Filtre MAF, PoolSeq (comptages de lecture poolés)
- 130 statistiques résumées (`sumstat.cpp`) portées en Python pur

**Côté séquences ADN** : le pipeline complet (modèle de substitution
`JK`/`K2P`/`HKY`/`TN`, placement des mutations via `msprime.sim_mutations`,
13 statistiques résumées, rejeu apparié DIYABC/msprime) est validé sur
`toy_example2_ms_dna`. Le déficit de variance observé sur les loci
mitochondriaux (`<M>`) — généalogie non partagée entre loci — a été
diagnostiqué et corrigé le 2026-09-02 (voir `notes/exploration.md`).

**MicroSat** : parsing du header complet ; le modèle de mutation GSM
(stepwise, via `msprime.TPM`) est simulé de bout en bout (généalogie +
mutation par locus, dispatch ploïdie/sexe `<A>/<X>/<Y>`, orchestration
`pipeline.py`/`reftable_loop.py` calquée sur le chemin ADN), la chaîne
de rejeu `_from_values` (mêmes tirages réels que DIYABC, pour une
comparaison appariée) est complète depuis le 2026-09-08, et les 11
statistiques résumées spécifiques (`NAL`/`HET`/`VAR`/`MGW`/`N2P`/`H2P`/
`V2P`/`FST`/`LIK`/`DAS`/`DM2`) sont implémentées et testées depuis le
2026-09-11 — voir "MicroSat GSM mutation model" et "MicroSat summary
statistics" dans `CLAUDE.md`. Ces 11 stats ont été comparées à un vrai
reftable DIYABC (`toy_example2_ms_dna`, 1000 particules) le 2026-09-14 :
9/11 collent bien, `FST` diverge d'un facteur ~1.8-2x — deux
hypothèses (canal SNI manquant, admixture) ont été rigoureusement
testées et réfutées, la cause reste ouverte (piste : décomposition
ANOVA `MSG`/`MSI`/`MSP`, la seule stat construite ainsi). Le canal de
mutation SNI est maintenant implémenté (2026-09-15, grille dense +
matrice de transition construite à la main) — voir "SNI mutation
channel" dans `CLAUDE.md` : ça confirme directement (pas seulement par
proxy) que SNI n'explique pas l'écart FST. L'implémentation initiale
introduisait un coût de performance réel (~34x sur `toy_example2_ms_dna`),
optimisé le 2026-09-16 (~1.5x seulement, en partageant `motif_size`
matrices `msprime.TPM` — une par classe de résidu — au lieu d'en
reconstruire une par état de la grille dense), `FST`/`DM2` inchangés
après l'optimisation.

**Mise à jour 2026-09-16/18** : la statistique `AML` (Choisy et al.
2004, coefficient d'admixture à 3 populations) est maintenant
implémentée (bissection sur `cal_Aml3p`) et validée à 100% contre un
vrai reftable DIYABC sur `toy_example1_ms_modified` (4 populations
réelles) — toutes les colonnes AML concordent, seules les colonnes FST
divergent encore. L'investigation FST a été approfondie avec ce nouveau
dataset à 4 populations : le déficit ne touche pas `<A>` et `<M>`
également (KS≈0.4 pour `<A>`, ≈0.75 pour `<M>` — quasiment le double),
et la formule a été revérifiée une 5ᵉ fois (correcte, y compris son
comportement négatif légitime à faible différenciation, vérifié par un
test synthétique exécuté). Piste actuelle, pas encore confirmée : un
déficit de structure de population dans la généalogie `<M>` partagée
elle-même, pas dans le calcul de la statistique.

Le pipeline a aussi été refactorisé pour ne lire `header.txt`/`.snp`/
`.mss` qu'**une seule fois par run** (pas une fois par particule) via
trois dataclasses contexte (`SnpReplayContext`/`MicrosatReplayContext`/
`DnaReplayContext`, voir `bridge/header_dataclasses.py`), transmises à
travers `ProcessPoolExecutor` — couvre les trois familles de données
(SNP/MicroSat/ADN) et leurs deux chemins chacune (tirage aléatoire et
rejeu `_from_values`).

Deux chantiers explicitement identifiés et non commencés : un dataset
(`toy_example1_ms`, non modifié) échantillonne une seule population
biologique à 4 temps différents plutôt que 4 populations distinctes —
non supporté par l'architecture actuelle (`demography_builder.py` ne
construit qu'une population msprime là où DIYABC en attend 4) ; et un
mécanisme DIYABC de bascule vers un modèle de substitution simplifié en
dessous d'un certain seuil de séquences, mentionné mais pas encore
localisé dans le C++.

Voir `CLAUDE.md` pour l'architecture détaillée et l'historique complet
des investigations, `notes/exploration.md` pour le journal de recherche
brut (citations de code source, diagnostics, bugs trouvés/corrigés).

## Structure
- `reference/` — fichiers générés par le DIYABC historique (NE PAS MODIFIER)
- `bridge/` — pipeline Python (`header.txt` → `msprime` → stats → `reftable.bin`) :
  - `scenario_parser.py` / `header_dataclasses.py` — parsing des scénarios
    (`header_dataclasses.py` héberge aussi les dataclasses contexte
    `SnpReplayContext`/`MicrosatReplayContext`/`DnaReplayContext`)
  - `prior_parser.py` / `parameter_sampling.py` — priors et tirage sous contraintes
  - `loci_parser.py` — parsing de la description des loci (formats condensé et détaillé)
  - `demography_builder.py` — construction de la `Demography` msprime
  - `observed_data.py` — lecture des fichiers observés (`.snp`/`.mss`), mapping population
  - `ancestry_simulation.py` — coalescence + mutation (SNP, séquences ADN et MicroSat)
  - `summary_statistics.py` — statistiques résumées (SNP, PoolSeq, ADN, MicroSat)
  - `stats_group_parser.py` — filtrage des colonnes de stats réellement demandées
  - `snp_writer.py` / `statobs_parser.py` — écriture/lecture au format DIYABC (chemin de validation croisée, plus le chemin par défaut)
  - `pipeline.py` — orchestration de haut niveau
  - `reftable_loop.py` — boucle multi-particules, écriture `reftable.bin`
- `tests/` — suite de tests automatisés (`pytest tests/ -v`)
- `notes/` — journal de recherche (découvertes, bugs, décisions)
- `scripts/` — scripts d'investigation ad hoc (jetables, non testés)

## Environnement
```bash
conda activate diyabc_msprime   # Python 3.11, msprime, tskit, numpy, scipy
pytest tests/ -v
```

## Documentation
- `CLAUDE.md` — référence complète : architecture, historique daté de
  chaque étape/bug/investigation, limitations connues.
- `notes/exploration.md` — journal brut de recherche sur le code source
  DIYABC (`particuleC.cpp`, `history.cpp`, `sumstat.cpp`, `data.cpp`,
  `header.cpp`).
