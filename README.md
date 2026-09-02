# POC — Remplacement du simulateur DIYABC par msprime

## Objectif
Démontrer la faisabilité de remplacer le simulateur génétique de DIYABC
(C++, `particuleC.cpp::dosimulpart`) par `msprime` : un pipeline
`header.txt` → `msprime.Demography` → coalescence+mutation → statistiques
résumées, entièrement en Python, produisant un `reftable.bin`
structurellement et statistiquement équivalent à celui du vrai DIYABC.

## État actuel (2026-09-02)

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

**MicroSat** : parsing du header (format détaillé, priors de groupe)
complet ; aucune simulation (modèle de mutation stepwise, statistiques
`NAL`/`HET`...) n'est encore implémentée.

Voir `CLAUDE.md` pour l'architecture détaillée et l'historique complet
des investigations, `notes/exploration.md` pour le journal de recherche
brut (citations de code source, diagnostics, bugs trouvés/corrigés).

## Structure
- `reference/` — fichiers générés par le DIYABC historique (NE PAS MODIFIER)
- `bridge/` — pipeline Python (`header.txt` → `msprime` → stats → `reftable.bin`) :
  - `scenario_parser.py` / `header_dataclasses.py` — parsing des scénarios
  - `prior_parser.py` / `parameter_sampling.py` — priors et tirage sous contraintes
  - `loci_parser.py` — parsing de la description des loci (formats condensé et détaillé)
  - `demography_builder.py` — construction de la `Demography` msprime
  - `observed_data.py` — lecture des fichiers observés (`.snp`/`.mss`), mapping population
  - `ancestry_simulation.py` — coalescence + mutation (SNP et séquences ADN)
  - `summary_statistics.py` — statistiques résumées (SNP, PoolSeq, ADN)
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
