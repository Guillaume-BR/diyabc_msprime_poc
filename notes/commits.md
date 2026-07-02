# Git Commits

```
7f88b9e bug fixed by ruff
a510cbc nettoyage validate_stats.py
4f92d82 ruff style
298388b test ruff
f51007f test ruff precommit
724e72a POC validé avec temps réduit
cf6b76e Optimisation vectorielle numpy : _fst_wc, NEI, et algorithme de Hudson (tables tskit) -- 1000 particules en 8.6min vs 20min, plus rapide que DIYABC -g 1 -t 8
a0b3ae3 Pipeline 100% Python complet et validé : 1000 particules × 5000 loci en ~20min (8 cœurs), vs 2m48s DIYABC avec -g 50 -t 8
3a97c4a nettoyage
4f10273 Pipeline 100% Python : compute_summary_statistics branchée sur summary_statistics.py (130 stats validées, sans subprocess ni fichier .snp intermédiaire) -- 6.9s/particule sur 5000 loci
7fc5cdd summary_statistics.py : toutes les 130 stats SNP implémentées et validées par comparaison au binaire general (ML1-3, HW, HB, FST1-4, NEI, AML, F3, F4)
2cc2053 CORRECTIF PERFORMANCE CRITIQUE : -g 50 -> -g 1 (facteur 29x), -g est la taille de batch interne, sans rapport avec notre architecture par particule
47ff325 Bug fixed with ruff
065b5df Bug fixed with ruff
5ece139 VALIDATION COMPLETE : reftable.bin produit par notre pipeline Python lu avec succès par readReftable.R (référence indépendante)
05a8759 reftable_loop.py : write_reftable_bin, écrivain du format binaire avec filtrage des paramètres constants
8aedfb5 reftable_loop.py : boucle parallélisée (ProcessPoolExecutor) produisant nrec particules, corrige seed=0 rejetée par msprime
be2b1a0 compute_summary_statistics : orchestration complète (simulation -> .snp -> binaire C++ -> 130 stats), validé empiriquement et testé
ad8e5e0 Validation empirique : architecture .snp simulé -> binaire C++ -> statobsRF.txt fonctionne (testé manuellement sur 10 loci)
c6c0780 snp_writer.py : écriture de fichier .snp DIYABC depuis des génotypes simulés, validé sur cas minimal
66ba5db bug fixed : all scenarios pass
fbebfa2 ajout de SplitEvent
d84a716 class SplitEvent créée
8700760 ancestry_simulation.py : remplace le modèle à taux fixe par l'algorithme de Hudson (une mutation par locus), génotypes regroupés par population
ef9593c  fix : nb_attemps
0cd2186 pipeline.py : run_poc_for_directory (point d'entrée par chemin de dossier) ; clarification : écriture reftable.bin dépend des stats résumées (étape suivante)
6145811 Point d'étape : pipeline complet jusqu'à la simulation mutée, reste les statistiques résumées (sumstat.cpp à explorer)
4fd14e8 ancestry_simulation.py : mutate_independent_loci avec modèle binaire et graines dérivées par locus, validé sur scénario 1 human
5853b90 ancestry_simulation.py : simulate_independent_loci + build_samples_argument, simulation msprime complète validée sur scénario 1 human (10 loci, 240 lignées)
4db681e observed_data.py : mapping indice de population -> nom réel, documenté et testé sur human
a55f67d [200~observed_data.py : comptage des échantillons par population depuis le fichier .snp, validé sur human (4 pops x 30 ind)~
4115466 pipeline.py : orchestration complère header.txt -> Demography msprime testée de bout en but sur le scenaio 1
4b66551 parameter_sampling : tirage avec contraintes d'ordre, validé sur les 21 priors / 4 contraintes de human (déterminisme confirmé)
35f22b7 demography_builder : evaluate_expression + build_demography, validés sur scénario 1 human
48cd8aa prior_parser + parameter sampling : tirage de valeurs avec contraintes
712ee76 scenario_parser : parsing fonctionnel de sample/varNE/merge
b913861 structure initiale du POC + fichier de ref human

```