# Git Commits

```
35d3164 refatble_loop.py prêt pour les microsats
3066a16 fix docstings
fbada23 summaray_stats + pipeline + tests_dédiés
a89720b fichier mss nécessaire pour certains tests
c50e37c ajout de fonctions d'enchainement dans ancestry.py + test dédié
740b9e6 rajout de tests de foncitons sur ancestry_simulation.py
ba5371c avancée sur les microsat
90d655e observed_data.py + tests
8386166 observed_microsat + tests
068fee1 fix configuration.py
f890587 refactor contante dans configuration.py
fdba21c dispatch XY pour les séquences ADN.
b35ab51 scripts de rejeu
f0c9ed7 observed_data pour les microsat
2cfab8d MAJ Claude.md
c6fe273 suppression msprime_cpp
c6780b3 MAJ readme
dc0a7c9 bug G3 M fixed
1e558c7 Revert "écart sur G3 mitochondial fixé"
bf6d358 écart sur G3 mitochondial fixé
9bfaedc reftable_loop et test_reftable_loop
b414f75 claude.md mis à jour sur les docstring en mode google
8462288 docstrings de snp_writer mis ua format google
d209905  docstring au foramt google de reftable_loop.py et réorganisation
af974d4  docstring au foramt google de reftable_loop.py et réorganisation
c72d9d1 docstrings de pipeline.py au format google
04a1a4a summary_statistics.py docstrings au format google
3757b76 Docsstrings ancestry_simulation.py en format google + header_dataclasses.py
918f2a4 modification de noms de variables
5b923ea docstrings stats_group_parser format google
b6bf3f6 Docstring de parameter_sampling en format Google
4f0869e docstrings parameter_sampling mis au format Google
fdff51e docstrings de demo_graphy_builder passé en mode Google
a3d9c21 docstrings de demo_graphy_builder passé en mode Google
3fca24d docstrings de observed_data.py en mode Google
329e771 docstring statobs_parser.py
3bd9085 changement d'un nom d'attribut total_loci -> loci_count_by_heritage
44bc413 docsstring revue pour priorparser et scenario_parser
c090ac6 changment dans la classe Prior où bounds devient min,max,mean,sdshape
f9c0409 rename scenario_types en header_dataclasses et docstrings modèle google
4ca688f ancestry_simulation.py
e5d2eef ajout de Claude readme
0bf6931 fichiers de tests pipeline et summary_statistics
5006bc3 ajout des scripts de rejeux pour les séquences dna
447eebd rejeu du fichier entier des paramètres fait : reftable_loop.py
96e2fbb reftable_loop.py mAJ avec _run_sigle_particle_dna et from_value ajoutés
df27638 mise à jour de pipeline.py avec ajout de compute_stats dna et la version from_values
12ecb22 simulation des mutation ADN
9185862 Récupération par locus des valuers de simulation diyabc ancestry_simualtion.py + test
b449fe6 Add parse_real_reftable_params_with_group_priors for DNA replay (paired historical + group prior values)
d167f5f Rename DNA summary stats to compute_* convention, add compute_all_statistics_dna, fix multi-group column naming, generate real DIYABC reftable for toy_example2_ms_dna
d2a17c6 Claude et resume_stat
400c574 HST codé et testé
0df1b6b MP2 codé et testé
b03b958 NS2 implémenté avec test
2537a46 NH2 implémenté avec test
e7f65a5 Stats DTA, PSS, MNS, VNS avec tests implémentées
cf27e96 Stats MDP et VDP avec tests implémentées
16b47c2 Add DNA sequence summary statistics (NSS, NHA) and fix population_layout naming bug
8d3e217 Pipeline de mutation ADN complet : matQ, RateMap, placement via msprime.sim_mutations
9348b93 cablage afin d'obtenir les matrices de transition par locus faisant parti d'un groupe avec modèle  mutationnel
604cb8a modèle de mutation et matrice de transition implémentés
e873338 parsing des fichiers mss pour obtenir la fréqeunce des bases
91cd3de tirage des valeurs pour les groupes de priors dans les fichiers microsat
b12fa1f test no draw until te1
f5f7a32 fonctions et tests du parsing du header pour les microsat/dna-sequences
9dd3fab relecture du projet et refactorisations mineures
8e6b60b with_mrc_filter : pool de généalogies partagé entre tous les loci (au lieu d'un pool par-locus)
466705e with_maf_filter : batch_size scale avec num_loci (max(20, num_loci//4))
929d407 Aligne ruff-pre-commit (v0.8.0 -> v0.15.1) sur la version de l'env conda diyabc_msprime
02bcf23 Supprime les tests MAF dépendant de reference/toy_example3_scenario1 (dossier retiré)
c114ab0 maf_filter mrc_filter batching pour améliorer la vitesse d'exécution
4c340f6 Amélioration perf poolseq
7585e35 Corrige 2 bugs PoolSeq causant une divergence massive (130/130 stats) sur reftable complet
ffd0bb7 Corrige observed_reads : purge MRC manquante faussait les 130 stats PoolSeq
ab3ff13 branchement de la version poolseq
e0ba46b Formules stats et tests correspondants pour POOLSEQ
f040071 Ajoute simulate_poolseq_reads_with_mrc_filter et corrige le lint ruff
26815ce Simulation pour poolseq avec prise en charge du mrc et test ad hoc
c638154 Codage de simulate_pooseq_reads
a3e1525 parsing des lignes de fichiers snp POOLSEQ aevc test
9210198 gitignore
9f2a71b observed_data.py modifié pour fichier SNP POOLSEQ et tests correspondants
f53c21e Remplace le gain de perf extrapole par le chiffre reel mesure sur 1000 particules
9b96f69 Deduplique build_samples_argument : plus de double scan du .snp par particule
965f31e Etend le cache population/samples au chemin with_maf_filter (maf>0)
e6176be Optimise simulate_snp_genotypes : cache population/samples par particule au lieu de par locus
9137151 Met a jour CLAUDE.md pour refleter la portee actuelle du POC et documente observed_data.py
c347b8c Documente l'investigation de l'ecart de performance DIYABC/msprime et petits ajustements
a4b8b5e Limite les threads BLAS/OMP par worker pour éviter la contention CPU en parallélisation
e285947 Ajoute le filtre MAF pour les loci SNP (with_maf_filter / with_maf_filter_shared_ancestry)
3503c41 Petits ajustements de commentaires et docstrings
12736c7 Simplifie loci_parser.py et clarifie l'organisation de 3 modules bridge/
ea53c00 Corrige la correlation seed scenario/parametres, ajoute le rejeu des priors reels
64f87bd num_loci accepte None pour utiliser les vrais comptes du header.txt
1ef3bde Ignore tmp/ et temp/ (scratch de dev) au lieu de les suivre
94a1a18 Support des loci <X>/<Y>/<M> dans le pipeline msprime, en plus de <A>
dbcfc02 initially_active=True sur les populations + petits correctifs
ee1a59b loci_parser gere le format multi-types condense
27b16ea ajoute rewrite_real_reftable_txt pour lire le reftable DIYABC brut
13ed4d0 corrige le rejeu/ecriture multi-scenario (cases vides DIYABC)
3368b1e corrige la doc obsolete sur le support multi-scenario
9a71503 documente la cause racine du bug de ligne de fin de header
30d5a15 ajoute un script pour rejouer les tirages reels de DIYABC
4767fa9 valide les stats sur human_modif_scenario1
96f05e5 rejeu des tirages reels DIYABC cote msprime
1df049c test DRAW UNTIL ts>ta avec N fixes
95e8094 Ajoute un test contrôlé bivarié N2×ts au notebook d'anomalie HWm
e5722fd Ajoute un test contrôlé bivarié N2×N3 au notebook d'anomalie HWm
f4a6f2e Ajoute un notebook documentant l'anomalie de corrélation N2/N3 <-> HWm_2/HWm_3
7276a64 Ajoute la comparaison DIYABC/msprime à paramètres démographiques fixes
aedc77c Parse la section 'group summary statistics' de header.txt
dd45c19 Corrige un décalage de colonnes dans write_reftable_txt (cases vides)
c3ecdfc Arrondit les priors N (taille) et T (temps) à l'entier après tirage
15d000b Ajoute simulate_reference_directory : point d'entrée par sous-dossier de test
3a05868 Supprime les dossiers de travail vides créés par particule dans reftable_loop.py
ae74bca Tirage du scénario pondéré par son poids + support multi-scénario du reftable
43e3dbc Implémente l'admixture (SplitEvent) dans build_demography
29fffbc Regroupe les scripts d'investigation ad hoc de la racine dans scripts/
ccc1adb Modification de reftable_loop.py : ajout de la fonction write_txt pour avoir aussi les sorties en txt afin de comparer avec le first_record de diyabc
1a6996a script rapport auto
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