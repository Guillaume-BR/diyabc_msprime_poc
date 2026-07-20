"""
Rejoue, particule par particule, les tirages de priors RÉELLEMENT
effectués par DIYABC dans un reftable existant (ex:
first_records_of_the_reference_table_0.txt) -- au lieu d'en tirer de
nouveaux indépendamment côté msprime (voir reftable_loop.
run_reftable_simulation pour ça).

Permet une comparaison appariée DIYABC/msprime : chaque particule
msprime utilise EXACTEMENT le même (N1,N2,N3,ta,ts,...) que la
particule DIYABC de même rang -- tout écart entre les deux résulte
alors uniquement du moteur de simulation, jamais d'un tirage de prior
différent (voir notebook/correlation_N2_N3_HWm_anomaly.ipynb, section
12, pour l'analyse qui a motivé ce script à l'origine).
"""

import shutil
from pathlib import Path
from time import time

from bridge.loci_parser import parse_loci_description
from bridge.prior_parser import parse_priors
from bridge.reftable_loop import replay_reftable_simulation, write_reftable_txt
from bridge.scenario_parser import parse_header_scenarios

REFERENCE_DIR = Path("reference/human_modif_scenario1_5000loci")
REAL_REFTABLE_PATH = REFERENCE_DIR / "first_records_of_the_reference_table_0.txt"
# None = utilise les vrais comptes par type déclarés dans header.txt
# (500 <A> / 50 <X> / 50 <M> / 50 <Y> pour ce dataset), pas un override.
NUM_LOCI = None
STATS_FILTER = "HEADER"
OUTPUT_PATH = REFERENCE_DIR / "reftable_msprime_replay.txt"

start_time = time()
print("Début du lecture du headerRF.txt")
header_text = (REFERENCE_DIR / "headerRF.txt").read_text()

priors, _constraints = parse_priors(header_text)
# TOUS les scénarios candidats, pas un seul : le reftable réel mélange
# les 3 scénarios de toy_example5 (chaque ligne garde le scenario_index
# RÉELLEMENT tiré par DIYABC, lu par parse_real_reftable_params) -- ne
# filtrer qu'un seul scenario ferait planter le rejeu des lignes tirées
# sur un autre scénario (KeyError dans _kept_param_names_by_scenario).
scenarios = parse_header_scenarios(header_text)
t1 = time()
print(
    f"Lecture du headerRF.txt terminée en {t1 - start_time:.2f}s : {len(priors)} priors, {len(scenarios)} scénarios"
)

# read_header_text (bridge/pipeline.py) préfère header.txt à headerRF.txt
# si les deux existent dans le dossier -- or c'est headerRF.txt qui
# déclare le vocabulaire de statistiques ML1p/HWm/... compatible avec
# compute_all_statistics (header.txt utilise souvent un vocabulaire
# obsolète, ex: HP0/HM1 sur human_modif_scenario1). On simule donc dans
# une copie de travail SANS header.txt, pour forcer la lecture de
# headerRF.txt -- même piège que documenté dans
# notebook/correlation_N2_N3_HWm_anomaly.ipynb, section 10.
snp_filename = header_text.splitlines()[0].strip()
WORK_DIR = Path("./tmp/replay_diyabc_priors") / REFERENCE_DIR.name
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)
shutil.copy(REFERENCE_DIR / "headerRF.txt", WORK_DIR / "headerRF.txt")
shutil.copy(REFERENCE_DIR / snp_filename, WORK_DIR / snp_filename)

t2 = time()
print(f"Copie du répertoire de travail terminée en {t2 - t1:.2f}s : {WORK_DIR}")


header_text = (WORK_DIR / "headerRF.txt").read_text()
total_loci = parse_loci_description(header_text).total_loci.values()
total_loci = sum(total_loci)

print(
    f"Début de la simulation des {total_loci} loci pour rejouer les tirages de priors DIYABC..."
)
results = replay_reftable_simulation(
    reference_directory=WORK_DIR,
    priors=priors,
    scenarios=scenarios,
    real_reftable_path=REAL_REFTABLE_PATH,
    num_loci=NUM_LOCI,
    stats_filter=STATS_FILTER,
    max_workers=16,
)

t3 = time()
print(
    f"Rejeu des {len(results)} particules de {total_loci} loci terminé en {t3 - t2:.2f}s"
)

write_reftable_txt(results, priors, scenarios, OUTPUT_PATH)
t4 = time()
print(f"Écriture du reftable terminé en {t4 - t3:.2f}s")
print(f"{len(results)} particules rejouées avec les tirages réels de DIYABC.")
print(f"Écrit dans {OUTPUT_PATH}")
