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
12, pour l'analyse qui a motivé ce script).
"""

import shutil
from pathlib import Path

from bridge.prior_parser import parse_priors
from bridge.reftable_loop import replay_reftable_simulation, write_reftable_txt
from bridge.scenario_parser import parse_header_scenarios

REFERENCE_DIR = Path("reference/human_modif_scenario1")
REAL_REFTABLE_PATH = REFERENCE_DIR / "first_records_of_the_reference_table_0.txt"
SCENARIO_INDEX = 1
NUM_LOCI = 100
STATS_FILTER = "HEADER"
OUTPUT_PATH = REFERENCE_DIR / "reftable_msprime_replay.txt"

header_text = (REFERENCE_DIR / "headerRF.txt").read_text()

priors, _constraints = parse_priors(header_text)
all_scenarios = parse_header_scenarios(header_text)
scenario = [s for s in all_scenarios if s.index == SCENARIO_INDEX]

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

results = replay_reftable_simulation(
    reference_directory=WORK_DIR,
    priors=priors,
    scenarios=scenario,
    real_reftable_path=REAL_REFTABLE_PATH,
    num_loci=NUM_LOCI,
    stats_filter=STATS_FILTER,
    max_workers=16,
)

write_reftable_txt(results, priors, scenario, OUTPUT_PATH)
print(f"{len(results)} particules rejouées avec les tirages réels de DIYABC.")
print(f"Écrit dans {OUTPUT_PATH}")
