"""
Identique à replay_diyabc_priors_dna.py mais pointant sur
reference/toy_example2_ms_dna_50loci (50 loci <A> + 50 loci <M> au lieu
de 5+5) -- test de l'hypothèse "bruit d'échantillonnage lié au faible
nombre de loci" pour l'écart statistique observé sur MNS/VNS/DTA/VPD
(voir notes/exploration.md, précédent SNP du 2026-07-17).
"""

import shutil
from pathlib import Path
from time import time

from bridge.loci_parser import parse_loci_description
from bridge.prior_parser import parse_group_priors, parse_priors
from bridge.reftable_loop import (
    group_prior_column_names,
    replay_reftable_simulation_dna,
    write_reftable_txt,
)
from bridge.scenario_parser import parse_header_scenarios

REFERENCE_DIR = Path("reference/toy_example2_ms_dna_50loci")
REAL_REFTABLE_PATH = REFERENCE_DIR / "first_records_of_the_reference_table_0.txt"
OUTPUT_PATH = REFERENCE_DIR / "reftable_msprime_replay.txt"

start_time = time()
print("Début du lecture du headerRF.txt")
header_text = (REFERENCE_DIR / "headerRF.txt").read_text()

priors, _constraints = parse_priors(header_text)
group_priors = parse_group_priors(header_text)
group_priors_names = group_prior_column_names(header_text)

# TOUS les scénarios candidats, pas un seul : le reftable réel mélange
scenarios = parse_header_scenarios(header_text)
t1 = time()
print(
    f"Lecture du headerRF.txt terminée en {t1 - start_time:.2f}s : {len(priors)} priors, {len(group_priors)} group priors, {len(scenarios)} scénarios"
)

mss_filename = header_text.splitlines()[0].strip()
WORK_DIR = Path("./tmp/replay_diyabc_priors") / REFERENCE_DIR.name
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)
shutil.copy(REFERENCE_DIR / "headerRF.txt", WORK_DIR / "headerRF.txt")
shutil.copy(REFERENCE_DIR / mss_filename, WORK_DIR / mss_filename)

t2 = time()
print(f"Copie du répertoire de travail terminée en {t2 - t1:.2f}s : {WORK_DIR}")


header_text = (WORK_DIR / "headerRF.txt").read_text()
total_loci = len(parse_loci_description(header_text))

print(
    f"Début de la simulation des {total_loci} loci pour rejouer les tirages de priors DIYABC..."
)
results = replay_reftable_simulation_dna(
    reference_directory=WORK_DIR,
    priors=priors,
    group_priors_names=group_priors_names,
    scenarios=scenarios,
    real_reftable_path=REAL_REFTABLE_PATH,
    stats_filter="ALL",
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
