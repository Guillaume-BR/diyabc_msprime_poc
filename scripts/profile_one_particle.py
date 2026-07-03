import os
import shutil
import subprocess
import time
from pathlib import Path

from bridge.loci_parser import rewrite_loci_count
from bridge.pipeline import run_poc_for_directory
from bridge.snp_writer import write_snp_file

GENERAL_BINARY_PATH = os.environ["DIYABC_GENERAL_PATH"]
WORK_DIR = Path("./tmp/profile_one")
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

t0 = time.time()
genotypes_per_locus, values = run_poc_for_directory(
    "reference/human", scenario_index=1, num_loci=5000, seed=1
)
genotypes_list = list(genotypes_per_locus)
t1 = time.time()
print(f"Simulation msprime (5000 loci) : {t1 - t0:.1f}s")

with open("reference/human/header.txt") as f:
    header_text = f.read()
snp_filename = header_text.splitlines()[0].strip()
write_snp_file(genotypes_list, WORK_DIR / snp_filename)
t2 = time.time()
print(f"Écriture .snp : {t2 - t1:.1f}s")

adapted = rewrite_loci_count(header_text, 5000)
(WORK_DIR / "header.txt").write_text(adapted)
shutil.copy("reference/human/RNG_state_0000.bin", WORK_DIR / "RNG_state_0000.bin")
t3 = time.time()
print(f"Préparation fichiers : {t3 - t2:.1f}s")

subprocess.run(
    [
        GENERAL_BINARY_PATH,
        "-p",
        "./",
        "-R",
        "ALL",
        "-r",
        "1",
        "-g",
        "1",
        "-m",
        "-t",
        "1",
    ],
    cwd=WORK_DIR,
    check=True,
    capture_output=True,
)
t4 = time.time()
print(f"Appel binaire general : {t4 - t3:.1f}s")

print(f"\nTOTAL : {t4 - t0:.1f}s")
