# validate_stats.py
import os
import shutil
import subprocess
from pathlib import Path

from bridge.loci_parser import rewrite_loci_count
from bridge.pipeline import run_poc_for_directory
from bridge.snp_writer import write_snp_file
from bridge.statobs_parser import parse_statobs
from bridge.summary_statistics import compute_all_statistics

GENERAL = os.environ["DIYABC_GENERAL_PATH"]
WORK = Path("./tmp/validate_stats")
if WORK.exists():
    shutil.rmtree(WORK)
WORK.mkdir(parents=True)

# 1. Simuler un petit cas (50 loci pour avoir des stats stables)
num_loci = 50
genotypes_per_locus, _ = run_poc_for_directory(
    "reference/human", scenario_index=1, num_loci=num_loci, seed=42
)
genotypes_list = list(genotypes_per_locus)

# 2. Calculer les stats avec NOS formules Python
pop_names = ["pop1", "pop2", "pop3", "pop4"]
our_stats = compute_all_statistics(genotypes_list, pop_names)

# 3. Passer les MÊMES données au binaire general
with open("reference/human/header.txt") as f:
    header_text = f.read()
snp_filename = header_text.splitlines()[0].strip()
write_snp_file(genotypes_list, WORK / snp_filename)
adapted = rewrite_loci_count(header_text, num_loci)
(WORK / "header.txt").write_text(adapted)
shutil.copy("reference/human/RNG_state_0000.bin", WORK / "RNG_state_0000.bin")
subprocess.run(
    [GENERAL, "-p", "./", "-R", "ALL", "-r", "1", "-g", "1", "-m", "-t", "1"],
    cwd=WORK,
    check=True,
    capture_output=True,
)
general_stats = parse_statobs((WORK / "statobsRF.txt").read_text())

# 4. Comparer terme à terme
print(
    f"\n{'Statistique':<20} {'Nos formules':>15} {'general':>15} {'écart relatif':>15}"
)
print("-" * 70)
for key in sorted(our_stats):
    if key in general_stats:
        ours = our_stats[key]
        theirs = general_stats[key]
        if theirs != 0:
            rel_err = abs(ours - theirs) / abs(theirs)
        else:
            rel_err = abs(ours - theirs)
        ok = "✓" if rel_err < 0.01 else "✗"
        print(f"{ok} {key:<20} {ours:>15.6f} {theirs:>15.6f} {rel_err:>15.2e}")
    else:
        print(f"? {key:<20} {our_stats[key]:>15.6f} {'(absent)':>15}")
