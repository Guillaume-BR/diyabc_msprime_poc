from bridge.pipeline import run_poc_for_directory
from bridge.snp_writer import write_snp_file

genotypes_per_locus, values = run_poc_for_directory(
    "reference/human",
    scenario_index=1,
    num_loci=10,
    seed=42,
)

write_snp_file(list(genotypes_per_locus), "./tmp/test_python_stats/human_snp_all22chr_maf5.snp")

print("Fichier .snp écrit. Valeurs de paramètres tirées :")
for name, value in values.items():
    print(f"  {name}: {value:.2f}")
