from bridge.pipeline import compute_summary_statistics

summary_statistics, values = compute_summary_statistics(
    reference_directory="reference/human",
    scenario_index=1,
    num_loci=10,
    seed=42,
    general_binary_path="/home/bernardr/Documents/Github/diyabc/build/src-JMC-C++/general",
    work_directory="./tmp/test_python_stats",
    stats_filter="ALL",
)

print("Fichier .snp écrit. Valeurs de paramètres tirées :")
