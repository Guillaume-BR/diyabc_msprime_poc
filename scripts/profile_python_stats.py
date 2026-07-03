import time

from bridge.pipeline import compute_summary_statistics

t0 = time.time()
stats, values = compute_summary_statistics(
    reference_directory="reference/human",
    scenario_index=1,
    num_loci=5000,
    seed=1,
)
t1 = time.time()
print(f"Temps total (5000 loci, 130 stats Python) : {t1 - t0:.1f}s")
print(f"Nombre de stats calculées : {len(stats)}")
print(
    f"Estimation pour 1000 particules (séquentiel) : {(t1 - t0) * 1000 / 60:.1f} minutes"
)
