"""
Pour mesurer le temps de calcul de la simulation de 1000 particules de 5000 loci chacune, et l'écriture du fichier reftable.bin.
"""

import time

from bridge.reftable_loop import simulate_reference_directory

t0 = time.time()
results = simulate_reference_directory(
    test_directory="reference/toy_example5_1000loci",
    num_loci=None,
    nrec=1000,
    stats_filter="HEADER",
    max_workers=16,
)
t1 = time.time()

print(
    f"Simulation pour {len(results)} particules : {t1 - t0:.1f}s ({(t1 - t0) / 60:.1f} min)"
)
