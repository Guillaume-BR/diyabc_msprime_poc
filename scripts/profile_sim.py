import random
import time

# Ajoute à profile_sim.py, après le bloc Hudson existant :
import numpy as np

from bridge.ancestry_simulation import (
    _draw_single_mutation_node,
    build_samples_argument,
    simulate_independent_loci,
)
from bridge.pipeline import build_random_demography_for_scenario_index


def draw_fast(tree, ts, rng):
    node_times = ts.tables.nodes.time
    nodes = np.fromiter((u for u in tree.nodes() if u != tree.root), dtype=np.int64)
    parents = np.fromiter((tree.parent(u) for u in nodes), dtype=np.int64)
    lengths = node_times[parents] - node_times[nodes]
    total = lengths.sum()
    target = rng.uniform(0, total)
    idx = np.searchsorted(np.cumsum(lengths), target)
    if idx >= len(nodes):
        idx = len(nodes) - 1
    return int(nodes[idx])


with open("reference/human/header.txt") as f:
    header_text = f.read()
demography, _ = build_random_demography_for_scenario_index(header_text, 1, 1)
samples = build_samples_argument("reference/human/human_snp_all22chr_maf5.snp")

# 1. Temps de sim_ancestry seul (coalescence, sans mutation)
t0 = time.time()
tree_sequences = list(simulate_independent_loci(demography, samples, 5000, 1))
t1 = time.time()
print(f"sim_ancestry (5000 loci, coalescence) : {t1 - t0:.2f}s")

# 2. Temps de l'algo de Hudson (tirage mutation + génotypes)
rng = random.Random(1)
t2 = time.time()
for ts in tree_sequences:
    tree = ts.first()
    node = _draw_single_mutation_node(tree, rng)
    derived = set(tree.samples(node))
t3 = time.time()

print(f"Algorithme de Hudson (5000 loci)      : {t3 - t2:.2f}s")

rng2 = random.Random(1)
t4 = time.time()
for ts in tree_sequences:
    tree = ts.first()
    node = draw_fast(tree, ts, rng2)
    derived = set(tree.samples(node))
t5 = time.time()
print(f"Hudson vectorisé (tirage)             : {t5 - t4:.2f}s")

rng3 = random.Random(1)
t6 = time.time()
for ts in tree_sequences:
    tree = ts.first()
    node = _draw_single_mutation_node(tree, rng3)
    # SANS tree.samples()
t7 = time.time()
print(f"Hudson SANS tree.samples()            : {t7 - t6:.2f}s")

# Et juste tree.samples() seul
rng4 = random.Random(1)
t8 = time.time()
for ts in tree_sequences:
    tree = ts.first()
    _ = set(tree.samples(tree.root))  # tous les samples
t9 = time.time()
print(f"tree.samples() seul                   : {t9 - t8:.2f}s")

# Ajoute à profile_sim.py :
t10 = time.time()
count = 0
for ts in tree_sequences:
    tree = ts.first()
    for _u in tree.nodes():
        count += 1
t11 = time.time()
print(f"Itération tree.nodes() seule          : {t11 - t10:.2f}s ({count} noeuds)")

# Et le coût de tree.first() seul (matérialiser l'arbre)
t12 = time.time()
for ts in tree_sequences:
    tree = ts.first()
t13 = time.time()
print(f"tree.first() seul                     : {t13 - t12:.2f}s")

# Ajoute à profile_sim.py :
t14 = time.time()
for ts in tree_sequences:
    tree = ts.first()
    lengths = [tree.branch_length(u) for u in tree.nodes() if u != tree.root]
t15 = time.time()
print(f"branch_length() dans la boucle        : {t15 - t14:.2f}s")


def draw_vec(ts, rng):
    edges = ts.tables.edges
    nt = ts.tables.nodes.time
    lengths = nt[edges.parent] - nt[edges.child]
    total = lengths.sum()
    target = rng.uniform(0, total)
    idx = np.searchsorted(np.cumsum(lengths), target)
    if idx >= len(edges.child):
        idx = len(edges.child) - 1
    return int(edges.child[idx])


rng5 = random.Random(1)
t16 = time.time()
for ts in tree_sequences:
    node = draw_vec(ts, rng5)
    derived = set(ts.first().samples(node))
t17 = time.time()
print(f"Hudson tables vectorisé (complet)     : {t17 - t16:.2f}s")
