# compare_reftable_te4.py
"""
Compare le reftable réel DIYABC (first_records_of_the_reference_table_0.txt)
et son rejeu msprime (reftable_msprime_replay.txt) sur toy_example4 --
comparaison appariée (mêmes tirages de priors des deux côtés, voir
scripts/replay_diyabc_priors.py pour comment reftable_msprime_replay.txt
est généré).

IMPORTANT : le fichier brut DIYABC écrit des cases VIDES (pas "nan")
pour les paramètres non utilisés par le scénario de chaque ligne --
un pandas.read_csv(sep=r"\s+") naïf dessus décale silencieusement les
colonnes (voir bridge.reftable_loop.parse_real_reftable_params). On
repasse donc d'abord par rewrite_real_reftable_txt pour obtenir une
version propre (cases vides -> "nan"), directement comparable à
reftable_msprime_replay.txt (même convention d'écriture, voir
bridge.reftable_loop.write_reftable_txt).
"""

import pandas as pd
from scipy import stats as ss

from bridge.prior_parser import parse_priors
from bridge.reftable_loop import rewrite_real_reftable_txt
from bridge.scenario_parser import parse_header_scenarios

REFERENCE_DIR = "reference/toy_example4"
REAL_REFTABLE_RAW = f"{REFERENCE_DIR}/first_records_of_the_reference_table_0.txt"
REAL_REFTABLE_CLEAN = f"{REFERENCE_DIR}/first_records_clean.txt"
MSPRIME_REFTABLE = f"{REFERENCE_DIR}/reftable_msprime_replay.txt"

with open(f"{REFERENCE_DIR}/headerRF.txt") as f:
    header_text = f.read()
priors, _constraints = parse_priors(header_text)
scenarios = parse_header_scenarios(header_text)

rewrite_real_reftable_txt(REAL_REFTABLE_RAW, REAL_REFTABLE_CLEAN, priors, scenarios)

diyabc = pd.read_csv(REAL_REFTABLE_CLEAN, sep=r"\s+")
msprime = pd.read_csv(MSPRIME_REFTABLE, sep=r"\s+")
print(f"diyabc: {diyabc.shape}, msprime: {msprime.shape}")

param_names = [p.name for p in priors]
param_cols = [c for c in param_names if c in diyabc.columns]
stat_cols = [c for c in diyabc.columns if c not in param_cols and c != "scenario"]

# 1. Verifie l'appariement des priors et des scenarios (doit etre ~0 ecart)
print("\n--- Verification de l'appariement des priors ---")
for c in param_cols:
    d = (diyabc[c] - msprime[c]).abs()
    print(f"{c:8s} ecart max = {d.max():.6f}")

print("\nRepartition des scenarios (diyabc):")
print(diyabc["scenario"].value_counts().sort_index())
print("Repartition des scenarios (msprime):")
print(msprime["scenario"].value_counts().sort_index())

# 2. Compare les statistiques
print(f"\n--- Comparaison des {len(stat_cols)} statistiques ---")
rows = []
for c in stat_cols:
    d = diyabc[c].dropna()
    m = msprime[c].dropna()
    if len(d) == 0 or len(m) == 0:
        continue
    mean_d, mean_m = d.mean(), m.mean()
    rel_diff = abs(mean_m - mean_d) / abs(mean_d) if mean_d != 0 else float("nan")
    ks_stat, ks_p = ss.ks_2samp(d, m)
    rows.append((c, mean_d, mean_m, rel_diff, ks_p))

df = pd.DataFrame(
    rows, columns=["stat", "mean_diyabc", "mean_msprime", "rel_diff", "ks_p"]
)
n_sig = (df["ks_p"] < 0.05).sum()
print(f"{n_sig} / {len(df)} statistiques avec KS p<0.05")
print()
print(df.sort_values("rel_diff", ascending=False).head(30).to_string(index=False))
