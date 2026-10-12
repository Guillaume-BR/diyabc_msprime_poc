import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from bridge import pipeline
from bridge.header_dataclasses import MicrosatReplayContext
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    allele_bounds_per_locus,
    observed_count_sample,
    observed_microsatellites,
    parse_sex_ratio,
)

WORK_DIR = Path(__file__).parent

REFERENCE_DIR = Path(
    sys.argv[1] if len(sys.argv) > 1 else "reference/toy_example1_ms_modified"
)
N_PARTICLES = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SCENARIO_INDEX = 1

elapsed: dict[str, float] = defaultdict(float)
calls: dict[str, int] = defaultdict(int)


header_text = pipeline.read_header_text(REFERENCE_DIR)
mss_path = REFERENCE_DIR / header_text.splitlines()[0].strip()
list_loci = parse_loci_description(header_text)
microsat_observed = observed_microsatellites(mss_path, list_loci)
context = MicrosatReplayContext(
    header_text=header_text,
    mss_path=mss_path,
    list_loci=list_loci,
    microsat_observed=microsat_observed,
    bounds_per_locus=allele_bounds_per_locus(microsat_observed, list_loci),
    samples_default=observed_count_sample(mss_path),
    sex_ratio=parse_sex_ratio(mss_path),
)


results = pipeline.compute_summary_statistics_microsat(
    context, SCENARIO_INDEX, seed=1000
)
list_aml = [k for k in results[0] if k.startswith("AML")]
result_aml = {k: results[0][k] for k in list_aml}
df_result_aml = pd.DataFrame([result_aml])
df_result_aml.to_csv(WORK_DIR / f"result_aml_{N_PARTICLES + 1}.csv", index=False)
print(result_aml)
