from pathlib import Path

from bridge import pipeline
from bridge.ancestry_simulation import build_local_param_microsat_per_locus

REFERENCE_DIR = Path("reference/toy_example1_ms_modified")

header_text = pipeline.read_header_text(REFERENCE_DIR)

params, values = build_local_param_microsat_per_locus(header_text, seed=42)
print(len({Pgeom for _, Pgeom, _ in params.values()}))
