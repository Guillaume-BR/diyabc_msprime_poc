"""
Contrôle d'identité avant/après une optimisation : mêmes particules MicroSat
(mêmes seeds), mêmes statistiques, EXACTEMENT.

Une optimisation ne doit rien changer aux sorties à seed égal (CLAUDE.md,
« Every optimization touches code validated against DIYABC »). Ce script
enregistre les statistiques de N particules dans un JSON (les floats y sont
écrits avec leur repr, donc sans perte), puis les compare à celles d'un autre
état du code.

Usage, depuis la racine du dépôt :

    # 1. Référence, sur le code AVANT modification
    git stash
    python3 -m scripts.compare_particles_microsat save ref.json
    git stash pop

    # 2. Comparaison, sur le code APRÈS modification
    python3 -m scripts.compare_particles_microsat compare ref.json

Arguments facultatifs après le fichier : dossier_reference (défaut
reference/toy_example1_ms_modified) et nombre de particules (défaut 30).
Les seeds sont 1..N. La comparaison relance les mêmes seeds sur le même
dossier : le dossier et N sont relus du fichier de référence.

Piège : la dichotomie d'AML bascule sur le SIGNE d'une pente ; un écart de
1e-13 sur une vraisemblance ne change presque jamais une sortie, mais
peut le faire. D'où N grand plutôt qu'un test à tolérance.
"""

import json
import sys
from pathlib import Path

from bridge import pipeline
from bridge.header_dataclasses import MicrosatReplayContext
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import (
    allele_bounds_per_locus,
    observed_count_sample,
    observed_microsatellites,
    parse_sex_ratio,
)

SCENARIO_INDEX = 1


def build_context(reference_dir: Path) -> MicrosatReplayContext:
    """Construit le contexte MicroSat d'un dossier de référence.

    Args:
        reference_dir: Dossier contenant le header et le fichier .mss.

    Returns:
        Le MicrosatReplayContext correspondant.
    """
    header_text = pipeline.read_header_text(reference_dir)
    mss_path = reference_dir / header_text.splitlines()[0].strip()
    list_loci = parse_loci_description(header_text)
    observed = observed_microsatellites(mss_path, list_loci)
    return MicrosatReplayContext(
        header_text=header_text,
        mss_path=mss_path,
        list_loci=list_loci,
        microsat_observed=observed,
        bounds_per_locus=allele_bounds_per_locus(observed, list_loci),
        samples_default=observed_count_sample(mss_path),
        sex_ratio=parse_sex_ratio(mss_path),
    )


def run_particles(reference_dir: Path, n_particles: int) -> dict[str, dict[str, float]]:
    """Simule les particules de seeds 1..n_particles.

    Args:
        reference_dir: Dossier de référence.
        n_particles: Nombre de particules.

    Returns:
        {seed: {nom_stat: valeur}}, clés en str pour le JSON.
    """
    context = build_context(reference_dir)
    results = {}
    for seed in range(1, n_particles + 1):
        stats, _, _ = pipeline.compute_summary_statistics_microsat(
            context, SCENARIO_INDEX, seed=seed
        )
        results[str(seed)] = {k: float(v) for k, v in stats.items()}
        print(f"  particule {seed}/{n_particles}", end="\r", flush=True)
    print()
    return results


def main() -> None:
    """Point d'entrée : sous-commandes `save` et `compare`."""
    if len(sys.argv) < 3 or sys.argv[1] not in ("save", "compare"):
        sys.exit(__doc__)
    mode, ref_path = sys.argv[1], Path(sys.argv[2])

    if mode == "save":
        reference_dir = Path(
            sys.argv[3] if len(sys.argv) > 3 else "reference/toy_example1_ms_modified"
        )
        n_particles = int(sys.argv[4]) if len(sys.argv) > 4 else 30
        payload = {
            "reference_dir": str(reference_dir),
            "n_particles": n_particles,
            "results": run_particles(reference_dir, n_particles),
        }
        ref_path.write_text(json.dumps(payload))
        print(f"Référence enregistrée : {ref_path} ({n_particles} particules)")
        return

    payload = json.loads(ref_path.read_text())
    new = run_particles(Path(payload["reference_dir"]), payload["n_particles"])
    n_diff_stats, particles_with_diff = 0, set()
    worst = (0.0, None)
    for seed, ref_stats in payload["results"].items():
        if ref_stats.keys() != new[seed].keys():
            sys.exit(f"Seed {seed} : les colonnes diffèrent.")
        for name, ref_value in ref_stats.items():
            value = new[seed][name]
            if value != ref_value and not (value != value and ref_value != ref_value):
                n_diff_stats += 1
                particles_with_diff.add(seed)
                gap = abs(value - ref_value)
                if gap > worst[0]:
                    worst = (gap, (seed, name, ref_value, value))
    n_cells = sum(len(s) for s in payload["results"].values())
    print(
        f"{payload['n_particles']} particules, {n_cells} valeurs comparées : "
        f"{n_diff_stats} différentes sur {len(particles_with_diff)} particule(s)."
    )
    if n_diff_stats:
        gap, (seed, name, ref_value, value) = worst
        print(
            f"Pire écart : seed {seed}, {name} : {ref_value!r} -> {value!r} ({gap:.3e})"
        )
        sys.exit(1)
    print("IDENTIQUE, valeur par valeur.")


if __name__ == "__main__":
    main()
