#!/usr/bin/env python3
"""
Stage 99 — Synthetic Data Generator.

Creates realistic synthetic CSV and metadata files in ``data/raw/`` to enable
end-to-end pipeline testing without depending on external data sources.

Generated files:
- ``data/raw/samples_A.csv`` — 200-row synthetic measurement dataset
- ``data/raw/samples_B.csv`` — 150-row synthetic measurement dataset
- ``data/raw/metadata.yaml`` — provenance record for the generated data

Follows the MetaInformAnt Thin Orchestration Pattern:
- Data-generation and provenance-metadata logic lives in the
  ``template_bioinformatics_project`` library; this script only handles
  argparse and delegation.

Usage::

    uv run scripts/99_create_synthetic_data.py
    uv run scripts/99_create_synthetic_data.py --n-samples 500 --seed 42
    uv run scripts/99_create_synthetic_data.py --config config/default.yaml
"""

import sys
import time
import argparse
from pathlib import Path

# Path bootstrap: make the project's src/ library importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from template_bioinformatics_project import synthetic  # noqa: E402
from template_bioinformatics_project.common import load_optional_config  # noqa: E402


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate synthetic raw data for pipeline testing.",
    )
    parser.add_argument("--config", default="config/default.yaml",
                        help="Config YAML (used to resolve data/raw/ path)")
    parser.add_argument("--n-samples", type=int, default=200,
                        help="Rows in the first dataset (default: 200)")
    parser.add_argument("--n-features", type=int, default=8,
                        help="Number of numeric feature columns (default: 8)")
    parser.add_argument("--seed", type=int, default=2026,
                        help="Random seed for reproducibility (default: 2026)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing synthetic files")
    args = parser.parse_args()

    config = load_optional_config(args.config)
    raw_dir = synthetic.get_raw_dir(config)

    t0 = time.perf_counter()
    synthetic.run(
        raw_dir,
        n_samples=args.n_samples,
        n_features=args.n_features,
        seed=args.seed,
        force=args.force,
    )

    print(f"\n✅ Synthetic data generated in {time.perf_counter() - t0:.2f} s")
    print(f"   Raw data dir: {raw_dir.resolve()}")


if __name__ == "__main__":
    main()
