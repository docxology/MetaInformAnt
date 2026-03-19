#!/usr/bin/env python3
"""
Stage 99 — Synthetic Data Generator.

Creates realistic synthetic CSV and metadata files in ``data/raw/`` to enable
end-to-end pipeline testing without depending on external data sources.

Generated files:
- ``data/raw/samples_A.csv`` — 200-row synthetic measurement dataset
- ``data/raw/samples_B.csv`` — 150-row synthetic measurement dataset
- ``data/raw/metadata.yaml`` — provenance record for the generated data

Usage::

    uv run scripts/99_create_synthetic_data.py
    uv run scripts/99_create_synthetic_data.py --n-samples 500 --seed 42
    uv run scripts/99_create_synthetic_data.py --config config/default.yaml
"""

import sys
import time
import random
import argparse
import datetime
from pathlib import Path

import yaml
import numpy as np


# ── Config ─────────────────────────────────────────────────────────────────────

def load_config(config_path: str | Path) -> dict:
    config_path = Path(config_path)
    if not config_path.is_file():
        return {}
    with config_path.open() as fh:
        return yaml.safe_load(fh) or {}


def get_raw_dir(config: dict) -> Path:
    return Path(config.get("paths", {}).get("data_raw", "data/raw/"))


# ── Synthetic generation ───────────────────────────────────────────────────────

def generate_sample_dataframe(
    n_rows: int,
    n_features: int,
    rng: np.random.Generator,
    missing_fraction: float = 0.05,
) -> "list[dict]":
    """
    Generate ``n_rows`` synthetic observations with ``n_features`` numeric
    features plus a categorical ``group`` column.

    A small fraction of values is randomly set to NaN to simulate real data.
    """
    groups = ["control", "treatment_A", "treatment_B"]

    rows = []
    for i in range(n_rows):
        row: dict = {"sample_id": f"S{i:04d}", "group": rng.choice(groups)}
        for j in range(n_features):
            val = float(rng.normal(loc=j * 0.5, scale=1.0 + j * 0.1))
            if rng.random() < missing_fraction:
                val = float("nan")
            row[f"feature_{j + 1:02d}"] = val
        rows.append(row)
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    """Write a list of row dicts to a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    headers = list(rows[0].keys())
    lines = [",".join(headers)]
    for row in rows:
        fields = []
        for h in headers:
            v = row.get(h)
            if v is None or (isinstance(v, float) and v != v):  # NaN check
                fields.append("")
            else:
                fields.append(str(v))
        lines.append(",".join(fields))
    path.write_text("\n".join(lines) + "\n")


def write_metadata(raw_dir: Path, datasets: list[dict]) -> None:
    """Write a provenance YAML file describing the generated data."""
    metadata = {
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "generator": "scripts/99_create_synthetic_data.py",
        "purpose": "Synthetic data for end-to-end pipeline testing",
        "datasets": datasets,
        "notes": (
            "This data is entirely synthetic. "
            "Do not use for scientific conclusions."
        ),
    }
    meta_path = raw_dir / "metadata.yaml"
    with meta_path.open("w") as fh:
        yaml.dump(metadata, fh, sort_keys=False)
    print(f"  metadata → {meta_path}")


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

    config = load_config(args.config)
    raw_dir = get_raw_dir(config)

    # ── Idempotency ────────────────────────────────────────────────────────────
    existing = list(raw_dir.glob("samples_*.csv"))
    if existing and not args.force:
        print(f"Synthetic data already present ({len(existing)} file(s)); use --force to regenerate.")
        sys.exit(0)

    rng = np.random.default_rng(args.seed)
    t0 = time.perf_counter()

    datasets_meta = []

    # Dataset A
    rows_a = generate_sample_dataframe(args.n_samples, args.n_features, rng)
    path_a = raw_dir / "samples_A.csv"
    write_csv(rows_a, path_a)
    print(f"  samples_A → {path_a}  ({len(rows_a)} rows × {args.n_features + 2} cols)")
    datasets_meta.append({"file": "samples_A.csv", "n_rows": len(rows_a), "n_features": args.n_features})

    # Dataset B (smaller, different seed offset)
    n_b = max(50, args.n_samples // 2)
    rows_b = generate_sample_dataframe(n_b, args.n_features, rng)
    path_b = raw_dir / "samples_B.csv"
    write_csv(rows_b, path_b)
    print(f"  samples_B → {path_b}  ({len(rows_b)} rows × {args.n_features + 2} cols)")
    datasets_meta.append({"file": "samples_B.csv", "n_rows": len(rows_b), "n_features": args.n_features})

    write_metadata(raw_dir, datasets_meta)

    print(f"\n✅ Synthetic data generated in {time.perf_counter() - t0:.2f} s")
    print(f"   Raw data dir: {raw_dir.resolve()}")


if __name__ == "__main__":
    main()
