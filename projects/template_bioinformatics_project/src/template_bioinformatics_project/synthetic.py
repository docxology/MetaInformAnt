"""Stage 99 logic — synthetic data generation and provenance metadata."""

from __future__ import annotations

import datetime
import sys
from datetime import UTC
from pathlib import Path

import numpy as np
import yaml


def get_raw_dir(config: dict) -> Path:
    return Path(config.get("paths", {}).get("data_raw", "data/raw/"))


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
        "generated_at": datetime.datetime.now(UTC).isoformat(),
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


def run(
    raw_dir: Path,
    n_samples: int,
    n_features: int,
    seed: int,
    force: bool = False,
) -> None:
    """Generate the synthetic datasets and provenance metadata in ``raw_dir``."""
    # ── Idempotency ────────────────────────────────────────────────────────────
    existing = list(raw_dir.glob("samples_*.csv"))
    if existing and not force:
        print(f"Synthetic data already present ({len(existing)} file(s)); use --force to regenerate.")
        sys.exit(0)

    rng = np.random.default_rng(seed)

    datasets_meta = []

    # Dataset A
    rows_a = generate_sample_dataframe(n_samples, n_features, rng)
    path_a = raw_dir / "samples_A.csv"
    write_csv(rows_a, path_a)
    print(f"  samples_A → {path_a}  ({len(rows_a)} rows × {n_features + 2} cols)")
    datasets_meta.append({"file": "samples_A.csv", "n_rows": len(rows_a), "n_features": n_features})

    # Dataset B (smaller, different seed offset)
    n_b = max(50, n_samples // 2)
    rows_b = generate_sample_dataframe(n_b, n_features, rng)
    path_b = raw_dir / "samples_B.csv"
    write_csv(rows_b, path_b)
    print(f"  samples_B → {path_b}  ({len(rows_b)} rows × {n_features + 2} cols)")
    datasets_meta.append({"file": "samples_B.csv", "n_rows": len(rows_b), "n_features": n_features})

    write_metadata(raw_dir, datasets_meta)
