"""Stage 1 logic — ingest, filter, and normalise raw data."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd


def process_data(config: dict, logger: logging.Logger, force: bool = False) -> None:
    """
    Ingest raw CSV files from ``data/raw/`` and write processed output.

    Steps:
    1. Discover raw CSV files.
    2. Concatenate them into a unified DataFrame.
    3. Apply configurable filtering (missing-value threshold, count threshold).
    4. Optionally z-score normalise numeric columns.
    5. Write ``data/processed/processed_data.csv``.
    """
    raw_dir = Path(config["paths"]["data_raw"])
    processed_dir = Path(config["paths"]["data_processed"])
    output_file = processed_dir / "processed_data.csv"

    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    # ── Idempotency guard ──────────────────────────────────────────────────────
    if output_file.exists() and not force:
        logger.info("Output already exists, skipping (use --force to reprocess): %s", output_file)
        return

    # ── Discover input files ───────────────────────────────────────────────────
    raw_files = sorted(raw_dir.glob("*.csv"))
    if not raw_files:
        logger.warning("No CSV files found in %s — nothing to process.", raw_dir)
        return

    logger.info("Discovered %d raw file(s) in %s", len(raw_files), raw_dir)

    # ── Load ───────────────────────────────────────────────────────────────────
    frames = []
    for path in raw_files:
        df = pd.read_csv(path)
        df["_source_file"] = path.name
        frames.append(df)
        logger.debug("Loaded %s  (%d rows × %d cols)", path.name, *df.shape)

    data = pd.concat(frames, ignore_index=True)
    logger.info("Combined shape: %d rows × %d columns", *data.shape)

    # ── Filter: missing values ─────────────────────────────────────────────────
    max_missing = config["processing"]["missing_fraction_max"]
    numeric_cols = data.select_dtypes(include="number").columns.tolist()
    before = len(numeric_cols)
    missing_frac = data[numeric_cols].isnull().mean()
    keep = missing_frac[missing_frac <= max_missing].index.tolist()
    data = data[keep + [c for c in data.columns if c not in numeric_cols]]
    logger.info(
        "Missing-value filter (threshold=%.2f): kept %d / %d numeric columns",
        max_missing, len(keep), before,
    )

    # ── Filter: minimum sample count ──────────────────────────────────────────
    min_count = config["processing"]["min_sample_count"]
    initial_rows = len(data)
    data = data.dropna(thresh=min_count)
    logger.info(
        "Row-count filter (min_sample_count=%d): kept %d / %d rows",
        min_count, len(data), initial_rows,
    )

    # ── Normalise ─────────────────────────────────────────────────────────────
    if config["processing"]["normalize"]:
        numeric_now = data.select_dtypes(include="number").columns.tolist()
        for col in numeric_now:
            std = data[col].std()
            if std > 0:
                data[col] = (data[col] - data[col].mean()) / std
        logger.info("Applied z-score normalisation to %d numeric columns", len(numeric_now))

    # ── Write output ───────────────────────────────────────────────────────────
    data.to_csv(output_file, index=False)
    logger.info("Wrote processed data → %s  (%d rows × %d cols)", output_file, *data.shape)
