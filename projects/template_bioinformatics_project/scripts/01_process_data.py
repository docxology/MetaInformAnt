#!/usr/bin/env python3
"""
Stage 1 — Data Processing.

Reads raw data files from ``data/raw/``, applies configurable filtering and
normalisation, and writes cleaned outputs to ``data/processed/``.

Follows the MetaInformAnt Thin Orchestration Pattern:
- All paths from ``config/default.yaml`` — no hardcoded strings.
- Idempotent: skips work when all expected outputs already exist.
- Structured logging to ``logs/01_process_data.log``.

Usage::

    uv run scripts/01_process_data.py --config config/default.yaml
    uv run scripts/01_process_data.py --config config/default.yaml --force
"""

import sys
import time
import logging
import argparse
from pathlib import Path

import yaml
import numpy as np
import pandas as pd


# ── Logging ────────────────────────────────────────────────────────────────────

def setup_logging(config: dict) -> logging.Logger:
    """Configure dual file + console logging from config."""
    log_dir = Path(config["paths"]["logs"])
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "01_process_data.log"

    logging.basicConfig(
        level=getattr(logging, config["logging"]["level"], logging.INFO),
        format=config["logging"]["format"],
        handlers=[
            logging.FileHandler(log_file, mode="a"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("stage_01")


# ── Config ─────────────────────────────────────────────────────────────────────

def load_config(config_path: str | Path) -> dict:
    """Load and return the YAML configuration file."""
    config_path = Path(config_path)
    if not config_path.is_file():
        print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open() as fh:
        return yaml.safe_load(fh)


# ── Processing Logic ────────────────────────────────────────────────────────────

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


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 1 — Process raw data into cleaned, normalised form.",
    )
    parser.add_argument("--config", default="config/default.yaml",
                        help="Path to config YAML (default: config/default.yaml)")
    parser.add_argument("--force", action="store_true",
                        help="Reprocess even if output already exists")
    args = parser.parse_args()

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 1: Data Processing ── config=%s", args.config)

    try:
        process_data(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    elapsed = time.perf_counter() - t0
    logger.info("Stage 1 complete in %.2f s", elapsed)


if __name__ == "__main__":
    main()
