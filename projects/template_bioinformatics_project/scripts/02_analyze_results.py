#!/usr/bin/env python3
"""
Stage 2 — Downstream Statistical Analysis.

Reads processed data from ``data/processed/processed_data.csv``, computes
summary statistics, and optionally runs PCA or correlation analysis as
configured.  Outputs summary tables to ``results/tables/``.

Follows the MetaInformAnt Thin Orchestration Pattern:
- All paths from ``config/default.yaml``.
- Idempotent: skips if summary table already present.
- Structured logging to ``logs/02_analyze_results.log``.

Usage::

    uv run scripts/02_analyze_results.py --config config/default.yaml
    uv run scripts/02_analyze_results.py --config config/default.yaml --force
"""

import sys
import time
import json
import logging
import argparse
from pathlib import Path

import yaml
import numpy as np
import pandas as pd


# ── Logging ────────────────────────────────────────────────────────────────────

def setup_logging(config: dict) -> logging.Logger:
    log_dir = Path(config["paths"]["logs"])
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "02_analyze_results.log"
    logging.basicConfig(
        level=getattr(logging, config["logging"]["level"], logging.INFO),
        format=config["logging"]["format"],
        handlers=[
            logging.FileHandler(log_file, mode="a"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("stage_02")


# ── Config ─────────────────────────────────────────────────────────────────────

def load_config(config_path: str | Path) -> dict:
    config_path = Path(config_path)
    if not config_path.is_file():
        print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open() as fh:
        return yaml.safe_load(fh)


# ── Analysis Functions ─────────────────────────────────────────────────────────

def compute_summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Return descriptive statistics for all numeric columns."""
    numeric = df.select_dtypes(include="number")
    stats = numeric.describe().T
    stats["cv"] = stats["std"] / stats["mean"].abs()  # coefficient of variation
    return stats


def compute_correlation_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Return pairwise Pearson correlation matrix for numeric columns."""
    return df.select_dtypes(include="number").corr()


def compute_pca_summary(df: pd.DataFrame, n_components: int) -> dict:
    """
    Compute PCA via SVD (no sklearn dependency) and return explained variance.

    Returns
    -------
    dict with keys: ``n_components``, ``explained_variance_ratio``, ``total_variance_explained``
    """
    numeric = df.select_dtypes(include="number").dropna(axis=1)
    if numeric.shape[1] < 2:
        return {"error": "Insufficient numeric columns for PCA"}

    X = numeric.values
    X = X - X.mean(axis=0)
    _, s, _ = np.linalg.svd(X, full_matrices=False)
    variance = s ** 2 / (len(X) - 1)
    total = variance.sum()
    n = min(n_components, len(variance))
    ratio = (variance[:n] / total).tolist()
    return {
        "n_components": n,
        "explained_variance_ratio": ratio,
        "total_variance_explained": float(sum(ratio)),
    }


# ── Orchestration ──────────────────────────────────────────────────────────────

def run_analysis(config: dict, logger: logging.Logger, force: bool = False) -> None:
    """Load processed data and run the configured analysis method."""
    processed_dir = Path(config["paths"]["data_processed"])
    tables_dir = Path(config["paths"]["results_tables"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    input_file = processed_dir / "processed_data.csv"
    output_summary = tables_dir / "summary_statistics.csv"
    output_meta = tables_dir / "analysis_metadata.json"

    # ── Idempotency ────────────────────────────────────────────────────────────
    if output_summary.exists() and not force:
        logger.info("Summary already exists, skipping (--force to rerun): %s", output_summary)
        return

    if not input_file.exists():
        logger.error("Input not found — run Stage 1 first: %s", input_file)
        sys.exit(1)

    data = pd.read_csv(input_file)
    logger.info("Loaded processed data: %d rows × %d cols", *data.shape)

    method = config["analysis"]["method"]
    logger.info("Analysis method: %s", method)

    # ── Summary statistics (always computed) ───────────────────────────────────
    summary = compute_summary_statistics(data)
    summary.to_csv(output_summary)
    logger.info("Summary statistics → %s", output_summary)

    metadata: dict = {
        "method": method,
        "n_rows": len(data),
        "n_numeric_cols": int(data.select_dtypes(include="number").shape[1]),
    }

    # ── Optional correlation ───────────────────────────────────────────────────
    if method in ("correlation", "summary"):
        corr = compute_correlation_matrix(data)
        corr_path = tables_dir / "correlation_matrix.csv"
        corr.to_csv(corr_path)
        logger.info("Correlation matrix → %s", corr_path)
        metadata["correlation_computed"] = True

    # ── Optional PCA ──────────────────────────────────────────────────────────
    if method == "pca":
        n_components = config["analysis"]["n_components"]
        pca_result = compute_pca_summary(data, n_components)
        pca_path = tables_dir / "pca_summary.json"
        with pca_path.open("w") as fh:
            json.dump(pca_result, fh, indent=2)
        logger.info("PCA summary → %s  (%.1f%% variance explained)",
                    pca_path, pca_result.get("total_variance_explained", 0) * 100)
        metadata["pca"] = pca_result

    # ── Write metadata ─────────────────────────────────────────────────────────
    with output_meta.open("w") as fh:
        json.dump(metadata, fh, indent=2)
    logger.info("Analysis metadata → %s", output_meta)


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 2 — Downstream statistical analysis of processed data.",
    )
    parser.add_argument("--config", default="config/default.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Rerun even if outputs already exist")
    args = parser.parse_args()

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 2: Results Analysis ── config=%s", args.config)

    try:
        run_analysis(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    logger.info("Stage 2 complete in %.2f s", time.perf_counter() - t0)


if __name__ == "__main__":
    main()
