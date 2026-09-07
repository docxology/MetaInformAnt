"""Stage 2 logic — summary statistics, correlation, and PCA analysis."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd


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
