#!/usr/bin/env python3
"""
Stage 3 — Visualisation.

Reads processed data from ``data/processed/processed_data.csv`` and summary
statistics from ``results/tables/`` to generate publication-quality figures
saved to ``results/figures/``.

Figures produced:
- ``distribution_grid.{fmt}`` — per-column histogram grid
- ``correlation_heatmap.{fmt}`` — heatmap of pairwise correlations (if available)

Follows the MetaInformAnt Thin Orchestration Pattern:
- All paths and visual parameters come from ``config/default.yaml``.
- Idempotent: skips if all expected figures already exist.
- Structured logging to ``logs/03_visualize.log``.

Usage::

    uv run scripts/03_visualize.py --config config/default.yaml
    uv run scripts/03_visualize.py --config config/default.yaml --force
"""

import sys
import time
import logging
import argparse
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns


# ── Logging ────────────────────────────────────────────────────────────────────

def setup_logging(config: dict) -> logging.Logger:
    log_dir = Path(config["paths"]["logs"])
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "03_visualize.log"
    logging.basicConfig(
        level=getattr(logging, config["logging"]["level"], logging.INFO),
        format=config["logging"]["format"],
        handlers=[
            logging.FileHandler(log_file, mode="a"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("stage_03")


# ── Config ─────────────────────────────────────────────────────────────────────

def load_config(config_path: str | Path) -> dict:
    config_path = Path(config_path)
    if not config_path.is_file():
        print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open() as fh:
        return yaml.safe_load(fh)


# ── Plot Functions ─────────────────────────────────────────────────────────────

def plot_distribution_grid(
    data: pd.DataFrame,
    output_path: Path,
    vcfg: dict,
) -> None:
    """
    Create a grid of histograms, one per numeric column.

    Parameters
    ----------
    data : DataFrame
    output_path : destination file path
    vcfg : visualization config sub-dict
    """
    numeric = data.select_dtypes(include="number").drop(columns=["_source_file"], errors="ignore")
    cols = numeric.columns.tolist()
    if not cols:
        return

    ncols = min(4, len(cols))
    nrows = (len(cols) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(vcfg["figure_width"] * ncols / 4, vcfg["figure_height"] * nrows / 2),
        squeeze=False,
    )
    palette = sns.color_palette(vcfg["color_palette"], len(cols))

    for idx, col in enumerate(cols):
        r, c = divmod(idx, ncols)
        ax = axes[r][c]
        ax.hist(numeric[col].dropna(), bins=30, color=palette[idx], edgecolor="white", linewidth=0.4)
        ax.set_title(col, fontsize=9, pad=4)
        ax.set_xlabel("Value", fontsize=7)
        ax.set_ylabel("Count", fontsize=7)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for idx in range(len(cols), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle("Feature Distributions", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, dpi=vcfg["dpi"], bbox_inches="tight")
    plt.close(fig)


def plot_correlation_heatmap(
    corr_path: Path,
    output_path: Path,
    vcfg: dict,
) -> None:
    """Create a seaborn heatmap of the pre-computed correlation matrix."""
    corr = pd.read_csv(corr_path, index_col=0)
    if corr.empty:
        return

    fig, ax = plt.subplots(figsize=(vcfg["figure_width"], vcfg["figure_height"]))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(
        corr,
        mask=mask,
        cmap="coolwarm",
        center=0,
        vmin=-1, vmax=1,
        annot=corr.shape[0] <= 12,
        fmt=".2f",
        linewidths=0.3,
        ax=ax,
    )
    ax.set_title("Feature Correlation Matrix", fontsize=11)
    fig.tight_layout()
    fig.savefig(output_path, dpi=vcfg["dpi"], bbox_inches="tight")
    plt.close(fig)


# ── Orchestration ──────────────────────────────────────────────────────────────

def run_visualize(config: dict, logger: logging.Logger, force: bool = False) -> None:
    """Generate all configured figures from processed data and summary tables."""
    figures_dir = Path(config["paths"]["results_figures"])
    tables_dir = Path(config["paths"]["results_tables"])
    processed_file = Path(config["paths"]["data_processed"]) / "processed_data.csv"
    vcfg = config["visualization"]
    fmt = vcfg["file_format"]

    figures_dir.mkdir(parents=True, exist_ok=True)

    dist_path = figures_dir / f"distribution_grid.{fmt}"
    corr_path = figures_dir / f"correlation_heatmap.{fmt}"
    corr_table = tables_dir / "correlation_matrix.csv"

    expected = [dist_path]
    if corr_table.exists():
        expected.append(corr_path)

    if all(p.exists() for p in expected) and not force:
        logger.info("All figures already exist, skipping (--force to regenerate)")
        return

    if not processed_file.exists():
        logger.error("Processed data not found — run Stage 1 first: %s", processed_file)
        sys.exit(1)

    data = pd.read_csv(processed_file)
    logger.info("Loaded processed data: %d rows × %d cols", *data.shape)

    # Distribution grid
    plot_distribution_grid(data, dist_path, vcfg)
    logger.info("Distribution grid → %s", dist_path)

    # Correlation heatmap (only if correlation matrix was generated in Stage 2)
    if corr_table.exists():
        plot_correlation_heatmap(corr_table, corr_path, vcfg)
        logger.info("Correlation heatmap → %s", corr_path)
    else:
        logger.info("No correlation matrix found; skipping heatmap (run Stage 2 first)")


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 3 — Generate visualisation figures.",
    )
    parser.add_argument("--config", default="config/default.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Regenerate figures even if they already exist")
    args = parser.parse_args()

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 3: Visualisation ── config=%s", args.config)

    try:
        run_visualize(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    logger.info("Stage 3 complete in %.2f s", time.perf_counter() - t0)


if __name__ == "__main__":
    main()
