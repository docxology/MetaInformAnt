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
- All figure-construction logic lives in the ``template_bioinformatics_project`` library.

Usage::

    uv run scripts/03_visualize.py --config config/default.yaml
    uv run scripts/03_visualize.py --config config/default.yaml --force
"""

import sys
import time
import logging
import argparse
from pathlib import Path

# Path bootstrap: make the project's src/ library importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from template_bioinformatics_project import visualization  # noqa: E402


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


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 3 — Generate visualisation figures.",
    )
    parser.add_argument("--config", default="config/default.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Regenerate figures even if they already exist")
    args = parser.parse_args()

    from template_bioinformatics_project.common import load_config

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 3: Visualisation ── config=%s", args.config)

    try:
        visualization.run_visualize(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    logger.info("Stage 3 complete in %.2f s", time.perf_counter() - t0)


if __name__ == "__main__":
    main()
