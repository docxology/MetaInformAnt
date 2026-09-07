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
- All analysis logic lives in the ``template_bioinformatics_project`` library.

Usage::

    uv run scripts/02_analyze_results.py --config config/default.yaml
    uv run scripts/02_analyze_results.py --config config/default.yaml --force
"""

import sys
import time
import logging
import argparse
from pathlib import Path

# Path bootstrap: make the project's src/ library importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from template_bioinformatics_project import analysis  # noqa: E402


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


# ── Entry Point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stage 2 — Downstream statistical analysis of processed data.",
    )
    parser.add_argument("--config", default="config/default.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Rerun even if outputs already exist")
    args = parser.parse_args()

    from template_bioinformatics_project.common import load_config

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 2: Results Analysis ── config=%s", args.config)

    try:
        analysis.run_analysis(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    logger.info("Stage 2 complete in %.2f s", time.perf_counter() - t0)


if __name__ == "__main__":
    main()
