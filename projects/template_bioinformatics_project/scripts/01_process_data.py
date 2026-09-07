#!/usr/bin/env python3
"""
Stage 1 — Data Processing.

Reads raw data files from ``data/raw/``, applies configurable filtering and
normalisation, and writes cleaned outputs to ``data/processed/``.

Follows the MetaInformAnt Thin Orchestration Pattern:
- All paths from ``config/default.yaml`` — no hardcoded strings.
- Idempotent: skips work when all expected outputs already exist.
- Structured logging to ``logs/01_process_data.log``.
- All processing logic lives in the ``template_bioinformatics_project`` library.

Usage::

    uv run scripts/01_process_data.py --config config/default.yaml
    uv run scripts/01_process_data.py --config config/default.yaml --force
"""

import sys
import time
import logging
import argparse
from pathlib import Path

# Path bootstrap: make the project's src/ library importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from template_bioinformatics_project import processing  # noqa: E402


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

    from template_bioinformatics_project.common import load_config

    config = load_config(args.config)
    logger = setup_logging(config)

    t0 = time.perf_counter()
    logger.info("── Stage 1: Data Processing ── config=%s", args.config)

    try:
        processing.process_data(config, logger, force=args.force)
    except Exception as exc:
        logger.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)

    elapsed = time.perf_counter() - t0
    logger.info("Stage 1 complete in %.2f s", elapsed)


if __name__ == "__main__":
    main()
