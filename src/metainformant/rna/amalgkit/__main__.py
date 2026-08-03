"""Run the current streaming Amalgkit producer for one configured species.

The producer owns metadata, acquisition, integration, and quantification.
Merge, within-species filtering, finalization, and sanity are performed by the
project's lock-owned downstream checkpoint runner after the producer stops.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from metainformant.core.utils import logging

REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
if not DEFAULT_CONFIG_DIR.is_dir():
    DEFAULT_CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"

logger = logging.get_logger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="metainformant.rna.amalgkit",
        description="Run the current streaming Amalgkit producer for one species",
    )
    parser.add_argument("--config", "-c", type=Path, required=True, help="Path to a current species YAML configuration")
    parser.add_argument("--data-root", type=Path, help="External Amalgkit data root")
    parser.add_argument("--max-gb", type=float, default=float(os.environ.get("PIPELINE_MAX_GB", "50")))
    parser.add_argument("--workers", type=int, default=int(os.environ.get("PIPELINE_WORKERS", "4")))
    parser.add_argument("--threads", type=int, default=int(os.environ.get("PIPELINE_THREADS", "6")))
    parser.add_argument("--quant-slots", type=int, default=None)
    parser.add_argument("--fastq-threads", type=int, default=None)
    parser.add_argument("--fastq-slots", type=int, default=None)
    parser.add_argument("--compression-threads", type=int, default=None)
    parser.add_argument("--validation-slots", type=int, default=None)
    parser.add_argument("--max-in-flight", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true", help="Resolve the producer inputs without executing")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Resolve one current config and run its streaming producer."""

    args = _parser().parse_args(argv)
    config_path = args.config.expanduser().resolve()
    if not config_path.is_file():
        logger.error("Configuration file not found: %s", config_path)
        return 2
    if not config_path.name.startswith("amalgkit_") or config_path.name in {
        "amalgkit_template.yaml",
        "amalgkit_test.yaml",
        "amalgkit_cross_species.yaml",
    }:
        logger.error("Configuration is not a runnable species config: %s", config_path)
        return 2

    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())

    from metainformant.rna.engine.species import configured_data_root, species_name_from_config
    from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

    config_dir = config_path.parent
    data_root = configured_data_root()
    species_config = config_path.name
    species = species_name_from_config(species_config)

    print(f"Species: {species}")
    print(f"Config directory: {config_dir}")
    print(f"Data root: {data_root}")
    print(f"Max sample size: {args.max_gb} GB | workers: {args.workers} | threads: {args.threads}")
    if args.dry_run:
        return 0

    orchestrator = StreamingPipelineOrchestrator(
        config_dir=config_dir,
        log_dir=data_root / "logs",
        db_path=data_root / "pipeline_progress.db",
    )
    orchestrator.run_all(
        [species_config],
        args.max_gb,
        args.workers,
        args.threads,
        quant_slots=args.quant_slots,
        fasterq_threads=args.fastq_threads,
        fasterq_slots=args.fastq_slots,
        compression_threads=args.compression_threads,
        validation_slots=args.validation_slots,
        max_in_flight=args.max_in_flight,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
