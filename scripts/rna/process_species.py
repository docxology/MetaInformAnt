#!/usr/bin/env python3
"""Run the current Amalgkit 0.16.60 workflow for one configured species.

This is a convenience wrapper around the same orchestrator used by
``run_all_species.py``.  It intentionally accepts only current workflow
parameters so single-species runs cannot drift onto a separate downloader,
index contract, or quantification implementation.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.cohort import (  # noqa: E402
    CohortError,
    add_producer_cohort_arguments,
    default_max_sample_gb,
    producer_cohort_from_args,
    select_cohort_species,
    verify_cohort_configs,
)
from metainformant.rna.engine.species import configured_data_root, discover_species_config_names  # noqa: E402
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
if not DEFAULT_CONFIG_DIR.is_dir():
    DEFAULT_CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, help="Configured species identifier, e.g. apis_mellifera")
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--data-root", type=Path, help="External Amalgkit data root")
    parser.add_argument(
        "--max-gb",
        type=float,
        default=default_max_sample_gb(),
        help="Max sample size in GB (default: AMALGKIT_MAX_SAMPLE_GB, then PIPELINE_MAX_GB, then 1000)",
    )
    parser.add_argument("--workers", type=int, default=int(os.environ.get("PIPELINE_WORKERS", 4)))
    parser.add_argument("--threads", type=int, default=int(os.environ.get("PIPELINE_THREADS", 6)))
    parser.add_argument(
        "--fastq-threads",
        type=int,
        default=None,
        help="Threads for fasterq-dump fallback (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--fastq-slots",
        type=int,
        default=None,
        help="Maximum concurrent fasterq-dump fallbacks (default: up to four; limits external-disk contention)",
    )
    parser.add_argument(
        "--compression-threads",
        type=int,
        default=None,
        help="Threads per pigz compression process (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--validation-slots",
        type=int,
        default=None,
        help="Concurrent full local FASTQ validations (default: up to four; limits volume contention)",
    )
    parser.add_argument(
        "--max-in-flight",
        type=int,
        default=None,
        help="Maximum submitted sample tasks (default: two times active workers)",
    )
    parser.add_argument("--dry-run", action="store_true", help="List the selected species without executing")
    add_producer_cohort_arguments(parser)
    args = parser.parse_args()
    args.parser = parser
    return args


def main() -> int:
    args = parse_args()
    config_dir = args.config_dir.expanduser().resolve()
    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root()

    try:
        cohort = producer_cohort_from_args(args)
        if cohort is not None:
            (species,) = select_cohort_species(cohort, [args.species])
            verify_cohort_configs(cohort, config_dir, [species])
            wanted = cohort.entry(species).config_name
    except CohortError as error:
        args.parser.error(str(error))
    if cohort is None:
        config_names = discover_species_config_names(config_dir)
        wanted = f"amalgkit_{args.species}.yaml"
        if wanted not in config_names:
            raise SystemExit(f"no runnable configuration found for {args.species!r} under {config_dir}")

    print(f"Species: {args.species}")
    print(f"Config directory: {config_dir}")
    print(f"Data root: {data_root}")
    if cohort is not None:
        print(f"Frozen cohort: {cohort.cohort_id} (definition sha256 {cohort.definition_sha256})")
    else:
        print(f"Metadata refresh (NCBI): {'allowed' if args.refresh_metadata else 'refused'}")
    print(f"Max sample size: {args.max_gb} GB | workers: {args.workers} | threads: {args.threads}")
    if args.dry_run:
        return 0

    orchestrator = StreamingPipelineOrchestrator(
        config_dir=config_dir,
        log_dir=data_root / "logs",
        db_path=data_root / "pipeline_progress.db",
        cohort=cohort,
        allow_metadata_refresh=args.refresh_metadata,
    )
    orchestrator.run_all(
        [wanted],
        args.max_gb,
        args.workers,
        args.threads,
        fasterq_threads=args.fastq_threads,
        fasterq_slots=args.fastq_slots,
        compression_threads=args.compression_threads,
        validation_slots=args.validation_slots,
        max_in_flight=args.max_in_flight,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
