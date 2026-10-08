"""Run an immutable Amalgkit acquisition task manifest on this host."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from metainformant.rna.engine.acquisition_manifest import load_snapshot, load_task_selection, sha256_file
from metainformant.rna.engine.acquisition_snapshot import stage_worker_inputs
from metainformant.rna.engine.acquisition_worker import run_manifest


def build_parser() -> argparse.ArgumentParser:
    """Build the worker parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=18)
    parser.add_argument("--threads", type=int, default=18)
    parser.add_argument("--fastq-threads", type=int, default=1)
    parser.add_argument("--compression-threads", type=int, default=1)
    parser.add_argument("--compression-level", type=int, choices=range(1, 10), default=1)
    parser.add_argument("--ena-file-workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--validation-slots", type=int, default=4)
    parser.add_argument("--quant-slots", type=int)
    parser.add_argument("--fastq-slots", type=int)
    parser.add_argument("--max-in-flight", type=int)
    parser.add_argument("--task-selection", type=Path)
    parser.add_argument(
        "--stage-inputs",
        action="store_true",
        help="Copy frozen metadata/index inputs without replacing conflicting local files",
    )
    parser.add_argument(
        "--max-raw-bytes",
        type=int,
        help="Bound aggregate in-flight raw byte reservations; requires positive size evidence",
    )
    parser.add_argument("--durable-bucket")
    parser.add_argument("--durable-cohort")
    parser.add_argument("--profile", help="Optional AWS credential profile for durable output reuse/publication")
    parser.add_argument("--region")
    parser.add_argument(
        "--reclaim-raw-after-quant", action="store_true", help="Opt in to existing provenance-gated raw reclamation"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the cloud task worker."""

    args = build_parser().parse_args(argv)
    os.environ["AMALGKIT_PIPELINE_ENA_FILE_WORKERS"] = str(args.ena_file_workers)
    os.environ["AMALGKIT_PIPELINE_COMPRESSION_LEVEL"] = str(args.compression_level)
    if args.max_raw_bytes is not None:
        if args.max_raw_bytes <= 0:
            raise ValueError("max-raw-bytes must be positive")
        os.environ["AMALGKIT_CLOUD_MAX_RAW_BYTES"] = str(args.max_raw_bytes)
    for value, name in (
        (args.durable_bucket, "AMALGKIT_DURABLE_BUCKET"),
        (args.durable_cohort, "AMALGKIT_DURABLE_COHORT"),
        (args.profile, "AWS_PROFILE"),
        (args.region, "AWS_DEFAULT_REGION"),
    ):
        if value:
            os.environ[name] = value
    if args.reclaim_raw_after_quant:
        os.environ["AMALGKIT_RECLAIM_RAW_AFTER_QUANT"] = "yes"
    if args.stage_inputs:
        snapshot, tasks = load_snapshot(args.manifest)
        _, selected = load_task_selection(
            args.task_selection, snapshot, tasks, snapshot_sha256=sha256_file(args.manifest.with_name("snapshot.json"))
        )
        stage_worker_inputs(args.manifest, args.data_root, frozenset(task["species"] for task in selected))
    summary = run_manifest(
        manifest_path=args.manifest.expanduser().resolve(),
        data_root=args.data_root.expanduser().resolve(),
        config_dir=args.config_dir.expanduser().resolve(),
        workers=args.workers,
        threads=args.threads,
        fastq_threads=args.fastq_threads,
        compression_threads=args.compression_threads,
        validation_slots=args.validation_slots,
        quant_slots=args.quant_slots,
        fasterq_slots=args.fastq_slots,
        max_in_flight=args.max_in_flight,
        task_selection=args.task_selection.expanduser().resolve() if args.task_selection else None,
    )
    print(
        json.dumps(
            {"task_count": summary["task_count"], "counts": summary["counts"]},
            indent=2,
            sort_keys=True,
        )
    )
    return 1 if any(summary["counts"].get(key, 0) for key in ("failed", "unresolved")) else 0
