"""Record and inspect sample exclusions in the RNA-seq progress database.

Exclusions are the durable record backing campaign-level cohort decisions:

- ``permanent_drop`` — the accession can never produce a current
  quantification (for example 16S amplicon submissions in an mRNA cohort).
  The streaming orchestrator skips these during task building.
- ``re_download`` — the local transfer is unusable but the accession itself
  is valid; it is marked for a fresh ENA transfer and stays eligible.

Usage::

    python -m metainformant.rna.engine.exclusions record \\
        --species apis_mellifera --tsv exclusions.tsv --reason-code permanent_drop
    python -m metainformant.rna.engine.exclusions list --species apis_mellifera
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from metainformant.core.utils.logging import get_logger
from metainformant.rna.engine.progress_db import EXCLUSION_REASON_CODES, ProgressDB

logger = get_logger(__name__)


def _resolve_db_path(data_root: Path | None, db_path: Path | None) -> Path:
    """Resolve the progress database path from explicit args or the data root."""
    if db_path is not None:
        return db_path
    if data_root is not None:
        return Path(data_root) / "pipeline_progress.db"
    from metainformant.rna.engine.species import configured_data_root

    return configured_data_root() / "pipeline_progress.db"


def _load_tsv_entries(tsv_path: Path) -> list[dict[str, str]]:
    """Load exclusion entries from a TSV with an ``srr_id`` column."""
    entries: list[dict[str, str]] = []
    with open(tsv_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            srr_id = (row.get("srr_id") or "").strip()
            if not srr_id:
                continue
            entries.append({"srr_id": srr_id})
    if not entries:
        raise ValueError(f"No srr_id column values found in {tsv_path}")
    return entries


def _parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--data-root",
        type=Path,
        help="External Amalgkit data root (default: configured data root)",
    )
    common.add_argument(
        "--db",
        type=Path,
        help="Explicit progress database path (overrides --data-root)",
    )

    parser = argparse.ArgumentParser(
        prog="metainformant.rna.engine.exclusions",
        description="Record and inspect sample exclusions in the progress database",
        parents=[common],
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    record = subparsers.add_parser(
        "record",
        parents=[common],
        help="Record exclusions from a TSV or explicit SRR list",
    )
    record.add_argument("--species", required=True)
    record.add_argument("--tsv", type=Path, help="TSV file with an srr_id column")
    record.add_argument("--srr", action="append", default=[], help="Accession; repeatable")
    record.add_argument(
        "--reason-code",
        choices=sorted(EXCLUSION_REASON_CODES),
        default="permanent_drop",
    )
    record.add_argument("--reason-detail", default=None)
    record.add_argument("--recorded-by", default=None)

    listing = subparsers.add_parser("list", parents=[common], help="List recorded exclusions")
    listing.add_argument("--species", default=None)

    remove = subparsers.add_parser("remove", parents=[common], help="Delete one recorded exclusion")
    remove.add_argument("--species", required=True)
    remove.add_argument("--srr", required=True)

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for exclusion recording and inspection."""
    args = _parser().parse_args(argv)
    db = ProgressDB(_resolve_db_path(args.data_root, args.db))
    try:
        if args.command == "record":
            if bool(args.tsv) == bool(args.srr):
                logger.error("Provide exactly one of --tsv or --srr")
                return 2
            entries = _load_tsv_entries(args.tsv) if args.tsv else [{"srr_id": srr} for srr in args.srr]
            for entry in entries:
                entry["reason_code"] = args.reason_code
                entry["reason_detail"] = args.reason_detail
                entry["recorded_by"] = args.recorded_by
            recorded = db.record_exclusions(args.species, entries)
            print(f"recorded {recorded} exclusions for {args.species} ({args.reason_code})")
        elif args.command == "list":
            rows = db.get_exclusions(args.species)
            for row in rows:
                print(
                    "\t".join(
                        [
                            row["species"],
                            row["srr_id"],
                            row["reason_code"],
                            row["reason_detail"] or "",
                            row["recorded_at"],
                        ]
                    )
                )
            print(f"total: {len(rows)}", file=sys.stderr)
        elif args.command == "remove":
            if db.remove_exclusion(args.species, args.srr):
                print(f"removed exclusion {args.species}/{args.srr}")
            else:
                print(f"no exclusion recorded for {args.species}/{args.srr}")
                return 1
    finally:
        db.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
