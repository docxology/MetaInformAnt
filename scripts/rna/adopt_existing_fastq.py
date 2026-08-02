#!/usr/bin/env python3
"""Adopt already-downloaded raw runs into an active data-root workspace.

The operation is non-destructive with respect to sample content: it moves
whole accession directories within the selected external volume and never
overwrites an existing target.  The default is a dry run; ``--execute`` is
required for the move and a JSON manifest records every decision.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

_RUN_ACCESSION = re.compile(r"(?:SRR|ERR|DRR)\d+$")
_RUN_COLUMNS = ("run_accession", "run", "run_id", "sra_run", "run_id")


def _selected_runs(metadata_path: Path) -> set[str]:
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        column = next((name for name in _RUN_COLUMNS if name in fieldnames), None)
        if column is None:
            raise ValueError(f"metadata has no run-accession column: {metadata_path}")
        return {
            str(row[column]).strip()
            for row in reader
            if row.get(column) and _RUN_ACCESSION.fullmatch(str(row[column]).strip())
        }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True, help="Existing raw-run directory")
    parser.add_argument("--target-root", type=Path, required=True, help="Active raw-run directory")
    parser.add_argument("--metadata", type=Path, required=True, help="Current selected metadata TSV")
    parser.add_argument("--data-root", type=Path, required=True, help="Data root used for the audit manifest")
    parser.add_argument("--execute", action="store_true", help="Move directories; default is dry-run")
    parser.add_argument("--manifest", type=Path, help="JSON audit path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    source = args.source_root.expanduser().resolve()
    target = args.target_root.expanduser().resolve()
    metadata = args.metadata.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    if not source.is_dir():
        raise SystemExit(f"source root is not a directory: {source}")
    selected = _selected_runs(metadata)
    target.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    for candidate in sorted(source.iterdir()):
        if not candidate.is_dir() or not _RUN_ACCESSION.fullmatch(candidate.name):
            continue
        record: dict[str, object] = {"run_accession": candidate.name, "source": str(candidate)}
        if candidate.name not in selected:
            record["action"] = "skipped_not_in_selected_metadata"
        elif (target / candidate.name).exists():
            record["action"] = "skipped_target_exists"
            record["target"] = str(target / candidate.name)
        elif args.execute:
            destination = target / candidate.name
            shutil.move(str(candidate), str(destination))
            record["action"] = "moved"
            record["target"] = str(destination)
        else:
            record["action"] = "would_move"
            record["target"] = str(target / candidate.name)
        records.append(record)

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = (args.manifest or data_root / "results" / f"raw_adoption_{timestamp}.json").expanduser()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for record in records:
        action = str(record["action"])
        counts[action] = counts.get(action, 0) + 1
    payload = {
        "schema": "metainformant.rna.raw_adoption.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_root": str(source),
        "target_root": str(target),
        "metadata": str(metadata),
        "data_root": str(data_root),
        "mode": "execute" if args.execute else "dry-run",
        "selected_run_count": len(selected),
        "counts": counts,
        "records": records,
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(manifest_path), "mode": payload["mode"], **counts}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
