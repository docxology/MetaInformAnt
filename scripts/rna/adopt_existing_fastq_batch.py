#!/usr/bin/env python3
"""Adopt selected pre-existing FASTQ run directories for every current species.

The mounted data root may contain pre-existing ``<species>/fastq/getfastq`` trees and
the current runner consumes ``<active-work>/getfastq``.  This command moves
only accession directories present in current selected metadata, never
overwrites an existing target, and writes one auditable JSON manifest.

The default is a dry run.  Pass ``--execute`` only after reviewing the counts.
The operation is intended for a single mounted data volume; no raw reads are
deleted, and a rerun is idempotent because moved directories are no longer in
the source and existing targets are skipped.
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
_RUN_COLUMNS = ("run_accession", "run", "run_id", "sra_run")
_NON_SPECIES_MARKERS = ("template", "test", "cross_species")


def _selected_runs(metadata_path: Path) -> set[str]:
    """Read concrete run accessions from a selected metadata TSV."""

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


def _species_work_dir(data_root: Path, species: str) -> Path:
    """Resolve the same active aliases used by the streaming orchestrator."""

    candidates = [species]
    if species == "apis_mellifera":
        candidates.extend(("apis_mellifera_all", "amellifera"))
    elif species == "pogonomyrmex_barbatus":
        candidates.append("pbarbatus")

    for candidate in candidates:
        work_dir = data_root / candidate / "work"
        if work_dir.is_dir() and any((work_dir / marker).exists() for marker in ("metadata", "index", "quant", "merge")):
            return work_dir
    for candidate in candidates:
        if (data_root / candidate).is_dir():
            return data_root / candidate / "work"
    return data_root / species / "work"


def _species_from_config(path: Path) -> str:
    return path.stem.removeprefix("amalgkit_")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--species", action="append", help="Limit to one or more species identifiers")
    parser.add_argument("--execute", action="store_true", help="Move selected directories; default is dry-run")
    parser.add_argument("--manifest", type=Path, help="JSON audit path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    data_root = args.data_root.expanduser().resolve()
    config_dir = args.config_dir.expanduser().resolve()
    if not data_root.is_dir():
        raise SystemExit(f"data root is not a directory: {data_root}")
    config_paths = [
        path
        for path in sorted(config_dir.glob("amalgkit_*.yaml"))
        if not any(marker in _species_from_config(path).lower() for marker in _NON_SPECIES_MARKERS)
    ]
    if args.species:
        requested = set(args.species)
        config_paths = [path for path in config_paths if _species_from_config(path) in requested]
    if not config_paths:
        raise SystemExit(f"no species configurations found under {config_dir}")

    records: list[dict[str, object]] = []
    for config_path in config_paths:
        species = _species_from_config(config_path)
        active_work = _species_work_dir(data_root, species)
        source = data_root / species / "fastq" / "getfastq"
        target = active_work / "getfastq"
        metadata = active_work / "metadata" / "metadata_selected.tsv"
        if not metadata.is_file():
            records.append(
                {
                    "species": species,
                    "action": "skipped_no_selected_metadata",
                    "source_root": str(source),
                    "target_root": str(target),
                    "metadata": str(metadata),
                }
            )
            continue
        if not source.is_dir():
            records.append(
                {
                    "species": species,
                    "action": "skipped_no_preexisting_source",
                    "source_root": str(source),
                    "target_root": str(target),
                    "metadata": str(metadata),
                }
            )
            continue

        selected = _selected_runs(metadata)
        target.mkdir(parents=True, exist_ok=True)
        species_counts: dict[str, int] = {}
        for candidate in sorted(source.iterdir()):
            if not candidate.is_dir() or not _RUN_ACCESSION.fullmatch(candidate.name):
                continue
            record: dict[str, object] = {
                "species": species,
                "run_accession": candidate.name,
                "source": str(candidate),
                "target_root": str(target),
            }
            if candidate.name not in selected:
                action = "skipped_not_in_selected_metadata"
            elif (target / candidate.name).exists():
                record["target"] = str(target / candidate.name)
                action = "skipped_target_exists"
            elif args.execute:
                destination = target / candidate.name
                shutil.move(str(candidate), str(destination))
                record["target"] = str(destination)
                action = "moved"
            else:
                record["target"] = str(target / candidate.name)
                action = "would_move"
            record["action"] = action
            species_counts[action] = species_counts.get(action, 0) + 1
            records.append(record)
        if not species_counts:
            records.append(
                {
                    "species": species,
                    "action": "source_has_no_accession_directories",
                    "source_root": str(source),
                    "target_root": str(target),
                    "metadata": str(metadata),
                }
            )

    counts: dict[str, int] = {}
    for record in records:
        action = str(record["action"])
        counts[action] = counts.get(action, 0) + 1
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = (args.manifest or data_root / "results" / f"raw_adoption_batch_{timestamp}.json").expanduser()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "metainformant.rna.raw_adoption_batch.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "config_dir": str(config_dir),
        "data_root": str(data_root),
        "mode": "execute" if args.execute else "dry-run",
        "species_count": len(config_paths),
        "counts": counts,
        "records": records,
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(manifest_path), "mode": payload["mode"], **counts}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
