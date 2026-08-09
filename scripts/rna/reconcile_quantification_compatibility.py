#!/usr/bin/env python3
"""Classify and safely restore existing quantification outputs.

The default is a dry run. --apply updates only the progress database;
sample outputs, raw inputs, and archive directories are never modified.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from metainformant.rna.core.sample_utils import find_quantification_file  # noqa: E402
from metainformant.rna.engine.progress_db import ProgressDB  # noqa: E402
from metainformant.rna.engine.provenance import (  # noqa: E402
    QUANT_STATUS_CURRENT,
    QUANT_STATUS_VERSION_DRIFT,
    classify_quantification,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--species", action="append", help="Restrict to one or more species")
    parser.add_argument("--apply", action="store_true", help="Apply database classifications; default is dry-run")
    parser.add_argument(
        "--requantification-policy",
        choices=("preserve", "version-drift", "all"),
        default="preserve",
        help="Queue compatible outputs for explicit rebuild instead of preserving them",
    )
    parser.add_argument("--manifest", type=Path, help="JSON audit manifest path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    data_root = args.data_root.expanduser().resolve()
    db = ProgressDB(data_root / "pipeline_progress.db")
    wanted = set(args.species or [])
    records: list[dict[str, object]] = []
    work_dirs: dict[str, Path] = {}

    for candidate in sorted(data_root.glob("*/work")):
        species = candidate.parent.name.removesuffix("_all")
        if wanted and species not in wanted:
            continue
        if (candidate / "quant").is_dir():
            work_dirs.setdefault(species, candidate)

    for species, work_dir in sorted(work_dirs.items()):
        quant_root = work_dir / "quant"
        sample_dirs = sorted(path for path in quant_root.iterdir() if path.is_dir() and not path.name.startswith("."))
        for sample_dir in sample_dirs:
            srr_id = sample_dir.name
            quant_file = find_quantification_file(sample_dir, srr_id)
            if quant_file is None or not quant_file.is_file():
                continue
            audit = classify_quantification(sample_dir, srr_id)
            status = str(audit["status"])
            if status in {QUANT_STATUS_CURRENT, QUANT_STATUS_VERSION_DRIFT}:
                if args.requantification_policy == "all":
                    target_state = "pending"
                elif args.requantification_policy == "version-drift" and status == QUANT_STATUS_VERSION_DRIFT:
                    target_state = "pending"
                else:
                    target_state = "quantified"
            else:
                target_state = "quarantined"
            record = {
                "species": species,
                "srr_id": srr_id,
                "status": status,
                "target_state": target_state,
                "reason": audit["reason"],
                "contract_id": audit.get("contract_id"),
                "observed_amalgkit_version": audit.get("observed_amalgkit_version"),
                "output": str(quant_file),
            }
            records.append(record)
            if args.apply:
                db.init_species(species, [srr_id])
                db.record_quantification_audit(
                    species,
                    srr_id,
                    status=status,
                    reason=str(audit["reason"]),
                    contract_id=audit.get("contract_id"),
                    observed_amalgkit_version=audit.get("observed_amalgkit_version"),
                    observed_release_tag=audit.get("observed_release_tag"),
                    observed_source_revision=audit.get("observed_source_revision"),
                    state=target_state,
                )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest = args.manifest or data_root / "results" / f"quantification_compatibility_{timestamp}.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    summary: dict[str, int] = {}
    for record in records:
        key = f"{record['status']}->{record['target_state']}"
        summary[key] = summary.get(key, 0) + 1
    payload = {
        "schema": "metainformant.rna.quantification_compatibility.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "mode": "apply" if args.apply else "dry-run",
        "requantification_policy": args.requantification_policy,
        "summary": summary,
        "records": records,
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(manifest), "mode": payload["mode"], "summary": summary}, sort_keys=True))
    db.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
