#!/usr/bin/env python3
"""Reclaim raw inputs only for exact current-contract quantifications.

This is intentionally separate from broad cleanup commands. It never uses a
stale progress-table state or a readable table without current provenance as evidence;
the per-sample 0.16.32 provenance sidecar and a non-empty abundance output are
required.  The default is a dry run.  Pass ``--execute`` to remove the listed
FASTQ/SRA inputs and write a JSON audit manifest.
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
from metainformant.rna.engine.provenance import is_current_quantification  # noqa: E402
from metainformant.rna.engine.raw_cleanup import reclaim_sample_raw_inputs  # noqa: E402
from metainformant.rna.engine.species import (  # noqa: E402
    configured_data_root,
    discover_species_config_names,
)
from metainformant.rna.engine.streaming_orchestrator import _species_work_dir  # noqa: E402


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        help="Selected Amalgkit data root (default: AMALGKIT_DATA_ROOT or output/amalgkit)",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=REPO_ROOT / "projects/hymenoptera_amalgkit/config/amalgkit",
    )
    parser.add_argument("--species", action="append", help="Restrict to one or more species identifiers")
    parser.add_argument("--execute", action="store_true", help="Actually remove raw inputs; default is dry-run")
    parser.add_argument("--manifest", type=Path, help="JSON audit path (default: <data-root>/results/raw_reclamation_*.json)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.data_root:
        import os

        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root().expanduser().resolve()
    config_dir = args.config_dir.expanduser().resolve()
    wanted = set(args.species or [])
    configs = discover_species_config_names(config_dir)
    if wanted:
        configs = [name for name in configs if name.removeprefix("amalgkit_").removesuffix(".yaml") in wanted]
    if not configs:
        raise SystemExit(f"no runnable species configurations found under {config_dir}")

    records: list[dict[str, object]] = []
    totals = {"species": 0, "current_quantified": 0, "files": 0, "bytes": 0}
    for config_name in configs:
        species = config_name.removeprefix("amalgkit_").removesuffix(".yaml")
        work_dir = _species_work_dir(species)
        quant_root = work_dir / "quant"
        if not quant_root.is_dir():
            continue
        totals["species"] += 1
        for sample_dir in sorted(quant_root.iterdir()):
            if not sample_dir.is_dir():
                continue
            run_accession = sample_dir.name
            quant_file = find_quantification_file(sample_dir, run_accession)
            if quant_file is None or not quant_file.is_file() or quant_file.stat().st_size <= 100:
                continue
            if not is_current_quantification(sample_dir, run_accession):
                continue
            result = reclaim_sample_raw_inputs(work_dir, run_accession, dry_run=not args.execute)
            record = {"species": species, **result}
            records.append(record)
            totals["current_quantified"] += 1
            totals["files"] += int(result["files_deleted"])
            totals["bytes"] += int(result["bytes_freed"])

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = (args.manifest or data_root / "results" / f"raw_reclamation_{timestamp}.json").expanduser()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "metainformant.rna.raw_reclamation.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "config_dir": str(config_dir),
        "mode": "execute" if args.execute else "dry-run",
        "totals": totals,
        "records": records,
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(manifest_path), **totals, "mode": payload["mode"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
