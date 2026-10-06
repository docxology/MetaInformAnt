"""CLI for live, read-only cloud/local Hymenoptera status snapshots."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

from metainformant.rna.engine.campaign_status import load_inventory, markdown_tables, reconcile
from metainformant.rna.engine.campaign_status_io import collect_cloud, file_coverage, read_database


def main() -> int:
    """Write timestamped Markdown, JSON and per-sample transfer/status TSVs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--local-root", type=Path, required=True)
    parser.add_argument("--diagnostic-root", type=Path)
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--profile", default="dev-agent")
    parser.add_argument("--region", default="us-east-2")
    parser.add_argument("--output-dir", type=Path, default=Path("output/hymenoptera_status"))
    parser.add_argument("--no-worker-probe", action="store_true", help="Assigned worker stages become unknown; never assumed pending")
    parser.add_argument("--stages-only", action="store_true", help="Fast DB/worker status only; omit file coverage and transfer claims")
    args = parser.parse_args()
    started = datetime.now(UTC).isoformat()
    inventory_bytes = (args.campaign_root / "inventory.json").read_bytes()
    inventory = load_inventory(inventory_bytes)
    inventory.task_ids()
    cloud = collect_cloud(args.campaign_root, args.bucket, args.profile, args.region, probe_workers=not args.no_worker_probe)
    local = read_database(args.local_root / "pipeline_progress.db")
    present, partial = frozenset(), frozenset()
    if not args.stages_only:
        present, partial = file_coverage(inventory, args.local_root)
    diagnostic = frozenset()
    if args.diagnostic_root and not args.stages_only:
        diagnostic, _ = file_coverage(inventory, args.diagnostic_root)
    report = reconcile(inventory, cloud.locked, cloud.assigned, cloud.worker, local, present, partial, diagnostic)
    finished = datetime.now(UTC).isoformat()
    directory = args.output_dir / datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    payload = {
        "schema": "metainformant.rna.campaign_status.v1",
        "started_at": started, "finished_at": finished,
        "inventory_sha256": hashlib.sha256(inventory_bytes).hexdigest(),
        "local_root": str(args.local_root.resolve()),
        "diagnostic_root": str(args.diagnostic_root.resolve()) if args.diagnostic_root else None,
        "cloud_observation": asdict(cloud), "report": asdict(report),
        "file_coverage": "NOT_OBSERVED" if args.stages_only else "NONEMPTY_FILE_PRESENCE_ONLY",
    }
    if args.stages_only:
        payload["report"] = {
            "rows": [{"species": row.species, "eligible": row.eligible, "cloud": row.cloud, "local": row.local} for row in report.rows],
            "totals": {"eligible": report.totals.eligible, "cloud": report.totals.cloud, "local": report.totals.local},
            "local_outside_inventory": report.local_outside_inventory,
        }
    # Sets are protocol representations only; deterministic JSON lists preserve observations.
    (directory / "status.json").write_text(json.dumps(payload, indent=2, default=lambda value: sorted(value)) + "\n")
    notes = (
        f"# Hymenoptera cloud/local status\n\nObservation window: {started} to {finished}.\n\n"
        "Cloud stages are mutually exclusive. `locked` means a reference-bound S3 receipt exists; "
        "receipt/blob content is not revalidated by this status command. Live worker DB observations "
        "cover assigned tasks only; `unassigned` means eligible, unlocked and not on a running worker. "
        "Failures are current worker observations, not permanent exclusions. Unknown telemetry remains unknown.\n\n"
        "Local stage counts are recorded SQLite states and may be stale; they do not prove active local processes. "
        "File coverage is separate: all four expected quant/provenance files must exist and be nonempty. "
        "`present_locked` is local file presence intersected with S3 receipts, not checksum-verified transfer. "
        "`transfer_gap` is S3-locked samples lacking complete files at the canonical local root. "
        "Diagnostic copies are an overlapping subset in a separate root. Stage columns sum to Total; "
        "the four coverage columns also sum to Total; transfer_gap and diagnostic_present are overlapping marginals.\n\n"
        f"Running cloud instances: {len(cloud.instances)}. Local DB rows outside inventory: {report.local_outside_inventory}.\n\n"
    )
    if cloud.diagnostics:
        notes += "Telemetry diagnostics:\n\n" + "\n".join(f"- {d}" for d in cloud.diagnostics) + "\n\n"
    if args.stages_only:
        notes += "File coverage and transfer candidates were NOT OBSERVED (--stages-only). No file/transfer counts are emitted.\n\n"
    markdown = notes + markdown_tables(report, include_coverage=not args.stages_only)
    (directory / "status.md").write_text(markdown)
    with (directory / "samples.tsv").open("w", newline="") as handle:
        fields = ("task_id", "species", "cloud", "local") if args.stages_only else tuple(asdict(report.samples[0]))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(asdict(sample) for sample in report.samples)
    if not args.stages_only:
        with (directory / "transfer_candidates.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("task_id", "canonical_local_coverage", "diagnostic_present"))
            writer.writerows((s.task_id, s.coverage, s.diagnostic_present) for s in report.samples if s.transfer_gap)
    print(markdown)
    print(f"Saved snapshot: {directory}")
    return 0
