"""Generic Amalgkit acquisition planning, execution and scenario estimation CLI."""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from metainformant.rna.engine.acquisition_allocation import allocate_tasks, write_allocation
from metainformant.rna.engine.acquisition_estimates import ThroughputEvidence, LaneCosts, estimate_lane, combine_estimates
from metainformant.rna.engine.acquisition_snapshot import create_campaign_manifest


def main(argv: list[str] | None = None) -> int:
    """Delegate worker/controller actions; plans and estimates never launch compute."""
    args_list = list(sys.argv[1:] if argv is None else argv)
    if args_list and args_list[0] in {"local", "worker"}:
        from metainformant.rna.engine.acquisition_worker_cli import main as worker_main
        return worker_main(args_list[1:])
    if args_list and args_list[0] == "aws":
        from metainformant.rna.engine.aws_completion import main as aws_main
        if not any(x == "--config-dir" or x.startswith("--config-dir=") for x in args_list) and not {"--help", "-h"}.intersection(args_list):
            raise ValueError("generic AWS acquisition requires an explicit --config-dir")
        if not any(x == "--priority-species" or x.startswith("--priority-species=") for x in args_list):
            args_list.extend(("--priority-species", ""))
        return aws_main(args_list[1:])
    parser = argparse.ArgumentParser(description=__doc__, epilog="Execution: local/worker --help or aws --help. Configure every worker's stage limits explicitly.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    freeze = subparsers.add_parser("freeze", help="Freeze selected metadata plus ENA-discovered runs for any configured species set")
    freeze.add_argument("--data-root", type=Path, required=True)
    freeze.add_argument("--config-dir", type=Path, required=True)
    freeze.add_argument("--output-dir", type=Path, required=True)
    quote = subparsers.add_parser("quote-aws", help="Read current regional compute and gp3 prices without launching workers")
    quote.add_argument("--instance-type", required=True)
    quote.add_argument("--disk-gib", type=int, required=True)
    quote.add_argument("--region", default="us-east-2")
    quote.add_argument("--profile")
    quote.add_argument("--hourly-floor", type=float, default=0.55)
    plan = subparsers.add_parser("plan", help="Freeze a disjoint local/AWS allocation without launching workers")
    source = plan.add_mutually_exclusive_group(required=True)
    source.add_argument("--campaign-root", type=Path)
    source.add_argument("--manifest", type=Path)
    plan.add_argument("--backend", choices=("local", "aws", "hybrid"), required=True)
    plan.add_argument("--local-fraction", type=float, default=0.5)
    plan.add_argument("--cloud-status", type=Path, help="An existing campaign_status JSON snapshot; its observation time is retained")
    plan.add_argument("--output-dir", type=Path, required=True)
    estimate = subparsers.add_parser("estimate", help="Model elapsed time and gross cost at explicit capacity/rate assumptions")
    estimate.add_argument("--allocation", type=Path, required=True)
    estimate.add_argument("--rates", type=Path, required=True, help="JSON lane evidence and cost assumptions; see generic acquisition guide")
    estimate.add_argument("--local-units", type=int, default=4)
    estimate.add_argument("--aws-units", type=int, default=6)
    estimate.add_argument("--spent-usd", type=float, default=0)
    estimate.add_argument("--reserved-usd", type=float, default=0)
    estimate.add_argument("--ceiling-usd", type=float)
    estimate.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(args_list)
    if args.command == "quote-aws":
        from metainformant.rna.engine.acquisition_pricing import quote_aws_worker
        print(json.dumps(asdict(quote_aws_worker(region=args.region, instance_type=args.instance_type, disk_gib=args.disk_gib,
                                                profile=args.profile, hourly_floor=args.hourly_floor)), indent=2))
        return 0
    if args.command == "freeze":
        from metainformant.rna.engine.completion_inventory import freeze_inventory
        frozen = freeze_inventory(args.data_root, args.config_dir, args.output_dir, expected_species_count=None)
        print(json.dumps({"species": frozen["species_count"], "eligible": frozen["task_count"], "inventory": str(args.output_dir / "inventory.json")}, indent=2))
        return 0
    if args.command == "plan":
        manifest = create_campaign_manifest(args.campaign_root) if args.campaign_root else args.manifest
        completed, reserved = frozenset(), frozenset()
        source_observation = None
        if args.cloud_status:
            status = json.loads(args.cloud_status.read_text())
            cloud = status.get("cloud_observation", status.get("cloud"))
            if not isinstance(cloud, dict):
                raise ValueError("cloud-status JSON lacks cloud observations")
            completed, reserved = frozenset(cloud["locked"]), frozenset(cloud["assigned"])
            source_observation = cloud["observed_at"]
        allocation = allocate_tasks(manifest, backend=args.backend, completed=completed, reserved=reserved, local_fraction=args.local_fraction)
        target = write_allocation(allocation, args.output_dir)
        print(json.dumps({"allocation": str(target), "eligible": allocation.eligible,
                          "completed": len(allocation.completed_task_ids), "reserved": len(allocation.reserved_task_ids),
                          "local_pending": len(allocation.local_task_ids), "aws_pending": len(allocation.aws_task_ids),
                          "source_observed_at": source_observation}, indent=2))
        return 0
    allocation = json.loads(args.allocation.read_text())
    if allocation.get("schema") != "metainformant.rna.acquisition_allocation.v1":
        raise ValueError("unsupported acquisition allocation for estimate")
    assumptions = json.loads(args.rates.read_text())
    lanes = {}
    for lane, units in (("local", args.local_units), ("aws", args.aws_units)):
        samples = len(allocation[f"{lane}_task_ids"])
        if not samples:
            continue
        evidence = ThroughputEvidence(**assumptions[lane]["evidence"])
        costs = LaneCosts(**assumptions[lane]["costs"])
        lanes[lane] = estimate_lane(samples, units, evidence, costs)
    if not lanes:
        raise ValueError("allocation has no pending lane work; reserved tasks are not a completion estimate")
    total = combine_estimates(tuple(lanes.values()), spent_usd=args.spent_usd,
                              reserved_usd=args.reserved_usd, ceiling_usd=args.ceiling_usd)
    result = {"schema": "metainformant.rna.acquisition_estimate.v1", "lanes": {k: asdict(v) for k, v in lanes.items()},
              "total": asdict(total), "rates_source": str(args.rates),
              "reserved_work_scope": "Active work is outside pending-lane estimates; do not double count its reservation.",
              "completion_scope": "Pending lane work only; final completion also depends on reserved tasks."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0
