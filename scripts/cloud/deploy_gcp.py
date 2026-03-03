#!/usr/bin/env python3
"""GCP Cloud Deployment CLI for the MetaInformAnt pipeline.

Thin orchestrator script for deploying, monitoring, and managing
a high-core GCP VM running the amalgkit RNA-seq pipeline.

Usage:
    python scripts/cloud/deploy_gcp.py deploy   --project MY_PROJECT
    python scripts/cloud/deploy_gcp.py status
    python scripts/cloud/deploy_gcp.py logs      [--lines 100]
    python scripts/cloud/deploy_gcp.py startup-log
    python scripts/cloud/deploy_gcp.py download  [--output output/amalgkit]
    python scripts/cloud/deploy_gcp.py stop
    python scripts/cloud/deploy_gcp.py start
    python scripts/cloud/deploy_gcp.py destroy

    # Dry run (prints gcloud command without executing):
    python scripts/cloud/deploy_gcp.py deploy --project MY_PROJECT --dry-run
"""
import argparse
import sys
from pathlib import Path

# Add src to Python path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.cloud.cloud_config import CloudConfig
from metainformant.cloud.gcp_deployer import GCPDeployer


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Deploy MetaInformAnt pipeline to GCP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── deploy ───────────────────────────────────────────────────────
    deploy = sub.add_parser("deploy", help="Create VM and start pipeline")
    deploy.add_argument("--project", required=True, help="GCP project ID")
    deploy.add_argument("--zone", default="us-central1-a", help="GCP zone")
    deploy.add_argument("--machine-type", default="n2-standard-32",
                        help="VM machine type (default: n2-standard-32)")
    deploy.add_argument("--disk-gb", type=int, default=1000,
                        help="Boot disk size in GB (default: 1000)")
    deploy.add_argument("--spot", action="store_true", default=True,
                        help="Use spot/preemptible pricing (default: true)")
    deploy.add_argument("--no-spot", action="store_false", dest="spot",
                        help="Use on-demand pricing")
    deploy.add_argument("--max-gb", type=float, default=20.0,
                        help="Max sample size in GB (default: 20.0)")
    deploy.add_argument("--workers", type=int, default=28,
                        help="Parallel workers (default: 28)")
    deploy.add_argument("--threads", type=int, default=32,
                        help="Total threads (default: 32)")
    deploy.add_argument("--name", default="metainformant-pipeline",
                        help="VM instance name")
    deploy.add_argument("--gcs-bucket", default="",
                        help="GCS bucket for periodic result sync")
    deploy.add_argument("--dry-run", action="store_true",
                        help="Print gcloud command without executing")

    # ── status ───────────────────────────────────────────────────────
    status = sub.add_parser("status", help="Check pipeline progress")
    _add_common(status)

    # ── logs ─────────────────────────────────────────────────────────
    logs = sub.add_parser("logs", help="Tail pipeline logs")
    logs.add_argument("--lines", type=int, default=50, help="Lines to show")
    _add_common(logs)

    # ── startup-log ──────────────────────────────────────────────────
    slog = sub.add_parser("startup-log", help="Check VM startup script log")
    _add_common(slog)

    # ── download ─────────────────────────────────────────────────────
    dl = sub.add_parser("download", help="Download results locally")
    dl.add_argument("--output", default="output/amalgkit",
                    help="Local output directory")
    _add_common(dl)

    # ── stop / start / destroy ───────────────────────────────────────
    for name, hlp in [("stop", "Stop VM (keep disk)"),
                       ("start", "Start stopped VM"),
                       ("destroy", "Delete VM and disk")]:
        p = sub.add_parser(name, help=hlp)
        _add_common(p)

    return parser


def _add_common(parser: argparse.ArgumentParser) -> None:
    """Add common args for commands that operate on an existing VM."""
    parser.add_argument("--project", default="", help="GCP project ID")
    parser.add_argument("--zone", default="us-central1-a", help="GCP zone")
    parser.add_argument("--name", default="metainformant-pipeline",
                        help="VM instance name")


def make_deployer(args: argparse.Namespace) -> GCPDeployer:
    """Build a deployer from parsed args."""
    if not GCPDeployer.gcloud_installed():
        print("❌ gcloud CLI not found. Install it with:")
        print("   bash scripts/cloud/install_gcloud.sh")
        sys.exit(1)

    cfg = CloudConfig(
        project=args.project,
        zone=args.zone,
        instance_name=args.name,
    )

    # Deploy command has extra fields
    if args.command == "deploy":
        cfg.machine_type = args.machine_type
        cfg.disk_size_gb = args.disk_gb
        cfg.spot = args.spot
        cfg.max_gb = args.max_gb
        cfg.workers = args.workers
        cfg.threads = args.threads
        cfg.gcs_bucket = args.gcs_bucket

    return GCPDeployer(cfg)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    deployer = make_deployer(args)

    if args.command == "deploy":
        if args.dry_run:
            result = deployer.create_vm(dry_run=True)
            print("DRY RUN — would execute:")
            print(f"  {result['command']}")
        else:
            deployer.full_deploy()

    elif args.command == "status":
        vm_status = deployer.get_vm_status()
        print(f"VM Status: {vm_status.get('status', 'UNKNOWN')}")
        if vm_status.get("status") == "RUNNING":
            print()
            print(deployer.get_pipeline_status())

    elif args.command == "logs":
        print(deployer.tail_logs(lines=args.lines))

    elif args.command == "startup-log":
        print(deployer.get_startup_log())

    elif args.command == "download":
        print("📥 Downloading results...")
        if deployer.download_results(local_dir=args.output):
            print(f"   ✓ Results saved to {args.output}/")
        else:
            print("   ❌ Download failed")

    elif args.command == "stop":
        if deployer.stop_vm():
            print("✓ VM stopped (disk preserved)")
        else:
            print("❌ Failed to stop VM")

    elif args.command == "start":
        if deployer.start_vm():
            print("✓ VM started")
        else:
            print("❌ Failed to start VM")

    elif args.command == "destroy":
        confirm = input(f"⚠ Delete VM '{args.name}' and ALL data? [y/N] ")
        if confirm.lower() == "y":
            if deployer.delete_vm():
                print("✓ VM and disks deleted")
            else:
                print("❌ Failed to delete VM")
        else:
            print("Cancelled.")


if __name__ == "__main__":
    main()
