#!/usr/bin/env python3
"""
main.py — Top-level CLI entry point for the template bioinformatics project.

Delegates to numbered pipeline scripts in scripts/ based on the --stage flag.
All heavy computation is performed inside metainformant.* library modules.

Usage:
    uv run main.py --stage all --config config/default.yaml
    uv run main.py --stage process
    uv run main.py --stage analyze
    uv run main.py --stage visualize
"""

import sys
import subprocess
import argparse
from pathlib import Path


STAGES = {
    "process": "scripts/01_process_data.py",
    "analyze": "scripts/02_analyze_results.py",
    "visualize": "scripts/03_visualize.py",
}

ORDER = ["process", "analyze", "visualize"]


def run_stage(stage_name: str, script: str, config: str) -> int:
    """Run a single pipeline stage via subprocess and return exit code."""
    print(f"\n{'='*60}")
    print(f"  Stage: {stage_name.upper():20s} → {script}")
    print(f"{'='*60}")
    result = subprocess.run(
        [sys.executable, script, "--config", config],
        cwd=Path(__file__).parent,
    )
    return result.returncode


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Template bioinformatics pipeline runner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Stages:\n"
            "  process    Stage 1 — ingest and process raw data\n"
            "  analyze    Stage 2 — downstream statistical analysis\n"
            "  visualize  Stage 3 — figure generation\n"
            "  all        Run all stages in order\n"
        ),
    )
    parser.add_argument(
        "--stage",
        choices=[*STAGES.keys(), "all"],
        default="all",
        help="Which stage to run (default: all)",
    )
    parser.add_argument(
        "--config",
        default="config/default.yaml",
        help="Path to the configuration YAML (default: config/default.yaml)",
    )
    args = parser.parse_args()

    stages_to_run = ORDER if args.stage == "all" else [args.stage]

    print(f"\n🚀 Starting pipeline  (stage={args.stage}, config={args.config})")

    failed = []
    for stage in stages_to_run:
        script = STAGES[stage]
        rc = run_stage(stage, script, args.config)
        if rc != 0:
            print(f"\n❌  Stage '{stage}' failed with exit code {rc}.", file=sys.stderr)
            failed.append(stage)
            break  # Stop on first failure

    if failed:
        print(f"\n❌  Pipeline aborted at stage: {failed[0]}", file=sys.stderr)
        sys.exit(1)

    print("\n✅  All requested stages completed successfully.")
    print("   📂 Raw data:     data/raw/")
    print("   📂 Processed:    data/processed/")
    print("   📂 Results:      results/")
    print("   📂 Logs:         logs/")


if __name__ == "__main__":
    main()
