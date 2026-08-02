#!/usr/bin/env python3
"""Validate the current Hymenoptera Amalgkit configuration and plan contract.

The validator is intentionally read-only. It checks every runnable species
configuration, confirms the fixed nine-stage plan, and verifies that the
current producer and status scripts are present. It never contacts ENA,
starts Amalgkit, or writes to a data root.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
if not DEFAULT_CONFIG_DIR.is_dir():
    DEFAULT_CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"

EXPECTED_STEPS = (
    "metadata",
    "select",
    "getfastq",
    "integrate",
    "quant",
    "merge",
    "wsfilter",
    "finalize",
    "sanity",
)
EXCLUSIONS = {"amalgkit_template.yaml", "amalgkit_test.yaml", "amalgkit_cross_species.yaml"}
CURRENT_SCRIPTS = (
    REPO_ROOT / "scripts" / "rna" / "run_all_species.py",
    REPO_ROOT / "scripts" / "rna" / "process_species.py",
    REPO_ROOT / "scripts" / "rna" / "check_pipeline_status.py",
)


def discover_species_configs(config_dir: Path) -> list[Path]:
    """Return sorted runnable species configurations from ``config_dir``."""

    return sorted(
        path
        for path in config_dir.glob("amalgkit_*.yaml")
        if path.name not in EXCLUSIONS and path.is_file()
    )


def load_config(path: Path) -> dict[str, Any]:
    """Load one YAML mapping without mutating it."""

    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("top-level YAML value must be a mapping")
    return value


def validate_config(path: Path) -> list[str]:
    """Return structural violations for one runnable species config."""

    try:
        config = load_config(path)
    except Exception as exc:
        return [f"cannot read YAML: {exc}"]

    errors: list[str] = []
    for key in ("work_dir", "threads", "species_list", "genome", "steps"):
        if key not in config:
            errors.append(f"missing top-level key: {key}")
    if not isinstance(config.get("threads"), int) or config.get("threads", 0) < 1:
        errors.append("threads must be a positive integer")
    if not isinstance(config.get("species_list"), list) or not config.get("species_list"):
        errors.append("species_list must be a non-empty list")
    if not isinstance(config.get("genome"), dict):
        errors.append("genome must be a mapping")
    if not isinstance(config.get("steps"), dict):
        errors.append("steps must be a mapping")
    return errors


def planned_steps(path: Path) -> tuple[str, ...]:
    """Resolve the engine plan for one config without executing any step."""

    sys.path.insert(0, str(REPO_ROOT / "src"))
    from metainformant.rna.engine.workflow import load_workflow_config, plan_workflow

    config = load_workflow_config(path)
    return tuple(name for name, _ in plan_workflow(config))


def parse_args() -> argparse.Namespace:
    """Parse read-only validator options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=DEFAULT_CONFIG_DIR,
        help=f"Current species-config directory (default: {DEFAULT_CONFIG_DIR})",
    )
    return parser.parse_args()


def main() -> int:
    """Run configuration, plan, and current-entrypoint validation."""

    args = parse_args()
    config_dir = args.config_dir.expanduser().resolve()
    configs = discover_species_configs(config_dir)
    failures: list[str] = []

    print(f"Configuration directory: {config_dir}")
    print(f"Runnable species configurations: {len(configs)}")
    if not configs:
        failures.append("no runnable species configurations found")

    for path in configs:
        errors = validate_config(path)
        if errors:
            failures.extend(f"{path.name}: {error}" for error in errors)
            print(f"FAIL {path.name}: {'; '.join(errors)}")
            continue
        try:
            actual_steps = planned_steps(path)
        except Exception as exc:
            failures.append(f"{path.name}: plan resolution failed: {exc}")
            print(f"FAIL {path.name}: plan resolution failed: {exc}")
            continue
        if actual_steps != EXPECTED_STEPS:
            failures.append(f"{path.name}: expected {EXPECTED_STEPS}, got {actual_steps}")
            print(f"FAIL {path.name}: plan {actual_steps}")
        else:
            print(f"PASS {path.name}: {' -> '.join(actual_steps)}")

    for script in CURRENT_SCRIPTS:
        if not script.is_file():
            failures.append(f"missing current script: {script}")

    print()
    if failures:
        print(f"VALIDATION FAILED ({len(failures)} issue(s))")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("VALIDATION PASSED: current configs, nine-stage plans, and entrypoints are ready")
    print()
    print("Read-only inspection:")
    print(
        "  uv run python scripts/rna/run_all_species.py "
        f"--config-dir {config_dir} --data-root \"$AMALGKIT_DATA_ROOT\" --dry-run"
    )
    print("Execution:")
    print("  bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
