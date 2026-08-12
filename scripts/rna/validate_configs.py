#!/usr/bin/env python3
"""Validate species YAMLs and selection rules against Amalgkit 0.16.59."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

from metainformant.rna.amalgkit import REQUIRED_AMALGKIT_VERSION, validate_amalgkit_version
from metainformant.rna.amalgkit.commands import CURRENT_AMALGKIT_STEPS, SUPPORTED_CLI_OPTIONS

EXCLUSIONS = {"amalgkit_template.yaml", "amalgkit_cross_species.yaml", "amalgkit_test.yaml"}
REQUIRED_TOP_LEVEL_KEYS = {"work_dir", "species_list", "genome", "steps"}
RETIRED_DIRECT_STEP_KEYS = {
    "config_dir",
    "keep_fastq",
    "mark_missing_rank",
    "max_sample",
    "pfd",
    "pfd_print",
    "priority",
}
ORCHESTRATION_ONLY_STEP_KEYS = {
    "backend",
    "batch_size",
    "bootstrap",
    "fragment_length",
    "fragment_sd",
    "jobs",
    "log_dir",
    "num_download_workers",
    "num_retries",
    "output_format",
    "progress_style",
    "progress_update_interval",
    "retry_delay",
    "show_progress",
    "validate_md5",
    "work_dir",
}


def load_yaml(path: Path) -> dict[str, Any]:
    """Load one YAML mapping."""

    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("top-level YAML value must be a mapping")
    return value


def _resolve_rule_path(config_path: Path, raw_path: str) -> Path:
    """Resolve a rule path relative to the repository root when possible."""

    candidate = Path(raw_path).expanduser()
    if candidate.is_absolute():
        return candidate
    repository_root = config_path.parents[2]
    return repository_root / candidate


def _validate_select_rules(config_path: Path, select: dict[str, Any]) -> list[str]:
    """Validate the explicit current Amalgkit selection-rule file."""

    errors: list[str] = []
    raw_path = select.get("select_rules_tsv")
    if not isinstance(raw_path, str) or not raw_path.strip():
        return ["select.select_rules_tsv is required by Amalgkit 0.16.59"]
    rule_path = _resolve_rule_path(config_path, raw_path)
    if not rule_path.is_file():
        return [f"select rule file does not exist: {raw_path}"]
    try:
        from amalgkit.select import read_select_config

        parsed = read_select_config(str(rule_path))
        required = {"min_nspots", "max_sample", "mark_missing_rank", "mark_redundant_biosamples"}
        missing = sorted(required - set(parsed["parameters"]))
        if missing:
            errors.append(f"select rule file is missing required parameters: {', '.join(missing)}")
    except Exception as exc:  # pragma: no cover - exact parser owns detailed diagnostics
        errors.append(f"invalid Amalgkit select rule file {raw_path}: {exc}")
    return errors


def validate_config(config_path: Path) -> list[str]:
    """Return violations for one current species configuration."""

    try:
        config = load_yaml(config_path)
    except Exception as exc:
        return [f"failed to load YAML: {exc}"]

    errors: list[str] = []
    errors.extend(f"missing required top-level key: {key}" for key in sorted(REQUIRED_TOP_LEVEL_KEYS - set(config)))

    genome = config.get("genome")
    if not isinstance(genome, dict):
        errors.append("genome section must be a mapping")
    elif not genome.get("accession") and not genome.get("files"):
        errors.append("genome section must define accession or files")

    steps = config.get("steps")
    if not isinstance(steps, dict):
        return errors + ["steps section must be a mapping"]

    for step, params in steps.items():
        if step not in CURRENT_AMALGKIT_STEPS:
            errors.append(f"unknown Amalgkit command in steps: {step}")
            continue
        if not isinstance(params, dict):
            errors.append(f"steps.{step} must be a mapping")
            continue
        if step == "select":
            errors.extend(f"{config_path.name}: {error}" for error in _validate_select_rules(config_path, params))
        supported = SUPPORTED_CLI_OPTIONS.get(step, frozenset())
        for key in params:
            normalized = str(key).replace("-", "_")
            if normalized in RETIRED_DIRECT_STEP_KEYS:
                errors.append(
                    f"steps.{step}.{key} is not a direct 0.16.59 option; encode selection policy in select_rules.tsv"
                )
            elif normalized not in supported and normalized not in ORCHESTRATION_ONLY_STEP_KEYS:
                errors.append(f"steps.{step}.{key} is not accepted by the 0.16.59 command registry")

    getfastq = steps.get("getfastq")
    if isinstance(getfastq, dict):
        source_flags = ("ncbi", "aws", "gcp", "ena", "ddbj")
        if source_flags and all(
            str(getfastq.get(flag, "yes")).lower() in {"no", "false", "0"} for flag in source_flags
        ):
            errors.append("all current getfastq download sources are disabled")

    return errors


def main() -> int:
    """Validate the full 27-species configuration inventory."""

    root_dir = Path(__file__).resolve().parents[2]
    default_config_dir = root_dir / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
    if not default_config_dir.is_dir():
        default_config_dir = root_dir / "config" / "amalgkit"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=default_config_dir,
        help=f"Directory containing current Amalgkit YAML files (default: {default_config_dir})",
    )
    args = parser.parse_args()
    config_dir = args.config_dir.expanduser().resolve()
    config_files = sorted(path for path in config_dir.glob("amalgkit_*.yaml") if path.name not in EXCLUSIONS)

    version_ok, version_message = validate_amalgkit_version(REQUIRED_AMALGKIT_VERSION, exact=True)
    if not version_ok:
        print(f"Amalgkit version: FAIL ({version_message})")
        return 1
    print(f"Amalgkit version: OK ({version_message})")
    print(f"Validating {len(config_files)} species configuration files...")

    failed = False
    for config_path in config_files:
        errors = validate_config(config_path)
        if errors:
            failed = True
            print(f"{config_path.name}: FAIL")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"{config_path.name}: OK")

    if failed:
        return 1
    if len(config_files) != 27:
        print(f"Expected 27 species configs, found {len(config_files)}")
        return 1
    print("All species configurations and current selection rules are valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
