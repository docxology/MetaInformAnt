from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping
from datetime import datetime
import json
import hashlib

from .amalgkit import AmalgkitParams
from .steps import STEP_RUNNERS
from ..core.config import load_mapping_from_file, apply_env_overrides
from ..core.io import ensure_directory, write_jsonl, dump_json


@dataclass(slots=True)
class AmalgkitWorkflowConfig:
    work_dir: Path
    threads: int = 4
    species_list: list[str] = field(default_factory=list)
    # Additional common parameters can be added as needed
    log_dir: Path | None = None
    manifest_path: Path | None = None
    per_step: dict[str, AmalgkitParams] = field(default_factory=dict)

    def to_common_params(self) -> AmalgkitParams:
        params: dict[str, Any] = {"threads": self.threads}
        if self.species_list:
            params["species-list"] = list(self.species_list)
        return params


def plan_workflow(config: AmalgkitWorkflowConfig) -> list[tuple[str, AmalgkitParams]]:
    """Return an ordered list of (subcommand, params) representing a full run.

    This does not execute anything; it allows dry inspection and TDD.
    """
    common = config.to_common_params()

    def merge_params(extra: Mapping[str, Any] | None = None) -> AmalgkitParams:
        if not extra:
            return dict(common)
        merged = dict(common)
        merged.update(extra)
        return merged

    steps: list[tuple[str, AmalgkitParams]] = []
    ordered = [
        "metadata",
        "integrate",
        "config",
        "select",
        "getfastq",
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]
    for name in ordered:
        extra = config.per_step.get(name, {})
        steps.append((name, merge_params(extra)))
    return steps


def plan_workflow_with_params(
    config: AmalgkitWorkflowConfig,
    step_params: dict[str, AmalgkitParams],
) -> list[tuple[str, AmalgkitParams]]:
    """Plan workflow while merging per-step params on top of common ones."""
    common = config.to_common_params()

    ordered = [
        "metadata",
        "integrate",
        "config",
        "select",
        "getfastq",
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]
    steps: list[tuple[str, AmalgkitParams]] = []
    for name in ordered:
        params = dict(common)
        params.update(step_params.get(name, {}))
        steps.append((name, params))
    return steps


def execute_workflow(config: AmalgkitWorkflowConfig, *, check: bool = False) -> list[int]:
    """Execute the full amalgkit workflow in order.

    Returns a list of return codes per step in order.
    """
    steps = plan_workflow(config)
    return_codes: list[int] = []
    manifest_records: list[dict[str, Any]] = []
    for subcommand, params in steps:
        runner = STEP_RUNNERS.get(subcommand)
        if runner is None:  # pragma: no cover - defensive
            raise KeyError(f"No runner registered for step: {subcommand}")
        result = runner(
            params,
            work_dir=config.work_dir,
            log_dir=(config.log_dir or (config.work_dir / "logs")),
            check=check,
        )
        return_codes.append(result.returncode)
        manifest_records.append(
            {
                "step": subcommand,
                "return_code": result.returncode,
                "stdout_bytes": len(result.stdout or ""),
                "stderr_bytes": len(result.stderr or ""),
                "timestamp_utc": datetime.utcnow().isoformat() + "Z",
            }
        )
        if check and result.returncode != 0:
            break
    # Write manifest
    manifest_path = config.manifest_path or (config.work_dir / "amalgkit.manifest.jsonl")
    ensure_directory(manifest_path.parent)
    write_jsonl(manifest_records, manifest_path)
    return return_codes


def load_workflow_config(config_file: str | Path) -> AmalgkitWorkflowConfig:
    """Load `AmalgkitWorkflowConfig` from a config file with env overrides.

    Expected top-level keys:
      - work_dir (str)
      - log_dir (str, optional)
      - threads (int)
      - species_list (list[str])
      - steps (mapping of step name -> params mapping)
    """
    raw = load_mapping_from_file(config_file)
    raw = apply_env_overrides(raw, prefix="AK")

    work_dir = Path(raw.get("work_dir", "output/amalgkit/work")).expanduser().resolve()
    log_dir_val = raw.get("log_dir")
    log_dir = Path(log_dir_val).expanduser().resolve() if isinstance(log_dir_val, str) else None
    threads = int(raw.get("threads", 4))
    species_list = list(raw.get("species_list", []))
    steps_map = raw.get("steps", {}) or {}
    if not isinstance(steps_map, dict):
        steps_map = {}

    return AmalgkitWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        species_list=species_list,
        per_step={str(k): dict(v) for k, v in steps_map.items() if isinstance(v, dict)},
    )


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
    "load_workflow_config",
]


