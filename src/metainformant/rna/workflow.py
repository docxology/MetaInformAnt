from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping
from datetime import datetime
import json
import hashlib

from .amalgkit import (
    AmalgkitParams,
    run_amalgkit,
)


@dataclass(slots=True)
class AmalgkitWorkflowConfig:
    work_dir: Path
    threads: int = 4
    species_list: list[str] = field(default_factory=list)
    # Additional common parameters can be added as needed
    log_dir: Path | None = None
    manifest_path: Path | None = None

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

    steps: list[tuple[str, AmalgkitParams]] = [
        ("metadata", merge_params({})),
        ("integrate", merge_params({})),
        ("config", merge_params({})),
        ("select", merge_params({})),
        ("getfastq", merge_params({})),
        ("quant", merge_params({})),
        ("merge", merge_params({})),
        ("cstmm", merge_params({})),
        ("curate", merge_params({})),
        ("csca", merge_params({})),
        ("sanity", merge_params({})),
    ]
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
        result = run_amalgkit(
            subcommand,
            params,
            work_dir=config.work_dir,
            log_dir=(config.log_dir or (config.work_dir / "logs")),
            step_name=subcommand,
            check=check,
            capture_output=True,
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
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as f:
        for rec in manifest_records:
            f.write(json.dumps(rec) + "\n")
    return return_codes


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
]


