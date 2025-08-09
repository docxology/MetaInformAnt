from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping
from datetime import datetime
import json
import hashlib

from .amalgkit import AmalgkitParams, build_amalgkit_command, check_cli_available
from . import steps as _steps_mod
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
    # Preflight: ensure amalgkit is available
    ok, help_or_msg = check_cli_available()
    if not ok:
        ensure_directory(config.work_dir)
        manifest_path = config.manifest_path or (config.work_dir / "amalgkit.manifest.jsonl")
        ensure_directory(manifest_path.parent)
        write_jsonl(
            [
                {
                    "step": "preflight",
                    "return_code": 127,
                    "stdout_bytes": 0,
                    "stderr_bytes": len(help_or_msg or ""),
                    "started_utc": datetime.utcnow().isoformat() + "Z",
                    "finished_utc": datetime.utcnow().isoformat() + "Z",
                    "duration_seconds": 0.0,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": {},
                    "command": "amalgkit -h",
                    "note": "amalgkit not available on PATH",
                }
            ],
            manifest_path,
        )
        return [127]
    steps = plan_workflow(config)
    return_codes: list[int] = []
    manifest_records: list[dict[str, Any]] = []
    for subcommand, params in steps:
        runner = _steps_mod.STEP_RUNNERS.get(subcommand)
        if runner is None:  # pragma: no cover - defensive
            raise KeyError(f"No runner registered for step: {subcommand}")
        start_ts = datetime.utcnow()
        result = runner(
            params,
            work_dir=config.work_dir,
            log_dir=(config.log_dir or (config.work_dir / "logs")),
            check=check,
        )
        end_ts = datetime.utcnow()
        duration_s = max(0.0, (end_ts - start_ts).total_seconds())
        return_codes.append(result.returncode)
        manifest_records.append(
            {
                "step": subcommand,
                "return_code": result.returncode,
                "stdout_bytes": len(result.stdout or ""),
                "stderr_bytes": len(result.stderr or ""),
                "started_utc": start_ts.isoformat() + "Z",
                "finished_utc": end_ts.isoformat() + "Z",
                "duration_seconds": duration_s,
                "work_dir": str(config.work_dir),
                "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                "params": dict(params),
                "command": " ".join(build_amalgkit_command(subcommand, params)),
            }
        )
        if check and result.returncode != 0:
            break
    # Write manifest
    manifest_path = config.manifest_path or (config.work_dir / "amalgkit.manifest.jsonl")
    ensure_directory(manifest_path.parent)
    write_jsonl(manifest_records, manifest_path)
    # Also write a compact JSON and Markdown report
    _write_run_reports(config, manifest_records)
    return return_codes


def _write_run_reports(config: AmalgkitWorkflowConfig, records: list[dict[str, Any]]) -> None:
    """Emit JSON and Markdown summaries next to the manifest.

    Files:
      - amalgkit.report.json
      - amalgkit.report.md
    """
    report_json = config.work_dir / "amalgkit.report.json"
    report_md = config.work_dir / "amalgkit.report.md"
    ensure_directory(config.work_dir)

    # JSON: include config summary and records
    summary = {
        "work_dir": str(config.work_dir),
        "log_dir": str(config.log_dir or (config.work_dir / "logs")),
        "threads": config.threads,
        "species_list": list(config.species_list),
        "num_steps": len(records),
        "return_codes": [r.get("return_code", -1) for r in records],
    }
    dump_json({"summary": summary, "records": records}, report_json, indent=2)

    # Markdown: brief human-readable view
    lines: list[str] = []
    lines.append(f"# Amalgkit Run Report\n")
    lines.append(f"Work dir: `{config.work_dir}`  ")
    lines.append(f"Logs: `{config.log_dir or (config.work_dir / 'logs')}`  ")
    lines.append(f"Threads: {config.threads}  ")
    if config.species_list:
        lines.append(f"Species: {', '.join(config.species_list)}  ")
    lines.append("")
    lines.append("| Step | Code | Duration (s) |")
    lines.append("|------|------|--------------|")
    for rec in records:
        lines.append(f"| {rec['step']} | {rec['return_code']} | {rec.get('duration_seconds', 0.0):.2f} |")
    report_md.write_text("\n".join(lines), encoding="utf-8")


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
    # Accept YAML list for species_list
    species_raw = raw.get("species_list", [])
    if isinstance(species_raw, list):
        species_list = [str(x) for x in species_raw]
    else:
        species_list = []
    steps_map = raw.get("steps", {}) or {}
    if not isinstance(steps_map, dict):
        steps_map = {}
    # Keep only known step names if provided
    allowed = {
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
    }
    steps_map = {k: v for k, v in steps_map.items() if str(k) in allowed and isinstance(v, dict)}

    return AmalgkitWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        species_list=species_list,
        per_step={str(k): dict(v) for k, v in steps_map.items()},
    )


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
    "load_workflow_config",
]


