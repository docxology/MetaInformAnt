from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping
from datetime import datetime
import json
import hashlib
import os

from .amalgkit import AmalgkitParams, build_amalgkit_command, check_cli_available, ensure_cli_available
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
    auto_install_amalgkit: bool = False
    genome: dict[str, Any] | None = None  # e.g., { accession: GCF_*, index_out: path, dest_dir: path, ftp_url: str, include: [...] }

    def to_common_params(self) -> AmalgkitParams:
        params: dict[str, Any] = {"threads": self.threads}
        if self.species_list:
            params["species-list"] = list(self.species_list)
        return params


def _build_default_search_string(species_list: list[str]) -> str | None:
    if not species_list:
        return None
    if len(species_list) == 1:
        sp = species_list[0].replace("_", " ")
        return f'"{sp}"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
    # Multiple species: join with OR
    parts = [f'"{sp.replace("_", " ")}"[Organism]' for sp in species_list]
    return f'({" OR ".join(parts)}) AND RNA-Seq[Strategy] AND Illumina[Platform]'


def _apply_step_defaults(config: AmalgkitWorkflowConfig) -> None:
    """Fill in sensible defaults for amalgkit steps when not provided in config."""
    ps = config.per_step
    # metadata defaults
    md = dict(ps.get("metadata", {}))
    md.setdefault("out_dir", str(config.work_dir))
    search = _build_default_search_string(config.species_list)
    if search and not md.get("search_string"):
        md["search_string"] = search
    email = os.environ.get("NCBI_EMAIL")
    if email and not md.get("entrez_email"):
        md["entrez_email"] = email
    ps["metadata"] = md

    # directories
    defaults = {
        "getfastq": {"out_dir": str(config.work_dir / "fastq")},
        "quant": {"out_dir": str(config.work_dir / "quant")},
        "merge": {"out": str(config.work_dir.parent / "merged" / "merged_abundance.tsv")},
        "cstmm": {"out_dir": str(config.work_dir / "cstmm")},
        "curate": {"out_dir": str(config.work_dir / "curate")},
        "csca": {"out_dir": str(config.work_dir / "csca")},
    }
    for step, d in defaults.items():
        cur = dict(ps.get(step, {}))
        for k, v in d.items():
            cur.setdefault(k, v)
        ps[step] = cur


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


def _write_manifest_records(path: Path, records: list[dict[str, Any]]) -> None:
    ensure_directory(path.parent)
    write_jsonl(records, path)


def _sanitize_params_for_subcommand(subcommand: str, params: Mapping[str, Any]) -> dict[str, Any]:
    """Drop unsupported flags for specific subcommands.

    Currently, `amalgkit metadata` accepts only: out_dir, redo, search_string, entrez_email.
    """
    if subcommand == "metadata":
        allowed = {"out_dir", "redo", "search_string", "entrez_email"}
        return {k: v for k, v in params.items() if k in allowed}
    return dict(params)


def execute_workflow(config: AmalgkitWorkflowConfig, *, check: bool = False) -> list[int]:
    """Execute the full amalgkit workflow in order.

    Returns a list of return codes per step in order.
    """
    ensure_directory(config.work_dir)
    manifest_path = config.manifest_path or (config.work_dir / "amalgkit.manifest.jsonl")

    manifest_records: list[dict[str, Any]] = []
    return_codes: list[int] = []

    # Ensure defaults are applied for required amalgkit params
    _apply_step_defaults(config)

    # 1) Optional genome preparation FIRST (so it's logged even if amalgkit missing)
    if config.genome:
        from ..dna.ncbi import download_genome_package_best_effort

        acc = str(config.genome.get("accession", ""))
        include = config.genome.get("include") or ["gff3", "rna", "cds", "protein", "genome", "seq-report"]
        ftp_url = config.genome.get("ftp_url")
        default_dest = config.work_dir.parent / "genome"
        dest_dir = Path(config.genome.get("dest_dir", default_dest)).expanduser().resolve()
        start_ts = datetime.utcnow()
        dl_rec = download_genome_package_best_effort(acc, dest_dir, include=include, ftp_url=ftp_url)
        end_ts = datetime.utcnow()
        genome_rec = {
            "step": "genome-prepare",
            "return_code": int(dl_rec.get("return_code", 0)),
            "stdout_bytes": len(dl_rec.get("stdout", "")),
            "stderr_bytes": len(dl_rec.get("stderr", "")),
            "started_utc": start_ts.isoformat() + "Z",
            "finished_utc": end_ts.isoformat() + "Z",
            "duration_seconds": max(0.0, (end_ts - start_ts).total_seconds()),
            "work_dir": str(config.work_dir),
            "log_dir": str(config.log_dir or (config.work_dir / "logs")),
            "params": {"accession": acc, "include": include, "dest_dir": str(dest_dir), "ftp_url": ftp_url or ""},
            "command": dl_rec.get("command", dl_rec.get("url", dl_rec.get("method", ""))),
            "extracted_dir": dl_rec.get("extracted_dir", ""),
            "zip_path": dl_rec.get("zip_path", ""),
            "method": dl_rec.get("method", ""),
        }
        manifest_records.append(genome_rec)
        return_codes.append(genome_rec["return_code"])  # record genome step code
        # Inject quant param if index/genome-dir not present
        quant_params = dict(config.per_step.get("quant", {}))
        if "index" not in quant_params and "genome-dir" not in quant_params and "genome_dir" not in quant_params:
            genome_dir = dl_rec.get("extracted_dir") or str(dest_dir)
            quant_params["genome_dir"] = genome_dir
            config.per_step["quant"] = quant_params

    # 2) Ensure amalgkit is available; optionally auto-install
    install_record: dict[str, Any] | None = None
    ok, help_or_msg = check_cli_available()
    if not ok and config.auto_install_amalgkit:
        ok, help_or_msg, install_record = ensure_cli_available(auto_install=True)
        if install_record is not None:
            manifest_records.append(
                {
                    "step": "amalgkit-install",
                    "return_code": install_record.get("return_code", -1),
                    "stdout_bytes": len(install_record.get("stdout", "")),
                    "stderr_bytes": len(install_record.get("stderr", "")),
                    "started_utc": datetime.utcnow().isoformat() + "Z",
                    "finished_utc": datetime.utcnow().isoformat() + "Z",
                    "duration_seconds": 0.0,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": {},
                    "command": install_record.get("command", ""),
                    "note": "attempted auto-install of amalgkit",
                }
            )
            return_codes.append(install_record.get("return_code", -1))

    if not ok:
        manifest_records.append(
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
        )
        return_codes.append(127)
        _write_manifest_records(manifest_path, manifest_records)
        return return_codes

    # 3) Run amalgkit steps
    steps = plan_workflow(config)
    for subcommand, params in steps:
        runner = _steps_mod.STEP_RUNNERS.get(subcommand)
        if runner is None:  # pragma: no cover - defensive
            raise KeyError(f"No runner registered for step: {subcommand}")
        filtered = _sanitize_params_for_subcommand(subcommand, params)
        start_ts = datetime.utcnow()
        result = runner(
            filtered,
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
                "params": dict(filtered),
                "command": " ".join(build_amalgkit_command(subcommand, filtered)),
            }
        )
        if check and result.returncode != 0:
            break

    _write_manifest_records(manifest_path, manifest_records)
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
      - auto_install_amalgkit (bool, optional)
      - genome (mapping, optional)
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

    auto_install_amalgkit = bool(raw.get("auto_install_amalgkit", False))
    genome_cfg = raw.get("genome") if isinstance(raw.get("genome"), dict) else None

    return AmalgkitWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        species_list=species_list,
        per_step={str(k): dict(v) for k, v in steps_map.items()},
        auto_install_amalgkit=auto_install_amalgkit,
        genome=genome_cfg,
    )


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
    "load_workflow_config",
]


