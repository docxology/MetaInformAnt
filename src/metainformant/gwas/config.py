"""GWAS workflow configuration management."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..core.config import apply_env_overrides, load_mapping_from_file


@dataclass
class GWASWorkflowConfig:
    """Configuration for GWAS workflow execution."""

    work_dir: Path
    threads: int = 8
    log_dir: Path | None = None
    genome: dict[str, Any] | None = None
    variants: dict[str, Any] = field(default_factory=dict)
    qc: dict[str, Any] = field(default_factory=dict)
    samples: dict[str, Any] = field(default_factory=dict)
    structure: dict[str, Any] = field(default_factory=dict)
    association: dict[str, Any] = field(default_factory=dict)
    correction: dict[str, Any] = field(default_factory=dict)
    output: dict[str, Any] = field(default_factory=dict)


def load_gwas_config(config_file: str | Path) -> GWASWorkflowConfig:
    """Load GWASWorkflowConfig from a config file with env overrides.

    Expected top-level keys:
      - work_dir (str)
      - log_dir (str, optional)
      - threads (int, default: 8)
      - genome (mapping, optional)
      - variants (mapping, optional)
      - qc (mapping, optional)
      - samples (mapping, optional)
      - structure (mapping, optional)
      - association (mapping, optional)
      - correction (mapping, optional)
      - output (mapping, optional)
    """
    raw = load_mapping_from_file(config_file)
    raw = apply_env_overrides(raw, prefix="GWAS")

    work_dir = Path(raw.get("work_dir", "output/gwas/work")).expanduser().resolve()
    log_dir_val = raw.get("log_dir")
    log_dir = Path(log_dir_val).expanduser().resolve() if isinstance(log_dir_val, str) else None
    threads = int(raw.get("threads", 8))

    genome_cfg = raw.get("genome") if isinstance(raw.get("genome"), dict) else None
    variants_cfg = raw.get("variants") if isinstance(raw.get("variants"), dict) else {}
    qc_cfg = raw.get("qc") if isinstance(raw.get("qc"), dict) else {}
    samples_cfg = raw.get("samples") if isinstance(raw.get("samples"), dict) else {}
    structure_cfg = raw.get("structure") if isinstance(raw.get("structure"), dict) else {}
    association_cfg = raw.get("association") if isinstance(raw.get("association"), dict) else {}
    correction_cfg = raw.get("correction") if isinstance(raw.get("correction"), dict) else {}
    output_cfg = raw.get("output") if isinstance(raw.get("output"), dict) else {}

    return GWASWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        genome=genome_cfg,
        variants=variants_cfg,
        qc=qc_cfg,
        samples=samples_cfg,
        structure=structure_cfg,
        association=association_cfg,
        correction=correction_cfg,
        output=output_cfg,
    )

