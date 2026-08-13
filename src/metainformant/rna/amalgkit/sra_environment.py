"""Campaign-local subprocess environments for SRA and FASTQ tooling."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from metainformant.core.io.sra_environment import (
    _build_sra_environment,
    _ensure_ncbi_settings,
    _resolve_data_root,
)


def resolve_data_root(data_root: str | Path | None = None) -> Path:
    """Resolve the campaign data root without consulting user-global state."""

    return _resolve_data_root(data_root)


def ensure_ncbi_settings(data_root: str | Path | None = None) -> Path:
    """Create the deterministic campaign-local NCBI settings file."""

    return _ensure_ncbi_settings(data_root)


def build_sra_environment(
    data_root: str | Path | None = None,
    *,
    base_environment: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Return an inherited environment scoped to one campaign data root."""

    return _build_sra_environment(data_root, base_environment=base_environment)
