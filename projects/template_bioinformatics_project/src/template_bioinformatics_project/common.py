"""Shared configuration helpers for the pipeline stages."""

from __future__ import annotations

import sys
from pathlib import Path

import yaml


def load_config(config_path: str | Path) -> dict:
    """Load and return the YAML configuration file (exits if missing)."""
    config_path = Path(config_path)
    if not config_path.is_file():
        print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open() as fh:
        return yaml.safe_load(fh)


def load_optional_config(config_path: str | Path) -> dict:
    """Load a YAML configuration file, returning ``{}`` when absent."""
    config_path = Path(config_path)
    if not config_path.is_file():
        return {}
    with config_path.open() as fh:
        return yaml.safe_load(fh) or {}
