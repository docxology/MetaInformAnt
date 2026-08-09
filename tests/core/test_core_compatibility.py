"""Tests for canonical core helper paths."""

from __future__ import annotations

from metainformant.core.io import paths as canonical_paths
from metainformant.core.utils import config as canonical_config


def test_core_config_exports_are_available_from_canonical_module() -> None:
    """Configuration helpers are available from the supported module."""
    assert callable(canonical_config.load_mapping_from_file)
    assert callable(canonical_config.apply_env_overrides)


def test_core_paths_exports_are_available_from_canonical_module() -> None:
    """Path helpers are available from the supported module."""
    assert callable(canonical_paths.get_project_root)
    assert callable(canonical_paths.expand_and_resolve)
