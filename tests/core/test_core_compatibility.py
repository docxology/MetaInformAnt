"""Compatibility import tests for migrated core helper paths."""

from __future__ import annotations

from metainformant.core.config import apply_env_overrides, load_mapping_from_file
from metainformant.core.io import paths as canonical_paths
from metainformant.core.paths import expand_and_resolve, get_project_root
from metainformant.core.utils import config as canonical_config


def test_core_config_compatibility_exports_canonical_functions() -> None:
    """The deprecated core.config path re-exports canonical config helpers."""
    assert load_mapping_from_file is canonical_config.load_mapping_from_file
    assert apply_env_overrides is canonical_config.apply_env_overrides


def test_core_paths_compatibility_exports_canonical_functions() -> None:
    """The deprecated core.paths path re-exports canonical path helpers."""
    assert get_project_root is canonical_paths.get_project_root
    assert expand_and_resolve is canonical_paths.expand_and_resolve
