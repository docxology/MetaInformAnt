"""Static package-boundary checks."""

from __future__ import annotations

from pathlib import Path

from scripts.quality.check_module_boundaries import collect_violations


def test_package_boundaries_allow_only_core_or_adapter_imports() -> None:
    """Domain modules should not import each other except through adapter surfaces."""
    repo_root = Path(__file__).resolve().parents[2]

    violations = collect_violations(repo_root)

    assert violations == []
