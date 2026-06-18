"""Regression tests for eQTL script entrypoints."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _import_script(script_path: Path, module_name: str) -> object:
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_eqtl_scripts_import_with_fresh_output_tree(tmp_path: Path, monkeypatch) -> None:
    """Script logging setup should create directories before opening log files."""
    repo_root = Path(__file__).resolve().parents[2]
    monkeypatch.chdir(tmp_path)

    _import_script(repo_root / "scripts/eqtl/run_eqtl_demo.py", "run_eqtl_demo_import_test")
    _import_script(repo_root / "scripts/eqtl/run_eqtl_real.py", "run_eqtl_real_import_test")

    assert (tmp_path / "output/eqtl/amellifera/logs").is_dir()
    assert (tmp_path / "output/eqtl/amellifera/results").is_dir()
    assert (tmp_path / "output/eqtl/amellifera/plots").is_dir()
