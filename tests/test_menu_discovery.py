"""Tests for script discovery module."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from metainformant.menu.core.discovery import (
    ScriptInfo,
    categorize_script,
    discover_scripts,
    extract_script_metadata,
    generate_menu_from_scripts,
)


class TestCategorizeScript:
    """Tests for script categorization."""

    def test_categorize_script_rna(self) -> None:
        """Test categorizing RNA script."""
        script_path = Path("scripts/rna/run_workflow.py")
        category = categorize_script(script_path)
        assert category == "rna"

    def test_categorize_script_gwas(self) -> None:
        """Test categorizing GWAS script."""
        script_path = Path("scripts/gwas/run_analysis.py")
        category = categorize_script(script_path)
        assert category == "gwas"

    def test_categorize_script_core(self) -> None:
        """Test categorizing core script."""
        script_path = Path("scripts/core/run_demo.py")
        category = categorize_script(script_path)
        assert category == "core"

    def test_categorize_script_archive(self) -> None:
        """Test categorizing archive script (should be 'other')."""
        script_path = Path("scripts/archive/old_script.py")
        category = categorize_script(script_path)
        assert category == "other"

    def test_categorize_script_other(self) -> None:
        """Test categorizing script without category."""
        script_path = Path("script.py")
        category = categorize_script(script_path)
        assert category == "other"


class TestExtractScriptMetadata:
    """Tests for script metadata extraction."""

    def test_extract_python_metadata(self, tmp_path: Path) -> None:
        """Test extracting metadata from Python script."""
        script_content = '''#!/usr/bin/env python3
"""Test script for RNA workflow.

This script runs RNA analysis workflows.
"""
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--threads", type=int, default=4, help="Thread count")
'''
        script_path = tmp_path / "test_script.py"
        script_path.write_text(script_content)

        metadata = extract_script_metadata(script_path)
        assert metadata.name == "test_script"
        assert metadata.script_type == "python"
        assert "RNA workflow" in metadata.description or "RNA analysis" in metadata.description

    def test_extract_bash_metadata(self, tmp_path: Path) -> None:
        """Test extracting metadata from bash script."""
        script_content = '''#!/bin/bash
# Test script for workflow execution
# This script runs a complete workflow

set -euo pipefail
'''
        script_path = tmp_path / "test_script.sh"
        script_path.write_text(script_content)

        metadata = extract_script_metadata(script_path)
        assert metadata.name == "test_script"
        assert metadata.script_type == "bash"
        assert len(metadata.description) > 0

    def test_extract_metadata_no_docstring(self, tmp_path: Path) -> None:
        """Test extracting metadata from script without docstring."""
        script_content = "print('hello')"
        script_path = tmp_path / "simple.py"
        script_path.write_text(script_content)

        metadata = extract_script_metadata(script_path)
        assert metadata.name == "simple"
        assert metadata.script_type == "python"
        # Should have default description
        assert len(metadata.description) > 0


class TestDiscoverScripts:
    """Tests for script discovery."""

    def test_discover_scripts_empty(self, tmp_path: Path) -> None:
        """Test discovering scripts in empty directory."""
        scripts = discover_scripts(tmp_path)
        assert isinstance(scripts, dict)
        assert len(scripts) == 0

    def test_discover_scripts_python(self, tmp_path: Path) -> None:
        """Test discovering Python scripts."""
        scripts_dir = tmp_path / "scripts" / "rna"
        scripts_dir.mkdir(parents=True)

        script_path = scripts_dir / "run_workflow.py"
        script_path.write_text('"""RNA workflow script."""\nprint("test")')

        scripts = discover_scripts(tmp_path)
        assert "rna" in scripts
        assert len(scripts["rna"]) == 1
        assert scripts["rna"][0].name == "run_workflow"

    def test_discover_scripts_bash(self, tmp_path: Path) -> None:
        """Test discovering bash scripts."""
        scripts_dir = tmp_path / "scripts" / "core"
        scripts_dir.mkdir(parents=True)

        script_path = scripts_dir / "setup.sh"
        script_path.write_text("#!/bin/bash\n# Setup script\necho 'test'")

        scripts = discover_scripts(tmp_path)
        assert "core" in scripts
        assert len(scripts["core"]) == 1
        assert scripts["core"][0].name == "setup"

    def test_discover_scripts_skip_test(self, tmp_path: Path) -> None:
        """Test that test scripts are skipped."""
        scripts_dir = tmp_path / "scripts" / "rna"
        scripts_dir.mkdir(parents=True)

        # Regular script
        script_path = scripts_dir / "run_workflow.py"
        script_path.write_text('"""RNA workflow."""')

        # Test script (should be skipped)
        test_path = scripts_dir / "test_workflow.py"
        test_path.write_text('"""Test script."""')

        scripts = discover_scripts(tmp_path)
        assert "rna" in scripts
        # Should only find non-test script
        script_names = [s.name for s in scripts["rna"]]
        assert "run_workflow" in script_names
        assert "test_workflow" not in script_names

    def test_discover_scripts_skip_cache(self, tmp_path: Path) -> None:
        """Test that __pycache__ directories are skipped."""
        scripts_dir = tmp_path / "scripts" / "rna"
        scripts_dir.mkdir(parents=True)

        cache_dir = scripts_dir / "__pycache__"
        cache_dir.mkdir()

        cache_file = cache_dir / "script.pyc"
        cache_file.write_bytes(b"test")

        scripts = discover_scripts(tmp_path)
        # Should not find anything in __pycache__
        if "rna" in scripts:
            assert len(scripts["rna"]) == 0


class TestGenerateMenuFromScripts:
    """Tests for menu generation from scripts."""

    def test_generate_menu_empty(self) -> None:
        """Test generating menu from empty script list."""
        menus = generate_menu_from_scripts({})
        assert "root" in menus
        assert len(menus["root"].items) == 0

    def test_generate_menu_single_category(self) -> None:
        """Test generating menu with single category."""
        scripts = {
            "rna": [
                ScriptInfo(
                    path=Path("scripts/rna/run_workflow.py"),
                    name="run_workflow",
                    description="Run RNA workflow",
                    category="rna",
                    script_type="python",
                )
            ]
        }
        menus = generate_menu_from_scripts(scripts)
        assert "root" in menus
        assert len(menus["root"].items) == 1
        assert "menu_rna" in menus
        assert len(menus["menu_rna"].items) == 1

    def test_generate_menu_multiple_categories(self) -> None:
        """Test generating menu with multiple categories."""
        scripts = {
            "rna": [
                ScriptInfo(
                    path=Path("scripts/rna/run_workflow.py"),
                    name="run_workflow",
                    description="Run RNA workflow",
                    category="rna",
                    script_type="python",
                )
            ],
            "gwas": [
                ScriptInfo(
                    path=Path("scripts/gwas/run_analysis.py"),
                    name="run_analysis",
                    description="Run GWAS analysis",
                    category="gwas",
                    script_type="python",
                )
            ],
        }
        menus = generate_menu_from_scripts(scripts)
        assert "root" in menus
        assert len(menus["root"].items) == 2
        assert "menu_rna" in menus
        assert "menu_gwas" in menus

    def test_generate_menu_script_items(self) -> None:
        """Test that script items are created correctly."""
        scripts = {
            "rna": [
                ScriptInfo(
                    path=Path("scripts/rna/run_workflow.py"),
                    name="run_workflow",
                    description="Run RNA workflow",
                    category="rna",
                    script_type="python",
                )
            ]
        }
        menus = generate_menu_from_scripts(scripts)
        menu_rna = menus["menu_rna"]
        assert len(menu_rna.items) == 1
        item = menu_rna.items[0]
        assert item.action.startswith("script:")
        assert "run_workflow" in item.label




