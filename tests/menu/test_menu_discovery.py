"""Tests for script discovery and current RNA menu entries."""

from __future__ import annotations

from pathlib import Path

from metainformant.menu.core.discovery import (
    ScriptInfo,
    categorize_script,
    discover_scripts,
    extract_script_metadata,
    generate_menu_from_scripts,
)


def test_categorize_script_rna() -> None:
    """Current RNA producer scripts are categorized under RNA."""

    assert categorize_script(Path("scripts/rna/run_all_species.py")) == "rna"


def test_categorize_script_other_domains() -> None:
    """Directory-based categorization remains domain independent."""

    assert categorize_script(Path("scripts/gwas/run_analysis.py")) == "gwas"
    assert categorize_script(Path("scripts/core/run_demo.py")) == "core"
    assert categorize_script(Path("script.py")) == "other"


def test_extract_python_metadata(tmp_path: Path) -> None:
    """Python script metadata is extracted from a real file."""

    script = tmp_path / "producer.py"
    script.write_text(
        '"""Current RNA producer.\n\nRuns the configured cohort.\n"""\n'
        "import argparse\n"
        "parser = argparse.ArgumentParser()\n",
        encoding="utf-8",
    )
    metadata = extract_script_metadata(script)
    assert metadata.name == "producer"
    assert metadata.script_type == "python"
    assert "Current RNA producer" in metadata.description


def test_discover_scripts_finds_current_python_scripts(tmp_path: Path) -> None:
    """Discovery returns current scripts and ignores test files."""

    scripts_dir = tmp_path / "scripts" / "rna"
    scripts_dir.mkdir(parents=True)
    (scripts_dir / "run_all_species.py").write_text('"""Producer."""\n', encoding="utf-8")
    (scripts_dir / "test_producer.py").write_text('"""Test."""\n', encoding="utf-8")

    scripts = discover_scripts(tmp_path)
    assert [item.name for item in scripts["rna"]] == ["run_all_species"]


def test_discover_scripts_skips_cache(tmp_path: Path) -> None:
    """Discovery ignores Python cache directories."""

    cache_dir = tmp_path / "scripts" / "rna" / "__pycache__"
    cache_dir.mkdir(parents=True)
    (cache_dir / "producer.pyc").write_bytes(b"not source")
    assert discover_scripts(tmp_path) == {}


def test_generate_menu_from_current_scripts() -> None:
    """Menu generation preserves the current producer path."""

    scripts = {
        "rna": [
            ScriptInfo(
                path=Path("scripts/rna/run_all_species.py"),
                name="run_all_species",
                description="Inspect and run the current cohort",
                category="rna",
                script_type="python",
            )
        ]
    }
    menus = generate_menu_from_scripts(scripts)
    item = menus["menu_rna"].items[0]
    assert item.action == "script:scripts/rna/run_all_species.py"
