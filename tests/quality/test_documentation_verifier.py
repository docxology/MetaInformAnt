"""Tests for documentation code verification helpers."""

from __future__ import annotations

from pathlib import Path

from scripts.verify_documentation_code import (
    CodeValidator,
    DocumentationParser,
    PythonSymbolIndex,
    ReportGenerator,
    Violation,
)


def test_documentation_parser_skips_historical_reports_by_default(tmp_path: Path) -> None:
    """Historical reports are retained in-tree but excluded from current verification."""
    docs_dir = tmp_path / "docs"
    docs_dir.mkdir()
    current = docs_dir / "current.md"
    current.write_text("# Current\n\n```python\nprint('ok')\n```\n")
    validation = docs_dir / "VALIDATION_REPORT.md"
    validation.write_text("# Active validation\n\n```python\nprint('current')\n```\n")
    snapshot = docs_dir / "snapshot.md"
    snapshot.write_text("# Snapshot\n\n> Historical snapshot: retained for provenance.\n")

    parser = DocumentationParser(docs_dir)
    assert set(parser.find_all_docs()) == {current, validation}

    include_parser = DocumentationParser(docs_dir, include_historical=True)
    assert set(include_parser.find_all_docs()) == {current, validation, snapshot}


def test_symbol_index_preserves_package_name_for_subpackage_root() -> None:
    """A subpackage src-dir should still index importable metainformant module names."""
    index = PythonSymbolIndex(Path("src/metainformant/core")).build_index()

    assert index.module_exists("metainformant.core.io")
    assert index.module_exists("metainformant.core.utils.config")


def test_optional_imports_are_skipped_unless_strict() -> None:
    """Known optional dependencies are non-fatal by default and strict on request."""
    index = PythonSymbolIndex(Path("src/metainformant/core")).build_index()
    pyproject = Path("pyproject.toml")

    default_validator = CodeValidator(index, pyproject)
    strict_validator = CodeValidator(index, pyproject, strict_optional_imports=True)

    assert default_validator._module_exists("h5py")
    assert default_validator._module_exists("pysam")
    assert not strict_validator._module_exists("definitely_missing_optional_package")
    assert default_validator._import_exists("pydantic", "BaseModel")


def test_parser_dedents_markdown_blocks_and_skips_multiline_inline(tmp_path: Path) -> None:
    """Indented Markdown examples should parse as blocks, not malformed inline code."""
    docs_dir = tmp_path / "docs"
    docs_dir.mkdir()
    doc = docs_dir / "example.md"
    doc.write_text(
        "\n".join(
            [
                "# Example",
                "",
                "1. Nested block:",
                "",
                "   ```python",
                "   import json",
                "   print(json.dumps({'ok': True}))",
                "   ```",
                "",
                "`not inline",
                "python code`",
            ]
        )
    )

    examples = DocumentationParser(docs_dir).extract_code_examples(doc)

    assert len(examples) == 1
    assert examples[0].language == "python"
    assert examples[0].code.startswith("import json")


def test_parser_relabels_non_script_python_fences(tmp_path: Path) -> None:
    """Console, notebook, Snakemake, and placeholder snippets are not Python scripts."""
    docs_dir = tmp_path / "docs"
    docs_dir.mkdir()
    doc = docs_dir / "examples.md"
    doc.write_text(
        "\n".join(
            [
                "```python",
                ">>> import metainformant",
                "```",
                "```python",
                "!uv pip install metainformant",
                "```",
                "```python",
                "rule run:",
                "    input: 'x'",
                "```",
                "```python",
                "from metainformant.rna import ...",
                "```",
            ]
        )
    )

    languages = [example.language for example in DocumentationParser(docs_dir).extract_code_examples(doc)]

    assert languages == ["pycon", "python-notebook", "snakemake", "python-snippet"]


def test_reexported_package_imports_are_valid() -> None:
    """Public package hubs that re-export symbols should validate as importable."""
    index = PythonSymbolIndex(Path("src/metainformant")).build_index()
    validator = CodeValidator(index, Path("pyproject.toml"))

    assert validator._import_exists("metainformant.rna.engine.workflow", "execute_workflow")
    assert validator._import_exists("metainformant.rna.engine.workflow", "WorkflowExecutionResult")
    assert validator._import_exists("metainformant.rna.amalgkit", "AMALGKIT_INSTALL_SPEC")


def test_report_generator_escapes_markdown_table_cells(tmp_path: Path) -> None:
    """Report rows should remain valid Markdown tables when issue text contains pipes."""
    report_path = tmp_path / "report.md"
    violation = Violation(
        doc_file=str(Path("docs/example.md").resolve()),
        line_number=7,
        code_example="from x import y | z",
        issue_type="ImportError",
        details="Cannot import y | z\nsecond line",
    )

    ReportGenerator([violation], Path.cwd()).generate_report(report_path)

    report = report_path.read_text()
    assert "Cannot import y \\| z second line" in report
    assert "`from x import y \\| z`" in report
