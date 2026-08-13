"""Regression tests for active documentation-link validation."""

from __future__ import annotations

from pathlib import Path

from scripts.rna._verify_utils import check_doc_links


def test_directory_fragment_resolves_through_readme(tmp_path: Path) -> None:
    """A link to a directory fragment must not make the validator read a directory."""
    (tmp_path / "README.md").write_text("# Home\n\n## Stable API\n")
    doc = tmp_path / "guide.md"
    doc.write_text("[API](#stable-api)\n")

    assert check_doc_links(doc, tmp_path) == []


def test_fragment_only_link_resolves_in_current_document(tmp_path: Path) -> None:
    """A fragment-only link targets the document containing the link."""
    doc = tmp_path / "guide.md"
    doc.write_text("## Local API\n\n[API](#local-api)\n")

    assert check_doc_links(doc, tmp_path) == []


def test_broken_relative_link_is_reported(tmp_path: Path) -> None:
    """Missing active documentation targets remain blocking findings."""
    doc = tmp_path / "guide.md"
    doc.write_text("[missing](missing.md)\n")

    issues = check_doc_links(doc, tmp_path)
    assert len(issues) == 1
    assert "missing.md" in issues[0]
