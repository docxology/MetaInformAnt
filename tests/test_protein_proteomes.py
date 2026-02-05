"""Tests for protein proteomes module."""

from __future__ import annotations

from pathlib import Path

from metainformant.protein.sequence.proteomes import read_taxon_ids


def test_read_taxon_ids_basic(tmp_path: Path) -> None:
    """Test basic taxon ID reading."""
    taxon_file = tmp_path / "taxon_ids.txt"
    taxon_file.write_text("9606\n10090\n10116\n")

    result = read_taxon_ids(taxon_file)
    # read_taxon_ids returns List[str]
    expected = ["9606", "10090", "10116"]

    assert result == expected


def test_read_taxon_ids_with_comments(tmp_path: Path) -> None:
    """Test reading taxon IDs with comments and empty lines."""
    taxon_file = tmp_path / "taxon_ids_comments.txt"
    taxon_file.write_text(
        """# Human
9606

# Mouse
10090

# Rat
10116
"""
    )

    result = read_taxon_ids(taxon_file)
    expected = ["9606", "10090", "10116"]

    assert result == expected


def test_read_taxon_ids_empty_file(tmp_path: Path) -> None:
    """Test reading empty taxon ID file."""
    taxon_file = tmp_path / "empty.txt"
    taxon_file.write_text("")

    result = read_taxon_ids(taxon_file)
    assert result == []


def test_read_taxon_ids_invalid_lines(tmp_path: Path) -> None:
    """Test reading file with invalid lines."""
    taxon_file = tmp_path / "invalid.txt"
    taxon_file.write_text(
        """9606
invalid_text
10090
also_invalid
10116
"""
    )

    result = read_taxon_ids(taxon_file)
    expected = ["9606", "10090", "10116"]  # Invalid lines should be skipped

    assert result == expected


def test_read_taxon_ids_whitespace_handling(tmp_path: Path) -> None:
    """Test handling of whitespace in taxon ID file."""
    taxon_file = tmp_path / "whitespace.txt"
    taxon_file.write_text("  9606  \n\t10090\t\n  \n10116\n")

    result = read_taxon_ids(taxon_file)
    expected = ["9606", "10090", "10116"]

    assert result == expected
