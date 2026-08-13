"""Adversarial tests for external NCBI archive extraction."""

from __future__ import annotations

import stat
import zipfile
from pathlib import Path

import pytest

from metainformant.gwas.data.download import _extract_zip_safely


def test_zip_extraction_confines_members_to_destination(tmp_path: Path) -> None:
    """Traversal members are rejected and cannot create files outside output."""
    archive = tmp_path / "malicious.zip"
    outside = tmp_path / "escaped.txt"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../../escaped.txt", "should not be written")

    with pytest.raises(ValueError, match="Unsafe ZIP member path"):
        _extract_zip_safely(archive, tmp_path / "extracted")

    assert not outside.exists()


def test_zip_extraction_rejects_symlink_members(tmp_path: Path) -> None:
    """A symlink entry cannot redirect a later extracted file."""
    archive = tmp_path / "symlink.zip"
    info = zipfile.ZipInfo("redirect")
    info.external_attr = (stat.S_IFLNK | 0o777) << 16
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(info, "/tmp/should-not-follow")

    with pytest.raises(ValueError, match="Symlink ZIP members"):
        _extract_zip_safely(archive, tmp_path / "extracted")


def test_zip_extraction_writes_valid_members(tmp_path: Path) -> None:
    """Valid nested files remain available after confined extraction."""
    archive = tmp_path / "valid.zip"
    destination = tmp_path / "extracted"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("nested/genome.fna", ">seq\nACGT\n")

    _extract_zip_safely(archive, destination)

    assert (destination / "nested" / "genome.fna").read_text() == ">seq\nACGT\n"
