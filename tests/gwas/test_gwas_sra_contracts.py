"""Regression tests for truthful SRA acquisition outcomes."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.data import sra_download


def test_sra_run_does_not_return_empty_directory_when_tools_are_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Unavailable conversion tools must remain an explicit failure."""
    monkeypatch.setattr(sra_download, "check_sra_tools_available", lambda: False)

    with pytest.raises(RuntimeError, match="SRA Toolkit is unavailable"):
        sra_download.download_sra_run("SRR123", tmp_path)
    assert not (tmp_path / "SRR123").exists()


def test_sra_run_rejects_path_traversal_and_invalid_threads(tmp_path: Path) -> None:
    """Run accessions are identifiers, never filesystem path fragments."""

    with pytest.raises(ValueError, match="Invalid SRA run accession"):
        sra_download.download_sra_run("../../outside", tmp_path)
    with pytest.raises(ValueError, match="threads must be at least 1"):
        sra_download.download_sra_run("SRR123", tmp_path, threads=0)
