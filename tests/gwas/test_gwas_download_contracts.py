"""Regression tests for GWAS download completion contracts."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.data import download


def test_reference_genome_download_does_not_return_empty_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Both failed acquisition strategies must remain visible as failure."""

    def fail(*args: object, **kwargs: object) -> Path:
        raise RuntimeError("offline")

    monkeypatch.setattr(download, "_download_from_ncbi_datasets", fail)
    monkeypatch.setattr(download, "_download_from_ftp", fail)
    with pytest.raises(RuntimeError, match="All genome download methods failed"):
        download.download_reference_genome("GCF_000001405.39", tmp_path)
    assert not (tmp_path / "GCF_000001405.39").exists()
