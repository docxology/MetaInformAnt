"""Regression tests for genome-download success contracts."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna.external import genomes


def test_best_effort_download_does_not_promote_empty_directory(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A failed download must be observable instead of looking complete."""

    def fail(*args: object, **kwargs: object) -> Path:
        raise RuntimeError("offline")

    monkeypatch.setattr(genomes, "download_genome_package", fail)
    with pytest.raises(RuntimeError, match="All genome download methods failed"):
        genomes.download_genome_package_best_effort("GCF_000001405.39", tmp_path)
    assert not (tmp_path / "GCF_000001405.39").exists()


def test_best_effort_download_uses_ftp_only_when_primary_fails(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """The optional FTP fallback remains available after a primary failure."""

    def fail(*args: object, **kwargs: object) -> Path:
        raise RuntimeError("primary offline")

    expected = tmp_path / "GCF_000001405.39"
    monkeypatch.setattr(genomes, "download_genome_package", fail)
    monkeypatch.setattr(genomes, "_download_from_ftp", lambda url, output: expected)
    assert (
        genomes.download_genome_package_best_effort("GCF_000001405.39", tmp_path, ftp_url="ftp://example") == expected
    )
