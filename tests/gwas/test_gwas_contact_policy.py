"""Tests for explicit NCBI contact policy in GWAS/SRA helpers."""

from __future__ import annotations

import pytest
import requests

from metainformant.core.ncbi import NCBIContactError
from metainformant.gwas.data import download, sra_download


@pytest.mark.parametrize(
    "operation",
    [
        lambda: download._get_project_runs("PRJNA1"),
        lambda: download.search_sra_for_organism("Apis mellifera"),
        lambda: sra_download._get_experiment_runs("SRX1"),
        lambda: sra_download._get_biosample_runs("SAMN1"),
        lambda: sra_download.prefetch_sra_metadata(["SRR1"]),
    ],
)
def test_gwas_ncbi_helpers_require_contact(monkeypatch: pytest.MonkeyPatch, operation) -> None:
    """GWAS NCBI helpers fail before networking when contact is unspecified."""
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(NCBIContactError):
        operation()


def test_gwas_ncbi_helpers_support_explicit_anonymous_mode(monkeypatch: pytest.MonkeyPatch) -> None:
    """Anonymous operation is available only through an explicit opt-in."""
    monkeypatch.delenv("NCBI_EMAIL", raising=False)

    def timeout(*args: object, **kwargs: object) -> None:
        raise requests.Timeout()

    monkeypatch.setattr(download.requests, "get", timeout)
    with pytest.warns(RuntimeWarning, match="anonymous NCBI"):
        result = download._get_project_runs("PRJNA1", allow_anonymous=True)
    assert result == []
