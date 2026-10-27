"""Regression tests: ``download_variant_data`` does not pretend to download genotype VCFs.

Defect: for ``data_type`` ``"vcf"`` (the default) and ``"both"`` the function wrote an empty ``vcf/`` directory
holding a README ("This would typically involve downloading from dbSNP ... For now, create a placeholder
structure"), listed that directory among the *downloaded files* of the study metadata and returned the study
directory as if data had been fetched. The GWAS Catalog REST API publishes association records, not study-level
genotype files, so there is nothing to download: the function now raises ``NotImplementedError`` for those types.
An unrecognised ``data_type`` used to be a silent no-op that still wrote the metadata file; it is a ``ValueError``.

Everything here is decided before the function touches the disk or the network, so the tests run offline. The
module also points the HTTP proxy variables at a closed local port while it runs: a regression that reached for
the service would then fail on a refused local connection instead of contacting any external host.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from metainformant.gwas.data.download import download_variant_data

ACCESSION = "GCST000001"
PROXY_VARIABLES = ("HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "http_proxy", "https_proxy", "all_proxy")


@pytest.fixture(autouse=True)
def no_external_network() -> Iterator[None]:
    """Route any HTTP request of the code under test to a closed local port (nothing leaves the machine)."""
    saved = {name: os.environ.get(name) for name in (*PROXY_VARIABLES, "NO_PROXY", "no_proxy")}
    for name in PROXY_VARIABLES:
        os.environ[name] = "http://127.0.0.1:9"
    os.environ["NO_PROXY"] = os.environ["no_proxy"] = ""
    try:
        yield
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


@pytest.mark.parametrize("data_type", ["vcf", "both"])
def test_genotype_vcf_data_is_not_implemented(tmp_path: Path, data_type: str) -> None:
    out = tmp_path / "variants"

    with pytest.raises(NotImplementedError, match="not implemented"):
        download_variant_data(ACCESSION, out, data_type=data_type)

    # Nothing was written: no study directory, no placeholder vcf/ directory, no README.
    assert not out.exists()


def test_the_default_data_type_is_the_unimplemented_vcf_type(tmp_path: Path) -> None:
    out = tmp_path / "variants"

    with pytest.raises(NotImplementedError):
        download_variant_data(ACCESSION, out)

    assert not out.exists()


def test_the_error_says_what_can_be_done_instead(tmp_path: Path) -> None:
    with pytest.raises(NotImplementedError) as excinfo:
        download_variant_data(ACCESSION, tmp_path / "variants", data_type="vcf")

    message = str(excinfo.value)
    assert "summary_stats" in message
    assert "download_variant_database" in message
    assert "genotype" in message


@pytest.mark.parametrize("data_type", ["bam", "", "VCF", "summary", None])
def test_an_unknown_data_type_is_a_value_error_not_a_silent_no_op(tmp_path: Path, data_type: Any) -> None:
    out = tmp_path / "variants"

    with pytest.raises(ValueError, match="Unknown data_type"):
        download_variant_data(ACCESSION, out, data_type=data_type)

    assert not out.exists()


@pytest.mark.parametrize("data_type", ["vcf", "both", "summary_stats", "bam"])
def test_an_invalid_accession_is_still_rejected_first(tmp_path: Path, data_type: str) -> None:
    out = tmp_path / "variants"

    with pytest.raises(ValueError, match="Invalid GWAS Catalog study accession"):
        download_variant_data("INVALID_STUDY", out, data_type=data_type)

    assert not out.exists()


def test_summary_stats_is_the_supported_type_and_reaches_the_service_layer(tmp_path: Path) -> None:
    """With the supported type the validation passes and the first thing attempted is the HTTP request.

    The proxy variables point at a closed local port, so the request fails locally with a requests error
    (not a NotImplementedError or ValueError): proof that the supported path is untouched by the new checks,
    obtained without contacting the service.
    """
    import requests

    with pytest.raises(requests.RequestException):
        download_variant_data(ACCESSION, tmp_path / "variants", data_type="summary_stats")
