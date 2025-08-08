from __future__ import annotations

import pytest

from metainformant.dna import ncbi
import importlib.util
import os


try:
    installed = importlib.util.find_spec("ncbi.datasets") is not None
except ModuleNotFoundError:
    installed = False


@pytest.mark.skipif(installed, reason="ncbi-datasets installed; expecting runtime error only when missing")
def test_ncbi_datasets_optional_dependency_errors():
    # Expect a clear error if ncbi-datasets-pylib is not installed
    with pytest.raises(RuntimeError):
        ncbi.get_accession_by_tax_id("9606")
    with pytest.raises(RuntimeError):
        ncbi.get_metadata_by_single_accession(["GCF_000001405.39"])
    with pytest.raises(RuntimeError):
        ncbi.download_genome_data_package(["GCF_000001405.39"], filename="/tmp/out.zip")


