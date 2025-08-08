from __future__ import annotations

import pytest

from metainformant.dna import ncbi


def test_ncbi_datasets_optional_dependency_errors():
    # Expect a clear error if ncbi-datasets-pylib is not installed
    with pytest.raises(RuntimeError):
        ncbi.get_accession_by_tax_id("9606")
    with pytest.raises(RuntimeError):
        ncbi.get_metadata_by_single_accession(["GCF_000001405.39"])
    with pytest.raises(RuntimeError):
        ncbi.download_genome_data_package(["GCF_000001405.39"], filename="/tmp/out.zip")


