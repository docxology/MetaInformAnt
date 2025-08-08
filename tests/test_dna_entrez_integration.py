from __future__ import annotations

import os
import pytest

from metainformant.dna import entrez


@pytest.mark.skipif(not os.environ.get("NCBI_EMAIL"), reason="NCBI email not provided")
def test_entrez_fetch_phiX_fasta_roundtrip():
    email = os.environ["NCBI_EMAIL"]
    # PhiX174 is common small reference
    rec = entrez.get_genome_from_ncbi("NC_001422.1", email=email)
    assert rec.id.startswith("NC_001422")
    assert len(str(rec.seq)) > 1000


