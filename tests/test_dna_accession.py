from metainformant.dna.genomes import is_valid_assembly_accession


def test_is_valid_assembly_accession_valid_cases():
    assert is_valid_assembly_accession("GCF_000001405.39")
    assert is_valid_assembly_accession("GCA_000001405")


def test_is_valid_assembly_accession_invalid_cases():
    for invalid in [
        "",
        "GCF_",
        "GCF_123",
        "GCF_000001405.XY",
        "GCX_000001405.1",
        "GCA_00000140",
    ]:
        assert not is_valid_assembly_accession(invalid)
