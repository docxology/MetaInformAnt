from metainformant.dna import restriction


def test_find_restriction_sites_simple() -> None:
    seq = "AAGCTTAAGCTT"  # HindIII site AAGCTT occurs twice
    sites = restriction.find_restriction_sites(seq, {"HindIII": "AAGCTT"})
    assert sites["HindIII"] == [0, 6]


def test_find_restriction_sites_with_ambiguity() -> None:
    seq = "GGCCGAGGCCGA"
    # Example motif with N ambiguity (match any) and R (A/G)
    sites = restriction.find_restriction_sites(seq, {"Fake": "GGNCR"})
    assert isinstance(sites["Fake"], list)
    assert all(isinstance(x, int) for x in sites["Fake"])  # positions list
