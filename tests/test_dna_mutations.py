from metainformant.dna import mutations


def test_point_mutations_and_hamming() -> None:
    s = "ACGTACGT"
    res = mutations.apply_point_mutations(s, {2: "T", 7: "A"})
    assert res == "ACTTACGA"
    assert mutations.hamming_distance(s, res) == 2


def test_random_mutate_reproducible() -> None:
    s = "AAAAAA"
    out1 = mutations.random_point_mutations(s, num_mutations=2, seed=42)
    out2 = mutations.random_point_mutations(s, num_mutations=2, seed=42)
    assert out1 == out2
