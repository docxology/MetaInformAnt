from metainformant.dna import composition


def test_gc_skew_and_cumulative() -> None:
    seq = "GGGCCCATAA"
    skew = composition.gc_skew(seq)
    assert abs(skew) < 1e-9  # equal G and C
    cums = composition.cumulative_gc_skew(seq)
    assert len(cums) == len(seq)


def test_melting_temp_basic() -> None:
    seq_short = "ATGC"
    # Wallace rule 2*(A+T) + 4*(G+C) = 2*2 + 4*2 = 12
    assert abs(composition.melting_temperature(seq_short) - 12.0) < 1e-9
    # Longer sequences can use a basic estimate too
    seq_long = "ATGC" * 10
    tm = composition.melting_temperature(seq_long)
    assert tm > 12.0


