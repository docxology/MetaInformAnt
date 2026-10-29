"""Regression test: ``call_structural_variants`` genotypes with the caller's ``min_mapq``.

The call keeps the reads with MAPQ >= ``min_mapq`` as evidence, but then genotyped with ``genotype_sv``'s default
``min_mapq=20``. With ``min_mapq`` below 20, reads with MAPQ from ``min_mapq`` to 19 supported the variant yet were
not counted as reference reads, so the allele fraction was too high: five supporting and five reference reads
(a heterozygous 0/1) genotyped as 1/1.
"""

from __future__ import annotations

from typing import Any

import pytest

from metainformant.structural_variants.detection.sv_calling import (
    InsertSizeStats,
    SVType,
    call_structural_variants,
    genotype_sv,
)

NORMAL_PAIR = {"insert_size": 400, "is_reverse": False, "mate_is_reverse": True}
STATS = InsertSizeStats(mean=400.0, std=50.0)


def _split_reads(n: int, mapq: int) -> list[dict[str, Any]]:
    """``n`` split reads clipped at chr1:5000+i whose clipped part aligns to chr1:8001 on the other strand (a DEL)."""
    return [
        {
            "name": f"split_{i}",
            "chrom": "chr1",
            "pos": 4950 + i,
            "mapq": mapq,
            "cigar": "50M100S",
            "sa_tag": "chr1,8001,-,100S50M,60,0;",
            "mate_chrom": "chr1",
            "mate_pos": 5300,
            "read_length": 150,
            **NORMAL_PAIR,
        }
        for i in range(n)
    ]


def _reference_reads(n: int, mapq: int, prefix: str = "ref") -> list[dict[str, Any]]:
    """``n`` reads that span the first breakpoint (chr1:5002 +- 100) without a split."""
    return [{"name": f"{prefix}_{i}", "chrom": "chr1", "pos": 4800, "mapq": mapq, "read_length": 400} for i in range(n)]


def _genotype(reads: list[dict[str, Any]], **kwargs: Any) -> str:
    calls = call_structural_variants(reads, min_support=3, insert_size_stats=STATS, **kwargs)
    assert [c.sv_type for c in calls] == [SVType.DEL]
    assert calls[0].evidence.split_reads == 5
    return calls[0].genotype


class TestCallerMinMapqReachesTheGenotype:
    @pytest.mark.parametrize(
        ("min_mapq", "read_mapq"),
        [(10, 15), (10, 10), (0, 5), (0, 0)],
        ids=["mapq-15-min-10", "mapq-equal-to-min-10", "mapq-5-min-0", "mapq-0-min-0"],
    )
    def test_reads_at_or_above_the_callers_threshold_count_as_reference_reads(
        self, min_mapq: int, read_mapq: int
    ) -> None:
        # 5 supporting split reads against 5 reference reads: allele fraction 0.5, heterozygous. With the default
        # threshold of 20 in the genotype the reference reads (MAPQ below 20) did not count: fraction 1.0, "1/1".
        reads = _split_reads(5, read_mapq) + _reference_reads(5, read_mapq)
        assert _genotype(reads, min_mapq=min_mapq) == "0/1"

    def test_the_genotype_follows_the_number_of_reference_reads_that_pass(self) -> None:
        # 5 supporting reads and k reference reads at MAPQ 12 (min_mapq 10): fraction 5 / (5 + k), by hand:
        # k=0 -> 1.0 (1/1); k=5 -> 0.5 (0/1); k=30 -> 5/35 = 0.143 (0/0).
        for n_reference, expected in [(0, "1/1"), (5, "0/1"), (30, "0/0")]:
            reads = _split_reads(5, 12) + _reference_reads(n_reference, 12)
            assert _genotype(reads, min_mapq=10) == expected

    def test_reads_below_the_callers_threshold_are_still_ignored(self) -> None:
        # min_mapq 30: five reference reads at MAPQ 60 count, ten at MAPQ 25 do not.
        reads = _split_reads(5, 60) + _reference_reads(5, 60) + _reference_reads(10, 25, prefix="low")
        assert _genotype(reads, min_mapq=30) == "0/1"

    def test_the_default_threshold_behaves_as_before(self) -> None:
        reads = _split_reads(5, 60) + _reference_reads(5, 60)
        assert _genotype(reads) == "0/1"  # min_mapq=20 either way
        assert _genotype(reads, min_mapq=20) == "0/1"

    def test_genotype_sv_called_directly_keeps_its_own_default(self) -> None:
        reads = _split_reads(5, 15) + _reference_reads(5, 15)
        [call] = call_structural_variants(reads, min_support=3, insert_size_stats=STATS, min_mapq=10)
        assert genotype_sv(call, reads) == "1/1"  # its default threshold of 20 drops the MAPQ-15 reference reads
        assert genotype_sv(call, reads, min_mapq=10) == "0/1"
