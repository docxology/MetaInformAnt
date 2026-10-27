"""Regression tests for ``call_structural_variants``, ``genotype_sv`` and ``plot_circos``.

Three defects are pinned here, each with real reads and real calls (no test doubles):

* A translocation joins loci on two chromosomes. It used to be reported with ``start`` and ``end`` taken
  from the two breakpoints and sorted as if they were one interval, so ``start`` could be the position on
  the *second* chromosome, and ``plot_circos`` then used ``end`` as the second chromosome's position. A
  call now carries ``chrom``/``start`` for the first breakpoint and ``chrom2``/``pos2`` for the second, and
  ``plot_circos`` draws the chord to ``pos2``.
* An insertion's two breakpoints coincide (split reads at one junction), so its breakpoint distance is
  always under 50 and the default ``min_sv_size=50`` dropped every insertion call. The inserted length is
  not measured by that evidence, so insertions are no longer filtered on a size that is not measured.
* ``genotype_sv`` scanned every read for every call. A ``ReadIndex`` answers the same question from a
  position index, so a call costs a binary search instead of a pass over all reads.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pytest

try:
    import matplotlib

    matplotlib.use("Agg")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from metainformant.structural_variants.detection.sv_calling import (
    InsertSizeStats,
    StructuralVariant,
    SVEvidence,
    SVType,
    call_structural_variants,
    genotype_sv,
)

NORMAL_PAIR = {"insert_size": 400, "is_reverse": False, "mate_is_reverse": True}
STATS = InsertSizeStats(mean=400.0, std=50.0)


def _tra_reads(sa_pos: int, n: int = 5, chrom: str = "chr1", sa_chrom: str = "chr2") -> list[dict[str, Any]]:
    """Split reads on ``chrom`` whose clipped part aligns to ``sa_chrom:sa_pos`` (1-based SA position)."""
    return [
        {
            "name": f"tra_{i}",
            "chrom": chrom,
            "pos": 4950 + i,  # right clip at 5000 + i
            "mapq": 60,
            "cigar": "50M100S",
            "sa_tag": f"{sa_chrom},{sa_pos},+,100S50M,60,0;",
            "mate_chrom": chrom,
            "mate_pos": 5300,
            "read_length": 150,
            **NORMAL_PAIR,
        }
        for i in range(n)
    ]


def _ins_reads(junction: int = 1000, n: int = 5, clip: int = 60) -> list[dict[str, Any]]:
    """Split reads that all clip at one junction and carry no supplementary alignment."""
    return [
        {
            "name": f"ins_{i}",
            "chrom": "chr1",
            "pos": junction - 40,
            "mapq": 60,
            "cigar": f"40M{clip}S",
            "mate_chrom": "chr1",
            "mate_pos": junction + 300,
            "read_length": 40 + clip,
            **NORMAL_PAIR,
        }
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# Translocations
# ---------------------------------------------------------------------------


class TestTranslocationRepresentation:
    @pytest.mark.parametrize("sa_pos", [3000, 9000], ids=["second-breakpoint-lower", "second-breakpoint-higher"])
    def test_each_breakpoint_stays_on_its_own_chromosome(self, sa_pos: int) -> None:
        calls = call_structural_variants(_tra_reads(sa_pos), min_support=3, insert_size_stats=STATS)

        assert [c.sv_type for c in calls] == [SVType.TRA]
        tra = calls[0]
        # First breakpoint: chr1, the median clip position (5000..5004).
        assert (tra.chrom, tra.start) == ("chr1", 5002)
        # Second breakpoint: chr2 at the 0-based SA position, whatever its value relative to the first.
        assert (tra.chrom2, tra.pos2) == ("chr2", sa_pos - 1)
        # No interval and no length: the end repeats the start.
        assert tra.end == tra.start
        assert tra.size == 0

    def test_evidence_keeps_both_breakpoints(self) -> None:
        tra = call_structural_variants(_tra_reads(3000), min_support=3, insert_size_stats=STATS)[0]
        assert (tra.evidence.chrom1, tra.evidence.breakpoint1) == ("chr1", 5002)
        assert (tra.evidence.chrom2, tra.evidence.breakpoint2) == ("chr2", 2999)

    def test_pos2_is_only_set_for_translocations(self) -> None:
        reads = [
            {
                "name": f"del_{i}",
                "chrom": "chr1",
                "pos": 4800 + i * 20,
                "mapq": 55,
                "cigar": "150M",
                "mate_chrom": "chr1",
                "mate_pos": 8200 + i * 20,
                "read_length": 150,
                **dict(NORMAL_PAIR, insert_size=3400),
            }
            for i in range(5)
        ]
        calls = call_structural_variants(reads, min_support=3, insert_size_stats=STATS)
        assert [c.sv_type for c in calls] == [SVType.DEL]
        assert calls[0].pos2 is None
        assert calls[0].chrom2 == ""
        assert calls[0].start < calls[0].end

    def test_genotype_uses_reads_at_the_first_breakpoint(self) -> None:
        # Reference reads span chr1:5002 (the first breakpoint); none span chr1:2999, which is where the
        # old call put ``start`` when the second breakpoint's coordinate was the smaller one.
        spanning = [
            {"name": f"ref_{i}", "chrom": "chr1", "pos": 4800, "mapq": 60, "read_length": 600} for i in range(20)
        ]
        calls = call_structural_variants(_tra_reads(3000) + spanning, min_support=3, insert_size_stats=STATS)
        tra = calls[0]
        assert tra.sv_type == SVType.TRA
        # 5 supporting reads against 20 reference reads: allele fraction 0.2 -> heterozygous.
        assert tra.genotype == "0/1"


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestCircosTranslocation:
    CHROMOSOMES = {"chr1": 20_000, "chr2": 20_000}

    @classmethod
    def _expected_angle(cls, chrom: str, pos: int) -> float:
        """Polar angle of ``chrom:pos`` in the documented layout: arcs proportional to length, 2% gaps."""
        gap = 2 * math.pi * 0.02
        available = 2 * math.pi * (1 - 0.02 * len(cls.CHROMOSOMES))
        total = sum(cls.CHROMOSOMES.values())
        arc_start = 0.0
        for name in sorted(cls.CHROMOSOMES):
            arc = available * cls.CHROMOSOMES[name] / total
            if name == chrom:
                return arc_start + pos / cls.CHROMOSOMES[name] * arc
            arc_start += arc + gap
        raise AssertionError(chrom)

    def _chords(self, variants: list[dict[str, Any]], tmp_path: Path) -> list[tuple[float, float]]:
        from metainformant.structural_variants.visualization.plots import plot_circos

        fig = plot_circos(variants, self.CHROMOSOMES, tmp_path / "circos.png")
        return [(float(line.get_xdata()[0]), float(line.get_xdata()[-1])) for line in fig.axes[0].lines]

    def test_chord_ends_at_pos2_of_a_called_translocation(self, tmp_path: Path) -> None:
        tra = call_structural_variants(_tra_reads(3000), min_support=3, insert_size_stats=STATS)[0]
        [(start_angle, end_angle)] = self._chords([asdict(tra)], tmp_path)
        assert start_angle == pytest.approx(self._expected_angle("chr1", 5002))
        assert end_angle == pytest.approx(self._expected_angle("chr2", 2999))

    def test_chord_to_pos2_not_end(self, tmp_path: Path) -> None:
        variant = {"chrom": "chr1", "start": 4000, "end": 4000, "sv_type": "TRA", "chrom2": "chr2", "pos2": 1000}
        [(start_angle, end_angle)] = self._chords([variant], tmp_path)
        assert start_angle == pytest.approx(self._expected_angle("chr1", 4000))
        assert end_angle == pytest.approx(self._expected_angle("chr2", 1000))

    def test_chord_of_a_dict_without_pos2_ends_at_end(self, tmp_path: Path) -> None:
        # The VCF-style layout (END holds the position on CHR2) keeps working.
        variant = {"chrom": "chr1", "start": 4000, "end": 1000, "sv_type": "TRA", "chrom2": "chr2"}
        [(start_angle, end_angle)] = self._chords([variant], tmp_path)
        assert start_angle == pytest.approx(self._expected_angle("chr1", 4000))
        assert end_angle == pytest.approx(self._expected_angle("chr2", 1000))

    def test_intra_chromosomal_chord_ignores_empty_chrom2_and_none_pos2(self, tmp_path: Path) -> None:
        variant = {"chrom": "chr1", "start": 1000, "end": 6000, "sv_type": "DEL", "chrom2": "", "pos2": None}
        [(start_angle, end_angle)] = self._chords([variant], tmp_path)
        assert start_angle == pytest.approx(self._expected_angle("chr1", 1000))
        assert end_angle == pytest.approx(self._expected_angle("chr1", 6000))


# ---------------------------------------------------------------------------
# Insertions
# ---------------------------------------------------------------------------


class TestInsertionsAreNotFilteredOnUnmeasuredSize:
    def test_insertion_survives_the_default_min_sv_size(self) -> None:
        calls = call_structural_variants(_ins_reads(), min_support=3, insert_size_stats=STATS)

        assert [c.sv_type for c in calls] == [SVType.INS]
        ins = calls[0]
        assert (ins.chrom, ins.start, ins.end) == ("chr1", 1000, 1000)
        # The inserted length is not measured by split-read evidence: 0 means "not measured".
        assert ins.size == 0
        assert ins.evidence.split_reads == 5

    @pytest.mark.parametrize("min_sv_size", [0, 50, 10_000])
    def test_insertion_is_independent_of_min_sv_size(self, min_sv_size: int) -> None:
        calls = call_structural_variants(_ins_reads(), min_support=3, insert_size_stats=STATS, min_sv_size=min_sv_size)
        assert [c.sv_type for c in calls] == [SVType.INS]

    def test_the_insertion_survives_the_size_filter_of_the_standard_pipeline(self) -> None:
        from metainformant.structural_variants.filtering.quality_filter import filter_by_size

        calls = call_structural_variants(_ins_reads(), min_support=3, insert_size_stats=STATS)
        records = [asdict(call) for call in calls]
        kept, _ = filter_by_size(records, min_size=50)
        assert [r["sv_type"] for r in kept] == [SVType.INS]

    def test_min_sv_size_still_filters_measurable_variants(self) -> None:
        reads = [
            {
                "name": f"del_{i}",
                "chrom": "chr1",
                "pos": 4800 + i * 20,
                "mapq": 55,
                "cigar": "150M",
                "mate_chrom": "chr1",
                "mate_pos": 8200 + i * 20,
                "read_length": 150,
                **dict(NORMAL_PAIR, insert_size=3400),
            }
            for i in range(5)
        ]
        kept = call_structural_variants(reads, min_support=3, insert_size_stats=STATS, min_sv_size=3000)
        dropped = call_structural_variants(reads, min_support=3, insert_size_stats=STATS, min_sv_size=5000)
        assert [c.sv_type for c in kept] == [SVType.DEL]
        assert 3000 <= kept[0].size < 5000
        assert dropped == []


# ---------------------------------------------------------------------------
# Genotyping against a read index
# ---------------------------------------------------------------------------


def _random_reads(rng: random.Random, n: int) -> list[dict[str, Any]]:
    reads: list[dict[str, Any]] = []
    for i in range(n):
        read: dict[str, Any] = {"chrom": rng.choice(["chr1", "chr2", "chr3"]), "pos": rng.randint(0, 20_000)}
        # Leave keys out at random so the defaults (read_length 150, mapq 0, name "") are covered.
        if rng.random() < 0.9:
            read["read_length"] = rng.choice([36, 100, 150, 250, 2_000])
        if rng.random() < 0.9:
            read["mapq"] = rng.choice([0, 5, 19, 20, 40, 60])
        if rng.random() < 0.8:
            read["name"] = f"r{i % 400}"
        reads.append(read)
    return reads


class TestReadIndexGenotyping:
    def test_index_gives_the_same_genotype_as_scanning_the_reads(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import ReadIndex

        rng = random.Random(20261027)
        reads = _random_reads(rng, 3_000)
        index = ReadIndex(reads)
        assert len(index) == len(reads)

        seen: set[str] = set()
        for _ in range(400):
            names = [f"r{rng.randint(0, 399)}" for _ in range(rng.randint(0, 6))]
            variant = StructuralVariant(
                chrom=rng.choice(["chr1", "chr2", "chr3", "chrX"]),
                start=rng.randint(0, 20_000),
                end=20_000,
                sv_type=SVType.DEL,
                evidence=SVEvidence(
                    split_reads=rng.randint(0, 8), discordant_pairs=rng.randint(0, 4), evidence_reads=names
                ),
            )
            min_mapq = rng.choice([0, 10, 20, 40])
            expected = genotype_sv(variant, reads, min_mapq=min_mapq)
            assert genotype_sv(variant, index, min_mapq=min_mapq) == expected
            seen.add(expected)
        # The comparison is not trivial: every genotype class occurs.
        assert {"0/0", "0/1", "1/1"} <= seen

    def test_count_spanning_matches_a_direct_count(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import ReadIndex

        rng = random.Random(7)
        reads = _random_reads(rng, 2_000)
        index = ReadIndex(reads)
        for _ in range(200):
            chrom = rng.choice(["chr1", "chr2", "chr3", "chr9"])
            window_start = rng.randint(0, 20_000)
            window_end = window_start + rng.randint(0, 400)
            min_mapq = rng.choice([0, 20, 40])
            direct = sum(
                1
                for read in reads
                if read.get("chrom", "") == chrom
                and read.get("mapq", 0) >= min_mapq
                and read.get("pos", 0) < window_start
                and read.get("pos", 0) + read.get("read_length", 150) > window_end
            )
            assert index.count_spanning(chrom, window_start, window_end, min_mapq) == direct

    def test_called_genotypes_are_unchanged(self) -> None:
        """Every call's genotype equals the one from scanning the reads that passed the MAPQ filter."""
        reads = _tra_reads(3000, n=6) + _ins_reads(2000, n=7)
        reads += [
            {"name": f"ref_{i}", "chrom": "chr1", "pos": 1850 + 10 * i, "mapq": 60, "read_length": 400}
            for i in range(12)
        ]
        calls = call_structural_variants(reads, min_support=3, insert_size_stats=STATS)
        assert {c.sv_type for c in calls} == {SVType.TRA, SVType.INS}
        assert len({c.genotype for c in calls}) == 2  # reference reads span only the insertion junction
        for call in calls:
            assert call.genotype == genotype_sv(call, [r for r in reads if r.get("mapq", 0) >= 20])


def _scaling_workload(n_background: int, n_calls: int) -> list[dict[str, Any]]:
    """``n_background`` ordinary reads plus ``n_calls`` independent deletions (split reads with an SA tag)."""
    rng = random.Random(11)
    chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"]
    reads: list[dict[str, Any]] = []
    for i in range(n_background):
        chrom = chroms[i % len(chroms)]
        pos = rng.randint(1_000, 4_000_000)
        reads.append(
            {
                "name": f"bg_{i}",
                "chrom": chrom,
                "pos": pos,
                "mapq": 60,
                "cigar": "150M",
                "mate_chrom": chrom,
                "mate_pos": pos + 300,
                "read_length": 150,
                **NORMAL_PAIR,
            }
        )
    for call in range(n_calls):
        chrom = chroms[call % len(chroms)]
        base = 50_000 + call * 15_000
        for k in range(6):
            reads.append(
                {
                    "name": f"sv{call}_r{k}",
                    "chrom": chrom,
                    "pos": base - 90 + k,
                    "mapq": 60,
                    "cigar": "90M60S",
                    "sa_tag": f"{chrom},{base + 2000},-,60S90M,60,0;",
                    "mate_chrom": chrom,
                    "mate_pos": base + 300,
                    "read_length": 150,
                    **NORMAL_PAIR,
                }
            )
    return reads


def _best_time(reads: list[dict[str, Any]], runs: int = 3) -> tuple[float, int]:
    best = float("inf")
    n_calls = 0
    for _ in range(runs):
        started = time.perf_counter()
        n_calls = len(call_structural_variants(reads, min_support=3, insert_size_stats=STATS))
        best = min(best, time.perf_counter() - started)
    return best, n_calls


def test_genotyping_cost_does_not_grow_with_the_number_of_calls() -> None:
    """With 100k reads, 200 calls must cost about what 10 calls cost.

    Scanning every read for every call made the run time grow linearly with the number of calls (about 90% of
    the time at 100 calls); the ratio below was about 9 then and is under 2 with the read index (the rest is
    the per-call evidence clustering). The two workloads share their 100k background reads, so the comparison
    does not depend on the machine's speed.
    """
    few_time, few_calls = _best_time(_scaling_workload(100_000, 10))
    many_time, many_calls = _best_time(_scaling_workload(100_000, 200))

    assert (few_calls, many_calls) == (10, 200)
    assert many_time / few_time < 4.0, f"10 calls: {few_time:.2f}s, 200 calls: {many_time:.2f}s"
