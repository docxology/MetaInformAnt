"""Regression tests for the read-phasing pipeline.

``haplotag_reads`` tagged every read HP=0 for blocks from ``phase_reads`` or
``build_phase_blocks`` because those blocks carried no ``variant_haplotypes``, and
it only ever looked at the first block a read overlapped. ``analysis.phasing.phase_reads``
documented ``aligned_pairs`` input and a ``min_variant_quality`` threshold but read
neither. All tests run the real functions on small synthetic datasets; there are no mocks.
"""

from __future__ import annotations

import logging
from typing import Any

import pytest

from metainformant.longread.analysis import phasing as analysis_phasing
from metainformant.longread.phasing import haplotyping

# ---------------------------------------------------------------------------
# phasing.haplotyping: phase_reads -> haplotag_reads
# ---------------------------------------------------------------------------


def _dict_variants(chrom: str, positions: list[int]) -> list[dict[str, Any]]:
    return [{"chrom": chrom, "position": p, "ref": "A", "alt": "G"} for p in positions]


def _dict_reads(
    prefix: str, chrom: str, positions: list[int], patterns: list[list[int]], copies: int = 10
) -> list[dict[str, Any]]:
    """Reads spanning ``positions`` that each carry one of the allele ``patterns`` (``copies`` reads per pattern)."""
    reads = []
    for pattern_index, pattern in enumerate(patterns):
        for copy in range(copies):
            reads.append(
                {
                    "read_id": f"{prefix}{pattern_index}_{copy}",
                    "chrom": chrom,
                    "start": min(positions) - 10,
                    "end": max(positions) + 10,
                    "alleles": [{"position": p, "allele": a} for p, a in zip(positions, pattern)],
                }
            )
    return reads


ANTI_PHASE = [[0, 1, 0, 1, 1, 0], [1, 0, 1, 0, 0, 1]]  # two haplotypes that differ at every variant


class TestHaplotagFromPhaseReadsBlocks:
    def test_phase_reads_blocks_carry_variant_haplotypes(self) -> None:
        positions = [100, 200, 300, 400, 500, 600]
        result = haplotyping.phase_reads(
            _dict_variants("chr1", positions), _dict_reads("r", "chr1", positions, ANTI_PHASE)
        )
        (block,) = result["phase_blocks"]
        haplotype_map = block["variant_haplotypes"]
        assert sorted(haplotype_map) == positions
        assert set(haplotype_map.values()) <= {0, 1}
        # the map is one of the two haplotypes of the data (the labelling of haplotype 1 is arbitrary)
        assert [haplotype_map[p] for p in positions] in ANTI_PHASE

    def test_haplotag_reproduces_the_read_assignments_of_phase_reads(self) -> None:
        positions = [100, 200, 300, 400, 500, 600]
        reads = _dict_reads("r", "chr1", positions, ANTI_PHASE)
        result = haplotyping.phase_reads(_dict_variants("chr1", positions), reads)
        assert result["n_phased_reads"] == len(reads)

        tagged = haplotyping.haplotag_reads(reads, result["phase_blocks"])

        assert [r["HP"] for r in tagged] == [result["haplotype_assignments"][r["read_id"]] for r in reads]
        assert {r["HP"] for r in tagged} == {1, 2}
        assert all(r["PS"] == 0 for r in tagged)
        # reads with the same allele pattern share a haplotype, the two patterns do not
        by_pattern = {
            pattern: {r["HP"] for r in tagged if r["read_id"].startswith(f"r{pattern}_")} for pattern in (0, 1)
        }
        assert by_pattern[0] != by_pattern[1] and all(len(v) == 1 for v in by_pattern.values())

    def test_every_read_is_tagged_when_blocks_interleave(self) -> None:
        # reads link variants 100/300/500 and variants 200/400/600 separately: two blocks whose spans overlap
        x_positions, y_positions = [100, 300, 500], [200, 400, 600]
        reads = _dict_reads("x", "chr1", x_positions, [[0, 0, 0], [1, 1, 1]]) + _dict_reads(
            "y", "chr1", y_positions, [[0, 1, 0], [1, 0, 1]]
        )
        variants = _dict_variants("chr1", sorted(x_positions + y_positions))
        result = haplotyping.phase_reads(variants, reads)
        assert len(result["phase_blocks"]) == 2

        tagged = haplotyping.haplotag_reads(reads, result["phase_blocks"])

        assert all(r["HP"] in (1, 2) for r in tagged)
        assert [r["HP"] for r in tagged] == [result["haplotype_assignments"][r["read_id"]] for r in reads]
        # PS is the index of the block that holds the read's variants
        block_of = {p: i for i, block in enumerate(result["phase_blocks"]) for p in block["variant_haplotypes"]}
        for read, out in zip(reads, tagged):
            assert out["PS"] == block_of[read["alleles"][0]["position"]]

    def test_blocks_on_two_chromosomes(self) -> None:
        positions = [100, 200, 300, 400, 500, 600]
        reads = _dict_reads("a", "chr1", positions, ANTI_PHASE) + _dict_reads("b", "chr2", positions, ANTI_PHASE)
        variants = _dict_variants("chr1", positions) + _dict_variants("chr2", positions)
        result = haplotyping.phase_reads(variants, reads)
        assert sorted(b["chrom"] for b in result["phase_blocks"]) == ["chr1", "chr2"]

        tagged = haplotyping.haplotag_reads(reads, result["phase_blocks"])

        assert all(r["HP"] in (1, 2) for r in tagged)
        assert [r["HP"] for r in tagged] == [result["haplotype_assignments"][r["read_id"]] for r in reads]
        for out in tagged:
            assert result["phase_blocks"][out["PS"]]["chrom"] == out["chrom"]


class TestHaplotagBlockSelection:
    def test_read_takes_the_block_that_holds_its_variants_not_the_first_overlapping_one(self) -> None:
        first = {"chrom": "chr1", "start": 100, "end": 500, "variant_haplotypes": {100: 0, 300: 0, 500: 0}}
        second = {"chrom": "chr1", "start": 200, "end": 600, "variant_haplotypes": {200: 0, 400: 1, 600: 0}}
        read = {
            "read_id": "r",
            "chrom": "chr1",
            "start": 150,
            "end": 650,
            "alleles": [{"position": 200, "allele": 1}, {"position": 400, "allele": 0}, {"position": 600, "allele": 1}],
        }
        (tagged,) = haplotyping.haplotag_reads([read], [first, second])
        assert (tagged["HP"], tagged["PS"]) == (2, 1)

    def test_the_block_with_the_most_votes_wins(self) -> None:
        weak = {"chrom": "chr1", "start": 100, "end": 400, "variant_haplotypes": {100: 0}}
        strong = {"chrom": "chr1", "start": 100, "end": 400, "variant_haplotypes": {200: 0, 300: 0}}
        read = {
            "read_id": "r",
            "chrom": "chr1",
            "start": 90,
            "end": 410,
            "alleles": [{"position": 100, "allele": 1}, {"position": 200, "allele": 0}, {"position": 300, "allele": 0}],
        }
        (tagged,) = haplotyping.haplotag_reads([read], [weak, strong])
        assert (tagged["HP"], tagged["PS"]) == (1, 1)

    def test_read_starting_on_the_last_variant_of_a_block_is_tagged(self) -> None:
        block = {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 1}}
        # the read covers position 300 (start 300, exclusive end 400) and carries the haplotype-1 allele there
        read = {"read_id": "r", "chrom": "chr1", "start": 300, "end": 400, "alleles": [{"position": 300, "allele": 1}]}
        (tagged,) = haplotyping.haplotag_reads([read], [block])
        assert (tagged["HP"], tagged["PS"]) == (1, 0)

    def test_read_ending_at_the_first_variant_is_not_tagged(self) -> None:
        block = {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 1}}
        read = {"read_id": "r", "chrom": "chr1", "start": 0, "end": 100, "alleles": [{"position": 100, "allele": 0}]}
        (tagged,) = haplotyping.haplotag_reads([read], [block])
        assert (tagged["HP"], tagged["PS"]) == (0, 0)

    def test_tied_or_missing_votes_leave_the_read_untagged(self) -> None:
        block = {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 0}}
        tie = {
            "read_id": "tie",
            "chrom": "chr1",
            "start": 90,
            "end": 310,
            "alleles": [{"position": 100, "allele": 0}, {"position": 300, "allele": 1}],
        }
        no_alleles = {"read_id": "none", "chrom": "chr1", "start": 90, "end": 310}
        elsewhere = {
            "read_id": "far",
            "chrom": "chr1",
            "start": 90,
            "end": 310,
            "alleles": [{"position": 5, "allele": 0}],
        }
        tagged = haplotyping.haplotag_reads([tie, no_alleles, elsewhere], [block])
        assert [(t["HP"], t["PS"]) for t in tagged] == [(0, 0)] * 3

    def test_blocks_without_variant_haplotypes_are_reported(self, caplog: pytest.LogCaptureFixture) -> None:
        blocks = [
            {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 1}},
            {"chrom": "chr1", "start": 1000, "end": 1300},
        ]
        read = {
            "read_id": "r",
            "chrom": "chr1",
            "start": 990,
            "end": 1310,
            "alleles": [{"position": 1000, "allele": 0}],
        }
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            (tagged,) = haplotyping.haplotag_reads([read], blocks)
        assert tagged["HP"] == 0
        assert any("1 of 2 phase blocks carry no variant_haplotypes" in r.getMessage() for r in caplog.records)

    def test_reads_are_copied_not_modified(self) -> None:
        block = {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 0}}
        read = {"read_id": "r", "chrom": "chr1", "start": 90, "end": 310, "alleles": [{"position": 100, "allele": 0}]}
        haplotyping.haplotag_reads([read], [block])
        assert "HP" not in read and "PS" not in read


class TestBuildPhaseBlocksHaplotypeMap:
    @staticmethod
    def _phased(haplotypes: list[Any], positions: list[int] | None = None, **extra: Any) -> list[dict[str, Any]]:
        positions = positions or [100 * (i + 1) for i in range(len(haplotypes))]
        return [{"chrom": "chr1", "position": p, "haplotype": h, **extra} for p, h in zip(positions, haplotypes)]

    def test_haplotype_names_the_haplotype_that_carries_the_alternate_allele(self) -> None:
        (block,) = haplotyping.build_phase_blocks(self._phased([1, 1, 2, 2]))
        # alt on haplotype 1 -> haplotype 1 carries allele 1; alt on haplotype 2 -> haplotype 1 carries allele 0
        assert block["variant_haplotypes"] == {100: 1, 200: 1, 300: 0, 400: 0}

    def test_explicit_allele_names_the_allele_carried_by_the_haplotype(self) -> None:
        phased = [
            {"chrom": "chr1", "position": 100, "haplotype": 1, "allele": 1},
            {"chrom": "chr1", "position": 200, "haplotype": 1, "allele": 0},
            {"chrom": "chr1", "position": 300, "haplotype": 2, "allele": 1},
            {"chrom": "chr1", "position": 400, "haplotype": 2, "allele": 0},
        ]
        (block,) = haplotyping.build_phase_blocks(phased)
        assert block["variant_haplotypes"] == {100: 1, 200: 0, 300: 0, 400: 1}

    def test_variants_without_a_usable_haplotype_or_allele_are_left_out(self) -> None:
        phased = [
            {"chrom": "chr1", "position": 100, "haplotype": 1},
            {"chrom": "chr1", "position": 200},
            {"chrom": "chr1", "position": 300, "haplotype": 3},
            {"chrom": "chr1", "position": 400, "haplotype": 2, "allele": 2},
            {"chrom": "chr1", "position": 500, "haplotype": 2},
        ]
        (block,) = haplotyping.build_phase_blocks(phased)
        assert block["variant_haplotypes"] == {100: 1, 500: 0}
        assert block["n_variants"] == 5  # grouping is unchanged

    def test_build_phase_blocks_output_haplotags_reads(self) -> None:
        phased = self._phased([1, 1, 2, 2])
        blocks = haplotyping.build_phase_blocks(phased)
        carries_alt_on_1 = [100, 200]
        carries_alt_on_2 = [300, 400]

        def read(read_id: str, alt_positions: list[int]) -> dict[str, Any]:
            return {
                "read_id": read_id,
                "chrom": "chr1",
                "start": 50,
                "end": 450,
                "alleles": [{"position": p, "allele": int(p in alt_positions)} for p in (100, 200, 300, 400)],
            }

        reads = [read("hap1", carries_alt_on_1), read("hap2", carries_alt_on_2)]
        tagged = haplotyping.haplotag_reads(reads, blocks)
        assert [(t["read_id"], t["HP"], t["PS"]) for t in tagged] == [("hap1", 1, 0), ("hap2", 2, 0)]

    def test_each_block_gets_its_own_map(self) -> None:
        phased = self._phased([1, 2], [100, 200]) + self._phased([2, 1], [5_000_000, 5_000_100])
        blocks = haplotyping.build_phase_blocks(phased, max_gap=1000)
        assert [b["variant_haplotypes"] for b in blocks] == [{100: 1, 200: 0}, {5_000_000: 0, 5_000_100: 1}]


# ---------------------------------------------------------------------------
# analysis.phasing.phase_reads: aligned_pairs and min_variant_quality
# ---------------------------------------------------------------------------

POSITIONS = [100, 120, 140, 160, 180]
HAPLOTYPES = [[0, 1, 0, 1, 0], [1, 0, 1, 0, 1]]  # anti-phase pairs of haplotypes


def _snv_variants(positions: list[int] = POSITIONS, **extra: Any) -> list[dict[str, Any]]:
    return [{"chromosome": "chr1", "position": p, "ref_allele": "A", "alt_allele": "G", **extra} for p in positions]


def _aligned_read(name: str, alleles: dict[int, int], start: int = 90, length: int = 110) -> dict[str, Any]:
    """A read with an ungapped alignment whose bases carry the given allele (0 = A, 1 = G) at each position."""
    bases = ["C"] * length
    for position, allele in alleles.items():
        bases[position - start] = "AG"[allele]
    return {
        "read_name": name,
        "reference_name": "chr1",
        "reference_start": start,
        "reference_end": start + length,
        "query_sequence": "".join(bases),
        "aligned_pairs": [(q, start + q) for q in range(length)],
    }


def _variant_read(
    name: str, alleles: dict[int, int], quality: int = 30, start: int = 90, end: int = 200
) -> dict[str, Any]:
    return {
        "read_name": name,
        "reference_name": "chr1",
        "reference_start": start,
        "reference_end": end,
        "variants": [{"position": p, "allele": a, "quality": quality} for p, a in alleles.items()],
    }


def _haplotype_alleles(pattern: list[int]) -> dict[int, int]:
    return dict(zip(POSITIONS, pattern))


class TestPhaseReadsFromAlignedPairs:
    def _alignment_reads(self) -> list[dict[str, Any]]:
        return [
            _aligned_read(f"h{h}_{copy}", _haplotype_alleles(pattern))
            for h, pattern in enumerate(HAPLOTYPES)
            for copy in range(8)
        ]

    def test_reads_with_only_aligned_pairs_are_phased(self) -> None:
        blocks, assignments = analysis_phasing.phase_reads(_snv_variants(), self._alignment_reads())

        assert len(blocks) == 1
        assert blocks[0].variants == POSITIONS
        assert len(assignments) == 16
        by_haplotype = {h: {assignments[f"h{h}_{copy}"] for copy in range(8)} for h in (0, 1)}
        assert all(len(v) == 1 for v in by_haplotype.values())
        assert by_haplotype[0] != by_haplotype[1]
        # the block's two haplotypes are the two anti-phase patterns of the data
        assert sorted([blocks[0].haplotype1, blocks[0].haplotype2]) == sorted(HAPLOTYPES)

    def test_alignment_and_variant_list_inputs_give_the_same_phasing(self) -> None:
        from_alignment = analysis_phasing.phase_reads(_snv_variants(), self._alignment_reads())
        explicit_reads = [
            _variant_read(f"h{h}_{copy}", _haplotype_alleles(pattern))
            for h, pattern in enumerate(HAPLOTYPES)
            for copy in range(8)
        ]
        from_variants = analysis_phasing.phase_reads(_snv_variants(), explicit_reads)
        assert from_alignment == from_variants

    def test_alleles_follow_the_aligned_pairs_not_the_read_offset(self) -> None:
        # soft-clipped head (query bases with no reference position) and a deletion shift query vs reference
        def gapped(name: str, pattern: list[int]) -> dict[str, Any]:
            read = _aligned_read(name, _haplotype_alleles(pattern), start=90, length=110)
            pairs: list[tuple[int | None, int | None]] = [(q, None) for q in range(5)]
            pairs += [(q + 5, 90 + q) for q in range(110)]
            sequence = "T" * 5 + read["query_sequence"]
            read.update(query_sequence=sequence, aligned_pairs=pairs)
            return read

        reads = [gapped(f"h{h}_{copy}", p) for h, p in enumerate(HAPLOTYPES) for copy in range(6)]
        blocks, assignments = analysis_phasing.phase_reads(_snv_variants(), reads)
        assert len(blocks) == 1 and blocks[0].variants == POSITIONS
        assert assignments["h0_0"] != assignments["h1_0"]

    def test_bases_that_match_neither_allele_are_not_observations(self) -> None:
        read = _aligned_read("r", _haplotype_alleles(HAPLOTYPES[0]))
        # a third base at the first two variants leaves three informative variants
        sequence = list(read["query_sequence"])
        sequence[100 - 90] = "T"
        sequence[120 - 90] = "T"
        read["query_sequence"] = "".join(sequence)
        reads = [dict(read, read_name=f"r{i}") for i in range(4)]
        blocks, assignments = analysis_phasing.phase_reads(_snv_variants(), reads)
        assert len(blocks) == 1 and blocks[0].variants == POSITIONS[2:]
        assert len(assignments) == 4

    def test_low_quality_bases_are_excluded_by_min_base_quality(self) -> None:
        reads = []
        for h, pattern in enumerate(HAPLOTYPES):
            for copy in range(6):
                read = _aligned_read(f"h{h}_{copy}", _haplotype_alleles(pattern))
                qualities = [30] * 110
                qualities[100 - 90] = 5  # first variant
                read["query_qualities"] = qualities
                reads.append(read)
        blocks, _ = analysis_phasing.phase_reads(_snv_variants(), reads, min_base_quality=10)
        assert blocks[0].variants == POSITIONS[1:]
        blocks, _ = analysis_phasing.phase_reads(_snv_variants(), reads, min_base_quality=0)
        assert blocks[0].variants == POSITIONS

    def test_explicit_variant_observations_take_precedence(self) -> None:
        # the alignment says haplotype pattern A while the explicit list says the opposite
        reads = []
        for h, pattern in enumerate(HAPLOTYPES):
            for copy in range(6):
                read = _aligned_read(f"h{h}_{copy}", _haplotype_alleles(HAPLOTYPES[0]))
                read["variants"] = [
                    {"position": p, "allele": a, "quality": 30} for p, a in _haplotype_alleles(pattern).items()
                ]
                reads.append(read)
        _, assignments = analysis_phasing.phase_reads(_snv_variants(), reads)
        assert assignments["h0_0"] != assignments["h1_0"]

    def test_reads_without_a_query_sequence_are_reported_and_skipped(self, caplog: pytest.LogCaptureFixture) -> None:
        reads = self._alignment_reads()
        for read in reads[:2]:
            read["query_sequence"] = ""
        with caplog.at_level(logging.WARNING, logger=analysis_phasing.logger.name):
            _, assignments = analysis_phasing.phase_reads(_snv_variants(), reads)
        assert len(assignments) == 14 and "h0_0" not in assignments
        assert any("2 reads carry aligned_pairs but no query_sequence" in r.getMessage() for r in caplog.records)

    def test_variants_that_are_not_single_nucleotide_changes_cannot_be_called_and_are_reported(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        variants = _snv_variants()
        variants[2] = {"chromosome": "chr1", "position": 140, "ref_allele": "A", "alt_allele": "GT"}
        with caplog.at_level(logging.WARNING, logger=analysis_phasing.logger.name):
            blocks, _ = analysis_phasing.phase_reads(variants, self._alignment_reads())
        assert blocks[0].variants == [100, 120, 160, 180]
        assert any("1 of 5 variants are not single-nucleotide" in r.getMessage() for r in caplog.records)


class TestMinVariantQuality:
    def _reads(self) -> list[dict[str, Any]]:
        return [
            _variant_read(f"h{h}_{copy}", _haplotype_alleles(pattern))
            for h, pattern in enumerate(HAPLOTYPES)
            for copy in range(8)
        ]

    def _variants_with_quality(self, qualities: list[float | None]) -> list[dict[str, Any]]:
        variants = _snv_variants()
        for variant, quality in zip(variants, qualities):
            if quality is not None:
                variant["quality"] = quality
        return variants

    def test_variants_below_the_threshold_are_left_out_of_the_phasing(self) -> None:
        variants = self._variants_with_quality([50.0, 5.0, 50.0, 19.9, 50.0])
        blocks, assignments = analysis_phasing.phase_reads(variants, self._reads())  # default minimum 20
        assert len(blocks) == 1
        assert blocks[0].variants == [100, 140, 180]
        assert blocks[0].num_variants == 3
        assert len(assignments) == 16

    def test_the_threshold_is_configurable_and_inclusive(self) -> None:
        variants = self._variants_with_quality([50.0, 5.0, 50.0, 20.0, 50.0])
        blocks, _ = analysis_phasing.phase_reads(variants, self._reads(), min_variant_quality=20)
        assert blocks[0].variants == [100, 140, 160, 180]
        blocks, _ = analysis_phasing.phase_reads(variants, self._reads(), min_variant_quality=0)
        assert blocks[0].variants == POSITIONS
        blocks, _ = analysis_phasing.phase_reads(variants, self._reads(), min_variant_quality=60)
        assert blocks == []

    def test_variants_without_a_quality_are_kept(self) -> None:
        blocks, _ = analysis_phasing.phase_reads(_snv_variants(), self._reads(), min_variant_quality=1000)
        assert blocks[0].variants == POSITIONS

    def test_variant_objects_carry_their_quality(self) -> None:
        variants = [
            analysis_phasing.Variant("chr1", p, "A", "G", quality=q)
            for p, q in zip(POSITIONS, [50.0, 5.0, 50.0, 50.0, 50.0])
        ]
        blocks, _ = analysis_phasing.phase_reads(variants, self._reads())
        assert blocks[0].variants == [100, 140, 160, 180]

    def test_reads_at_an_excluded_variant_do_not_phase_through_it(self) -> None:
        # the middle variant is excluded, so reads linking only through it can no longer connect the others
        variants = self._variants_with_quality([50.0, 50.0, 1.0, 50.0, 50.0])
        reads = [_variant_read(f"a{i}", {100: 0, 120: 1, 140: 0}) for i in range(6)]
        reads += [_variant_read(f"b{i}", {140: 0, 160: 1, 180: 0}) for i in range(6)]
        blocks, _ = analysis_phasing.phase_reads(variants, reads)
        assert [b.variants for b in blocks] == [[100, 120], [160, 180]]


# ---------------------------------------------------------------------------
# analysis.phasing.tag_reads_by_haplotype: the same block selection as haplotag_reads
# ---------------------------------------------------------------------------


def _phase_block(block_id: int, variants: list[int], haplotype1: list[int]) -> analysis_phasing.PhaseBlock:
    return analysis_phasing.PhaseBlock(
        block_id=block_id,
        chromosome="chr1",
        start=min(variants),
        end=max(variants),
        variants=variants,
        haplotype1=haplotype1,
        haplotype2=[1 - a for a in haplotype1],
        num_variants=len(variants),
    )


class TestTagReadsByHaplotypeBlockSelection:
    def test_every_read_of_interleaved_phase_reads_blocks_is_tagged(self) -> None:
        # two components whose spans overlap: the second one's reads used to stay HP=0
        x_positions, y_positions = [100, 300, 500], [200, 400, 600]
        reads = []
        for i in range(6):
            reads.append(_variant_read(f"x{i}", {p: i % 2 for p in x_positions}, start=90, end=520))
            reads.append(_variant_read(f"y{i}", {200: i % 2, 400: 1 - i % 2, 600: i % 2}, start=190, end=620))
        variants = _snv_variants(sorted(x_positions + y_positions))
        blocks, assignments = analysis_phasing.phase_reads(variants, reads)
        assert [b.variants for b in blocks] == [x_positions, y_positions]
        assert len(assignments) == 12

        tagged = analysis_phasing.tag_reads_by_haplotype(reads, blocks)

        assert [t["HP"] for t in tagged] == [assignments[r["read_name"]] for r in reads]
        assert {t["PS"] for t in tagged if t["read_name"].startswith("y")} == {blocks[1].block_id}
        assert {t["PS"] for t in tagged if t["read_name"].startswith("x")} == {blocks[0].block_id}

    def test_read_takes_the_block_with_most_of_its_variants(self) -> None:
        first = _phase_block(0, [100, 300, 500], [0, 0, 0])
        second = _phase_block(1, [200, 400, 600], [0, 1, 0])
        read = {
            "read_name": "r",
            "reference_name": "chr1",
            "reference_start": 150,
            "reference_end": 650,
            "variants": [
                {"position": 200, "allele": 1},
                {"position": 400, "allele": 0},
                {"position": 600, "allele": 1},
            ],
        }
        (tagged,) = analysis_phasing.tag_reads_by_haplotype([read], [first, second])
        assert (tagged["HP"], tagged["PS"]) == (2, 1)

    def test_read_starting_on_the_last_variant_of_a_block_is_tagged(self) -> None:
        block = _phase_block(0, [100, 300], [0, 1])
        read = {
            "read_name": "r",
            "reference_name": "chr1",
            "reference_start": 300,
            "reference_end": 400,
            "variants": [{"position": 300, "allele": 1}],
        }
        (tagged,) = analysis_phasing.tag_reads_by_haplotype([read], [block])
        assert (tagged["HP"], tagged["PS"]) == (1, 0)

    def test_tied_support_and_missing_variants_leave_the_read_untagged(self) -> None:
        block = _phase_block(0, [100, 300], [0, 0])
        tie = {
            "read_name": "tie",
            "reference_name": "chr1",
            "reference_start": 90,
            "reference_end": 310,
            "variants": [{"position": 100, "allele": 0}, {"position": 300, "allele": 1}],
        }
        bare = {"read_name": "bare", "reference_name": "chr1", "reference_start": 90, "reference_end": 310}
        tagged = analysis_phasing.tag_reads_by_haplotype([tie, bare], [block])
        assert [(t["HP"], t["PS"]) for t in tagged] == [(0, 0), (0, 0)]
