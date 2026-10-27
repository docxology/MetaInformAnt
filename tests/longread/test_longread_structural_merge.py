"""Regression tests for merged structural variant calls.

``_cluster_sv_calls`` rebuilt each merged call from a handful of fields and dropped
``mate_chromosome``, ``mate_position`` and ``inserted_sequence``; it also merged breakends
that pointed to different partners, and listed read names and sources through ``set``, whose
order of strings changes from one Python process to the next. Insertion calls never carried
the bases they insert at all. All tests run the real detector and clustering; there are no mocks.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from pathlib import Path
from typing import Any

import metainformant
from metainformant.longread.analysis.structural import (
    StructuralVariant,
    _cluster_sv_calls,
    detect_insertions,
    detect_sv_from_long_reads,
)
from metainformant.longread.io.bam import LongReadAlignment
from metainformant.longread.utils.summary import RunSummary, compare_run_summaries


def _dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _alignment(
    name: str, cigar: str, sequence: str = "", sa: str = "", start: int = 1000, **extra: Any
) -> LongReadAlignment:
    return LongReadAlignment(
        read_name=name,
        query_sequence=sequence,
        reference_name=extra.pop("chrom", "chr1"),
        reference_start=start,
        mapping_quality=60,
        cigar_string=cigar,
        tags={"SA": sa} if sa else {},
        query_length=len(sequence) or 4000,
        **extra,
    )


def _breakend_reads(prefix: str, mate_chrom: str, mate_positions: list[int]) -> list[LongReadAlignment]:
    """Reads aligned at chr1:1000-3000 whose supplementary alignment starts at each (1-based) mate position."""
    return [
        _alignment(f"{prefix}{i}", "2000M", sa=f"{mate_chrom},{mate + 1},+,2000M,60,0;")
        for i, mate in enumerate(mate_positions)
    ]


def _breakends(svs: list[StructuralVariant]) -> list[StructuralVariant]:
    return [sv for sv in svs if sv.sv_type == "BND"]


class TestBreakendMerging:
    def test_merged_breakend_keeps_the_mate_chromosome_and_median_mate_position(self) -> None:
        reads = _breakend_reads("r", "chr2", [4990, 5000, 5010])
        (bnd,) = _breakends(detect_sv_from_long_reads(reads, min_support=2))
        assert bnd.mate_chromosome == "chr2"
        assert bnd.mate_position == 5000
        assert bnd.supporting_reads == 3
        assert (bnd.chromosome, bnd.start) == ("chr1", 3000)

    def test_a_single_unmerged_breakend_keeps_its_mate_too(self) -> None:
        reads = _breakend_reads("r", "chrX", [7000, 7000])
        (bnd,) = _breakends(detect_sv_from_long_reads(reads, min_support=2))
        assert (bnd.mate_chromosome, bnd.mate_position) == ("chrX", 7000)

    def test_breakends_to_different_chromosomes_stay_separate(self) -> None:
        reads = _breakend_reads("a", "chr2", [5000, 5000, 5000]) + _breakend_reads("b", "chr3", [8000, 8000, 8000])
        bnds = _breakends(detect_sv_from_long_reads(reads, min_support=2))
        assert sorted((b.mate_chromosome, b.mate_position, b.supporting_reads) for b in bnds) == [
            ("chr2", 5000, 3),
            ("chr3", 8000, 3),
        ]

    def test_breakends_with_distant_mates_on_one_chromosome_stay_separate_even_when_interleaved(self) -> None:
        # same locus, two partners 85 kb apart, and the reads alternate between them
        mates = [5000, 90000, 5005, 90005, 4995, 89995]
        bnds = _breakends(detect_sv_from_long_reads(_breakend_reads("r", "chr2", mates), min_support=2))
        assert sorted((b.mate_position, b.supporting_reads) for b in bnds) == [(5000, 3), (90000, 3)]
        for bnd in bnds:
            assert {int(name[1:]) % 2 for name in bnd.read_names} == {0 if bnd.mate_position == 5000 else 1}

    def test_breakends_whose_mates_chain_within_the_distance_merge(self) -> None:
        mates = [5000, 5300, 5600, 5900]  # each within 500 bp of the next, 900 bp end to end
        (bnd,) = _breakends(detect_sv_from_long_reads(_breakend_reads("r", "chr2", mates), min_support=2))
        assert bnd.supporting_reads == 4
        assert bnd.mate_position == 5450  # median of the four, truncated to an integer

    def test_calls_of_other_types_have_no_mate(self) -> None:
        reads = [_alignment(f"r{i}", "1000M300D1000M") for i in range(3)]
        (deletion,) = detect_sv_from_long_reads(reads, min_support=2)
        assert (deletion.mate_chromosome, deletion.mate_position, deletion.inserted_sequence) == ("", 0, "")


class TestInsertionSequence:
    def test_insertion_call_carries_the_inserted_bases(self) -> None:
        rng = random.Random(1)
        left, inserted, right = _dna(rng, 1000), _dna(rng, 300), _dna(rng, 1000)
        reads = [_alignment(f"r{i}", "1000M300I1000M", left + inserted + right) for i in range(3)]
        (insertion,) = detect_insertions(reads, min_support=2)
        assert insertion.inserted_sequence == inserted
        assert insertion.size == 300
        assert insertion.start == 2000

    def test_offset_counts_soft_clips_but_not_hard_clips(self) -> None:
        rng = random.Random(2)
        sequence = _dna(rng, 60 + 1000 + 200 + 500)  # 40 hard-clipped bases are not in the stored sequence
        reads = [_alignment(f"r{i}", "40H60S1000M200I500M", sequence) for i in range(2)]
        (insertion,) = detect_insertions(reads, min_support=2)
        assert insertion.inserted_sequence == sequence[1060:1260]

    def test_sequence_is_taken_in_the_stored_orientation_of_a_reverse_read(self) -> None:
        rng = random.Random(3)
        sequence = _dna(rng, 800 + 120 + 800)
        reads = [_alignment(f"r{i}", "800M120I800M", sequence, is_reverse=True) for i in range(2)]
        (insertion,) = detect_insertions(reads, min_support=2)
        assert insertion.inserted_sequence == sequence[800:920]

    def test_two_insertions_in_one_read_get_their_own_bases(self) -> None:
        rng = random.Random(4)
        sequence = _dna(rng, 500 + 100 + 700 + 150 + 500)
        reads = [_alignment(f"r{i}", "500M100I700M150I500M", sequence) for i in range(2)]
        insertions = sorted(detect_insertions(reads, min_support=2), key=lambda sv: sv.start)
        assert [sv.inserted_sequence for sv in insertions] == [sequence[500:600], sequence[1300:1450]]

    def test_missing_or_short_sequence_leaves_it_empty_rather_than_truncated(self) -> None:
        for sequence in ("", "ACGT" * 100):
            reads = [_alignment(f"r{i}", "1000M300I1000M", sequence) for i in range(3)]
            (insertion,) = detect_insertions(reads, min_support=2)
            assert insertion.inserted_sequence == ""
            assert insertion.size == 300

    def test_merged_insertion_keeps_the_sequence_whose_length_is_closest_to_the_merged_size(self) -> None:
        rng = random.Random(5)
        flank = _dna(rng, 1000)
        inserted = {length: _dna(rng, length) for length in (298, 300, 303)}
        reads = [
            _alignment(f"r{length}", f"1000M{length}I1000M", flank + inserted[length] + flank) for length in inserted
        ]
        (insertion,) = detect_insertions(reads, min_support=2)
        assert insertion.size == 300  # median of 298, 300, 303
        assert insertion.inserted_sequence == inserted[300]
        assert insertion.supporting_reads == 3

    def test_cluster_keeps_fields_of_hand_built_calls(self) -> None:
        calls = [
            StructuralVariant("INS", "chr1", 100, 101, 200, 1, ["a"], 200.0, inserted_sequence="A" * 200),
            StructuralVariant("INS", "chr1", 110, 111, 202, 1, ["b"], 202.0, inserted_sequence="C" * 202),
            StructuralVariant("BND", "chr1", 900, 901, 0, 1, ["c"], 30.0, mate_chromosome="chr5", mate_position=70),
            StructuralVariant("BND", "chr1", 905, 906, 0, 1, ["d"], 31.0, mate_chromosome="chr5", mate_position=90),
        ]
        merged = {sv.sv_type: sv for sv in _cluster_sv_calls(calls, max_distance=500, min_support=2)}
        assert merged["INS"].inserted_sequence in ("A" * 200, "C" * 202)
        assert (merged["BND"].mate_chromosome, merged["BND"].mate_position) == ("chr5", 80)


class TestDeterministicOutput:
    def test_read_names_keep_the_order_in_which_reads_were_seen(self) -> None:
        names = [f"read_{i:03d}_{'xyz'[i % 3]}" for i in range(40)]
        reads = [_alignment(name, "1000M300D1000M") for name in names]
        (deletion,) = detect_sv_from_long_reads(reads, min_support=2)
        assert deletion.read_names == names

    def test_sources_are_listed_in_sorted_order(self) -> None:
        sources = ["split_read", "zulu", "cigar", "omega", "alpha", "mid", "beta", "delta"]
        calls = [
            StructuralVariant("DEL", "chr1", 100 + i, 400 + i, 300, 1, [f"r{i}"], 1.0, info={"source": source})
            for i, source in enumerate(sources)
        ]
        (merged,) = _cluster_sv_calls(calls)
        assert merged.info["sources"] == sorted(sources)

    def test_compare_run_summaries_lists_platforms_in_sorted_order(self) -> None:
        platforms = ["pacbio", "ont", "ultima", "element", "nanopore", "bgi"]
        summaries = [RunSummary(run_id=f"r{i}", sample_name=f"s{i}", platform=p) for i, p in enumerate(platforms)]
        assert compare_run_summaries(summaries)["platforms"] == sorted(platforms)

    def test_output_does_not_depend_on_the_hash_seed_of_the_process(self) -> None:
        # str hashes (and so the iteration order of a set of strings) differ per process unless PYTHONHASHSEED is set
        snippet = (
            "import json\n"
            "from metainformant.longread.analysis.structural import detect_sv_from_long_reads\n"
            "from metainformant.longread.utils.summary import RunSummary, compare_run_summaries\n"
            "reads = [{'read_name': f'read_{i:03d}', 'reference_name': 'chr1', 'reference_start': 1000,\n"
            "          'cigar_string': '1000M300D1000M', 'mapping_quality': 60, 'tags': {}, 'query_length': 2300}\n"
            "         for i in range(40)]\n"
            "svs = detect_sv_from_long_reads(reads, min_support=2)\n"
            "names = ['pacbio', 'ont', 'ultima', 'element', 'nanopore', 'bgi']\n"
            "platforms = compare_run_summaries([RunSummary(platform=p) for p in names])['platforms']\n"
            "print(json.dumps({'svs': [[sv.read_names, sv.info['sources']] for sv in svs], 'platforms': platforms}))\n"
        )
        # the children must import the same source tree as this process (the suite's conftest otherwise
        # prepends the repository's src to PYTHONPATH)
        source_root = str(Path(metainformant.__file__).resolve().parent.parent)
        outputs = []
        for seed in ("1", "2"):
            result = subprocess.run(
                [sys.executable, "-c", snippet],
                env={**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": source_root},
                capture_output=True,
                text=True,
                check=True,
                timeout=120,
            )
            outputs.append(json.loads(result.stdout))
        assert outputs[0] == outputs[1]
        assert outputs[0]["svs"][0][0] == [f"read_{i:03d}" for i in range(40)]
        assert outputs[0]["platforms"] == ["bgi", "element", "nanopore", "ont", "pacbio", "ultima"]


def test_distant_deletions_are_still_separate_calls() -> None:
    reads = [_alignment(f"a{i}", "1000M300D1000M", start=1000) for i in range(2)] + [
        _alignment(f"b{i}", "1000M300D1000M", start=1000 + 5000) for i in range(2)
    ]
    svs = [sv for sv in detect_sv_from_long_reads(reads, min_support=2) if sv.sv_type == "DEL"]
    assert len(svs) == 2 and all(sv.supporting_reads == 2 for sv in svs)
