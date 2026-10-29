"""MM/ML base-modification decoding against the SAM optional-tags specification (hts-specs ``SAMtags``).

Every expectation below is worked out by hand from the specification on a 12 bp read and was also
checked against the htslib reference implementation (htslib 1.24, ``samtools mpileup --output-mods`` on
the same records written as SAM): the modified bases land on the same SEQ indices with the same ML values.
The reads::

    sequenced read   A C G C G T A C C G T A      index 0..11
    C at 1 3 7 8     G at 2 4 9     A at 0 6 11     T at 5 10

    stored for FLAG 0x10 (the reverse complement)   T A C G G T A C G C G T
    stored index i is sequenced index 11 - i; sequenced C 1 3 7 8 sit at stored G 10 8 4 3

What the old decoder got wrong, and what these tests pin:

* it counted the letter in SEQ as stored, so a reverse-aligned read was decoded on the wrong strand of
  the sequence (the skip counts refer to the read as sequenced, FLAG 0x10 stores its reverse complement);
* ``C+mh`` was one code ``mh`` with one ML value per position, so the h probabilities were read as the m
  probabilities of later positions, shifting every following entry;
* the ``.`` flag made the entry unparseable: it was skipped without consuming ML, which shifted the
  probabilities of every entry after it;
* the flag field ``implicit`` was True for ``?``, the opposite of htslib (True when skipped bases are
  unmodified), and a missing ML tag turned into probability 0.0 instead of unknown;
* skip counts running past the end of the read were dropped without a word, and positions that were not
  aligned were reported as if the read coordinate were a genomic one.

On the ``-`` strand: the specification marks the modification as sitting on the opposite strand while the
counted letter is still the one written in the entry (htslib counts it in the read as sequenced, and
``samtools mpileup -M`` shows ``G[-m50]`` for ``G-m``). ``G-m`` therefore counts Gs and flags the 5mC opposite
them, and ``C-m`` counts Cs; the base that carries the modification is ``modified_base``.

All tests decode real hand-built tags and sequences; there are no mocks.
"""

from __future__ import annotations

import logging
import random
from array import array
from typing import Any

import pytest

from metainformant.longread.io import bam
from metainformant.longread.io.bam import LongReadAlignment, _parse_methylation_from_tags, extract_methylation_tags
from metainformant.longread.workflow.orchestrator_core import _extract_methylation_features

FORWARD = "ACGCGTACCGTA"  # the read as sequenced
REVERSE = "TACGGTACGCGT"  # its reverse complement: what SEQ holds for an alignment with FLAG 0x10
START = 1000


def _aln(
    sequence: str, mm: str, ml: Any, *, reverse: bool = False, pairs: bool = True, **tags: Any
) -> LongReadAlignment:
    all_tags: dict[str, Any] = {"MM": mm, **tags}
    if ml is not None:
        all_tags["ML"] = array("B", ml) if isinstance(ml, list) else ml
    return LongReadAlignment(
        read_name="r1",
        query_sequence=sequence,
        reference_name="chr1",
        reference_start=START,
        is_reverse=reverse,
        tags=all_tags,
        aligned_pairs=[(i, START + i) for i in range(len(sequence))] if pairs else [],
    )


def _decode(alignment: LongReadAlignment) -> dict[str, Any]:
    return extract_methylation_tags(alignment)


def _sites(entry: dict[str, Any]) -> list[tuple[int, int]]:
    """``(read position, ML value)`` of an entry's listed sites."""
    return [
        (position, round(probability * 255))
        for position, probability in zip(entry["read_positions"], entry["probabilities"])
    ]


class TestReadOrientation:
    """The skip counts refer to the read as sequenced; FLAG 0x10 stores its reverse complement."""

    def test_forward_read_counts_cs_in_seq_as_stored(self) -> None:
        # C+m,1,0: skip the C at 1, modify the C at 3 (200); skip none, modify the next C, 7 (100)
        entry = _decode(_aln(FORWARD, "C+m,1,0;", [200, 100]))["C+m"]
        assert _sites(entry) == [(3, 200), (7, 100)]
        assert [FORWARD[p] for p in entry["read_positions"]] == ["C", "C"]
        assert entry["genomic_positions"] == [(START + 3, 200 / 255), (START + 7, 100 / 255)]

    def test_reverse_aligned_read_is_reverse_complemented_before_counting(self) -> None:
        # same tags, same read as sequenced, stored reversed: sequenced 3 (200) -> stored 8, 7 (100) -> stored 4
        entry = _decode(_aln(REVERSE, "C+m,1,0;", [200, 100], reverse=True))["C+m"]
        assert _sites(entry) == [(4, 100), (8, 200)]  # ascending along the stored SEQ
        # the modified C of the read is the G of the stored sequence
        assert [REVERSE[p] for p in entry["read_positions"]] == ["G", "G"]
        assert entry["genomic_positions"] == [(START + 4, 100 / 255), (START + 8, 200 / 255)]

    def test_a_reverse_read_is_not_decoded_as_if_it_were_forward(self) -> None:
        reverse = _decode(_aln(REVERSE, "C+m,1,0;", [200, 100], reverse=True))["C+m"]
        forward_reading = _decode(_aln(REVERSE, "C+m,1,0;", [200, 100], reverse=False))["C+m"]
        assert reverse["read_positions"] == [4, 8]
        assert forward_reading["read_positions"] == [7, 9]  # the Cs of the stored sequence: the old, wrong answer
        assert [REVERSE[p] for p in forward_reading["read_positions"]] == ["C", "C"]

    def test_the_two_orientations_describe_the_same_read(self) -> None:
        forward = _decode(_aln(FORWARD, "C+m,0,1,0;", [11, 22, 33]))["C+m"]
        reverse = _decode(_aln(REVERSE, "C+m,0,1,0;", [11, 22, 33], reverse=True))["C+m"]
        length = len(FORWARD)
        assert sorted((length - 1 - p, v) for p, v in _sites(reverse)) == _sites(forward) == [(1, 11), (7, 22), (8, 33)]

    def test_lowercase_sequences_and_rna_t_are_counted(self) -> None:
        assert _sites(_decode(_aln(FORWARD.lower(), "C+m,1;", [9]))["C+m"]) == [(3, 9)]
        assert _sites(_decode(_aln(REVERSE.lower(), "C+m,1;", [9], reverse=True))["C+m"]) == [(8, 9)]

    def test_prematerialised_methylation_tags_decode_like_raw_tags(self) -> None:
        raw = _aln(REVERSE, "C+mh,1,0;", [1, 2, 3, 4], reverse=True)
        stored = _aln(REVERSE, "C+mh,1,0;", [1, 2, 3, 4], reverse=True)
        stored.methylation_tags = _parse_methylation_from_tags(stored.tags)  # what read_long_read_bam keeps
        assert _decode(raw) == _decode(stored)


class TestStrand:
    """``-`` marks the modification as sitting on the opposite strand; the counted letter is the one written."""

    def test_g_minus_counts_gs_and_flags_the_c_opposite(self) -> None:
        # G at 2 4 9: G-m,0,1 -> the G at 2 (50), skip the G at 4, the G at 9 (60)
        entry = _decode(_aln(FORWARD, "G-m,0,1;", [50, 60]))["G-m"]
        assert _sites(entry) == [(2, 50), (9, 60)]
        assert [FORWARD[p] for p in entry["read_positions"]] == ["G", "G"]
        assert (entry["base"], entry["strand"], entry["modified_base"]) == ("G", "-", "C")

    def test_c_minus_counts_cs_like_c_plus_does(self) -> None:
        # C at 1 3 7 8: C-m,0,1 -> the C at 1 (50), skip the C at 3, the C at 7 (60)
        minus = _decode(_aln(FORWARD, "C-m,0,1;", [50, 60]))["C-m"]
        plus = _decode(_aln(FORWARD, "C+m,0,1;", [50, 60]))["C+m"]
        assert _sites(minus) == _sites(plus) == [(1, 50), (7, 60)]
        assert (minus["strand"], minus["modified_base"]) == ("-", "G")
        assert (plus["strand"], plus["modified_base"]) == ("+", "C")

    def test_minus_strand_on_a_reverse_aligned_read(self) -> None:
        # C-m,0,1 on the read as sequenced: sequenced 1 (50) and 7 (60) -> stored 10 and 4
        entry = _decode(_aln(REVERSE, "C-m,0,1;", [50, 60], reverse=True))["C-m"]
        assert _sites(entry) == [(4, 60), (10, 50)]

    def test_plus_and_minus_entries_of_one_base_do_not_collide(self) -> None:
        result = _decode(_aln(FORWARD, "C+m,0;C-m,1;", [10, 20]))
        assert sorted(result) == ["C+m", "C-m"]
        assert _sites(result["C+m"]) == [(1, 10)]
        assert _sites(result["C-m"]) == [(3, 20)]

    @pytest.mark.parametrize(
        ("mm", "modified_base"),
        [
            ("C+m,0;", "C"),
            ("G+m,0;", "G"),
            ("C-m,0;", "G"),
            ("G-m,0;", "C"),
            ("A-a,0;", "T"),
            ("T-g,0;", "A"),
            ("N-n,0;", "N"),
        ],
    )
    def test_modified_base_is_the_base_for_plus_and_its_complement_for_minus(self, mm: str, modified_base: str) -> None:
        (modification,) = _parse_methylation_from_tags({"MM": mm})["modifications"]
        assert modification["modified_base"] == modified_base


class TestMultiCodeEntries:
    """``C+mh,1,1`` lists two codes at every position; ML interleaves their values in code order."""

    def test_two_codes_get_their_own_probabilities_at_shared_positions(self) -> None:
        # C at 1 3 7 8: skip one -> the C at 3: m=10 h=20; skip one -> skip the C at 7, the C at 8: m=30 h=40
        result = _decode(_aln(FORWARD, "C+mh,1,1;", [10, 20, 30, 40]))
        assert sorted(result) == ["C+h", "C+m"]
        assert _sites(result["C+m"]) == [(3, 10), (8, 30)]
        assert _sites(result["C+h"]) == [(3, 20), (8, 40)]

    def test_parse_gives_one_dictionary_per_code_with_the_deinterleaved_values(self) -> None:
        m, h = _parse_methylation_from_tags({"MM": "C+mh,1,1;", "ML": array("B", [10, 20, 30, 40])})["modifications"]
        assert (m["modification"], m["probabilities"], m["delta_positions"]) == ("m", [10, 30], [1, 1])
        assert (h["modification"], h["probabilities"], h["delta_positions"]) == ("h", [20, 40], [1, 1])

    def test_three_codes(self) -> None:
        values = [1, 2, 3, 4, 5, 6]  # positions 1 and 3 (C at 1 and the C at 7), codes m h f
        result = _decode(_aln(FORWARD, "C+mhf,0,1;", values))
        assert [_sites(result[f"C+{code}"]) for code in "mhf"] == [[(1, 1), (7, 4)], [(1, 2), (7, 5)], [(1, 3), (7, 6)]]

    def test_multi_code_on_a_reverse_read(self) -> None:
        # sequenced C 3 (m=10 h=20) and C 8 (m=30 h=40) -> stored 8 and 3, ascending: 3 then 8
        result = _decode(_aln(REVERSE, "C+mh,1,1;", [10, 20, 30, 40], reverse=True))
        assert _sites(result["C+m"]) == [(3, 30), (8, 10)]
        assert _sites(result["C+h"]) == [(3, 40), (8, 20)]

    def test_following_entries_are_not_shifted(self) -> None:
        # the multi-code entry takes 2 values per position; the A+a entry after it gets the third: 33
        result = _decode(_aln(FORWARD, "C+mh,0;A+a,1;", [11, 22, 33]))
        assert _sites(result["C+m"]) == [(1, 11)]
        assert _sites(result["C+h"]) == [(1, 22)]
        assert _sites(result["A+a"]) == [(6, 33)]  # skip the A at 0, modify the A at 6

    def test_entry_order_decides_which_values_belong_where(self) -> None:
        result = _decode(_aln(FORWARD, "A+a,1;C+mh,0;", [33, 11, 22]))
        assert (_sites(result["A+a"]), _sites(result["C+m"]), _sites(result["C+h"])) == (
            [(6, 33)],
            [(1, 11)],
            [(1, 22)],
        )

    def test_a_chebi_number_is_a_single_code(self) -> None:
        (modification,) = _parse_methylation_from_tags({"MM": "C+27551,1;", "ML": array("B", [7])})["modifications"]
        assert (modification["modification"], modification["probabilities"]) == ("27551", [7])
        assert list(_decode(_aln(FORWARD, "C+27551,1;", [7]))) == ["C+27551"]

    def test_the_old_one_code_reading_would_have_mixed_the_probabilities(self) -> None:
        # with ML [10 20 30 40] read as one value per position the 2nd position would have m=20, not 30
        m, _ = _parse_methylation_from_tags({"MM": "C+mh,1,1;", "ML": array("B", [10, 20, 30, 40])})["modifications"]
        assert m["probabilities"] != [10, 20]


class TestModeFlags:
    """``.``: skipped bases are unmodified; ``?``: unknown; no flag: the specification's default, ``.``."""

    @pytest.mark.parametrize(("mm", "implicit"), [("C+m,1;", True), ("C+m.,1;", True), ("C+m?,1;", False)])
    def test_implicit_means_the_skipped_bases_are_unmodified(self, mm: str, implicit: bool) -> None:
        (modification,) = _parse_methylation_from_tags({"MM": mm, "ML": array("B", [200])})["modifications"]
        assert modification["implicit"] is implicit
        assert modification["delta_positions"] == [1] and modification["probabilities"] == [200]

    @pytest.mark.parametrize("mm", ["C+m,1;", "C+m.,1;", "C+m?,1;"])
    def test_the_flag_does_not_move_the_listed_site(self, mm: str) -> None:
        entry = _decode(_aln(FORWARD, mm, [200]))["C+m"]
        assert _sites(entry) == [(3, 200)]

    def test_the_skipped_bases_are_reported_with_the_mode(self) -> None:
        # C at 1 3 7 8; only the C at 3 is listed; the other three are skipped (unmodified or unknown by the flag)
        dot = _decode(_aln(FORWARD, "C+m.,1;", [200]))["C+m"]
        unknown = _decode(_aln(FORWARD, "C+m?,1;", [200]))["C+m"]
        assert dot["skipped_read_positions"] == unknown["skipped_read_positions"] == [1, 7, 8]
        assert dot["skipped_genomic_positions"] == [START + 1, START + 7, START + 8]
        assert (dot["implicit"], unknown["implicit"]) == (True, False)

    def test_skipped_bases_are_counted_after_the_last_listed_site_too(self) -> None:
        entry = _decode(_aln(FORWARD, "C+m,0;", [5]))["C+m"]
        assert entry["read_positions"] == [1] and entry["skipped_read_positions"] == [3, 7, 8]

    def test_skipped_bases_of_a_reverse_read_are_in_stored_coordinates(self) -> None:
        entry = _decode(_aln(REVERSE, "C+m.,1;", [200], reverse=True))["C+m"]
        # sequenced C 1 3 7 8 sit at stored 10 8 4 3; the listed one is sequenced 3 = stored 8
        assert entry["read_positions"] == [8] and entry["skipped_read_positions"] == [3, 4, 10]

    def test_the_dot_flag_no_longer_shifts_the_probabilities_of_later_entries(self) -> None:
        # the old pattern did not match "C+m.", skipped the entry without consuming its ML value,
        # and gave the A+a entry the 200 meant for C+m
        result = _decode(_aln(FORWARD, "C+m.,1;A+a,0;", [200, 90]))
        assert _sites(result["C+m"]) == [(3, 200)]
        assert _sites(result["A+a"]) == [(0, 90)]

    def test_n_counts_every_base(self) -> None:
        # N+m,0,2,0: base 0, skip two (1, 2), base 3, then base 4
        assert _sites(_decode(_aln(FORWARD, "N+m,0,2,0;", [53, 63, 73]))["N+m"]) == [(0, 53), (3, 63), (4, 73)]

    def test_an_entry_without_positions_is_an_empty_entry_that_consumes_no_ml(self) -> None:
        result = _decode(_aln(FORWARD, "C+m.;A+a,0;", [90]))
        assert result["C+m"]["read_positions"] == [] and result["C+m"]["skipped_read_positions"] == [1, 3, 7, 8]
        assert _sites(result["A+a"]) == [(0, 90)]


class TestUnknownAndMalformedInput:
    def test_mm_without_ml_leaves_the_probabilities_unknown_not_zero(self) -> None:
        entry = _decode(_aln(FORWARD, "C+m,1,0;", None))["C+m"]
        assert entry["read_positions"] == [3, 7]
        assert entry["probabilities"] == []
        assert entry["genomic_positions"] == [(START + 3, None), (START + 7, None)]  # not 0.0

    @pytest.mark.parametrize("ml", [[200], [200, 100, 50]], ids=["too-short", "too-long"])
    def test_an_ml_tag_of_the_wrong_length_is_not_used_and_says_so(
        self, ml: list[int], caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger=bam.logger.name):
            entry = _decode(_aln(FORWARD, "C+m,1,0;", ml))["C+m"]
        assert entry["probabilities"] == [] and entry["read_positions"] == [3, 7]
        assert any("ML tag has" in r.getMessage() and "MM tag lists 2" in r.getMessage() for r in caplog.records)

    @pytest.mark.parametrize(
        "mm", ["C+;", "X+m,1;", "C*m,1;", "C+m,a;", "C+m,-1;", "C+m27551,1;", "C+mm?.,1;", "+m,1;"]
    )
    def test_a_malformed_mm_tag_gives_nothing_and_says_so(self, mm: str, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING, logger=bam.logger.name):
            assert _parse_methylation_from_tags({"MM": mm, "ML": array("B", [1])}) == {}
        assert any("malformed MM entry" in r.getMessage() for r in caplog.records)

    def test_one_malformed_entry_voids_the_tag_instead_of_shifting_the_others(self) -> None:
        assert _parse_methylation_from_tags({"MM": "C+m,1;X+m,2;A+a,0;", "ML": array("B", [1, 2, 3])}) == {}

    def test_skips_that_run_past_the_end_of_the_read_are_reported(self, caplog: pytest.LogCaptureFixture) -> None:
        # four Cs: the third skip count (5) lands beyond the last one
        with caplog.at_level(logging.WARNING, logger=bam.logger.name):
            entry = _decode(_aln(FORWARD, "C+m,0,1,5;", [10, 20, 30]))["C+m"]
        assert _sites(entry) == [(1, 10), (7, 20)]  # the probabilities stay with their sites
        assert entry["n_unplaced"] == 1
        assert any("lists 3 positions but only 2 could be placed" in r.getMessage() for r in caplog.records)

    def test_a_read_without_sequence_places_nothing_and_says_so(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING, logger=bam.logger.name):
            entry = _decode(_aln("", "C+m,1,0;", [200, 100]))["C+m"]
        assert (entry["read_positions"], entry["probabilities"], entry["n_unplaced"]) == ([], [], 2)
        assert any("no sequence" in r.getMessage() for r in caplog.records)

    def test_positions_without_a_reference_coordinate_are_none_not_read_coordinates(self) -> None:
        # 2S6M: stored positions 0 and 1 are soft-clipped; the rest align to START..
        sequence = "CCACGTAC"  # C at 0 1 3 7
        alignment = _aln(sequence, "C+m,0,0,0,0;", [1, 2, 3, 4], pairs=False)
        alignment.aligned_pairs = bam._aligned_pairs_from_cigar([(4, 2), (0, 6)], START)
        entry = _decode(alignment)["C+m"]
        assert entry["read_positions"] == [0, 1, 3, 7]
        assert [position for position, _ in entry["genomic_positions"]] == [None, None, START + 1, START + 5]

    def test_an_unaligned_read_keeps_read_positions_and_has_no_reference_positions(self) -> None:
        entry = _decode(_aln(FORWARD, "C+m,1,0;", [200, 100], pairs=False))["C+m"]
        assert entry["read_positions"] == [3, 7]
        assert [position for position, _ in entry["genomic_positions"]] == [None, None]
        assert entry["probabilities"] == pytest.approx([200 / 255, 100 / 255])

    def test_the_draft_tag_names_are_read_when_the_final_ones_are_absent(self) -> None:
        alignment = LongReadAlignment(
            read_name="r",
            query_sequence=FORWARD,
            tags={"Mm": "C+m,1;", "Ml": array("B", [200])},
            aligned_pairs=[(i, i) for i in range(len(FORWARD))],
        )
        assert _sites(_decode(alignment)["C+m"]) == [(3, 200)]

    def test_mm_as_bytes(self) -> None:
        (modification,) = _parse_methylation_from_tags({"MM": b"C+m,1;", "ML": array("B", [200])})["modifications"]
        assert modification["delta_positions"] == [1]

    def test_rna_u_counts_the_ts_of_seq(self) -> None:
        # T at 5 and 10: U+17802,1 skips the T at 5 and modifies the T at 10
        assert _sites(_decode(_aln(FORWARD, "U+17802,1;", [77]))["U+17802"]) == [(10, 77)]


class TestAlignedPairsFromCigar:
    def test_pairs_follow_the_cigar_like_pysam_get_aligned_pairs(self) -> None:
        # 2S 3M 1I 2M 2D 1M: soft clips and insertions have no reference position, deletions no query position
        pairs = bam._aligned_pairs_from_cigar([(4, 2), (0, 3), (1, 1), (0, 2), (2, 2), (0, 1)], 100)
        assert pairs == [
            (0, None),
            (1, None),
            (2, 100),
            (3, 101),
            (4, 102),
            (5, None),
            (6, 103),
            (7, 104),
            (None, 105),
            (None, 106),
            (8, 107),
        ]

    def test_hard_clips_padding_and_skipped_regions(self) -> None:
        assert bam._aligned_pairs_from_cigar([(5, 7), (0, 1), (3, 2), (6, 1), (7, 1), (8, 1)], 10) == [
            (0, 10),
            (None, 11),
            (None, 12),
            (1, 13),
            (2, 14),
        ]


class TestWorkflowMethylationFeatures:
    """``_extract_methylation_features`` decodes with the same rules and no longer invents positions or signal."""

    @staticmethod
    def _read(sequence: str, mm: str, ml: list[int], *, reverse: bool = False, cigar: str = "12M") -> dict[str, Any]:
        return {
            "read_id": "r1",
            "sequence": sequence,
            "reference_name": "chr1",
            "reference_start": START,
            "cigar_string": cigar,
            "is_reverse": reverse,
            "tags": {"MM": mm, "ML": array("B", ml)},
        }

    def test_forward_read_features_sit_on_the_modified_cs(self) -> None:
        features = _extract_methylation_features([self._read(FORWARD, "C+m,1,0;", [200, 100])], "5mC")
        assert [(f["position"], f["strand"]) for f in features] == [(START + 3, "+"), (START + 7, "+")]
        assert [f["probability"] for f in features] == pytest.approx([200 / 255, 100 / 255])
        assert all(f["chromosome"] == "chr1" and f["modification_type"] == "5mC" for f in features)

    def test_the_old_position_arithmetic_is_gone(self) -> None:
        # reference_start + cumulative skips + index put these two sites on the adjacent positions 1001 and 1002
        # whatever the read held; they are the Cs at read indices 3 and 7, 1003 and 1007
        features = _extract_methylation_features([self._read(FORWARD, "C+m,1,0;", [200, 100])], "5mC")
        assert [f["position"] for f in features] != [START + 1, START + 2]
        assert [f["position"] for f in features] == [START + 3, START + 7]

    def test_no_signal_features_are_made_up(self) -> None:
        feature, _ = _extract_methylation_features([self._read(FORWARD, "C+m,1,0;", [200, 100])], "5mC")
        assert not {"context_shift", "dwell_time", "current_std"} & set(feature)
        assert feature["context"] == ""  # the tags do not say CpG

    def test_reverse_aligned_read_flips_the_reference_strand(self) -> None:
        # the sequenced C is the stored G, and the sequenced strand is the reference minus strand
        features = _extract_methylation_features([self._read(REVERSE, "C+m,1,0;", [200, 100], reverse=True)], "5mC")
        assert [(f["position"], f["strand"]) for f in features] == [(START + 4, "-"), (START + 8, "-")]
        assert [f["probability"] for f in features] == pytest.approx([100 / 255, 200 / 255])

    def test_g_minus_is_5mc_on_the_opposite_strand(self) -> None:
        forward = _extract_methylation_features([self._read(FORWARD, "G-m,0,1;", [50, 60])], "5mC")
        assert [(f["position"], f["strand"]) for f in forward] == [(START + 2, "-"), (START + 9, "-")]
        reverse = _extract_methylation_features([self._read(REVERSE, "G-m,0,1;", [50, 60], reverse=True)], "5mC")
        # G at sequenced 2 and 9 -> stored 9 and 2; the opposite strand of a reverse read is the reference plus strand
        assert [(f["position"], f["strand"]) for f in reverse] == [(START + 2, "+"), (START + 9, "+")]

    def test_modification_type_selects_the_entries(self) -> None:
        read = self._read(FORWARD, "C+m,1;A+a,1;C+h,0;", [200, 90, 33])
        assert [f["position"] for f in _extract_methylation_features([read], "5mC")] == [START + 3]
        assert [f["position"] for f in _extract_methylation_features([read], "6mA")] == [START + 6]
        assert _extract_methylation_features([read], "other") == []

    def test_multi_code_entry_gives_5mc_features_from_the_m_values_only(self) -> None:
        features = _extract_methylation_features([self._read(FORWARD, "C+mh,1,1;", [10, 20, 30, 40])], "5mC")
        assert [f["probability"] for f in features] == pytest.approx([10 / 255, 30 / 255])

    def test_cigar_gaps_are_respected(self) -> None:
        # 4M2D8M: the read base at index 4 sits at reference START + 6
        sequence = "ACGCCGTACCGT"  # C at 1 3 4 7 8
        features = _extract_methylation_features([self._read(sequence, "C+m,2;", [150], cigar="4M2D8M")], "5mC")
        assert [f["position"] for f in features] == [START + 6]  # third C = read index 4

    def test_sites_in_a_soft_clip_have_no_reference_position_and_are_reported(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        sequence = "CCACGTACGTAA"  # C at 0 1 3 7
        read = self._read(sequence, "C+m,0,0,0,0;", [1, 2, 3, 4], cigar="2S10M")
        with caplog.at_level(logging.WARNING, logger="metainformant.longread.workflow.orchestrator_core"):
            features = _extract_methylation_features([read], "5mC")
        assert [f["position"] for f in features] == [START + 1, START + 5]  # read indices 3 and 7
        assert any("2 tagged sites left out" in r.getMessage() for r in caplog.records)

    def test_a_read_without_sequence_or_cigar_gives_no_made_up_positions(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        bare = {
            "tags": {"MM": "C+m,0,1;", "ML": array("B", [200, 100])},
            "reference_name": "chr1",
            "reference_start": 5000,
        }
        with caplog.at_level(logging.WARNING, logger="metainformant.longread.workflow.orchestrator_core"):
            assert _extract_methylation_features([bare], "5mC") == []
        assert caplog.records

    def test_mm_without_ml_gives_no_features(self) -> None:
        read = self._read(FORWARD, "C+m,1,0;", [])
        read["tags"] = {"MM": "C+m,1,0;"}
        assert _extract_methylation_features([read], "5mC") == []

    def test_precomputed_features_are_passed_through(self) -> None:
        feature = {"position": 7, "modification_type": "5mC", "context_shift": 3.0}
        other = {"position": 9, "modification_type": "6mA"}
        read = {"methylation_features": [feature, other]}
        assert _extract_methylation_features([read], "5mC") == [feature]

    def test_features_feed_call_5mc_with_the_basecaller_probability(self) -> None:
        from metainformant.longread.analysis.modified_bases import call_5mc, call_6ma

        features = _extract_methylation_features([self._read(FORWARD, "C+m,1,0;", [255, 0])], "5mC")
        calls = call_5mc(features)
        assert [(c.position, c.strand, c.probability, c.modification_type) for c in calls] == [
            (START + 3, "+", 1.0, "5mC"),
            (START + 7, "+", 0.0, "5mC"),
        ]
        six = call_6ma(_extract_methylation_features([self._read(FORWARD, "A+a,0;", [128])], "6mA"))
        assert [(c.position, c.probability, c.modification_type) for c in six] == [(START, 128 / 255, "6mA")]

    def test_call_5mc_rejects_a_probability_outside_zero_one(self) -> None:
        from metainformant.longread.analysis.modified_bases import call_5mc

        for bad in (1.5, -0.1, float("nan")):
            with pytest.raises(ValueError, match="between 0 and 1"):
                call_5mc([{"position": 1, "probability": bad}])


class TestAgainstAnIndependentDecoder:
    """Random reads, tags and orientations decoded by the package and by a plain, literal reading of the specification.

    The reference decoder below is written the long way (build the sequenced read, walk it while counting,
    map back) so a slip in the package's index arithmetic cannot hide behind the same slip here. The 2000
    random records of the htslib comparison described in the module docstring gave identical sites.
    """

    COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    @classmethod
    def _decode_literally(cls, seq: str, reverse: bool, entries: list[tuple[str, str, str, list[int]]], ml: list[int]):
        sequenced = seq
        if reverse:
            sequenced = "".join(cls.COMPLEMENT[c] for c in reversed(seq))
        expected: dict[str, list[tuple[int, int]]] = {}
        cursor = 0
        for letter, strand, codes, skips in entries:
            positions = []
            index = 0  # walk the sequenced read, counting the letter
            for skip in skips:
                remaining = skip
                while True:
                    if letter == "N" or sequenced[index] == letter:
                        if remaining == 0:
                            break
                        remaining -= 1
                    index += 1
                positions.append(index)
                index += 1
            for code_number, code in enumerate(codes):
                sites = []
                for number, position in enumerate(positions):
                    value = ml[cursor + number * len(codes) + code_number]
                    stored = len(seq) - 1 - position if reverse else position
                    sites.append((stored, value))
                expected[f"{letter}{strand}{code}"] = sorted(sites)
            cursor += len(positions) * len(codes)
        return expected

    def test_random_records_decode_like_the_literal_reading(self) -> None:
        rng = random.Random(20261028)
        compared = 0
        for _ in range(300):
            length = rng.randint(15, 70)
            sequenced = "".join(rng.choice("ACGT") for _ in range(length))
            reverse = rng.random() < 0.5
            stored = "".join(self.COMPLEMENT[c] for c in reversed(sequenced)) if reverse else sequenced
            entries, parts, ml, used = [], [], [], set()
            for _ in range(rng.randint(1, 3)):
                letter, strand = rng.choice("ACGTN"), rng.choice("+-")
                codes = rng.choice(["m", "h", "mh", "a", "mhf"])
                if (letter, strand) in used:
                    continue
                used.add((letter, strand))
                counted = [i for i, c in enumerate(sequenced) if letter == "N" or c == letter]
                if not counted:
                    continue
                chosen = sorted(rng.sample(range(len(counted)), rng.randint(0, min(len(counted), 6))))
                skips = [index - previous - 1 for previous, index in zip([-1] + chosen[:-1], chosen)]
                flag = rng.choice(["", ".", "?"])
                entries.append((letter, strand, codes, skips))
                parts.append(f"{letter}{strand}{codes}{flag}" + "".join(f",{s}" for s in skips) + ";")
                ml.extend(rng.randint(0, 255) for _ in range(len(skips) * len(codes)))
            if not entries:
                continue
            alignment = _aln(stored, "".join(parts), ml or None, reverse=reverse)
            decoded = {key: _sites(value) for key, value in _decode(alignment).items()}
            assert decoded == self._decode_literally(stored, reverse, entries, ml), "".join(parts)
            compared += 1
        assert compared > 250
