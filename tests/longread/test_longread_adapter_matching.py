"""Regression tests for adapter detection, trimming and chimeric read splitting.

The default adapter table used to hold two 8 bp entries that match random
sequence, ``detect_adapters`` stepped through reads in strides of half the adapter
length, and matching had no indel model. Adapter-free reads lost about 200 bp to
``trim_adapters`` and were split by ``split_chimeric_reads``, while real adapters
at unaligned offsets were missed.

These tests pin the repaired behaviour on seeded random reads and on reads with
real adapters spliced in at random offsets, and check the false-positive bound
that justifies the thresholds against an exhaustive count and against random
sequence. Everything runs on real sequences; there are no mocks.
"""

from __future__ import annotations

import itertools
import logging
import math
import random

import pytest

from metainformant.core.sequence import reverse_complement
from metainformant.longread.quality import filtering
from metainformant.longread.quality.filtering import (
    ALL_KNOWN_ADAPTERS,
    ONT_ADAPTERS,
    PACBIO_ADAPTERS,
    detect_adapters,
    split_chimeric_reads,
    trim_adapters,
)

TOP = ONT_ADAPTERS["ONT_LSK109_top"]  # 28 bp
SMRTBELL = PACBIO_ADAPTERS["PacBio_SMRTbell"]  # 45 bp
SHORT_ADAPTER = "AGATCGGAAGAGC"  # 13 bp, below the 15 bp minimum at the default rate


def _dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _other_base(base: str) -> str:
    return "ACGT"["ACGT".index(base) - 3]  # a deterministic base different from ``base``


def _substitute(sequence: str, index: int) -> str:
    return sequence[:index] + _other_base(sequence[index]) + sequence[index + 1 :]


def _delete(sequence: str, index: int) -> str:
    return sequence[:index] + sequence[index + 1 :]


def _insert(sequence: str, index: int) -> str:
    # the inserted base differs from both neighbours so the edit is not absorbed by a homopolymer
    neighbours = {sequence[index - 1] if index > 0 else "", sequence[index] if index < len(sequence) else ""}
    base = next(b for b in "ACGT" if b not in neighbours)
    return sequence[:index] + base + sequence[index:]


def _random_edits(sequence: str, edits: int, rng: random.Random) -> str:
    """Apply ``edits`` random substitutions, insertions and deletions."""
    chars = list(sequence)
    for _ in range(edits):
        kind = rng.choice("sid")
        if kind == "s":
            i = rng.randrange(len(chars))
            chars[i] = rng.choice([b for b in "ACGT" if b != chars[i]])
        elif kind == "d":
            del chars[rng.randrange(len(chars))]
        else:
            chars.insert(rng.randrange(len(chars) + 1), rng.choice("ACGT"))
    return "".join(chars)


def _levenshtein(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, start=1):
        current = [i]
        for j, char_b in enumerate(b, start=1):
            current.append(min(previous[j - 1] + (char_a != char_b), previous[j] + 1, current[j - 1] + 1))
        previous = current
    return previous[-1]


class TestDefaultAdapterTable:
    def test_no_default_adapter_is_shorter_than_the_minimum_length(self) -> None:
        assert min(len(sequence) for sequence in ALL_KNOWN_ADAPTERS.values()) >= 15
        assert filtering.minimum_adapter_length() == 15

    def test_eight_bp_motifs_are_not_default_adapters(self) -> None:
        # CAGCACCT and AGGTGCTG (reverse complements of each other) matched random sequence
        for motif in ("CAGCACCT", "AGGTGCTG"):
            assert motif not in ALL_KNOWN_ADAPTERS.values()
        assert "ONT_LSK110_top" not in ALL_KNOWN_ADAPTERS
        assert "ONT_LSK110_bottom" not in ALL_KNOWN_ADAPTERS

    def test_edit_budgets_of_the_default_table(self) -> None:
        budgets = {name: filtering.adapter_error_budget(len(seq), 0.75) for name, seq in ALL_KNOWN_ADAPTERS.items()}
        assert budgets == {
            "ONT_LSK109_top": 3,
            "ONT_LSK109_bottom": 2,
            "ONT_RNA_adapter": 2,
            "ONT_rapid_top": 6,
            "ONT_rapid_barcode": 4,
            "PacBio_SMRTbell": 5,
            "PacBio_barcoded": 3,
        }
        # an 8 bp adapter cannot be searched within the default false-positive rate, not even exactly
        assert filtering.adapter_error_budget(8, 0.75) == -1
        # identity 1.0 allows no edits at all
        assert filtering.adapter_error_budget(28, 1.0) == 0


class TestAdapterFreeReads:
    """Random reads carry no adapter, so nothing may be trimmed, reported or split."""

    @staticmethod
    def _random_reads(count: int = 100, length: int = 3000, seed: int = 1234) -> list[dict[str, str]]:
        rng = random.Random(seed)
        return [{"read_id": f"r{i}", "sequence": _dna(rng, length)} for i in range(count)]

    def test_trim_adapters_leaves_random_reads_untouched(self) -> None:
        reads = self._random_reads()
        trimmed = trim_adapters(reads)
        assert [t.sequence for t in trimmed] == [r["sequence"] for r in reads]
        assert all(t.metadata["trimmed_start"] == 0 and t.metadata["trimmed_end"] == 0 for t in trimmed)

    def test_detect_adapters_finds_nothing_in_random_reads(self) -> None:
        reads = self._random_reads()
        assert [detect_adapters(r["sequence"]) for r in reads] == [[] for _ in reads]

    def test_split_chimeric_reads_keeps_random_reads_whole(self) -> None:
        reads = self._random_reads()
        result = split_chimeric_reads(reads)
        assert [r.read_id for r in result] == [r["read_id"] for r in reads]
        assert [r.sequence for r in result] == [r["sequence"] for r in reads]
        assert not any(r.metadata["chimeric"] for r in result)

    @pytest.mark.parametrize(
        "sequence",
        ["A" * 3000, "T" * 3000, "AC" * 1500, "GGCT" * 750, "AATGTACTTCGT" * 250],
        ids=["polyA", "polyT", "dinucleotide", "tetranucleotide", "adapter-fragment-repeat"],
    )
    def test_low_complexity_reads_are_not_trimmed_or_split(self, sequence: str) -> None:
        assert [t.sequence for t in trim_adapters([{"read_id": "x", "sequence": sequence}])] == [sequence]
        assert detect_adapters(sequence) == []

    def test_false_positives_stay_within_the_stated_bound_on_a_larger_sample(self) -> None:
        reads = self._random_reads(count=200, length=5000, seed=99)
        assert sum(len(detect_adapters(r["sequence"])) for r in reads) == 0


class TestAdaptersAtAnyOffset:
    @pytest.mark.parametrize("offset", range(0, 60))
    def test_adapter_is_found_at_every_offset(self, offset: int) -> None:
        # the old half-adapter stride (14 bp) only saw offsets that happened to line up
        rng = random.Random(offset)
        read = _dna(rng, offset) + TOP + _dna(rng, 100)
        matches = detect_adapters(read, {"top": TOP})
        assert len(matches) == 1
        assert (matches[0].position, matches[0].end_position) == (offset, offset + len(TOP))
        assert matches[0].adapter_name == "top"
        assert matches[0].identity == 1.0
        assert matches[0].edit_distance == 0

    @pytest.mark.parametrize("name", sorted(ALL_KNOWN_ADAPTERS))
    def test_both_orientations_at_random_offsets(self, name: str) -> None:
        adapter = ALL_KNOWN_ADAPTERS[name]
        rng = random.Random(f"orientation-{name}")
        for _ in range(25):
            left, right = _dna(rng, rng.randint(0, 900)), _dna(rng, rng.randint(0, 900))
            for oriented, label in ((adapter, name), (reverse_complement(adapter), name + "_rc")):
                matches = detect_adapters(left + oriented + right, {name: adapter})
                assert [(m.adapter_name, m.position, m.end_position) for m in matches] == [
                    (label, len(left), len(left) + len(oriented))
                ]

    def test_internal_adapter_is_labelled_internal(self) -> None:
        rng = random.Random(5)
        read = _dna(rng, 1500) + TOP + _dna(rng, 1500)
        (match,) = detect_adapters(read, {"top": TOP})
        assert match.location == "internal"
        assert detect_adapters(TOP + _dna(rng, 1500), {"top": TOP})[0].location == "start"
        assert detect_adapters(_dna(rng, 1500) + TOP, {"top": TOP})[0].location == "end"

    def test_two_adapters_in_one_read_are_both_reported(self) -> None:
        rng = random.Random(8)
        first, middle = _dna(rng, 700), _dna(rng, 900)
        read = first + TOP + middle + reverse_complement(TOP) + _dna(rng, 400)
        matches = detect_adapters(read, {"top": TOP})
        assert [(m.adapter_name, m.position) for m in matches] == [
            ("top", len(first)),
            ("top_rc", len(first) + len(TOP) + len(middle)),
        ]

    def test_chimeric_read_is_split_at_adapters_with_unaligned_offsets(self) -> None:
        rng = random.Random(31)
        reads, expected = [], {}
        for index in range(30):
            first, second = _dna(rng, rng.randint(600, 2500)), _dna(rng, rng.randint(600, 2500))
            reads.append({"read_id": f"c{index}", "sequence": first + TOP + second})
            expected[f"c{index}"] = (first, second)
        result = split_chimeric_reads(reads, {"top": TOP}, min_fragment_length=500)
        assert len(result) == 60
        for index in range(30):
            fragments = [r for r in result if r.metadata.get("parent_read") == f"c{index}"]
            assert [f.sequence for f in fragments] == list(expected[f"c{index}"])
            assert all(f.metadata["chimeric"] for f in fragments)

    def test_default_table_splits_a_read_at_an_internal_adapter(self) -> None:
        rng = random.Random(77)
        first, second = _dna(rng, 1800), _dna(rng, 2100)
        result = split_chimeric_reads([{"read_id": "c", "sequence": first + SMRTBELL + second}])
        assert [r.sequence for r in result] == [first, second]
        assert [r.read_id for r in result] == ["c_fragment_0", "c_fragment_1"]

    def test_trim_adapters_removes_adapters_at_unaligned_offsets_from_both_ends(self) -> None:
        rng = random.Random(14)
        for _ in range(30):
            insert = _dna(rng, 2000)
            head, tail = _dna(rng, rng.randint(0, 60)), _dna(rng, rng.randint(0, 60))
            read = head + TOP + insert + reverse_complement(TOP) + tail
            (trimmed,) = trim_adapters([{"read_id": "r", "sequence": read}], {"top": TOP})
            assert trimmed.sequence == insert
            assert trimmed.metadata["trimmed_start"] == len(head) + len(TOP)
            assert trimmed.metadata["trimmed_end"] == len(reverse_complement(TOP)) + len(tail)

    def test_default_table_trims_a_smrtbell_adapter_at_the_start(self) -> None:
        rng = random.Random(3)
        insert = _dna(rng, 2500)
        (trimmed,) = trim_adapters([{"read_id": "r", "sequence": _dna(rng, 17) + SMRTBELL + insert}])
        assert trimmed.sequence == insert


class TestErrorModel:
    @pytest.mark.parametrize(
        ("variant", "edits"),
        [
            (_substitute(TOP, 0), 1),
            (_substitute(TOP, 13), 1),
            (_substitute(TOP, 27), 1),
            (_delete(TOP, 10), 1),
            (_insert(TOP, 20), 1),
            (_insert(_substitute(TOP, 3), 18), 2),
            (_delete(_substitute(_insert(TOP, 6), 15), 22), 3),
        ],
        ids=["sub-first", "sub-middle", "sub-last", "deletion", "insertion", "sub+ins", "ins+sub+del"],
    )
    def test_adapters_within_the_edit_budget_are_found(self, variant: str, edits: int) -> None:
        rng = random.Random(edits)
        left, right = _dna(rng, 333), _dna(rng, 333)
        matches = detect_adapters(left + variant + right, {"top": TOP})
        assert len(matches) == 1
        match = matches[0]
        assert match.edit_distance == edits == _levenshtein(TOP, variant)
        assert match.identity == pytest.approx(1 - edits / len(TOP))
        assert abs(match.position - len(left)) <= edits
        assert abs(match.end_position - (len(left) + len(variant))) <= edits

    def test_budget_edits_at_random_positions_are_always_found(self) -> None:
        rng = random.Random(2026)
        budget = filtering.adapter_error_budget(len(TOP), 0.75)
        assert budget == 3
        for _ in range(300):
            variant = _random_edits(TOP, rng.randint(0, budget), rng)
            left, right = _dna(rng, rng.randint(0, 500)), _dna(rng, rng.randint(0, 500))
            matches = detect_adapters(left + variant + right, {"top": TOP})
            assert len(matches) == 1
            assert matches[0].position < len(left) + len(variant) and matches[0].end_position > len(left)
            assert matches[0].edit_distance <= budget

    def test_every_single_edit_variant_is_found(self) -> None:
        # exhaustive: every substitution, deletion and insertion of a 16 bp adapter, including both ends
        adapter = "ACGTTGCAAGCTTGAC"
        rate = 1e-7
        assert filtering.adapter_error_budget(len(adapter), 0.75, rate) == 1
        variants = {adapter}
        for i in range(len(adapter)):
            variants.add(_substitute(adapter, i))
            variants.add(_delete(adapter, i))
        for i in range(len(adapter) + 1):
            variants.update(adapter[:i] + base + adapter[i:] for base in "ACGT")
        assert len(variants) >= 80  # identical variants (edits inside a homopolymer) collapse
        rng = random.Random(16)
        for variant in sorted(variants):
            left, right = _dna(rng, 150), _dna(rng, 150)
            matches = detect_adapters(left + variant + right, {"a": adapter}, max_false_positive_rate=rate)
            assert any(
                m.position <= len(left) + len(variant) and m.end_position >= len(left) and m.edit_distance <= 1
                for m in matches
            ), variant

    def test_placements_beyond_the_edit_budget_are_not_reported(self) -> None:
        # eight substitutions leave 71 % identity, which the old 0.70 threshold accepted
        variant = TOP
        for i in [0, 4, 8, 12, 16, 20, 24, 27]:
            variant = _substitute(variant, i)
        assert _levenshtein(TOP, variant) == 8
        rng = random.Random(4)
        read = _dna(rng, 14 * 20) + variant + _dna(rng, 400)
        assert detect_adapters(read, {"top": TOP}) == []
        (trimmed,) = trim_adapters([{"read_id": "r", "sequence": variant + _dna(rng, 800)}], {"top": TOP})
        assert trimmed.sequence.startswith(variant)

    def test_min_identity_caps_the_tolerated_edits(self) -> None:
        variant = _substitute(TOP, 9)
        rng = random.Random(6)
        read = _dna(rng, 300) + variant + _dna(rng, 300)
        assert len(detect_adapters(read, {"top": TOP}, min_identity=0.70)) == 1
        assert detect_adapters(read, {"top": TOP}, min_identity=1.0) == []

    def test_match_reports_the_false_positive_bound_of_its_budget(self) -> None:
        rng = random.Random(10)
        read = _dna(rng, 300) + _substitute(TOP, 5) + _dna(rng, 300)
        (match,) = detect_adapters(read, {"top": TOP})
        assert match.false_positive_bound == filtering.adapter_false_positive_bound(len(TOP), 3)
        assert match.false_positive_bound <= filtering.DEFAULT_MAX_FALSE_POSITIVE_RATE


class TestFalsePositiveBound:
    @pytest.mark.parametrize("length", [1, 5, 8, 15, 28])
    def test_exact_match_bound_is_four_to_the_minus_length(self, length: int) -> None:
        assert filtering.adapter_false_positive_bound(length, 0) == pytest.approx(4.0**-length, rel=1e-12)

    @pytest.mark.parametrize(
        ("adapter", "errors", "slack"),
        [("ACGTAC", 0, 1.0), ("ACGTAC", 1, 1.1), ("ACGTAC", 2, 1.5), ("AAAAAA", 1, 2.0)],
    )
    def test_bound_dominates_an_exhaustive_count_and_is_not_vacuous(
        self, adapter: str, errors: int, slack: float
    ) -> None:
        # exact probability mass of every string within ``errors`` edits, summed over the spans a placement can have
        exact = 0.0
        for span in range(max(1, len(adapter) - errors), len(adapter) + errors + 1):
            within = sum(
                1 for chars in itertools.product("ACGT", repeat=span) if _levenshtein(adapter, "".join(chars)) <= errors
            )
            exact += within / 4**span
        bound = filtering.adapter_false_positive_bound(len(adapter), errors)
        assert exact <= bound <= slack * exact

    def test_bound_grows_with_the_edit_budget_and_shrinks_with_length(self) -> None:
        bounds = [filtering.adapter_false_positive_bound(22, k) for k in range(5)]
        assert bounds == sorted(bounds)
        assert [filtering.adapter_false_positive_bound(n, 2) for n in (20, 24, 28)] == sorted(
            [filtering.adapter_false_positive_bound(n, 2) for n in (20, 24, 28)], reverse=True
        )

    def test_minimum_adapter_length_follows_the_rate(self) -> None:
        assert filtering.minimum_adapter_length(1e-9) == 15
        assert filtering.minimum_adapter_length(1e-6) == 10
        shortest = filtering.minimum_adapter_length(1e-12)
        assert 4.0**-shortest <= 1e-12 < 4.0 ** -(shortest - 1)

    @pytest.mark.parametrize(
        ("adapter", "rate", "budget"),
        [("ACGTTGCAAG", 1e-6, 0), ("ACGTTGCAAGCTTG", 1e-6, 1)],
        ids=["10bp-exact", "14bp-one-edit"],
    )
    def test_spurious_hits_on_random_sequence_stay_within_the_bound(
        self, adapter: str, rate: float, budget: int
    ) -> None:
        assert filtering.adapter_error_budget(len(adapter), 0.75, rate) == budget
        bound = filtering.adapter_false_positive_bound(len(adapter), budget)
        assert bound <= rate
        rng = random.Random(f"fp-{adapter}")
        offsets = 0
        hits = 0
        for _ in range(400):
            read = _dna(rng, 5000)
            offsets += len(read)
            hits += len(detect_adapters(read, {"a": adapter}, max_false_positive_rate=rate))
        expected_at_most = 2 * offsets * bound  # both orientations
        assert hits <= expected_at_most + 6 * math.sqrt(expected_at_most) + 3

    def test_short_adapter_is_searched_exactly_and_warns(self, caplog: pytest.LogCaptureFixture) -> None:
        rng = random.Random(13)
        exact = _dna(rng, 400) + SHORT_ADAPTER + _dna(rng, 400)
        one_substitution = _dna(rng, 400) + _substitute(SHORT_ADAPTER, 6) + _dna(rng, 400)
        with caplog.at_level(logging.WARNING, logger=filtering.logger.name):
            (match,) = detect_adapters(exact, {"truseq": SHORT_ADAPTER})
            assert detect_adapters(one_substitution, {"truseq": SHORT_ADAPTER}) == []
        assert (match.position, match.edit_distance) == (400, 0)
        assert match.false_positive_bound == pytest.approx(4.0**-13)
        messages = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        assert messages and all("truseq" in m and "13 bp" in m and "15 bp" in m for m in messages)

    def test_eight_bp_adapter_given_explicitly_is_exact_only_and_warns(self, caplog: pytest.LogCaptureFixture) -> None:
        motif = "CAGCACCT"
        rng = random.Random(12)
        read = _dna(rng, 300) + _substitute(motif, 3) + _dna(rng, 300)
        with caplog.at_level(logging.WARNING, logger=filtering.logger.name):
            assert detect_adapters(read, {"flank": motif}, min_identity=0.75) == []
        assert any("flank" in r.getMessage() and "8 bp" in r.getMessage() for r in caplog.records)

    @pytest.mark.parametrize("bad_identity", [-0.1, 1.1, float("nan")])
    def test_invalid_min_identity_raises(self, bad_identity: float) -> None:
        with pytest.raises(ValueError, match="min_identity"):
            detect_adapters("ACGT" * 50, min_identity=bad_identity)
        with pytest.raises(ValueError, match="min_identity"):
            trim_adapters([{"read_id": "r", "sequence": "ACGT" * 50}], min_identity=bad_identity)
        with pytest.raises(ValueError, match="min_identity"):
            split_chimeric_reads([{"read_id": "r", "sequence": "ACGT" * 50}], min_adapter_identity=bad_identity)

    @pytest.mark.parametrize("bad_rate", [0.0, 1.0, -1e-9, float("nan")])
    def test_invalid_false_positive_rate_raises(self, bad_rate: float) -> None:
        with pytest.raises(ValueError, match="max_false_positive_rate"):
            detect_adapters("ACGT" * 50, max_false_positive_rate=bad_rate)
        with pytest.raises(ValueError, match="max_false_positive_rate"):
            filtering.minimum_adapter_length(bad_rate)
        with pytest.raises(ValueError, match="max_false_positive_rate"):
            filtering.adapter_error_budget(28, 0.75, bad_rate)

    def test_invalid_bound_arguments_raise(self) -> None:
        with pytest.raises(ValueError, match="adapter_length"):
            filtering.adapter_false_positive_bound(0, 0)
        with pytest.raises(ValueError, match="max_errors"):
            filtering.adapter_false_positive_bound(10, -1)
        with pytest.raises(ValueError, match="adapter_length"):
            filtering.adapter_error_budget(0)
