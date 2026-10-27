"""Regression tests for the repeat and palindrome scans behind ``analyze_sequence_information``.

``_find_repeats`` compared every substring with the whole sequence and rescanned its growing result list for
duplicates, and ``_find_palindromes`` tested every substring, so the cost grew with the fourth power (repeats) and
the third power (palindromes) of the input length: a 500 bp sequence took minutes while ``methods=["entropy"]``
took under a millisecond. Both are now a pruned k-mer hash index and a centre expansion with the same results.

The expected records below were worked out by hand from the definitions in the docstrings. The brute-force
references transcribe those definitions directly, so every comparison checks the fast scans against the original
semantics instead of against themselves. All tests use real sequences and real implementations
(real-implementation policy).
"""

from __future__ import annotations

import random
import time
from typing import Any

import pytest

from metainformant.information.metrics.analysis.analysis import (
    _find_palindromes,
    _find_repeats,
    _infer_sequence_type,
    _is_palindrome,
    analyze_sequence_information,
)

Record = dict[str, Any]


def _reference_repeats(sequence: str, min_length: int, max_length: int | None = None) -> list[Record]:
    """Brute-force definition: every (length, start) whose substring occurs twice without overlapping."""
    longest = len(sequence) // 2 if max_length is None else min(max_length, len(sequence) // 2)
    records: list[Record] = []
    for length in range(min_length, longest + 1):
        for start in range(len(sequence) - length + 1):
            pattern = sequence[start : start + length]
            count = sequence.count(pattern)  # str.count counts non-overlapping occurrences
            if count > 1:
                records.append({"pattern": pattern, "length": length, "count": count, "start": start})
    return records


def _reference_palindromes(sequence: str, min_length: int, max_length: int | None = None) -> list[Record]:
    """Brute-force definition: every substring of at least ``min_length`` characters that is a palindrome."""
    records: list[Record] = []
    for start in range(len(sequence) - min_length + 1):
        longest = len(sequence) - start if max_length is None else min(max_length, len(sequence) - start)
        for length in range(min_length, longest + 1):
            substring = sequence[start : start + length]
            if _is_palindrome(substring):
                records.append({"sequence": substring, "start": start, "length": length})
    return records


# --- hand-worked expectations -------------------------------------------------------------------------------------


def test_repeats_report_every_start_of_a_repeated_substring() -> None:
    assert _find_repeats("ACGTACGT") == [
        {"pattern": "ACG", "length": 3, "count": 2, "start": 0},
        {"pattern": "CGT", "length": 3, "count": 2, "start": 1},
        {"pattern": "ACG", "length": 3, "count": 2, "start": 4},
        {"pattern": "CGT", "length": 3, "count": 2, "start": 5},
        {"pattern": "ACGT", "length": 4, "count": 2, "start": 0},
        {"pattern": "ACGT", "length": 4, "count": 2, "start": 4},
    ]


def test_repeat_occurrences_must_not_overlap() -> None:
    # "ABA" starts at 0 and 2 of "ABABAB", but the two copies share a character, so str.count sees one occurrence
    assert _find_repeats("ABABAB") == []
    # in "aaaaaa" the copies at 0 and 3 are disjoint, and every start of "aaa" is reported with the same count
    assert _find_repeats("aaaaaa") == [{"pattern": "aaa", "length": 3, "count": 2, "start": s} for s in range(4)]


def test_repeats_are_case_sensitive() -> None:
    records = _find_repeats("ACGacgACGacg")

    # "ACG" occurs at 0 and 6 and "acg" at 3 and 9; a case-insensitive scan would count four of each
    assert {record["count"] for record in records if record["pattern"] in ("ACG", "acg")} == {2}
    assert {record["pattern"] for record in records if record["start"] == 0} == {"ACG", "ACGa", "ACGac", "ACGacg"}


def test_repeats_honour_min_length() -> None:
    assert _find_repeats("ACGTACGT", min_length=5) == []
    assert {record["length"] for record in _find_repeats("ACGTACGT", min_length=4)} == {4}
    assert {record["length"] for record in _find_repeats("ACGTACGT", min_length=1)} == {1, 2, 3, 4}


def test_palindromes_are_reverse_complements_for_nucleotides() -> None:
    # the EcoRI site GAATTC and the AATT inside it; ACCA reads the same both ways but is not its reverse complement
    assert _find_palindromes("GAATTC") == [
        {"sequence": "GAATTC", "start": 0, "length": 6},
        {"sequence": "AATT", "start": 1, "length": 4},
    ]
    assert _find_palindromes("ACCA") == []
    assert _find_palindromes("aaTT") == [{"sequence": "aaTT", "start": 0, "length": 4}]


def test_palindromes_nest_and_are_ordered_by_start_then_length() -> None:
    assert _find_palindromes("GGATCCGGATCC") == [
        {"sequence": "GGATCC", "start": 0, "length": 6},
        {"sequence": "GGATCCGGATCC", "start": 0, "length": 12},
        {"sequence": "GATC", "start": 1, "length": 4},
        {"sequence": "GATCCGGATC", "start": 1, "length": 10},
        {"sequence": "ATCCGGAT", "start": 2, "length": 8},
        {"sequence": "TCCGGA", "start": 3, "length": 6},
        {"sequence": "CCGG", "start": 4, "length": 4},
        {"sequence": "GGATCC", "start": 6, "length": 6},
        {"sequence": "GATC", "start": 7, "length": 4},
    ]


def test_palindromes_are_plain_for_proteins() -> None:
    assert _find_palindromes("MKKMLLQQL") == [
        {"sequence": "MKKM", "start": 0, "length": 4},
        {"sequence": "LQQL", "start": 5, "length": 4},
    ]


def test_palindrome_kind_is_chosen_per_substring() -> None:
    """The kind depends on each substring's own characters, so a palindrome need not contain smaller ones.

    "XAAX" is a plain palindrome, but the "AA" inside it reads as DNA and is not its own reverse complement.
    """
    assert _find_palindromes("XAAX", min_length=2) == [{"sequence": "XAAX", "start": 0, "length": 4}]


def test_palindrome_predicate_is_unchanged() -> None:
    assert _is_palindrome("ACGT") and _is_palindrome("GAATTC")
    assert not _is_palindrome("ACCA") and not _is_palindrome("AA")
    assert _is_palindrome("MKKM") and _is_palindrome("mkkm") and _is_palindrome("XAAX")


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("ACGT", "DNA"),
        ("acgt", "DNA"),
        ("ACGU", "RNA"),
        ("ACDEFGHIK", "protein"),
        ("XZ", "unknown"),
        ("ACGTX", "unknown"),  # four nucleotide letters and four amino-acid letters tie
        ("ACGUX", "RNA"),  # U is a nucleotide letter but not an amino-acid letter
        ("", "DNA"),
    ],
)
def test_sequence_type_inference_is_unchanged(sequence: str, expected: str) -> None:
    # the palindrome scan reuses this rule to decide, per substring, which kind of palindrome applies
    assert _infer_sequence_type(sequence) == expected


@pytest.mark.parametrize("sequence", ["", "A", "ACG", "ACGT", "ACGTA"])
def test_inputs_too_short_to_repeat_give_no_repeats(sequence: str) -> None:
    # two non-overlapping copies of a 3-character pattern need 6 characters
    assert _find_repeats(sequence) == []


@pytest.mark.parametrize("sequence", ["", "A", "ACG"])
def test_inputs_shorter_than_min_length_give_no_palindromes(sequence: str) -> None:
    assert _find_palindromes(sequence) == []


@pytest.mark.parametrize("min_length", [0, -1])
def test_min_length_below_one_is_rejected(min_length: int) -> None:
    # a zero-length pattern is meaningless; the old scans silently reported empty-string "repeats" and "palindromes"
    with pytest.raises(ValueError, match="min_length"):
        _find_repeats("ACGTACGT", min_length=min_length)
    with pytest.raises(ValueError, match="min_length"):
        _find_palindromes("ACGTACGT", min_length=min_length)


# --- identical results to the brute-force definitions ---------------------------------------------------------------

ALPHABETS = [
    "ACGT",  # DNA
    "ACGU",  # RNA
    "ACGTU",  # mixed DNA/RNA letters
    "ACGTN",  # ambiguity code
    "acgt",  # lower case
    "ACgt",  # mixed case
    "AC",  # tiny alphabet with many overlapping occurrences
    "A",  # homopolymer
    "ACGTX",  # a letter that is neither nucleotide nor amino acid
    "AUX",  # U without any amino-acid-only letter, next to a foreign letter
    "ACDEFGHIKLMNPQRSTVWY",  # protein
    "ACGT-",  # alignment gap
    "ACGTß",  # "ß".upper() == "SS", so the upper-cased text is longer than the sequence
]


def _sequences(alphabet: str, seed: int) -> list[str]:
    rng = random.Random(seed)
    sequences = ["".join(rng.choice(alphabet) for _ in range(length)) for length in range(0, 41, 2)]
    unit = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 5)))
    sequences.append((unit * 20)[:38])  # tandem repeat
    half = "".join(rng.choice(alphabet) for _ in range(15))
    sequences.append(half + half[::-1])  # mirror image
    sequences.append(half + rng.choice(alphabet) + half[::-1])  # mirror image around a centre character
    return sequences


@pytest.mark.parametrize("alphabet", ALPHABETS)
@pytest.mark.parametrize("min_length", [1, 3, 4])
def test_scans_match_the_brute_force_definitions(alphabet: str, min_length: int) -> None:
    for sequence in _sequences(alphabet, seed=len(alphabet) * 31 + min_length):
        assert _find_repeats(sequence, min_length) == _reference_repeats(sequence, min_length), sequence
        assert _find_palindromes(sequence, min_length) == _reference_palindromes(sequence, min_length), sequence


def test_default_arguments_match_the_brute_force_definitions() -> None:
    rng = random.Random(11)
    for _ in range(40):
        sequence = "".join(rng.choice("ACGT") for _ in range(rng.randint(10, 36)))
        assert _find_repeats(sequence) == _reference_repeats(sequence, 3)
        assert _find_palindromes(sequence) == _reference_palindromes(sequence, 4)


@pytest.mark.parametrize("unit", ["A", "AC", "ACGT", "ATCG", "GAATTC", "N", "MK", "UX"])
def test_low_complexity_sequences_match_the_brute_force_definitions(unit: str) -> None:
    sequence = (unit * 300)[:240]  # homopolymers, microsatellites and tandem repeats: huge numbers of overlapping hits

    assert _find_palindromes(sequence) == _reference_palindromes(sequence, 4)
    assert _find_repeats(sequence) == _reference_repeats(sequence, 3)


# --- run time -------------------------------------------------------------------------------------------------------


def _planted_dna(length: int, seed: int) -> str:
    """Random DNA with one 40-mer planted at four evenly spaced, non-overlapping places."""
    rng = random.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(length)]
    motif = [rng.choice("ACGT") for _ in range(40)]
    for start in range(100, length - 100, length // 4):
        bases[start : start + 40] = motif
    return "".join(bases)


@pytest.mark.timeout(120)  # the old scans needed hours for this input; the guard turns that into a failure
def test_two_kilobase_dna_is_analysed_in_seconds() -> None:
    sequence = _planted_dna(2000, seed=7)

    started = time.perf_counter()
    result = analyze_sequence_information(sequence)
    elapsed = time.perf_counter() - started

    assert elapsed < 5.0, f"analysing a 2 kb sequence took {elapsed:.1f} s"
    repeats = result["methods"]["patterns"]["repeats"]
    palindromes = result["methods"]["patterns"]["palindromes"]

    # the planted 40-mer occurs four times
    assert {"pattern": sequence[100:140], "length": 40, "count": 4, "start": 100} in repeats
    # records are ordered by length and then start, and every record is true of the sequence
    keys = [(record["length"], record["start"]) for record in repeats]
    assert keys == sorted(set(keys))
    for record in repeats:
        assert sequence[record["start"] : record["start"] + record["length"]] == record["pattern"]
        assert sequence.count(record["pattern"]) == record["count"] > 1
    assert all(_is_palindrome(record["sequence"]) for record in palindromes)
    # nothing is missing among the short repeats and palindromes, which the brute-force definitions can still check
    assert [record for record in repeats if record["length"] <= 8] == _reference_repeats(sequence, 3, max_length=8)
    assert [r for r in palindromes if r["length"] <= 12] == _reference_palindromes(sequence, 4, max_length=12)


@pytest.mark.timeout(120)
def test_two_kilobase_protein_is_analysed_in_seconds() -> None:
    rng = random.Random(3)
    sequence = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(2000))

    started = time.perf_counter()
    result = analyze_sequence_information(sequence)
    elapsed = time.perf_counter() - started

    assert elapsed < 5.0, f"analysing a 2 kb protein took {elapsed:.1f} s"
    repeats = result["methods"]["patterns"]["repeats"]
    palindromes = result["methods"]["patterns"]["palindromes"]
    assert [record for record in repeats if record["length"] <= 5] == _reference_repeats(sequence, 3, max_length=5)
    assert [r for r in palindromes if r["length"] <= 12] == _reference_palindromes(sequence, 4, max_length=12)


@pytest.mark.timeout(120)
def test_highly_repetitive_input_is_not_quadratic_in_its_own_output() -> None:
    # 600 bp of a 4-bp unit gives over 100,000 repeat records; the old duplicate check rescanned all of the records
    # found so far for each new one
    sequence = "ACGT" * 150

    started = time.perf_counter()
    repeats = _find_repeats(sequence)
    elapsed = time.perf_counter() - started

    assert elapsed < 10.0, f"scanning a 600 bp tandem repeat took {elapsed:.1f} s"
    assert len(repeats) > 100_000
    by_key = {(record["length"], record["start"]): record for record in repeats}
    assert len(by_key) == len(repeats)
    # spot-check both directions: a (length, start) is reported exactly when its substring repeats
    rng = random.Random(5)
    for _ in range(3000):
        length = rng.randint(3, 300)
        start = rng.randint(0, len(sequence) - length)
        pattern = sequence[start : start + length]
        count = sequence.count(pattern)
        expected = {"pattern": pattern, "length": length, "count": count, "start": start} if count > 1 else None
        assert by_key.get((length, start)) == expected


@pytest.mark.timeout(120)
@pytest.mark.parametrize("unit", ["A", "AC"])
def test_a_two_kilobase_microsatellite_is_scanned_for_palindromes_in_seconds(unit: str) -> None:
    # almost every substring is a plain palindrome here, which the old scan tested one by one (cubic time), but
    # nucleotide-like text only counts when it is its own reverse complement, and a repeat of A or AC never is
    sequence = (unit * 2000)[:2000]

    started = time.perf_counter()
    palindromes = _find_palindromes(sequence)
    elapsed = time.perf_counter() - started

    assert elapsed < 5.0, f"scanning a 2 kb microsatellite took {elapsed:.1f} s"
    assert palindromes == []
