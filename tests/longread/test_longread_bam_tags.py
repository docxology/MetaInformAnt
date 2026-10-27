"""Regression tests for MM/ML base-modification tag decoding in ``longread.io.bam``.

pysam returns the B-type ML tag as ``array('B')``. ``_parse_methylation_from_tags`` only
recognised lists, tuples and bytes, so ``int(array)`` raised ``TypeError`` for every read of a
real modified-base BAM, and ``tags.get("ML") or ...`` raised ``ValueError`` for numpy arrays.
``extract_methylation_tags`` keyed every entry ``X+code``, so ``C-m`` replaced ``C+m``.

pysam is not installed in this environment, so the tests build the tag dictionaries a pysam
read would give (``array('B')`` values) and run the real parsing on them.
"""

from __future__ import annotations

from array import array
from typing import Any

import numpy as np
import pytest

from metainformant.longread.io.bam import LongReadAlignment, _parse_methylation_from_tags, extract_methylation_tags
from metainformant.longread.workflow.orchestrator_core import _extract_methylation_features

SEQUENCE = "ACGTACGTACGTACGT"  # C at 1, 5, 9, 13; G at 2, 6, 10, 14; A at 0, 4, 8, 12


def _alignment(tags: dict[str, Any], parse_at_read_time: bool = False) -> LongReadAlignment:
    return LongReadAlignment(
        read_name="r1",
        query_sequence=SEQUENCE,
        reference_start=5000,
        reference_end=5016,
        tags=tags,
        methylation_tags=_parse_methylation_from_tags(tags) if parse_at_read_time else {},
        aligned_pairs=[(i, 5000 + i) for i in range(len(SEQUENCE))],
    )


ML_FORMS = {
    "array('B')": lambda values: array("B", values),
    "numpy uint8": lambda values: np.array(values, dtype=np.uint8),
    "numpy int64": lambda values: np.array(values, dtype=np.int64),
    "list": list,
    "tuple": tuple,
    "bytes": bytes,
    "bytearray": bytearray,
}


class TestMlTagTypes:
    @pytest.mark.parametrize("form", sorted(ML_FORMS))
    def test_every_array_type_is_decoded(self, form: str) -> None:
        tags = {"MM": "C+m,0,1;", "ML": ML_FORMS[form]([200, 100])}
        parsed = _parse_methylation_from_tags(tags)
        (modification,) = parsed["modifications"]
        assert modification["probabilities"] == [200, 100]
        assert all(type(p) is int for p in modification["probabilities"])
        assert modification["delta_positions"] == [0, 1]

    def test_values_are_shared_out_across_modifications_in_mm_order(self) -> None:
        tags = {"MM": "C+m,0,1;A+a,2;", "ML": array("B", [250, 10, 128])}
        first, second = _parse_methylation_from_tags(tags)["modifications"]
        assert (first["base"], first["probabilities"]) == ("C", [250, 10])
        assert (second["base"], second["probabilities"]) == ("A", [128])

    def test_legacy_lowercase_tag_names_with_an_array(self) -> None:
        parsed = _parse_methylation_from_tags({"Mm": "C+m,0;", "Ml": array("B", [7])})
        assert parsed["modifications"][0]["probabilities"] == [7]

    def test_single_value_forms_still_work(self) -> None:
        assert _parse_methylation_from_tags({"MM": "C+m,0;", "ML": 200})["modifications"][0]["probabilities"] == [200]
        assert _parse_methylation_from_tags({"MM": "C+m,0;", "ML": "200"})["modifications"][0]["probabilities"] == [200]

    def test_missing_or_empty_ml_gives_no_probabilities(self) -> None:
        for ml in (None, array("B"), np.array([], dtype=np.uint8), []):
            tags: dict[str, Any] = {"MM": "C+m,0,1;"}
            if ml is not None:
                tags["ML"] = ml
            assert _parse_methylation_from_tags(tags)["modifications"][0]["probabilities"] == []

    def test_no_mm_tag_gives_nothing(self) -> None:
        assert _parse_methylation_from_tags({"ML": array("B", [1, 2])}) == {}
        assert _parse_methylation_from_tags({}) == {}


class TestExtractWithPysamStyleTags:
    def test_extract_methylation_tags_from_array_tags(self) -> None:
        tags = {"MM": "C+m,0,1;", "ML": array("B", [200, 100])}
        for parse_at_read_time in (False, True):  # raw tags, and the dict read_long_read_bam stores
            result = extract_methylation_tags(_alignment(tags, parse_at_read_time))
            assert list(result) == ["C+m"]
            assert result["C+m"]["read_positions"] == [1, 9]
            assert result["C+m"]["probabilities"] == pytest.approx([200 / 255, 100 / 255])
            assert result["C+m"]["genomic_positions"] == [
                (5001, pytest.approx(200 / 255)),
                (5009, pytest.approx(100 / 255)),
            ]

    def test_extract_methylation_tags_from_numpy_tags(self) -> None:
        tags = {"MM": "A+a,0,0;", "ML": np.array([255, 0], dtype=np.uint8)}
        result = extract_methylation_tags(_alignment(tags))
        assert result["A+a"]["read_positions"] == [0, 4]
        assert result["A+a"]["probabilities"] == pytest.approx([1.0, 0.0])

    def test_methylation_features_are_built_from_array_tags(self) -> None:
        # the workflow path that feeds call_5mc from MM/ML tags
        reads = [
            {
                "tags": {"MM": "C+m,0,1;", "ML": array("B", [200, 100])},
                "reference_name": "chr1",
                "reference_start": 5000,
            }
        ]
        features = _extract_methylation_features(reads, "5mC")
        assert [f["context_shift"] for f in features] == pytest.approx([200 / 255 * 5.0, 100 / 255 * 5.0])
        assert all(f["chromosome"] == "chr1" for f in features)


class TestStrandInModificationKey:
    def test_plus_and_minus_strand_entries_do_not_overwrite_each_other(self) -> None:
        tags = {"MM": "C+m,0,1;C-m,2;", "ML": array("B", [200, 100, 50])}
        result = extract_methylation_tags(_alignment(tags))
        assert sorted(result) == ["C+m", "C-m"]
        assert result["C+m"]["read_positions"] == [1, 9]
        assert result["C+m"]["probabilities"] == pytest.approx([200 / 255, 100 / 255])
        assert result["C-m"]["probabilities"] == pytest.approx([50 / 255])

    def test_entries_differing_only_in_strand_keep_their_own_probabilities_in_either_order(self) -> None:
        for mm, ml, plus, minus in (
            ("C-m,0;C+m,1;", [10, 20], 20, 10),
            ("C+m,1;C-m,0;", [20, 10], 20, 10),
        ):
            result = extract_methylation_tags(_alignment({"MM": mm, "ML": array("B", ml)}))
            assert result["C+m"]["probabilities"] == pytest.approx([plus / 255])
            assert result["C-m"]["probabilities"] == pytest.approx([minus / 255])

    def test_different_codes_and_bases_are_still_separate_entries(self) -> None:
        tags = {"MM": "C+m,0;C+h,0;A+a,0;", "ML": array("B", [1, 2, 3])}
        assert sorted(extract_methylation_tags(_alignment(tags))) == ["A+a", "C+h", "C+m"]

    def test_read_with_only_plus_strand_entries_keeps_the_documented_key(self) -> None:
        # the key format for the common CpG case is unchanged
        assert list(extract_methylation_tags(_alignment({"MM": "C+m,0;", "ML": [9]}))) == ["C+m"]
