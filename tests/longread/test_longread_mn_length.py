"""MN consistency follows SAMtags: changed SEQ lengths invalidate MM/ML offsets."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.longread.io.bam import LongReadAlignment, _parse_methylation_from_tags, extract_methylation_tags


def alignment(*, mn: Any = 4, mode: str = "raw", reverse: bool = False) -> LongReadAlignment:
    tags = {"MM": "C+m,0;", "ML": [255], "MN": mn}
    record = LongReadAlignment(
        read_name="mn-read",
        query_sequence="GGGG" if reverse else "CCCC",
        query_length=100,
        is_reverse=reverse,
        tags=tags,
        aligned_pairs=[(i, i) for i in range(4)],
    )
    if mode != "raw":
        record.methylation_tags = _parse_methylation_from_tags(tags)
        if mode == "parsed_only":
            record.tags = {}
    return record


@pytest.mark.parametrize("mode", ["raw", "parsed", "parsed_only"])
@pytest.mark.parametrize("reverse", [False, True])
def test_matching_mn_uses_stored_seq_length_and_keeps_the_call(mode: str, reverse: bool) -> None:
    record = alignment(mode=mode, reverse=reverse)
    decoded = extract_methylation_tags(record)["C+m"]
    assert decoded["read_positions"] == ([3] if reverse else [0])
    assert decoded["probabilities"] == [1.0]


@pytest.mark.parametrize("mode", ["raw", "parsed", "parsed_only"])
@pytest.mark.parametrize("mn", [0, 3, 5])
@pytest.mark.parametrize("reverse", [False, True])
def test_stale_mn_refuses_calls_before_either_decoding_path(mode: str, mn: int, reverse: bool) -> None:
    with pytest.raises(ValueError, match=r"MN.*SEQ length"):
        extract_methylation_tags(alignment(mn=mn, mode=mode, reverse=reverse))


@pytest.mark.parametrize("mode", ["raw", "parsed", "parsed_only"])
@pytest.mark.parametrize("mn", [True, 4.0, "4", None, -1])
def test_mn_must_be_a_nonnegative_integer(mode: str, mn: Any) -> None:
    with pytest.raises(ValueError, match="MN"):
        extract_methylation_tags(alignment(mn=mn, mode=mode))


@pytest.mark.parametrize("mode", ["raw", "parsed", "parsed_only"])
def test_absent_mn_remains_compatible(mode: str) -> None:
    record = alignment()
    del record.tags["MN"]
    if mode != "raw":
        record.methylation_tags = _parse_methylation_from_tags(record.tags)
        if mode == "parsed_only":
            record.tags = {}
    assert extract_methylation_tags(record)["C+m"]["probabilities"] == [1.0]


def test_integer_scalar_mn_is_accepted() -> None:
    assert extract_methylation_tags(alignment(mn=np.int64(4)))["C+m"]["read_positions"] == [0]


def test_conflicting_cached_mn_cannot_be_hidden_by_current_raw_tag() -> None:
    record = alignment()
    record.methylation_tags = _parse_methylation_from_tags({**record.tags, "MN": 5})
    with pytest.raises(ValueError, match=r"MN.*SEQ length"):
        extract_methylation_tags(record)


@pytest.mark.external_tool
def test_native_samtools_accepts_matching_and_absent_mn_but_rejects_stale_mn(tmp_path: Path) -> None:
    samtools = shutil.which("samtools")
    if samtools is None:
        pytest.skip("samtools is required for the native MN control")
    for mn in (None, 4, 5):
        sam = tmp_path / f"mn-{mn}.sam"
        suffix = "" if mn is None else f"\tMN:i:{mn}"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
            "r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\tMM:Z:C+m,0;\tML:B:C,255" + suffix + "\n"
        )
        result = subprocess.run(
            [samtools, "mpileup", "--output-mods", str(sam)], text=True, capture_output=True, timeout=20
        )
        (tmp_path / f"mn-{mn}.stdout").write_text(result.stdout)
        (tmp_path / f"mn-{mn}.stderr").write_text(result.stderr)
        (tmp_path / f"mn-{mn}.json").write_text(
            json.dumps({"argv": result.args, "exit_code": result.returncode}, indent=2) + "\n"
        )
        if mn == 5:
            assert result.returncode != 0
            assert "MM/MN data length is incompatible with SEQ length" in result.stderr
        else:
            assert result.returncode == 0, result.stderr
            assert "[+m255]" in result.stdout
