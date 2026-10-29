"""Loaders never invent a quality or a volume for a record that does not give one.

``read_transcripts`` gave every transcript without a quality value Q20 (exactly the default threshold, so it passed
the default filter and failed any stricter one), and the MERFISH parsers gave spots without a quality column a
perfect 1.0 and cells without a volume 0.0. An unknown value is NaN, and a quality filter that cannot be applied
because the file has no quality column is an error. Real files written to a temporary directory; no test doubles.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from metainformant.spatial.io.merfish import load_transcript_spots, parse_cell_metadata
from metainformant.spatial.io.xenium import read_transcripts


def _write(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> Path:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _xenium_rows(qualities: list[str | None]) -> list[dict[str, str]]:
    rows = []
    for i, quality in enumerate(qualities):
        row = {"transcript_id": f"t{i}", "feature_name": f"G{i}", "x_location": str(i), "y_location": "0"}
        if quality is not None:
            row["qv"] = quality
        rows.append(row)
    return rows


class TestXeniumTranscriptQuality:
    BASE = ["transcript_id", "feature_name", "x_location", "y_location"]

    def test_a_file_without_a_quality_column_cannot_be_quality_filtered(self, tmp_path: Path) -> None:
        path = _write(tmp_path / "transcripts.csv", self.BASE, _xenium_rows([None, None]))

        with pytest.raises(ValueError, match="no quality column"):
            read_transcripts(path)  # the default threshold is 20

    def test_without_a_quality_column_and_without_a_threshold_the_quality_is_nan(self, tmp_path: Path) -> None:
        path = _write(tmp_path / "transcripts.csv", self.BASE, _xenium_rows([None, None]))

        transcripts = read_transcripts(path, min_quality=0)

        assert [t.gene for t in transcripts] == ["G0", "G1"]
        assert all(math.isnan(t.quality_value) for t in transcripts)  # not a made-up Q20

    def test_an_unknown_quality_does_not_pass_a_threshold(self, tmp_path: Path) -> None:
        rows = _xenium_rows(["30", "", "15", "25"])
        path = _write(tmp_path / "transcripts.csv", [*self.BASE, "qv"], rows)

        transcripts = read_transcripts(path, min_quality=20.0)

        # t1 has an empty quality (it used to become Q20 and pass), t2 is below the threshold
        assert [t.gene for t in transcripts] == ["G0", "G3"]
        assert [t.quality_value for t in transcripts] == [30.0, 25.0]

    def test_an_unknown_quality_is_kept_as_nan_when_nothing_is_filtered(self, tmp_path: Path) -> None:
        path = _write(tmp_path / "transcripts.csv", [*self.BASE, "qv"], _xenium_rows(["30", "", "15"]))

        transcripts = read_transcripts(path, min_quality=0)

        assert len(transcripts) == 3
        assert transcripts[0].quality_value == 30.0 and transcripts[2].quality_value == 15.0
        assert math.isnan(transcripts[1].quality_value)

    def test_a_stricter_threshold_filters_more(self, tmp_path: Path) -> None:
        path = _write(tmp_path / "transcripts.csv", [*self.BASE, "qv"], _xenium_rows(["10", "20", "30", "40"]))

        assert len(read_transcripts(path, min_quality=20.0)) == 3
        assert len(read_transcripts(path, min_quality=30.0)) == 2


class TestMerfishUnknownValues:
    def test_spots_without_a_quality_column_have_nan_quality(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path / "spots.csv",
            ["gene", "x", "y"],
            [{"gene": "A", "x": "1", "y": "2"}, {"gene": "B", "x": "3", "y": "4"}],
        )

        spots = load_transcript_spots(path)

        assert len(spots) == 2 and all(math.isnan(spot.quality) for spot in spots)

    def test_a_spot_with_an_empty_quality_is_nan_and_the_others_keep_theirs(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path / "spots.csv",
            ["gene", "x", "y", "quality"],
            [
                {"gene": "A", "x": "1", "y": "2", "quality": "0.4"},
                {"gene": "B", "x": "3", "y": "4", "quality": ""},
            ],
        )

        spots = load_transcript_spots(path)

        assert spots[0].quality == pytest.approx(0.4) and math.isnan(spots[1].quality)

    def test_cells_without_a_volume_have_nan_volume(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path / "cells.csv",
            ["cell_id", "center_x", "center_y", "volume"],
            [
                {"cell_id": "c1", "center_x": "1", "center_y": "2", "volume": "150.5"},
                {"cell_id": "c2", "center_x": "3", "center_y": "4", "volume": ""},
            ],
        )

        cells = parse_cell_metadata(path)

        assert cells[0].volume == pytest.approx(150.5) and math.isnan(cells[1].volume)

    def test_a_file_without_a_volume_column_gives_nan_volumes(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path / "cells.csv",
            ["cell_id", "center_x", "center_y"],
            [{"cell_id": "c1", "center_x": "1", "center_y": "2"}],
        )

        assert math.isnan(parse_cell_metadata(path)[0].volume)
