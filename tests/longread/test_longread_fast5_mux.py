"""Regression tests for the mux of reads parsed from FAST5 files.

``_parse_fast5_read_group`` set ``mux = 0`` for every read of a multi-read FAST5 and never read
the ``start_mux`` attribute that ONT stores on each read's ``Raw`` group. The fixtures here are real
HDF5 files written with h5py in the layouts ONT uses (multi-read ``read_<id>`` groups and the
legacy single-read ``Raw/Reads/Read_<n>`` layout) and parsed with the real reader.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.longread.io.fast5 import get_read_metadata, read_fast5

h5py = pytest.importorskip("h5py")  # FAST5 files are HDF5; not every CI lane installs it

READS = [
    {"read_id": "0f3c1a2e-0000-4000-8000-000000000001", "mux": 1, "channel": 17, "start": 1000, "duration": 5},
    {"read_id": "0f3c1a2e-0000-4000-8000-000000000002", "mux": 2, "channel": 17, "start": 2000, "duration": 6},
    {"read_id": "0f3c1a2e-0000-4000-8000-000000000003", "mux": 3, "channel": 311, "start": 3000, "duration": 7},
    {"read_id": "0f3c1a2e-0000-4000-8000-000000000004", "mux": 4, "channel": 512, "start": 4000, "duration": 8},
]


def _write_multi_read(path: Path, reads: list[dict[str, Any]], with_mux: bool = True) -> Path:
    """Write a multi-read FAST5 the way ONT does: attributes on Raw, channel_id and tracking_id."""
    with h5py.File(path, "w") as f5:
        for read in reads:
            group = f5.create_group(f"read_{read['read_id']}")
            raw = group.create_group("Raw")
            raw.attrs["read_id"] = read["read_id"]
            raw.attrs["start_time"] = np.uint64(read["start"])
            raw.attrs["duration"] = np.uint32(read["duration"])
            if with_mux:
                raw.attrs["start_mux"] = np.uint8(read["mux"])
            raw.create_dataset("Signal", data=np.arange(read["duration"], dtype=np.int16))
            channel = group.create_group("channel_id")
            channel.attrs["channel_number"] = str(read["channel"])  # ONT stores the channel as a string
            channel.attrs["sampling_rate"] = np.float64(4000.0)
            channel.attrs["digitisation"] = np.float64(8192.0)
            channel.attrs["offset"] = np.float64(-5.0)
            channel.attrs["range"] = np.float64(1467.0)
            group.create_group("tracking_id").attrs["run_id"] = "run-1234"
    return path


def _write_single_read(path: Path, read_id: str, mux: int) -> Path:
    """Write a legacy single-read FAST5 (Raw/Reads/Read_<n>)."""
    with h5py.File(path, "w") as f5:
        raw_read = f5.create_group("Raw/Reads/Read_42")
        raw_read.attrs["read_id"] = read_id
        raw_read.attrs["start_mux"] = np.uint8(mux)
        raw_read.attrs["start_time"] = np.uint64(900)
        raw_read.attrs["duration"] = np.uint32(4)
        raw_read.create_dataset("Signal", data=np.arange(4, dtype=np.int16))
        channel = f5.create_group("UniqueGlobalKey/channel_id")
        channel.attrs["channel_number"] = "88"
        channel.attrs["sampling_rate"] = np.float64(4000.0)
        f5.create_group("UniqueGlobalKey/tracking_id").attrs["run_id"] = "run-9"
    return path


class TestMultiReadMux:
    def test_each_read_gets_the_mux_stored_in_the_file(self, tmp_path: Path) -> None:
        reads = read_fast5(_write_multi_read(tmp_path / "multi.fast5", READS))
        by_id = {r.read_id: r for r in reads}
        assert {rid: r.mux for rid, r in by_id.items()} == {spec["read_id"]: spec["mux"] for spec in READS}
        assert {r.mux for r in reads} == {1, 2, 3, 4}

    def test_other_fields_are_unaffected(self, tmp_path: Path) -> None:
        reads = {r.read_id: r for r in read_fast5(_write_multi_read(tmp_path / "multi.fast5", READS))}
        for spec in READS:
            read = reads[spec["read_id"]]
            assert (read.channel_id, read.start_time, read.duration, read.run_id) == (
                spec["channel"],
                spec["start"],
                spec["duration"],
                "run-1234",
            )
            assert read.metadata["source_format"] == "fast5_multi"
            assert read.signal is not None and len(read.signal) == spec["duration"]

    def test_metadata_reports_the_mux(self, tmp_path: Path) -> None:
        reads = read_fast5(_write_multi_read(tmp_path / "multi.fast5", READS))
        assert sorted(get_read_metadata(r)["mux"] for r in reads) == [1, 2, 3, 4]

    def test_mux_is_a_plain_int(self, tmp_path: Path) -> None:
        (read,) = read_fast5(_write_multi_read(tmp_path / "one.fast5", READS[2:3]))
        assert type(read.mux) is int and read.mux == 3

    def test_a_file_without_start_mux_gives_zero_and_still_parses(self, tmp_path: Path) -> None:
        reads = read_fast5(_write_multi_read(tmp_path / "nomux.fast5", READS, with_mux=False))
        assert [r.mux for r in reads] == [0, 0, 0, 0]
        assert sorted(r.channel_id for r in reads) == [17, 17, 311, 512]

    @pytest.mark.parametrize("mux", [1, 2, 3, 4])
    def test_every_mux_value_round_trips(self, tmp_path: Path, mux: int) -> None:
        spec = dict(READS[0], mux=mux)
        (read,) = read_fast5(_write_multi_read(tmp_path / f"mux{mux}.fast5", [spec]))
        assert read.mux == mux


class TestSingleReadMux:
    def test_legacy_single_read_layout_keeps_its_mux(self, tmp_path: Path) -> None:
        # this layout already worked; the multi-read reader now agrees with it
        (read,) = read_fast5(_write_single_read(tmp_path / "single.fast5", "read-abc", mux=3))
        assert (read.read_id, read.mux, read.channel_id, read.metadata["source_format"]) == (
            "read-abc",
            3,
            88,
            "fast5_single",
        )
