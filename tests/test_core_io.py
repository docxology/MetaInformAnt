from __future__ import annotations

import json
from pathlib import Path

from metainformant.core import io as core_io


def test_ensure_directory_creates(tmp_path: Path) -> None:
    p = tmp_path / "nested" / "dir"
    assert not p.exists()
    made = core_io.ensure_directory(p)
    assert made.exists() and made.is_dir()


def test_json_roundtrip(tmp_path: Path) -> None:
    data = {"a": 1, "b": [1, 2, 3]}
    path = tmp_path / "data.json"
    core_io.dump_json(data, path)
    loaded = core_io.load_json(path)
    assert loaded == data


def test_json_gz_roundtrip(tmp_path: Path) -> None:
    data = {"x": "y"}
    path = tmp_path / "data.json.gz"
    core_io.dump_json(data, path)
    loaded = core_io.load_json(path)
    assert loaded == data


def test_jsonl_roundtrip(tmp_path: Path) -> None:
    rows = [{"i": i} for i in range(5)]
    path = tmp_path / "rows.jsonl"
    core_io.write_jsonl(rows, path)
    read_back = list(core_io.read_jsonl(path))
    assert read_back == rows


def test_tsv_roundtrip(tmp_path: Path) -> None:
    rows = [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}]
    path = tmp_path / "rows.tsv"
    core_io.write_delimited(rows, path, delimiter="\t")
    read_back = list(core_io.read_delimited(path, delimiter="\t"))
    assert read_back == rows


def test_open_text_auto_handles_gz(tmp_path: Path) -> None:
    txt_gz = tmp_path / "hello.txt.gz"
    with core_io.open_text_auto(txt_gz, mode="wt") as fh:
        fh.write("hello world\n")
    with core_io.open_text_auto(txt_gz, mode="rt") as fh:
        content = fh.read()
    assert content.strip() == "hello world"
