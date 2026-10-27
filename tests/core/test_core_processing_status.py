"""``download_and_process_data`` must not report processing it never ran.

The ``processing`` section of a config was handled by a placeholder that stored ``{"status": "completed"}`` for every
entry without looking at it, so a workflow that did no processing at all reported all of it as completed (and
``run_config_based_workflow`` reported success). Nothing in the package implements the ``type`` values such a section
names (there is no registry to look one up in), so the entries are now recorded as ``not_executed`` with the reason,
and are listed in ``errors`` so that the workflow is not reported as successful.

Downloads use a real local ``file://`` source, so no network is involved; no test doubles (real-implementation policy).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import pytest

from metainformant.core import io
from metainformant.core.execution.workflow import (
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
)

PROCESSING = {
    "analysis": {"type": "basic_analysis", "parameters": {"param1": "value1"}},
    "normalization": {"type": "normalize", "method": "log2"},
}


def test_processing_entries_are_recorded_as_not_executed(tmp_path: Path) -> None:
    results = download_and_process_data({"processing": PROCESSING}, tmp_path / "out")

    assert set(results["processing"]) == set(PROCESSING)
    for name, config in PROCESSING.items():
        entry = results["processing"][name]
        assert entry["status"] == "not_executed"  # was "completed" without anything having run
        assert repr(config["type"]) in entry["reason"]
        assert entry["config"] == config  # the entry's configuration is kept for the caller
    assert "completed" not in {entry["status"] for entry in results["processing"].values()}


def test_every_unexecuted_entry_is_listed_as_an_error(tmp_path: Path) -> None:
    results = download_and_process_data({"processing": PROCESSING}, tmp_path / "out")

    assert len(results["errors"]) == len(PROCESSING)
    for name in PROCESSING:
        assert any(name in error and "not executed" in error for error in results["errors"])


def test_the_saved_results_file_carries_the_same_statuses(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    download_and_process_data({"processing": PROCESSING}, output_dir)

    saved = io.load_json(output_dir / "processing_results.json")

    assert {name: entry["status"] for name, entry in saved["processing"].items()} == {
        "analysis": "not_executed",
        "normalization": "not_executed",
    }
    assert len(saved["errors"]) == len(PROCESSING)


@pytest.mark.parametrize(
    ("entry", "expected_in_reason"),
    [
        ({"method": "log2"}, "no 'type'"),
        ({"type": ""}, "no 'type'"),
        ({"type": None}, "no 'type'"),
        ("a bare string", "no 'type'"),
        ({"type": "bwa_mem", "threads": 8}, "'bwa_mem'"),
    ],
)
def test_unusual_entries_are_reported_not_crashed_on(tmp_path: Path, entry: Any, expected_in_reason: str) -> None:
    results = download_and_process_data({"processing": {"step": entry}}, tmp_path / "out")

    assert results["processing"]["step"]["status"] == "not_executed"
    assert expected_in_reason in results["processing"]["step"]["reason"]
    assert results["processing"]["step"]["config"] == entry


@pytest.mark.parametrize("sample_type", ["basic", "scientific", "advanced"])
def test_the_package_sample_configs_do_not_claim_completed_processing(tmp_path: Path, sample_type: str) -> None:
    sample = tmp_path / f"{sample_type}.json"
    create_sample_config(sample, sample_type)
    config: Dict[str, Any] = io.load_json(sample)
    config["downloads"] = {}  # the samples point at example.com; only the processing section is under test

    results = download_and_process_data(config, tmp_path / "out")

    assert set(results["processing"]) == set(config["processing"])
    assert {entry["status"] for entry in results["processing"].values()} == {"not_executed"}


def test_downloads_are_unaffected_and_still_succeed(tmp_path: Path) -> None:
    source = tmp_path / "source.json"
    source.write_text('{"source": "local fixture"}\n', encoding="utf-8")
    config = {
        "downloads": {"data": {"url": source.as_uri(), "filename": "copy.json"}},
        "processing": {"analyze": {"type": "json_analysis"}},
    }

    results = download_and_process_data(config, tmp_path / "out")

    assert results["downloads"]["data"]["status"] == "success"
    assert (tmp_path / "out" / "copy.json").read_text(encoding="utf-8") == source.read_text(encoding="utf-8")
    assert results["processing"]["analyze"]["status"] == "not_executed"


def test_a_config_without_processing_has_nothing_to_report(tmp_path: Path) -> None:
    source = tmp_path / "source.json"
    source.write_text("{}\n", encoding="utf-8")
    config = {"downloads": {"data": {"url": source.as_uri(), "filename": "copy.json"}}}

    results = download_and_process_data(config, tmp_path / "out")

    assert results["processing"] == {}
    assert results["errors"] == []


def test_run_config_based_workflow_is_not_successful_when_processing_did_not_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source.json"
    source.write_text("{}\n", encoding="utf-8")
    config_file = tmp_path / "with_processing.json"
    io.dump_json(
        {
            "downloads": {"data": {"url": source.as_uri(), "filename": "copy.json"}},
            "processing": {"analyze": {"type": "json_analysis"}},
        },
        config_file,
    )
    monkeypatch.chdir(tmp_path)

    results = run_config_based_workflow(config_file)

    assert results["downloads"]["data"]["status"] == "success"
    assert results["processing"]["analyze"]["status"] == "not_executed"
    assert results["success"] is False


def test_run_config_based_workflow_still_succeeds_when_there_is_nothing_to_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source.json"
    source.write_text("{}\n", encoding="utf-8")
    config_file = tmp_path / "downloads_only.json"
    io.dump_json({"downloads": {"data": {"url": source.as_uri(), "filename": "copy.json"}}}, config_file)
    monkeypatch.chdir(tmp_path)

    results = run_config_based_workflow(config_file)

    assert results["success"] is True


def test_failed_download_is_an_error(tmp_path: Path) -> None:
    """An unreachable download must keep the workflow from reporting success."""
    config = {"downloads": {"missing": {"url": "file:///nonexistent/metainformant-test.dat"}}}
    results = download_and_process_data(config, output_dir=tmp_path)

    assert results["downloads"]["missing"]["status"] == "failed"
    assert any("missing" in error for error in results["errors"])
