"""Tests for the standalone MCP-adjacent amalgkit monitor."""

from __future__ import annotations

from pathlib import Path

from metainformant.mcp.tools import amalgkit_monitor


def test_parse_log_progress_reads_latest_progress_line(tmp_path: Path, monkeypatch) -> None:
    """Progress parsing uses the real filesystem and the module's documented log path."""
    log_path = tmp_path / "output" / "amalgkit" / "run_all_species_incremental.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text("starting\n[3/10] SRR000003 complete\ntrailing status\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)

    progress = amalgkit_monitor.parse_log_progress()

    assert progress["processed"] == 3
    assert progress["total"] == 10
    assert progress["percent"] == 30.0
    assert progress["last_line"] == "[3/10] SRR000003 complete"


def test_parse_log_progress_handles_missing_log(tmp_path: Path, monkeypatch) -> None:
    """Missing logs return a structured stopped-progress response."""
    monkeypatch.chdir(tmp_path)

    progress = amalgkit_monitor.parse_log_progress()

    assert progress == {"processed": 0, "total": 0, "last_line": "No log file found"}


def test_process_detection_ignores_monitor_and_matches_workflows() -> None:
    """Process matching distinguishes the monitor from real workflow commands."""
    assert not amalgkit_monitor.is_pipeline_cmdline(
        ["uv", "run", "python", "-m", "metainformant.mcp.tools.amalgkit_monitor"]
    )
    assert amalgkit_monitor.is_pipeline_cmdline(["python", "scripts/rna/run_workflow.py", "--status"])
    assert amalgkit_monitor.is_pipeline_cmdline(["amalgkit", "quant", "--threads", "4"])
