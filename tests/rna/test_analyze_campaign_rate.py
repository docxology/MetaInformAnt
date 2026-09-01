"""Tests for scripts/rna/analyze_campaign_rate.py (zero-mocks: real files, real parsing)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "rna" / "analyze_campaign_rate.py"
spec = importlib.util.spec_from_file_location("analyze_campaign_rate", SCRIPT)
analyze_campaign_rate = importlib.util.module_from_spec(spec)
sys.modules["analyze_campaign_rate"] = analyze_campaign_rate
spec.loader.exec_module(analyze_campaign_rate)

PREFIX = "metainformant.rna.engine.streaming_orchestrator - INFO - "


def _line(ts: str, msg: str) -> str:
    return f"2026-08-31 {ts},000 - {PREFIX}{msg}\n"


SYNTHETIC_LOG = (
    _line("18:10:02", "=== Phase 1 start ===")
    + _line("18:10:05", "[1/10] apis_mellifera SRR1111111: Done")
    + _line("19:15:00", "Downloaded SRR2222222_1.fastq.gz: 0.40 GB")
    + _line("19:20:00", "[2/10] apis_mellifera SRR2222222: Done")
    + _line("20:30:00", "[3/10] solenopsis_invicta SRR3333333: Failed")
    + _line("21:45:00", "[4/10] solenopsis_invicta SRR4444444: Done")
    + _line("22:10:02", "[5/10] apis_mellifera SRR5555555: Done")
)


@pytest.fixture()
def synthetic_log(tmp_path: Path) -> Path:
    log = tmp_path / "streaming_orchestrator_test.log"
    log.write_text(SYNTHETIC_LOG, encoding="utf-8")
    return log


def test_synthetic_log_task_counts(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    assert result["tasks"]["done"] == 4
    assert result["tasks"]["failed"] == 1
    assert result["tasks"]["finished_total"] == 5
    assert result["tasks"]["declared_total_from_log"] == 10


def test_synthetic_log_per_hour(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    per_hour = result["per_hour_counts"]
    assert per_hour["2026-08-31 18:00"] == 1
    assert per_hour["2026-08-31 19:00"] == 1
    assert per_hour["2026-08-31 20:00"] == 1
    assert per_hour["2026-08-31 21:00"] == 1
    assert per_hour["2026-08-31 22:00"] == 1


def test_synthetic_log_download_throughput(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    dl = result["downloads"]
    assert dl["files"] == 1
    assert dl["total_gb"] == 0.4
    # Window 18:10:02 -> 22:10:02 = 4h exactly.
    assert dl["gb_per_hour"] == pytest.approx(0.1, rel=1e-3)
    assert dl["mb_per_s_equivalent"] == pytest.approx(round(0.4 * 1024 / (4 * 3600), 2), rel=1e-2)


def test_synthetic_log_species_mix_trailing(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log, trailing_hours=1.0)
    mix = result["species_mix_trailing_hours"]["counts"]
    # Trailing 1h from 22:10 -> includes the 21:45 and 22:10 Done events.
    assert mix == {"solenopsis_invicta": 1, "apis_mellifera": 1}


def test_synthetic_log_cumulative_rate(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    # 5 finished / 4h window.
    assert result["cumulative_finished_rate_per_hour"] == pytest.approx(1.25, rel=1e-3)


def test_eta_estimate_with_explicit_total(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    eta = analyze_campaign_rate.eta_estimate(result, 10)
    assert eta is not None
    assert eta["remaining"] == 5
    assert eta["eta_hours"] == pytest.approx(4.0, rel=1e-3)
    assert "caveat" in eta


def test_eta_estimate_none_without_total(synthetic_log: Path) -> None:
    result = analyze_campaign_rate.analyze_log(synthetic_log)
    # declared_total_from_log=10 present, so ETA should use it.
    eta = analyze_campaign_rate.eta_estimate(result, None)
    assert eta is not None and eta["total"] == 10


def test_empty_log_returns_error(tmp_path: Path) -> None:
    log = tmp_path / "empty.log"
    log.write_text("", encoding="utf-8")
    result = analyze_campaign_rate.analyze_log(log)
    assert "error" in result


def test_cli_json_mode(synthetic_log: Path, capsys: pytest.CaptureFixture) -> None:
    rc = analyze_campaign_rate.main([str(synthetic_log), "--json", "--total", "10"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["tasks"]["finished_total"] == 5
    assert payload["eta"]["total"] == 10


def test_cli_missing_log_returns_error() -> None:
    rc = analyze_campaign_rate.main(["/nonexistent/nope.log"])
    assert rc == 2


LIVE_LOG = Path("/Volumes/external_drive/Data/amalgkit/logs/streaming_orchestrator_20260831_181002.log")


@pytest.mark.skipif(not LIVE_LOG.is_file(), reason="live campaign log not present")
def test_live_log_integration() -> None:
    result = analyze_campaign_rate.analyze_log(LIVE_LOG)
    assert result["tasks"]["finished_total"] > 0
    assert result["downloads"]["files"] > 0
    assert result["downloads"]["total_gb"] > 0
    assert result["cumulative_finished_rate_per_hour"] > 0
