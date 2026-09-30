"""Tests for the campaign environment preflight and durable error classes.

All tests exercise the real preflight probes against real directories and the
real progress-database classifier following the real-implementation policy.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from metainformant.rna.engine.progress_db import classify_sample_error
from metainformant.rna.engine.preflight import (
    PROBE_FILE_NAME,
    PreflightError,
    main,
    probe_data_root_writable,
    resolve_amalgkit_cli,
    run_campaign_preflight,
)


class TestDataRootProbe:
    """The data-root probe must exercise real create/rename/stat/unlink calls."""

    def test_writable_root_passes_without_leftovers(self, tmp_path: Path):
        resolved = probe_data_root_writable(tmp_path)

        assert resolved == tmp_path.resolve()
        assert not (tmp_path / PROBE_FILE_NAME).exists()
        assert not (tmp_path / f"{PROBE_FILE_NAME}.moved").exists()
        assert list(tmp_path.iterdir()) == []

    def test_read_only_root_raises_writable_error(self, tmp_path: Path):
        os.chmod(tmp_path, stat.S_IRUSR | stat.S_IXUSR)
        try:
            with pytest.raises(PreflightError, match="not writable"):
                probe_data_root_writable(tmp_path)
        finally:
            os.chmod(tmp_path, stat.S_IRWXU)
        assert list(tmp_path.iterdir()) == []

    def test_missing_root_raises_writable_error(self, tmp_path: Path):
        with pytest.raises(PreflightError, match="not writable"):
            probe_data_root_writable(tmp_path / "absent")


class TestAmalgkitResolution:
    """CLI resolution must mirror the producer's bare-command subprocess PATH."""

    def test_missing_cli_raises(self):
        with pytest.raises(PreflightError, match="amalgkit.*not found"):
            resolve_amalgkit_cli(search_path="")

    def test_resolved_from_explicit_search_path(self, tmp_path: Path):
        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        cli = fake_bin / "amalgkit"
        cli.write_text("#!/bin/sh\nexit 0\n")
        cli.chmod(0o755)

        resolved = resolve_amalgkit_cli(search_path=str(fake_bin))

        assert resolved == str(cli)


class TestCampaignPreflight:
    """The combined preflight reports every failed check, not only the first."""

    def test_happy_path_returns_resolved_facts(self, tmp_path: Path):
        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        cli = fake_bin / "amalgkit"
        cli.write_text("#!/bin/sh\nexit 0\n")
        cli.chmod(0o755)
        data_root = tmp_path / "root"
        data_root.mkdir()

        facts = run_campaign_preflight(data_root, search_path=str(fake_bin))

        assert facts["data_root"] == str(data_root.resolve())
        assert facts["amalgkit_cli"] == str(cli)

    def test_collects_all_failures(self, tmp_path: Path):
        read_only = tmp_path / "read_only"
        read_only.mkdir()
        os.chmod(read_only, stat.S_IRUSR | stat.S_IXUSR)
        try:
            with pytest.raises(PreflightError) as excinfo:
                run_campaign_preflight(read_only, search_path="")
        finally:
            os.chmod(read_only, stat.S_IRWXU)

        message = str(excinfo.value)
        assert "not writable" in message
        assert "amalgkit" in message

    def test_cli_exit_codes(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        data_root = tmp_path / "root"
        data_root.mkdir()

        assert main(["--data-root", str(data_root), "--search-path", ""]) == 1
        assert "PREFLIGHT FAILED" in capsys.readouterr().err

        fake_bin = tmp_path / "bin"
        fake_bin.mkdir()
        cli = fake_bin / "amalgkit"
        cli.write_text("#!/bin/sh\nexit 0\n")
        cli.chmod(0o755)

        assert main(["--data-root", str(data_root), "--search-path", str(fake_bin)]) == 0
        out = capsys.readouterr().out
        assert "preflight: OK" in out
        assert "amalgkit_cli" in out


class TestSampleErrorClasses:
    """Durable failure classes let the M-02 audit separate environment damage."""

    def test_environment_write_denied(self):
        error = "Unexpected sample task error: [Errno 1] Operation not permitted: 'a' -> 'b'"
        assert classify_sample_error(error) == "environment_write_denied"

    def test_timeout_classes(self):
        assert classify_sample_error("Quant timeout (>2h)") == "quantification_timeout"
        assert (
            classify_sample_error("fasterq-dump timeout for SRR1 (>2h)")
            == "extraction_timeout"
        )

    def test_environment_missing_tool(self):
        error = "Quant exception batch 2: [Errno 2] No such file or directory: 'amalgkit'"
        assert classify_sample_error(error) == "environment_missing_tool"

    def test_transfer_and_quantification_classes(self):
        assert classify_sample_error("Download Failed (all sources: ENA FTP/HTTP, NCBI)") == (
            "transfer_all_sources_failed"
        )
        assert classify_sample_error("Quantification Failed") == "quantification_failed"

    def test_unrecorded_and_unclassified(self):
        assert classify_sample_error(None) == "unrecorded"
        assert classify_sample_error("") == "unrecorded"
        assert classify_sample_error("Something novel happened") == "unclassified"


class TestOrchestratorPreflightWiring:
    """run_all must run the preflight before any discovery or scheduling work."""

    def test_preflight_failure_aborts_before_discovery(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        from metainformant.rna.engine import streaming_orchestrator

        def _fail(data_root: Path, **_kwargs: object) -> dict[str, str]:
            raise PreflightError(f"simulated broken environment: {data_root}")

        monkeypatch.setattr(streaming_orchestrator, "run_campaign_preflight", _fail)
        orchestrator = streaming_orchestrator.StreamingPipelineOrchestrator(
            config_dir=tmp_path / "configs",
            log_dir=tmp_path / "logs",
            db_path=tmp_path / "progress.db",
        )

        with pytest.raises(PreflightError, match="simulated broken environment"):
            orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)

        # Discovery never started: the invalid configuration was never parsed.
        assert orchestrator.db.get_total_counts() == {}

    def test_preflight_runs_with_configured_data_root(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        from metainformant.rna.engine import streaming_orchestrator

        observed: dict[str, Path] = {}

        def _record(data_root: Path, **_kwargs: object) -> dict[str, str]:
            observed["data_root"] = Path(data_root)
            return {"data_root": str(data_root), "amalgkit_cli": "fake"}

        monkeypatch.setattr(streaming_orchestrator, "run_campaign_preflight", _record)
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path / "data"))
        orchestrator = streaming_orchestrator.StreamingPipelineOrchestrator(
            config_dir=tmp_path / "configs",
            log_dir=tmp_path / "logs",
            db_path=tmp_path / "progress.db",
        )

        # A missing configuration completes discovery with zero tasks; the
        # recorded preflight root proves the preflight resolved the
        # configured data root before any of that work.
        orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)

        assert observed["data_root"] == (tmp_path / "data").resolve()
