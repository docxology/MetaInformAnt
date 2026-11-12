"""Tests for the campaign environment preflight and durable error classes.

All tests exercise the real preflight probes against real directories and the
real progress-database classifier following the real-implementation policy.
"""

from __future__ import annotations

import json
import os
import stat
import struct
from pathlib import Path

import pytest

from metainformant.rna.engine.preflight import (
    PROBE_FILE_NAME,
    PRODUCER_TOOLS,
    PreflightError,
    ReuseProbe,
    check_free_space,
    check_kallisto_compatibility,
    kallisto_index_version,
    main,
    probe_data_root_writable,
    resolve_amalgkit_cli,
    resolve_producer_tools,
    reuse_dry_run,
    run_campaign_preflight,
)
from metainformant.rna.engine.progress_db import classify_sample_error


def _fake_tools(bin_dir: Path, *, kallisto_version: str = "0.52.0", skip: tuple[str, ...] = ()) -> Path:
    """Real executables for every producer tool; ``kallisto version`` prints ``kallisto_version``."""

    bin_dir.mkdir(parents=True, exist_ok=True)
    for tool in PRODUCER_TOOLS:
        if tool in skip:
            continue
        script = bin_dir / tool
        if tool == "kallisto":
            script.write_text(
                f'#!/bin/sh\nif [ "$1" = version ]; then echo "kallisto, version {kallisto_version}"; fi\n'
            )
        else:
            script.write_text("#!/bin/sh\nexit 0\n")
        script.chmod(0o755)
    return bin_dir


def _run_info(quant_dir: Path, run: str, *, kallisto_version: str = "0.52.0", index_version: int = 13) -> Path:
    sample = quant_dir / run
    sample.mkdir(parents=True)
    path = sample / f"{run}_run_info.json"
    path.write_text(json.dumps({"kallisto_version": kallisto_version, "index_version": index_version}))
    return path


@pytest.fixture
def no_space_floors(monkeypatch: pytest.MonkeyPatch) -> None:
    """Tests that are not about free space must not depend on this host's disks."""

    monkeypatch.setenv("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0")
    monkeypatch.setenv("AMALGKIT_MIN_SYSTEM_FREE_GB", "0")


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

    def test_happy_path_returns_resolved_facts(self, tmp_path: Path, no_space_floors: None):
        fake_bin = _fake_tools(tmp_path / "bin")
        data_root = tmp_path / "root"
        data_root.mkdir()

        facts = run_campaign_preflight(data_root, search_path=str(fake_bin))

        assert facts["data_root"] == str(data_root.resolve())
        assert facts["amalgkit_cli"] == str(fake_bin / "amalgkit")
        assert facts["tool_fasterq_dump"] == str(fake_bin / "fasterq-dump")
        assert facts["tool_vdb_validate"] == str(fake_bin / "vdb-validate")
        assert float(facts["data_root_free_gib"]) > 0
        assert list(data_root.iterdir()) == []

    def test_missing_binary_is_named(self, tmp_path: Path, no_space_floors: None):
        fake_bin = _fake_tools(tmp_path / "bin", skip=("kallisto", "pigz"))
        data_root = tmp_path / "root"
        data_root.mkdir()

        with pytest.raises(PreflightError, match="not on the effective PATH: kallisto, pigz") as excinfo:
            run_campaign_preflight(data_root, search_path=str(fake_bin))
        assert "amalgkit" not in str(excinfo.value).split("PATH:")[1].splitlines()[0]
        with pytest.raises(PreflightError, match="fasterq-dump, prefetch, vdb-validate, pigz"):
            resolve_producer_tools(search_path=str(tmp_path / "empty"))

    def test_kallisto_must_match_the_recorded_outputs(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, no_space_floors: None
    ):
        data_root = tmp_path / "root"
        data_root.mkdir()
        quant = data_root / "species_a" / "work" / "quant"
        _run_info(quant, "SRR1")
        empty = data_root / "species_b" / "work" / "quant"
        empty.mkdir(parents=True)
        index = tmp_path / "a.idx"
        index.write_bytes(struct.pack("<Q", 13) + b"rest of a kallisto index")
        matching = _fake_tools(tmp_path / "bin_ok")
        facts = check_kallisto_compatibility(
            str(matching / "kallisto"), {"species_a": quant, "species_b": empty}, index_paths={"species_a": index}
        )
        assert facts["installed"] == "0.52.0"
        assert facts["checked"]["species_a"]["index_header_version"] == 13
        assert facts["without_outputs"] == ["species_b"]
        assert kallisto_index_version(index) == 13

        drifted = _fake_tools(tmp_path / "bin_new", kallisto_version="0.51.1")
        with pytest.raises(PreflightError, match="installed kallisto 0.51.1 differs from 0.52.0"):
            run_campaign_preflight(data_root, search_path=str(drifted), kallisto_quant_dirs={"species_a": quant})
        monkeypatch.setenv("AMALGKIT_ALLOW_KALLISTO_DRIFT", "1")
        facts = run_campaign_preflight(data_root, search_path=str(drifted), kallisto_quant_dirs={"species_a": quant})
        assert facts["kallisto_drift_allowed"] == "yes"
        monkeypatch.delenv("AMALGKIT_ALLOW_KALLISTO_DRIFT")

        index.write_bytes(struct.pack("<Q", 12) + b"an older index format")
        with pytest.raises(PreflightError, match="has format version 12, outputs recorded index_version 13"):
            run_campaign_preflight(
                data_root,
                search_path=str(matching),
                kallisto_quant_dirs={"species_a": quant},
                kallisto_index_paths={"species_a": index},
            )

    def test_free_space_floors(self, tmp_path: Path):
        data_root = tmp_path / "root"
        data_root.mkdir()
        facts = check_free_space(data_root, min_data_root_free_gb=0, min_system_free_gb=0, system_root=tmp_path)
        assert float(facts["data_root_free_gib"]) > 0 and float(facts["system_free_gib"]) > 0
        with pytest.raises(PreflightError, match="below the 1e\\+09 GiB floor \\(AMALGKIT_MIN_EXTERNAL_FREE_GB\\)"):
            check_free_space(data_root, min_data_root_free_gb=1e9, min_system_free_gb=0, system_root=tmp_path)
        with pytest.raises(PreflightError, match="AMALGKIT_MIN_SYSTEM_FREE_GB"):
            check_free_space(data_root, min_data_root_free_gb=0, min_system_free_gb=1e9, system_root=tmp_path)
        # A system root that does not exist (Linux, CI) is skipped.
        facts = check_free_space(data_root, min_data_root_free_gb=0, system_root=tmp_path / "absent")
        assert "system_free_gib" not in facts

    def test_reuse_dry_run_refuses_a_mass_quarantine(self, tmp_path: Path):
        quant = tmp_path / "quant"
        runs = [f"SRR{i:04d}" for i in range(200)]
        for run in runs:
            (quant / run).mkdir(parents=True)
            (quant / run / "verdict").write_text("current")
        calls: list[str] = []

        def classify(path: Path, run: str) -> dict[str, str]:
            calls.append(run)
            return {"status": (path / "verdict").read_text(), "reason": "fixture verdict file"}

        report = reuse_dry_run([ReuseProbe("species_a", quant, runs, classify)], fraction=0.01)
        assert report["species_a"] == {"recorded": 200, "sampled": 2, "would_quarantine": 0, "examples": []}
        assert len(calls) == 2
        again = reuse_dry_run([ReuseProbe("species_a", quant, runs, classify)], fraction=0.01)
        assert again == report and calls[:2] == calls[2:]  # the sample is deterministic

        for run in runs[:150]:
            (quant / run / "verdict").write_text("invalid")
        with pytest.raises(PreflightError, match="species_a: .* would be quarantined"):
            reuse_dry_run([ReuseProbe("species_a", quant, runs, classify)], fraction=0.05)
        assert reuse_dry_run([ReuseProbe("species_b", quant, [], classify)])["species_b"]["sampled"] == 0

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

    def test_cli_exit_codes(self, tmp_path: Path, capsys: pytest.CaptureFixture[str], no_space_floors: None):
        data_root = tmp_path / "root"
        data_root.mkdir()

        assert main(["--data-root", str(data_root), "--search-path", ""]) == 1
        assert "PREFLIGHT FAILED" in capsys.readouterr().err

        fake_bin = _fake_tools(tmp_path / "bin")

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
        assert classify_sample_error("fasterq-dump timeout for SRR1 (>2h)") == "extraction_timeout"

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

    def test_preflight_failure_aborts_before_discovery(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        from metainformant.rna.engine import streaming_orchestrator

        def _fail(data_root: Path, **_kwargs: object) -> dict[str, str]:
            raise PreflightError(f"simulated broken environment: {data_root}")

        monkeypatch.setattr(streaming_orchestrator, "run_campaign_preflight", _fail)
        (tmp_path / "data").mkdir()
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path / "data"))
        orchestrator = streaming_orchestrator.StreamingPipelineOrchestrator(
            config_dir=tmp_path / "configs",
            log_dir=tmp_path / "logs",
            db_path=tmp_path / "progress.db",
        )
        orchestrator.db.init_species("test_species", ["SRR1"])
        orchestrator.db.set_state("test_species", "SRR1", "quantifying")

        with pytest.raises(PreflightError, match="simulated broken environment") as excinfo:
            orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)
        assert streaming_orchestrator.producer_exit_code_for_exception(excinfo.value) == 69

        # Discovery never started (the invalid configuration was never parsed) and no startup reset ran.
        assert orchestrator.db.get_total_counts() == {"test_species": 1}
        assert orchestrator.db.get_state("test_species", "SRR1") == "quantifying"

    def test_real_preflight_refuses_a_host_without_tools(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, no_space_floors: None
    ):
        from metainformant.rna.engine import streaming_orchestrator

        data_root = tmp_path / "data"
        data_root.mkdir()
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
        monkeypatch.setenv("PATH", str(_fake_tools(tmp_path / "bin", skip=("kallisto",))))
        orchestrator = streaming_orchestrator.StreamingPipelineOrchestrator(
            config_dir=tmp_path / "configs",
            log_dir=tmp_path / "logs",
            db_path=tmp_path / "progress.db",
        )
        with pytest.raises(PreflightError, match="not on the effective PATH: kallisto"):
            orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)
        assert orchestrator.db.get_total_counts() == {}

    def test_preflight_runs_with_configured_data_root(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        from metainformant.rna.engine import streaming_orchestrator

        observed: dict[str, Path] = {}

        def _record(data_root: Path, **_kwargs: object) -> dict[str, str]:
            observed["data_root"] = Path(data_root)
            return {"data_root": str(data_root), "amalgkit_cli": "fake"}

        monkeypatch.setattr(streaming_orchestrator, "run_campaign_preflight", _record)
        (tmp_path / "data").mkdir()
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
