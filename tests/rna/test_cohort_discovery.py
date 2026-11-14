"""Producer discovery in frozen-cohort mode: frozen tasks only, no metadata query, no selection rewrite."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

import metainformant.rna.engine.streaming_orchestrator as orchestrator_module
from metainformant.rna.engine.cohort import (
    CohortError,
    cohort_task_table_name,
    install_cohort,
    load_frozen_cohort,
    quant_batch_run,
)
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
)
from tests.rna.frozen_cohort_support import (
    COHORT_ID,
    CONFIG_NAME,
    SPECIES,
    DataRoot,
    build_frozen_source,
    write_species_yaml,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _forbidden(*_args: object, **_kwargs: object) -> object:
    raise AssertionError("cohort discovery must not run Amalgkit metadata/select or any subprocess")


class Campaign:
    """A data root, the checked-in YAML and an installed frozen cohort whose digests match them."""

    def __init__(self, base: Path, monkeypatch: pytest.MonkeyPatch, *, forbid_subprocess: bool = True) -> None:
        self.base = base
        self.data = DataRoot(base / "data")
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(self.data.base))
        self.config_dir = base / "config"
        self.config_dir.mkdir()
        self.yaml = write_species_yaml(self.config_dir, self.data.index_dir)
        build_frozen_source(base / "frozen", config_bytes=self.yaml.read_bytes())
        install_cohort(base / "frozen", base / "cohort", cohort_id=COHORT_ID)
        self.cohort = load_frozen_cohort(base / "cohort")
        self.metadata_dir = self.data.work / "metadata"
        self.metadata_dir.mkdir()
        self.selected = self.metadata_dir / "metadata_selected.tsv"
        self.selected.write_bytes(b"run\tis_sampled\r\nSRR000\tyes\r\n")
        self.table = self.metadata_dir / cohort_task_table_name(COHORT_ID)
        if forbid_subprocess:
            # ``orchestrator_module.subprocess`` is the process-wide module, so this also guards helpers.
            monkeypatch.setattr(orchestrator_module.subprocess, "run", _forbidden)
            monkeypatch.setattr(StreamingPipelineOrchestrator, "run_tissue_normalization", _forbidden)

    def orchestrator(self, **kwargs: object) -> StreamingPipelineOrchestrator:
        kwargs.setdefault("cohort", self.cohort)
        return StreamingPipelineOrchestrator(
            self.config_dir, self.base / "logs", db_path=self.base / "progress.db", **kwargs
        )


@pytest.fixture
def campaign(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Campaign:
    return Campaign(tmp_path, monkeypatch)


def test_discovery_uses_only_frozen_tasks_and_never_rewrites_metadata(campaign: Campaign) -> None:
    campaign.data.write_sample("SRR101", config=campaign.yaml)  # reusable under the cohort contract
    campaign.data.write_sample("SRR102", config=campaign.data.config)  # other configuration bytes
    selected_before = (campaign.selected.read_bytes(), campaign.selected.stat().st_mtime_ns)
    orchestrator = campaign.orchestrator()
    orchestrator.db.init_species(SPECIES, ["SRR102"])
    orchestrator.db.set_state(SPECIES, "SRR102", "quantified")

    tasks = orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)

    assert [task["srr"] for task in tasks] == ["SRR103", "SRR105"]  # by total_bases
    positions = {"SRR103": 1, "SRR101": 2, "SRR102": 3, "SRR105": 4}  # frozen order minus SRR104
    for task in tasks:
        assert task["batch_idx"] == positions[task["srr"]]
        assert quant_batch_run(campaign.table, task["batch_idx"]) == task["srr"]
        assert task["metadata_path"] == campaign.table and task["cohort_id"] == COHORT_ID
        assert task["config_path"] == campaign.yaml
    assert tasks[0]["expected_paired"] is True and tasks[1]["expected_paired"] is False
    assert (campaign.selected.read_bytes(), campaign.selected.stat().st_mtime_ns) == selected_before
    assert not (campaign.data.work / ".metainformant_metadata_provenance.json").exists()
    assert sorted(path.name for path in campaign.metadata_dir.iterdir()) == sorted(
        ["metadata_selected.tsv", campaign.table.name]
    )
    db = orchestrator.db
    assert db.get_state(SPECIES, "SRR101") == "quantified"
    assert db.get_state(SPECIES, "SRR102") == "quarantined"
    assert db.get_state(SPECIES, "SRR103") == "pending"
    exclusions = db.get_exclusions(SPECIES)
    assert [(row["srr_id"], row["recorded_by"]) for row in exclusions] == [("SRR104", f"cohort:{COHORT_ID}")]
    assert db.get_state(SPECIES, "SRR104") is None

    # Resume: identical table bytes are not rewritten; the size budget defers large tasks.
    mtime = campaign.table.stat().st_mtime_ns
    resumed = orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=1000.0, threads=1)
    assert [task["srr"] for task in resumed] == ["SRR103"]
    assert campaign.table.stat().st_mtime_ns == mtime

    # An operator permanent_drop leaves the schedule and the task table; batches follow the table.
    db.record_exclusions(SPECIES, [{"srr_id": "SRR103", "reason_code": "permanent_drop"}])
    after_drop = orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)
    assert [(task["srr"], task["batch_idx"]) for task in after_drop] == [("SRR105", 3)]
    assert b"SRR103" not in campaign.table.read_bytes()
    assert quant_batch_run(campaign.table, 3) == "SRR105"
    orchestrator.db.close()


def test_discovery_refuses_changed_configuration_bytes(campaign: Campaign) -> None:
    campaign.yaml.write_bytes(campaign.yaml.read_bytes() + b"# edited\n")
    orchestrator = campaign.orchestrator()
    with pytest.raises(CohortError, match="differs from the frozen cohort"):
        orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=1000.0, threads=1)
    assert not campaign.table.exists()
    orchestrator.db.close()


def test_discovery_refuses_other_reference_index_bytes(campaign: Campaign) -> None:
    campaign.data.index.write_bytes(b"a rebuilt index\n")
    orchestrator = campaign.orchestrator()
    with pytest.raises(CohortError, match="differ from the cohort index"):
        orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=1000.0, threads=1)
    orchestrator.db.close()


def test_discovery_refuses_configurations_outside_the_cohort(campaign: Campaign) -> None:
    orchestrator = campaign.orchestrator()
    with pytest.raises(CohortError, match="not part of cohort"):
        orchestrator.discover_species_tasks("amalgkit_other_species.yaml", max_gb=1000.0, threads=1)
    orchestrator.db.close()


def test_progress_db_bound_to_another_cohort_is_refused(campaign: Campaign) -> None:
    with ProgressDB(campaign.base / "progress.db") as db:
        db.migrate_to_cohort(campaign.cohort)
        db._conn.execute("UPDATE campaign_meta SET value = 'other' WHERE key = 'cohort_id'")
        db._conn.commit()
    with pytest.raises(CohortError, match="bound to cohort other"):
        campaign.orchestrator()


def test_run_all_passes_the_task_table_to_every_sample(campaign: Campaign, monkeypatch: pytest.MonkeyPatch) -> None:
    orchestrator = campaign.orchestrator()
    calls: list[tuple[str, int, object]] = []

    def fake_process(srr: str, batch: int, *_args: object, **kwargs: object) -> dict[str, object]:
        calls.append((srr, batch, kwargs.get("metadata_path")))
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(
        orchestrator_module, "run_campaign_preflight", lambda *_a, **_k: {"data_root": "x", "amalgkit_cli": "x"}
    )
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert sorted(calls) == [
        ("SRR101", 2, campaign.table),
        ("SRR102", 3, campaign.table),
        ("SRR103", 1, campaign.table),
        ("SRR105", 4, campaign.table),
    ]
    orchestrator.db.close()


def test_quant_sample_refuses_a_batch_that_selects_another_run(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    orchestrator = campaign.orchestrator()
    orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)
    monkeypatch.setattr(orchestrator_module, "_run_command_in_process_group", _forbidden)
    ok, error = orchestrator.quant_sample(campaign.yaml, 1, SPECIES, 1, "SRR105", metadata_path=campaign.table)
    assert not ok and "Quant batch mismatch: row 1" in error and "'SRR103'" in error
    outside = campaign.base / "elsewhere.tsv"
    outside.write_bytes(campaign.table.read_bytes())
    ok, error = orchestrator.quant_sample(campaign.yaml, 1, SPECIES, 1, "SRR103", metadata_path=outside)
    assert not ok and "not a regular file in" in error
    orchestrator.db.close()


def test_non_cohort_discovery_refuses_ncbi_refresh_unless_allowed(campaign: Campaign) -> None:
    campaign.selected.unlink()
    orchestrator = campaign.orchestrator(cohort=None, allow_metadata_refresh=False)
    with pytest.raises(RuntimeError, match="metadata refresh is disabled"):
        orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=1000.0, threads=1)
    orchestrator.db.close()


# ---------------------------------------------------------------------------
# Producer scripts
# ---------------------------------------------------------------------------


def _run(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(REPO_ROOT / "src")
    environment.pop("AMALGKIT_MAX_SAMPLE_GB", None)
    environment.pop("PIPELINE_MAX_GB", None)
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "rna" / script), *args],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )


@pytest.mark.parametrize("script", ["run_all_species.py", "process_species.py"])
def test_scripts_resolve_species_from_the_cohort(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, script: str) -> None:
    campaign = Campaign(tmp_path, monkeypatch, forbid_subprocess=False)
    common = ["--config-dir", str(campaign.config_dir), "--data-root", str(campaign.data.base), "--dry-run"]
    species = ["--species", SPECIES]
    result = _run(script, *common, *species, "--cohort", str(campaign.base / "cohort"))
    assert result.returncode == 0, result.stderr
    assert f"Frozen cohort: {COHORT_ID}" in result.stdout
    assert "1000.0" in result.stdout  # the launcher's per-sample budget is the default
    assert not (campaign.data.base / "pipeline_progress.db").exists()
    conflict = _run(script, *common, *species, "--cohort", str(campaign.base / "cohort"), "--refresh-metadata")
    assert conflict.returncode == 2 and "--refresh-metadata cannot be combined" in conflict.stderr
    unknown = _run(script, *common, "--species", "other_species", "--cohort", str(campaign.base / "cohort"))
    assert unknown.returncode == 2 and "not in cohort" in unknown.stderr
    project_configs = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
    if project_configs.is_dir():
        non_cohort = _run(
            script,
            "--config-dir",
            str(project_configs),
            "--data-root",
            str(campaign.data.base),
            "--dry-run",
            "--species",
            "apis_mellifera",
        )
        assert non_cohort.returncode == 0, non_cohort.stderr
        assert "Metadata refresh (NCBI): refused" in non_cohort.stdout
