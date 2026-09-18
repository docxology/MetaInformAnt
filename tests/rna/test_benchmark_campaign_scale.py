"""Tests for scripts/rna/benchmark_campaign_scale.py (disposable fixtures, zero network).

Real-implementation policy: exercises the real ProgressDB, provenance sidecars,
and species discovery against tiny deterministic fixtures created in tmp_path.
The live data root is never opened.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

# Repo root: three levels up from tests/rna/, or the committed checkout when
# this file is executed outside the tree (workspace staging copies).
_REPO_ROOT = Path(os.environ.get("MAI_REPO_ROOT", Path(__file__).resolve().parents[2]))
SCRIPT = _REPO_ROOT / "scripts" / "rna" / "benchmark_campaign_scale.py"
spec = importlib.util.spec_from_file_location("benchmark_campaign_scale", SCRIPT)
benchmark_campaign_scale = importlib.util.module_from_spec(spec)
sys.modules["benchmark_campaign_scale"] = benchmark_campaign_scale
spec.loader.exec_module(benchmark_campaign_scale)

LIVE_DATA_ROOT = benchmark_campaign_scale.LIVE_DATA_ROOT


@pytest.fixture(scope="module")
def tiny_fixture(tmp_path_factory: pytest.TempPathFactory) -> benchmark_campaign_scale.CampaignFixture:
    """One small deterministic fixture (3 species x 40 rows) shared by tests."""

    target = tmp_path_factory.mktemp("bench_fixture")
    return benchmark_campaign_scale.build_campaign_fixture(
        target,
        120,
        species_count=3,
        seed=benchmark_campaign_scale.DEFAULT_SEED,
        provenance_sample_cap=2,
    )


class TestGuard:
    """The harness refuses any target at or beneath the live data root."""

    def test_refuses_live_data_root(self) -> None:
        with pytest.raises(SystemExit):
            benchmark_campaign_scale.guard_target_dir(LIVE_DATA_ROOT)

    def test_refuses_descendant_of_live_data_root(self) -> None:
        with pytest.raises(SystemExit):
            benchmark_campaign_scale.guard_target_dir(LIVE_DATA_ROOT / "work" / "bench")

    def test_allows_disposable_target(self, tmp_path: Path) -> None:
        resolved = benchmark_campaign_scale.guard_target_dir(tmp_path / "scratch")
        assert resolved == tmp_path / "scratch"

    def test_default_target_is_not_under_live_root(self) -> None:
        default_dir = _REPO_ROOT / "output" / "campaign_scale_benchmarks"
        resolved = benchmark_campaign_scale.guard_target_dir(default_dir)
        live = benchmark_campaign_scale._lexical_abspath(LIVE_DATA_ROOT)
        assert live not in resolved.parents and resolved != live


class TestFixture:
    """The synthetic fixture is deterministic and structurally complete."""

    def test_row_counts_and_species(self, tiny_fixture: benchmark_campaign_scale.CampaignFixture) -> None:
        total = sum(len(ids) for ids in tiny_fixture.cohort.values())
        assert total == 120
        assert sorted(tiny_fixture.cohort) == sorted(tiny_fixture.species)
        assert all(ids for ids in tiny_fixture.cohort.values())

    def test_state_partition_covers_cohort(self, tiny_fixture: benchmark_campaign_scale.CampaignFixture) -> None:
        for name, bucket in tiny_fixture.state_map.items():
            assert set(bucket) == {state for state, _ in benchmark_campaign_scale.STATE_WEIGHTS}
            assert sum(len(ids) for ids in bucket.values()) == len(tiny_fixture.cohort[name])

    def test_deterministic_given_same_seed(self, tmp_path: Path) -> None:
        first = benchmark_campaign_scale.build_campaign_fixture(
            tmp_path / "a", 40, species_count=2, seed=7, provenance_sample_cap=1
        )
        second = benchmark_campaign_scale.build_campaign_fixture(
            tmp_path / "b", 40, species_count=2, seed=7, provenance_sample_cap=1
        )
        assert first.cohort == second.cohort
        assert first.state_map == second.state_map

    def test_db_populated_and_queried(self, tiny_fixture: benchmark_campaign_scale.CampaignFixture) -> None:
        db = benchmark_campaign_scale.ProgressDB(db_path=tiny_fixture.db_path)
        try:
            counts = db.get_counts()
            assert sum(sum(states.values()) for states in counts.values()) == 120
            # Every species has a pending queue, so queue depth is meaningful.
            for name in tiny_fixture.species:
                assert db.get_samples(name, "pending")
        finally:
            db.close()

    def test_provenance_sidecars_classify_current(self, tiny_fixture: benchmark_campaign_scale.CampaignFixture) -> None:
        from metainformant.rna.engine.provenance import QUANT_STATUS_CURRENT, classify_quantification

        for name, dirs in tiny_fixture.sample_dirs.items():
            assert dirs, "provenance_sample_cap>0 must create sample dirs"
            for sample_dir in dirs:
                classification = classify_quantification(sample_dir, sample_dir.name, verify_content=True)
                assert classification["status"] == QUANT_STATUS_CURRENT, classification


class TestBenchmarks:
    """Benchmark functions return finite measurements over the tiny fixture."""

    def test_benchmarks_run(self, tiny_fixture: benchmark_campaign_scale.CampaignFixture) -> None:
        db = benchmark_campaign_scale.ProgressDB(db_path=tiny_fixture.db_path)
        try:
            discovery = benchmark_campaign_scale.bench_discovery(tiny_fixture, repeats=1)
            queue_depth = benchmark_campaign_scale.bench_queue_depth(tiny_fixture, db, repeats=1)
            retry = benchmark_campaign_scale.bench_retry_backoff(tiny_fixture, db, subset_size=4)
        finally:
            db.close()
        provenance_io = benchmark_campaign_scale.bench_provenance_io(tiny_fixture)

        assert discovery["cohort_rows"] == 120
        assert discovery["cohort_rows_per_second"] > 0
        assert queue_depth["total_task_rows"] == 120
        assert queue_depth["pending_queue_total"] >= 0
        assert retry["transition_pairs"] == 4
        assert retry["transition_seconds"] >= 0.0
        # Backdated downloading rows were reset by the stale scan.
        assert retry["stale_downloading_resets"] > 0
        assert provenance_io["sample_dirs"] == 6
        assert provenance_io["hashed_bytes"] > 0
        assert provenance_io["classification_status_counts"].get("current", 0) == 6


class TestReport:
    """Reports embed a machine-provenance block and render Markdown."""

    def test_machine_provenance_shape(self) -> None:
        provenance = benchmark_campaign_scale.machine_provenance(seed=1, scales=[10], target_dir=Path("/tmp/bench"))
        assert provenance["schema"] == benchmark_campaign_scale.BENCHMARK_SCHEMA
        assert provenance["git_commit"]
        assert provenance["seed"] == 1
        assert provenance["scales"] == [10]
        assert provenance["execution_context"] == "local disposable fixture (not a hosted run)"

    def test_safety_block_reports_zero_live_access(self, tmp_path: Path) -> None:
        result = benchmark_campaign_scale.run_scale_benchmark(
            40,
            target_dir=tmp_path / "scale_40",
            species_count=2,
            seed=benchmark_campaign_scale.DEFAULT_SEED,
            repeats=1,
            provenance_sample_cap=1,
            subset_size=4,
        )
        report = benchmark_campaign_scale.build_report([result], seed=7, scales=[40], target_dir=tmp_path)
        assert report["safety"]["live_data_root_read"] is False
        assert report["safety"]["live_db_opened"] is False
        assert report["scales"][0]["task_rows"] == 40

        markdown = benchmark_campaign_scale.render_markdown(report)
        assert "not a hosted run" in markdown or "local measurements only" in markdown
        assert "```json" in markdown
        provenance_json = json.loads(markdown.split("## Machine provenance")[1].split("```json")[1].split("```")[0])
        assert provenance_json["schema"] == benchmark_campaign_scale.BENCHMARK_SCHEMA
        assert "Budgets (advisory, descriptive)" in markdown
        assert "Metric definitions (telemetry)" in markdown


class TestMain:
    """The CLI end-to-end path works on a disposable target."""

    def test_main_small_scale(self, tmp_path: Path) -> None:
        target = tmp_path / "run"
        code = benchmark_campaign_scale.main(
            [
                "--scales",
                "60",
                "--species-count",
                "2",
                "--repeats",
                "1",
                "--target-dir",
                str(target),
            ]
        )
        assert code == 0
        # Fixture is deleted by default; reports remain.
        report_json = json.loads((target / "report.json").read_text(encoding="utf-8"))
        assert report_json["scales"][0]["task_rows"] == 60
        assert (target / "report.md").is_file()
        assert not list(target.glob("scale_60"))

    def test_main_keep_fixture(self, tmp_path: Path) -> None:
        target = tmp_path / "run_kept"
        code = benchmark_campaign_scale.main(
            [
                "--scales",
                "60",
                "--species-count",
                "2",
                "--repeats",
                "1",
                "--target-dir",
                str(target),
                "--keep-fixture",
            ]
        )
        assert code == 0
        assert list(target.glob("scale_60/pipeline_progress.db"))

    def test_main_refuses_live_root(self) -> None:
        with pytest.raises(SystemExit):
            benchmark_campaign_scale.main(["--scales", "60", "--target-dir", str(LIVE_DATA_ROOT / "bench")])
