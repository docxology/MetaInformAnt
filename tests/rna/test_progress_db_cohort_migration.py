"""Progress DB: non-destructive pruning, cohort migration, adoption and cohort reconcile."""

from __future__ import annotations

import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.cohort import (
    CohortError,
    install_cohort,
    load_frozen_cohort,
)
from metainformant.rna.engine.progress_db import (
    COHORT_MIGRATION_ERROR_PREFIX,
    RETRYABLE_ENVIRONMENT_CLASSES,
    RETRYABLE_RESET_CLASSES,
    ProgressDB,
    classify_sample_error,
)
from tests.rna.frozen_cohort_support import COHORT_ID, SPECIES, build_frozen_source

PRODUCER_LOCK_HOLDER = """
import sys, time
from metainformant.rna.engine.campaign_lock import acquire
lock = acquire(sys.argv[1], "producer")
print("held", flush=True)
time.sleep(600)
"""

EPERM = "Unexpected sample task error: [Errno 1] Operation not permitted: '/Volumes/external_drive/x'"
TRANSFER = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"
KALLISTO = "Quantification Failed: kallisto exited 1"


@pytest.fixture
def cohort(tmp_path: Path):
    build_frozen_source(tmp_path / "frozen")
    install_cohort(tmp_path / "frozen", tmp_path / "cohort", cohort_id=COHORT_ID)
    return load_frozen_cohort(tmp_path / "cohort")


def _rows(db_path: Path) -> dict[tuple[str, str], tuple[str, str | None]]:
    with sqlite3.connect(db_path) as connection:
        return {
            (species, srr): (state, error)
            for species, srr, state, error in connection.execute("SELECT species, srr_id, state, error FROM samples")
        }


def test_prune_keeps_every_row_and_reports_outside_ones(tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR101", "SRR102", "SRR900", "SRR901"])
    db.set_state(SPECIES, "SRR900", "quantified")
    db.set_state(SPECIES, "SRR901", "failed", error=KALLISTO)
    db.record_exclusions(SPECIES, [{"srr_id": "SRR901", "reason_code": "permanent_drop"}])
    assert db.rows_outside(SPECIES, ["SRR101", "SRR102"]) == {"quantified": 1, "failed": 1}
    assert db.prune_species(SPECIES, ["SRR101", "SRR102"]) == 2
    assert db.get_total_counts() == {SPECIES: 4}
    assert db.get_state(SPECIES, "SRR900") == "quantified"
    assert db.get_excluded_srr_ids(SPECIES) == {"SRR901"}
    db.close()


def test_error_classes_cover_host_failures() -> None:
    assert classify_sample_error(EPERM) == "environment_write_denied"
    assert classify_sample_error("[Errno 13] Permission denied: '/Volumes/x/getfastq'") == "environment_access_denied"
    assert classify_sample_error("[Errno 30] Read-only file system: '/Volumes/x'") == "environment_read_only"
    assert classify_sample_error("[Errno 28] No space left on device") == "environment_no_space"
    assert classify_sample_error("[Errno 2] No such file or directory: 'amalgkit'") == "environment_missing_tool"
    assert classify_sample_error(TRANSFER) == "transfer_all_sources_failed"
    assert RETRYABLE_ENVIRONMENT_CLASSES < RETRYABLE_RESET_CLASSES
    assert "transfer_all_sources_failed" not in RETRYABLE_ENVIRONMENT_CLASSES


def test_migration_adds_rows_resets_environment_failures_and_keeps_everything_else(cohort, tmp_path: Path) -> None:
    db_path = tmp_path / "p.db"
    db = ProgressDB(db_path)
    db.init_species(SPECIES, ["SRR101", "SRR102", "SRR103", "SRR104"])
    db.init_species("other_species", ["SRR777"])
    db.init_species(SPECIES, ["SRR900"])  # outside the cohort
    db.set_state(SPECIES, "SRR101", "quantified")
    db.set_state(SPECIES, "SRR102", "failed", error=EPERM)
    db.set_state(SPECIES, "SRR103", "failed", error=TRANSFER)
    db.set_state(SPECIES, "SRR104", "failed", error=KALLISTO)  # the cohort's permanent_drop
    db.set_state(SPECIES, "SRR900", "failed", error=EPERM)  # environmental, but not a cohort task
    db.record_exclusions(SPECIES, [{"srr_id": "SRR104", "reason_code": "permanent_drop", "recorded_by": "audit"}])
    before = _rows(db_path)

    summary = db.migrate_to_cohort(cohort, backup_path=tmp_path / "backup.db")
    after = _rows(db_path)
    assert summary["inserted_task_rows"] == 1  # SRR105
    assert summary["rows_after"] == summary["rows_before"] + 1 == len(after)
    assert after[(SPECIES, "SRR105")] == ("pending", None)
    assert after[(SPECIES, "SRR101")] == ("quantified", None)
    state, error = after[(SPECIES, "SRR102")]
    assert state == "pending" and error.startswith(f"{COHORT_MIGRATION_ERROR_PREFIX}[environment_write_denied]: ")
    assert EPERM in error
    assert after[(SPECIES, "SRR103")] == before[(SPECIES, "SRR103")]  # transfer failures stay failed by default
    assert after[(SPECIES, "SRR104")] == before[(SPECIES, "SRR104")]
    assert after[(SPECIES, "SRR900")] == before[(SPECIES, "SRR900")]
    assert after[("other_species", "SRR777")] == before[("other_species", "SRR777")]
    assert summary["reset_by_class"] == {"environment_write_denied": 1}
    assert summary["failed_kept_by_class"] == {
        "environment_write_denied": 1,
        "quantification_failed": 1,
        "transfer_all_sources_failed": 1,
    }
    assert summary["outside_cohort_rows"] == {"failed": 1, "pending": 1}
    assert summary["excluded_rows_by_state"] == {"failed": 1}
    assert summary["cohort_exclusions_inserted"] == 0  # the audit's exclusion is kept as recorded
    assert db.get_exclusions(SPECIES)[0]["recorded_by"] == "audit"
    assert db.cohort_binding() == {"cohort_id": COHORT_ID, "cohort_definition_sha256": cohort.definition_sha256}
    with sqlite3.connect(tmp_path / "backup.db") as backup:
        assert backup.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == summary["rows_before"]

    again = db.migrate_to_cohort(cohort, backup_path=tmp_path / "backup.db")
    assert again["inserted_task_rows"] == 0 and again["reset_by_class"] == {}
    assert _rows(db_path) == after
    db.close()


def test_migration_records_cohort_and_extra_exclusions_and_optional_transfer_reset(cohort, tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR103"])
    db.set_state(SPECIES, "SRR103", "failed", error=TRANSFER)
    summary = db.migrate_to_cohort(
        cohort,
        reset_classes=RETRYABLE_RESET_CLASSES,
        extra_exclusions=[{"species": SPECIES, "srr_id": "SRR105", "reason_detail": "kallisto exit 1 twice"}],
    )
    assert summary["cohort_exclusions_inserted"] == 1 and summary["extra_exclusions_inserted"] == 1
    exclusions = {row["srr_id"]: row for row in db.get_exclusions(SPECIES)}
    assert exclusions["SRR104"]["recorded_by"] == f"cohort:{COHORT_ID}"
    assert exclusions["SRR105"]["reason_code"] == "permanent_drop"
    assert db.get_state(SPECIES, "SRR103") == "pending"
    assert db.get_state(SPECIES, "SRR104") is None  # an exclusion is not turned into a task row
    with pytest.raises(CohortError, match="not part of cohort"):
        db.migrate_to_cohort(cohort, extra_exclusions=[{"species": SPECIES, "srr_id": "SRR999"}])
    db.close()


def test_migration_refuses_a_store_bound_to_another_cohort(cohort, tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    db.migrate_to_cohort(cohort)
    with sqlite3.connect(tmp_path / "p.db") as connection:
        connection.execute("UPDATE campaign_meta SET value = 'other-cohort' WHERE key = 'cohort_id'")
    with pytest.raises(CohortError, match="bound to cohort other-cohort"):
        db.migrate_to_cohort(cohort)
    db.close()


def test_a_migration_interrupted_after_its_backup_can_be_rerun(cohort, tmp_path: Path) -> None:
    """The first attempt wrote its backup and died before committing: the rerun keeps it and backs up anew."""

    (tmp_path / "backup.db").write_bytes(b"backup of an interrupted attempt")
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR101"])
    assert db.cohort_binding() == {}

    summary = db.migrate_to_cohort(cohort, backup_path=tmp_path / "backup.db")

    assert (tmp_path / "backup.db").read_bytes() == b"backup of an interrupted attempt"  # never overwritten
    fresh = Path(summary["backup"])
    assert fresh.parent == tmp_path and fresh.name.startswith("backup.db.") and fresh != tmp_path / "backup.db"
    with sqlite3.connect(fresh) as backup:
        assert backup.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 1
    assert db.cohort_binding()["cohort_id"] == COHORT_ID
    assert not [path.name for path in tmp_path.iterdir() if path.name.endswith(".tmp")]
    # Once bound, a rerun keeps the existing backup and changes nothing.
    assert db.migrate_to_cohort(cohort, backup_path=tmp_path / "backup.db")["backup"].startswith("kept existing")
    db.close()


def test_adopt_quantified_is_one_transaction_and_skips_permanent_drops(tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    db.record_exclusions(SPECIES, [{"srr_id": "SRR104"}])
    result = db.adopt_quantified(
        SPECIES,
        [
            {"srr_id": "SRR101", "status": "current", "reason": "imported", "contract_id": "c1"},
            {"srr_id": "SRR104", "status": "current", "reason": "imported", "contract_id": "c2"},
        ],
    )
    assert result == {"adopted": ["SRR101"], "excluded": ["SRR104"]}
    assert db.get_state(SPECIES, "SRR101") == "quantified" and db.get_state(SPECIES, "SRR104") is None
    assert db.get_quantification_audit_counts(SPECIES) == {SPECIES: {"current": 1}}
    with pytest.raises(ValueError, match="non-reusable"):
        db.adopt_quantified(SPECIES, [{"srr_id": "SRR102", "status": "invalid", "reason": "x"}])
    assert db.get_state(SPECIES, "SRR102") is None
    db.close()


def test_reconcile_cohort_promotes_quarantines_and_requeues(tmp_path: Path) -> None:
    quant = tmp_path / "quant"
    for name in ("SRR101", "SRR102", "SRR103", "SRR106", "SRR900"):
        (quant / name).mkdir(parents=True)
    verdicts = {
        "SRR101": {"status": "current", "reason": "ok", "contract_id": "a"},
        "SRR102": {"status": "invalid", "reason": "configuration differs from the cohort configuration"},
        "SRR103": {"status": "version_drift_compatible", "reason": "drift", "contract_id": "c"},
        "SRR106": {"status": "invalid", "reason": "quantification checksum mismatch"},
    }
    seen: list[str] = []

    def classify(path: Path, run: str) -> dict[str, str]:
        seen.append(run)
        assert path == quant / run
        return verdicts[run]

    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR101", "SRR102", "SRR103", "SRR105", "SRR106"])
    db.set_state(SPECIES, "SRR102", "quantified")
    db.set_state(SPECIES, "SRR105", "quantified")  # no output directory any more
    summary = db.reconcile_cohort(
        SPECIES,
        quant,
        ["SRR101", "SRR102", "SRR103", "SRR105", "SRR106"],
        classify,
        requantify_statuses={"version_drift_compatible"},
        workers=2,
    )
    assert sorted(seen) == ["SRR101", "SRR102", "SRR103", "SRR106"]  # SRR900 is not a schedulable task
    assert summary == {
        "scanned": 4,
        "quantified": 1,
        "requeued": 1,
        "quarantined": ["SRR102"],
        "vanished": ["SRR105"],
    }
    assert db.get_state(SPECIES, "SRR101") == "quantified"
    assert db.get_state(SPECIES, "SRR102") == "quarantined"
    assert db.get_state(SPECIES, "SRR103") == "pending"
    assert db.get_state(SPECIES, "SRR105") == "pending"
    assert db.get_state(SPECIES, "SRR106") == "pending"  # never claimed; the producer re-quantifies it
    assert db.get_state(SPECIES, "SRR900") is None
    db.close()


# ---------------------------------------------------------------------------
# Real progress DB copy (skipped unless the copied fixture is present)
# ---------------------------------------------------------------------------

REAL_DB_ENV = "RNA_E2_FIXTURE_DIR"


def test_real_canonical_db_copy_keeps_every_row_when_pruned(tmp_path: Path) -> None:
    value = os.environ.get(REAL_DB_ENV)
    if not value:
        pytest.skip(f"{REAL_DB_ENV} is not set")
    source = Path(value) / "candb" / "pipeline_progress.db"
    if not source.is_file():
        pytest.skip(f"{REAL_DB_ENV} does not hold candb/pipeline_progress.db")
    if str(source.resolve()).startswith("/Volumes/"):
        pytest.fail("real fixtures must be copies, never the data volume itself")
    copy = tmp_path / "pipeline_progress.db"
    shutil.copyfile(source, copy)
    db = ProgressDB(copy)
    before = db.get_counts()
    totals = db.get_total_counts()
    quantified = {species: set(db.get_samples(species, "quantified")) for species in before}
    exclusions = db.get_exclusions()
    for species in before:
        keep = sorted(quantified[species])[: len(quantified[species]) // 2]
        kept_outside = db.prune_species(species, keep)
        assert kept_outside == totals[species] - len(keep)
    assert db.get_counts() == before
    assert {species: set(db.get_samples(species, "quantified")) for species in before} == quantified
    assert db.get_exclusions() == exclusions
    assert sum(counts.get("quantified", 0) for counts in before.values()) == 10438
    db.close()


def test_migrate_db_cli_backs_up_and_refuses_while_a_producer_holds_the_root(
    cohort, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from metainformant.rna.engine.cohort import main

    root = tmp_path / "data"
    root.mkdir()
    db_path = root / "pipeline_progress.db"
    with ProgressDB(db_path) as db:
        db.init_species(SPECIES, ["SRR102"])
        db.set_state(SPECIES, "SRR102", "failed", error=EPERM)
    (root / "extra.tsv").write_text("species\taccession\treason_detail\n" f"{SPECIES}\tSRR105\tkallisto exit 1\n")
    args = ["migrate-db", "--cohort", str(tmp_path / "cohort"), "--db", str(db_path)]
    # A live producer holds the campaign writer lock of the data root (a separate process, like the real one).
    holder = subprocess.Popen(
        [sys.executable, "-c", PRODUCER_LOCK_HOLDER, str(root)], stdout=subprocess.PIPE, text=True
    )
    try:
        assert holder.stdout is not None and holder.stdout.readline().strip() == "held"
        assert main(args) == 2
        assert "locked by another writer" in capsys.readouterr().err
    finally:
        holder.kill()
        holder.wait(timeout=30)
        if holder.stdout is not None:
            holder.stdout.close()
    with ProgressDB(db_path) as db:
        assert db.get_state(SPECIES, "SRR102") == "failed"
    assert main([*args, "--extra-exclusions", str(root / "extra.tsv")]) == 0
    printed = capsys.readouterr().out
    assert '"inserted_task_rows": 3' in printed
    assert (root / f"pipeline_progress.db.pre-{COHORT_ID}.bak").is_file()
    with ProgressDB(db_path) as db:
        assert db.get_state(SPECIES, "SRR102") == "pending"
        assert db.get_excluded_srr_ids(SPECIES, reason_code="permanent_drop") == {"SRR104", "SRR105"}
