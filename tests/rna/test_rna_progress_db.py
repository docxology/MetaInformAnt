"""Tests for SQLite-backed progress tracking (progress_db).

All tests use real SQLite databases on disk following the real-implementation policy.
"""

from __future__ import annotations

import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.progress_db import ProgressDB


class TestProgressDB:
    """Test ProgressDB core operations."""

    def test_create_db(self, tmp_path: Path):
        """Database file is created on initialization."""
        db_path = tmp_path / "test.db"
        db = ProgressDB(db_path)
        assert db_path.exists()
        db.close()

    def test_init_species(self, tmp_path: Path):
        """init_species inserts samples as pending."""
        db = ProgressDB(tmp_path / "test.db")
        inserted = db.init_species("test_ant", ["SRR001", "SRR002", "SRR003"])
        assert inserted == 3

        counts = db.get_counts("test_ant")
        assert counts["test_ant"]["pending"] == 3
        db.close()

    def test_init_species_idempotent(self, tmp_path: Path):
        """Re-initializing a species does not duplicate or reset existing samples."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("test_ant", ["SRR001", "SRR002"])
        db.set_state("test_ant", "SRR001", "quantified")

        # Re-init with overlapping + new IDs
        inserted = db.init_species("test_ant", ["SRR001", "SRR002", "SRR003"])
        assert inserted == 1  # Only SRR003 is new

        # SRR001 should still be quantified, not reset to pending
        assert db.get_state("test_ant", "SRR001") == "quantified"
        assert db.get_state("test_ant", "SRR003") == "pending"
        db.close()

    def test_prune_species_removes_only_stale_progress_rows(self, tmp_path: Path):
        """A refreshed metadata cohort cannot inherit obsolete pending rows."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("test_ant", ["SRR001", "SRR002", "SRR003"])
        db.set_state("test_ant", "SRR001", "quantified")

        assert db.prune_species("test_ant", ["SRR001", "SRR003"]) == 1
        assert db.get_state("test_ant", "SRR001") == "quantified"
        assert db.get_state("test_ant", "SRR002") is None
        assert db.get_state("test_ant", "SRR003") == "pending"
        db.close()

    def test_state_transitions(self, tmp_path: Path):
        """Full lifecycle: pending → downloading → downloaded → quantifying → quantified."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR100"])

        assert db.get_state("ant", "SRR100") == "pending"

        db.set_state("ant", "SRR100", "downloading")
        assert db.get_state("ant", "SRR100") == "downloading"

        db.set_state("ant", "SRR100", "downloaded")
        assert db.get_state("ant", "SRR100") == "downloaded"

        db.set_state("ant", "SRR100", "quantifying")
        assert db.get_state("ant", "SRR100") == "quantifying"

        db.set_state("ant", "SRR100", "quantified")
        assert db.get_state("ant", "SRR100") == "quantified"
        db.close()

    def test_set_state_with_error(self, tmp_path: Path):
        """Failed state includes error message."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR100"])

        db.set_state("ant", "SRR100", "failed", error="Download timeout")
        assert db.get_state("ant", "SRR100") == "failed"

        failed = db.get_failed("ant")
        assert len(failed) == 1
        assert failed[0]["error"] == "Download timeout"
        db.close()

    def test_invalid_state_raises(self, tmp_path: Path):
        """Setting an invalid state raises ValueError."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR100"])

        with pytest.raises(ValueError, match="Invalid state"):
            db.set_state("ant", "SRR100", "bogus_state")
        db.close()

    def test_get_state_unknown_sample(self, tmp_path: Path):
        """get_state returns None for an unknown sample."""
        db = ProgressDB(tmp_path / "test.db")
        assert db.get_state("nonexistent", "SRR999") is None
        db.close()

    def test_get_counts_multi_species(self, tmp_path: Path):
        """get_counts returns counts grouped by species and state."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant_a", ["SRR1", "SRR2", "SRR3"])
        db.init_species("ant_b", ["SRR4", "SRR5"])

        db.set_state("ant_a", "SRR1", "quantified")
        db.set_state("ant_a", "SRR2", "downloading")
        db.set_state("ant_b", "SRR4", "failed", error="oops")

        counts = db.get_counts()
        assert counts["ant_a"]["quantified"] == 1
        assert counts["ant_a"]["downloading"] == 1
        assert counts["ant_a"]["pending"] == 1
        assert counts["ant_b"]["failed"] == 1
        assert counts["ant_b"]["pending"] == 1
        db.close()

    def test_get_total_counts(self, tmp_path: Path):
        """get_total_counts returns total per species."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant_a", ["SRR1", "SRR2", "SRR3"])
        db.init_species("ant_b", ["SRR4", "SRR5"])

        totals = db.get_total_counts()
        assert totals["ant_a"] == 3
        assert totals["ant_b"] == 2
        db.close()

    def test_get_samples_by_state(self, tmp_path: Path):
        """get_samples returns IDs matching a specific state."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1", "SRR2", "SRR3"])
        db.set_state("ant", "SRR1", "quantified")
        db.set_state("ant", "SRR2", "quantified")

        quant = db.get_samples("ant", "quantified")
        assert sorted(quant) == ["SRR1", "SRR2"]

        pending = db.get_samples("ant", "pending")
        assert pending == ["SRR3"]
        db.close()

    def test_bulk_set_state(self, tmp_path: Path):
        """bulk_set_state updates multiple samples at once."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1", "SRR2", "SRR3"])

        db.bulk_set_state("ant", ["SRR1", "SRR2"], "quantified")
        assert db.get_state("ant", "SRR1") == "quantified"
        assert db.get_state("ant", "SRR2") == "quantified"
        assert db.get_state("ant", "SRR3") == "pending"
        db.close()

    def test_get_failed(self, tmp_path: Path):
        """get_failed returns all failed samples with errors."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1", "SRR2"])
        db.set_state("ant", "SRR1", "failed", error="network timeout")
        db.set_state("ant", "SRR2", "failed", error="md5 mismatch")

        failed = db.get_failed()
        assert len(failed) == 2
        errors = {f["srr_id"]: f["error"] for f in failed}
        assert errors["SRR1"] == "network timeout"
        assert errors["SRR2"] == "md5 mismatch"
        db.close()


class TestProgressDBReconciliation:
    """Test filesystem reconciliation."""

    def test_reconcile_does_not_promote_legacy_abundance_files(self, tmp_path: Path):
        """Readable outputs without current provenance remain unverified."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1", "SRR2", "SRR3"])

        # Simulate quant output on disk
        quant_dir = tmp_path / "quant"
        (quant_dir / "SRR1").mkdir(parents=True)
        (quant_dir / "SRR1" / "SRR1_abundance.tsv").write_text("header\ndata\n")
        (quant_dir / "SRR2").mkdir(parents=True)
        (quant_dir / "SRR2" / "quant.sf").write_text("header\ndata\n")

        reconciled = db.reconcile("ant", quant_dir)
        assert reconciled == 0

        assert db.get_state("ant", "SRR1") == "pending"
        assert db.get_state("ant", "SRR2") == "pending"
        assert db.get_state("ant", "SRR3") == "pending"
        db.close()

    def test_reconcile_defaults_to_contract_only_not_full_digest(self, tmp_path: Path):
        """Resume reconciliation must not re-hash every quantification payload.

        Regression: the default resume scan re-digested the full quantified
        corpus inside Phase-1 discovery, stalling a large campaign's task
        submission for 15+ minutes per restart (post-discovery deadlock).
        The default scan is contract-only; ``verify_hashes=True`` still
        performs the full content audit and rejects tampered payloads.
        """

        from metainformant.rna.engine.provenance import write_quant_provenance

        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1"])

        quant_dir = tmp_path / "quant"
        sample_dir = quant_dir / "SRR1"
        sample_dir.mkdir(parents=True)
        abundance = sample_dir / "SRR1_abundance.tsv"
        abundance.write_text("target_id\ttpm\nTX1\t1\n", encoding="utf-8")
        config = tmp_path / "config.yaml"
        config.write_text("species_list: [Test_species]\n", encoding="utf-8")
        write_quant_provenance(
            sample_dir,
            species="test_species",
            run_accession="SRR1",
            config_path=config,
            command=["amalgkit", "quant", "--threads", "1"],
            quantification_file=abundance,
        )
        # Tamper with the payload after the sidecar recorded its digest.
        abundance.write_text("target_id\ttpm\nTX1\t999\n", encoding="utf-8")

        # Default resume scan: contract-verified, content not re-digested.
        reconciled = db.reconcile("ant", quant_dir)
        assert reconciled == 1
        assert db.get_state("ant", "SRR1") == "quantified"

        # Opt-in deep audit: content digest mismatch rejects the sample.
        db.set_state("ant", "SRR1", "pending")
        assert db.reconcile("ant", quant_dir, verify_hashes=True) == 0
        assert db.get_state("ant", "SRR1") == "pending"
        db.close()

    def test_reconcile_nonexistent_dir(self, tmp_path: Path):
        """reconcile returns 0 for nonexistent directory."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1"])
        reconciled = db.reconcile("ant", tmp_path / "nonexistent")
        assert reconciled == 0
        db.close()


class TestProgressDBConcurrency:
    """Test thread safety."""

    def test_concurrent_writes(self, tmp_path: Path):
        """Multiple threads can write to the DB without errors."""
        db = ProgressDB(tmp_path / "test.db")
        srr_ids = [f"SRR{i:04d}" for i in range(100)]
        db.init_species("ant", srr_ids)

        errors = []

        def worker(start: int, end: int):
            try:
                for i in range(start, end):
                    db.set_state("ant", f"SRR{i:04d}", "downloading")
                    db.set_state("ant", f"SRR{i:04d}", "downloaded")
                    db.set_state("ant", f"SRR{i:04d}", "quantified")
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=worker, args=(i * 25, (i + 1) * 25)) for i in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors, f"Thread errors: {errors}"

        counts = db.get_counts("ant")
        assert counts["ant"]["quantified"] == 100
        db.close()


class TestProgressDBPersistence:
    """Test that state survives DB close/reopen."""

    def test_state_persists_across_sessions(self, tmp_path: Path):
        """State is persisted on disk and reloaded correctly."""
        db_path = tmp_path / "test.db"

        # Session 1: create and modify
        db1 = ProgressDB(db_path)
        db1.init_species("ant", ["SRR1", "SRR2", "SRR3"])
        db1.set_state("ant", "SRR1", "quantified")
        db1.set_state("ant", "SRR2", "failed", error="timeout")
        db1.close()

        # Session 2: reopen and verify
        db2 = ProgressDB(db_path)
        assert db2.get_state("ant", "SRR1") == "quantified"
        assert db2.get_state("ant", "SRR2") == "failed"
        assert db2.get_state("ant", "SRR3") == "pending"

        failed = db2.get_failed()
        assert len(failed) == 1
        assert failed[0]["error"] == "timeout"
        db2.close()


class TestProgressDBRepr:
    """Test repr and context manager."""

    def test_repr(self, tmp_path: Path):
        """__repr__ shows useful summary."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1", "SRR2"])
        r = repr(db)
        assert "1 species" in r
        assert "2 samples" in r
        db.close()

    def test_context_manager(self, tmp_path: Path):
        """ProgressDB works as a context manager."""
        db_path = tmp_path / "test.db"
        with ProgressDB(db_path) as db:
            db.init_species("ant", ["SRR1"])
            assert db.get_state("ant", "SRR1") == "pending"
        # After exit, connection should be closed
        # Reopening should still work
        with ProgressDB(db_path) as db2:
            assert db2.get_state("ant", "SRR1") == "pending"


class TestSampleExclusions:
    """Permanent-drop and re-download exclusion recording (real SQLite)."""

    def test_record_and_read_exclusions(self, tmp_path: Path):
        """Exclusions round-trip with reason codes and details."""
        db = ProgressDB(tmp_path / "test.db")
        recorded = db.record_exclusions(
            "apis_mellifera",
            [
                {"srr_id": "SRR1036397", "reason_code": "permanent_drop", "reason_detail": "16S amplicon"},
                {"srr_id": "SRR7412516", "reason_detail": "27F primer at read start"},
            ],
        )
        assert recorded == 2
        rows = db.get_exclusions("apis_mellifera")
        assert {r["srr_id"] for r in rows} == {"SRR1036397", "SRR7412516"}
        by_id = {r["srr_id"]: r for r in rows}
        assert by_id["SRR1036397"]["reason_code"] == "permanent_drop"
        assert by_id["SRR7412516"]["reason_code"] == "permanent_drop"  # default
        assert db.get_exclusions()  # unfiltered query returns the rows too
        db.close()

    def test_excluded_srr_ids_filters_by_reason_code(self, tmp_path: Path):
        """permanent_drop blocks eligibility; re_download stays eligible."""
        db = ProgressDB(tmp_path / "test.db")
        db.record_exclusions(
            "apis_mellifera",
            [
                {"srr_id": "SRR1", "reason_code": "permanent_drop"},
                {"srr_id": "SRR2", "reason_code": "re_download"},
            ],
        )
        assert db.get_excluded_srr_ids("apis_mellifera", reason_code="permanent_drop") == {"SRR1"}
        assert db.get_excluded_srr_ids("apis_mellifera") == {"SRR1", "SRR2"}
        db.close()

    def test_record_exclusions_upserts(self, tmp_path: Path):
        """Re-recording the same accession updates reason fields, not duplicates."""
        db = ProgressDB(tmp_path / "test.db")
        db.record_exclusions("ant", [{"srr_id": "SRR9", "reason_code": "re_download"}])
        db.record_exclusions("ant", [{"srr_id": "SRR9", "reason_code": "permanent_drop", "reason_detail": "amplicon"}])
        rows = db.get_exclusions("ant")
        assert len(rows) == 1
        assert rows[0]["reason_code"] == "permanent_drop"
        db.close()

    def test_invalid_reason_code_rejected(self, tmp_path: Path):
        """Unknown reason codes fail loudly instead of recording garbage."""
        db = ProgressDB(tmp_path / "test.db")
        with pytest.raises(ValueError):
            db.record_exclusions("ant", [{"srr_id": "SRR1", "reason_code": "bogus"}])
        assert db.get_exclusions("ant") == []
        db.close()

    def test_remove_exclusion(self, tmp_path: Path):
        """Removed exclusions disappear; removing an absent row reports False."""
        db = ProgressDB(tmp_path / "test.db")
        db.record_exclusions("ant", [{"srr_id": "SRR1"}])
        assert db.remove_exclusion("ant", "SRR1") is True
        assert db.remove_exclusion("ant", "SRR1") is False
        assert db.get_exclusions("ant") == []
        db.close()

    def test_exclusions_independent_from_sample_states(self, tmp_path: Path):
        """Recording an exclusion never mutates the sample state machine."""
        db = ProgressDB(tmp_path / "test.db")
        db.init_species("ant", ["SRR1"])
        db.set_state("ant", "SRR1", "failed", error="kallisto quant failed")
        db.record_exclusions("ant", [{"srr_id": "SRR1", "reason_code": "permanent_drop"}])
        assert db.get_state("ant", "SRR1") == "failed"
        db.close()
