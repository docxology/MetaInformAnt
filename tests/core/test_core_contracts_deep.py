"""Deep behavior-level contract tests for stable metainformant.core APIs.

REAL IMPLEMENTATION - all tests exercise the actual implementations with
small deterministic fixtures; no network, no test doubles.

Contracts covered:
- core.io.atomic: atomic write interruption semantics (no partial file ever
  visible, no temp residue, KeyboardInterrupt-safe).
- core.io.checksums: digest round-trips against hashlib reference vectors,
  sidecar write/verify round-trip, refusal of unsupported input.
- core.io.paths + core.data.validation: path containment refusals
  (traversal, symlink escape) and security-oriented validation.
- core.utils.optional_deps: optional-dependency reporting message format,
  once-only issuance, suppression round-trip.
- core.execution.parallel: worker-limit enforcement, order preservation,
  error propagation, rate-limit lower bound.
- core.execution.workflow: DAG illegal-transition refusal (cycles, unknown
  deps), failure-skip propagation, step state-transition invariants.
"""

from __future__ import annotations

import hashlib
import logging
import threading
import time
from pathlib import Path

import pytest

from metainformant.core.data.validation import validate_path_within, validate_range, validate_schema
from metainformant.core.execution import parallel
from metainformant.core.execution.workflow import BaseWorkflowOrchestrator, WorkflowStep
from metainformant.core.io.atomic import atomic_replace, atomic_write, safe_write_bytes, safe_write_text
from metainformant.core.io.checksums import (
    compute_checksums_batch,
    compute_md5,
    compute_sha256,
    verify_checksum,
    verify_checksum_file,
    write_checksum_file,
)
from metainformant.core.io.paths import is_safe_path, is_within, sanitize_filename
from metainformant.core.utils import optional_deps
from metainformant.core.utils.errors import ValidationError

# ---------------------------------------------------------------------------
# core.io.atomic -- atomic write interruption semantics
# ---------------------------------------------------------------------------


class TestAtomicWriteInterruption:
    """A crashed write must never expose a partial or placeholder file."""

    def test_exception_leaves_no_target_and_no_temp_residue(self, tmp_path: Path) -> None:
        target = tmp_path / "out.txt"
        with pytest.raises(RuntimeError, match="crash"):
            with atomic_write(target) as handle:
                handle.write("half-written payload")
                raise RuntimeError("crash")

        assert not target.exists()
        assert list(tmp_path.iterdir()) == [], "temp file must be cleaned up on failure"

    def test_keyboard_interrupt_is_also_cleaned_up(self, tmp_path: Path) -> None:
        # atomic_write must handle BaseException (KeyboardInterrupt), not just Exception.
        target = tmp_path / "out.bin"
        with pytest.raises(KeyboardInterrupt):
            with atomic_write(target, mode="wb") as handle:
                handle.write(b"\x00\x01")
                raise KeyboardInterrupt

        assert not target.exists()
        assert list(tmp_path.iterdir()) == []

    def test_no_partial_file_visible_during_write(self, tmp_path: Path) -> None:
        target = tmp_path / "results.txt"
        observed: dict[str, object] = {}
        with atomic_write(target) as handle:
            handle.write("complete output")
            # While the write is in flight the target must not exist yet...
            observed["target_exists"] = target.exists()
            # ...and the staging temp file must live in the same directory
            # (rename-across-filesystems would not be atomic).
            observed["dir_entries"] = sorted(p.name for p in tmp_path.iterdir())

        assert observed["target_exists"] is False
        entries = observed["dir_entries"]
        assert isinstance(entries, list) and len(entries) == 1
        assert str(entries[0]).startswith(".results.txt.") and str(entries[0]).endswith(".tmp")
        assert target.read_text() == "complete output"

    def test_atomic_replace_refuses_missing_source_and_keeps_dst(self, tmp_path: Path) -> None:
        dst = tmp_path / "dst.txt"
        dst.write_text("previous good content")
        with pytest.raises(FileNotFoundError):
            atomic_replace(tmp_path / "does-not-exist.txt", dst)
        assert dst.read_text() == "previous good content"

    def test_safe_write_fully_replaces_previous_content(self, tmp_path: Path) -> None:
        target = tmp_path / "seq.fa"
        safe_write_text(target, "ACGTACGTACGT")
        safe_write_text(target, "TTTT")
        assert target.read_text() == "TTTT"
        safe_write_bytes(target, b"\xca\xfe\xba\xbe")
        assert target.read_bytes() == b"\xca\xfe\xba\xbe"


# ---------------------------------------------------------------------------
# core.io.checksums -- digest round-trips and integrity verification
# ---------------------------------------------------------------------------

_KNOWN_PAYLOAD = b"hello world"
_KNOWN_MD5 = "5eb63bbbe01eeed093cb22bb8f5acdc3"
_KNOWN_SHA256 = "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"


class TestChecksumRoundTrips:
    """Digests must match hashlib reference vectors and round-trip through sidecars."""

    def test_known_reference_vectors(self, tmp_path: Path) -> None:
        target = tmp_path / "payload.bin"
        target.write_bytes(_KNOWN_PAYLOAD)
        assert compute_md5(target) == _KNOWN_MD5
        assert compute_sha256(target) == _KNOWN_SHA256

    def test_chunked_reading_matches_reference_across_chunk_boundaries(self, tmp_path: Path) -> None:
        # 20000 bytes spans several 1024-byte chunks; chunk size must not affect digest.
        payload = bytes(range(256)) * 79  # 20224 bytes, deterministic pattern
        target = tmp_path / "big.bin"
        target.write_bytes(payload)
        assert compute_sha256(target) == hashlib.sha256(payload).hexdigest()
        assert compute_md5(target, chunk_size=1) == hashlib.md5(payload).hexdigest()

    def test_verify_checksum_is_case_insensitive_and_detects_corruption(self, tmp_path: Path) -> None:
        target = tmp_path / "data.txt"
        target.write_bytes(_KNOWN_PAYLOAD)
        assert verify_checksum(target, _KNOWN_SHA256.upper()) is True
        wrong = "0" * 64
        assert verify_checksum(target, wrong) is False

    def test_missing_sidecar_raises(self, tmp_path: Path) -> None:
        target = tmp_path / "orphan.txt"
        target.write_text("no sidecar here")
        with pytest.raises(FileNotFoundError, match="sidecar"):
            verify_checksum_file(target)

    def test_sidecar_round_trip_and_md5_fallback(self, tmp_path: Path) -> None:
        target = tmp_path / "sample.fa"
        target.write_text("ACGTACGT")

        sidecar = write_checksum_file(target)
        assert sidecar == tmp_path / "sample.fa.sha256"
        assert sidecar.read_text().split()[0] == compute_sha256(target)
        assert sidecar.read_text().split()[1] == "sample.fa"
        assert verify_checksum_file(target) is True

        # Tamper: content no longer matches the sidecar digest.
        target.write_text("ACGTACGA")
        assert verify_checksum_file(target) is False

        # md5-only sidecar is picked up by the fallback path.
        target.write_text("ACGTACGT")
        (tmp_path / "sample.fa.sha256").unlink()
        assert write_checksum_file(target, algorithm="md5").name == "sample.fa.md5"
        (tmp_path / "sample.fa.md5").unlink()
        (tmp_path / "sample.fa.md5").write_text(f"{hashlib.md5(b'ACGTACGT').hexdigest()}  sample.fa\n")
        assert verify_checksum_file(target) is True

    def test_unsupported_algorithm_refused(self, tmp_path: Path) -> None:
        target = tmp_path / "x.txt"
        target.write_text("x")
        # verify_checksum rejects algorithms outside the supported set.
        with pytest.raises(ValueError, match="Unsupported algorithm"):
            verify_checksum(target, "0" * 64, algorithm="crc32")

    def test_missing_file_refused(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            compute_sha256(tmp_path / "ghost.bin")

    def test_batch_skips_missing_files_without_failing(self, tmp_path: Path) -> None:
        a = tmp_path / "a.txt"
        b = tmp_path / "b.txt"
        a.write_bytes(b"alpha")
        b.write_bytes(b"beta")
        results = compute_checksums_batch([a, tmp_path / "missing.txt", b])
        assert set(results) == {str(a), str(b)}
        assert results[str(a)] == hashlib.sha256(b"alpha").hexdigest()
        assert results[str(b)] == hashlib.sha256(b"beta").hexdigest()


# ---------------------------------------------------------------------------
# core.io.paths + core.data.validation -- containment refusals
# ---------------------------------------------------------------------------


class TestPathContainment:
    """Traversal and symlink-escape must be refused at the containment boundary."""

    def test_is_within_accepts_nested_and_refuses_escape(self, tmp_path: Path) -> None:
        vault = tmp_path / "vault"
        vault.mkdir()
        (vault / "sub").mkdir()

        assert is_within(vault / "sub" / "file.txt", vault) is True
        assert is_within(vault, vault) is True  # the parent itself is within
        assert is_within(vault / ".." / "stolen.txt", vault) is False
        assert is_within(tmp_path / "elsewhere.txt", vault) is False

    def test_validate_path_within_returns_resolved_path_when_inside(self, tmp_path: Path) -> None:
        vault = tmp_path / "vault"
        vault.mkdir()
        resolved = validate_path_within(vault, vault / "sub" / "f.txt", name="user_path")
        assert resolved == (vault / "sub" / "f.txt").resolve()

    def test_validate_path_within_refuses_traversal(self, tmp_path: Path) -> None:
        vault = tmp_path / "vault"
        vault.mkdir()
        with pytest.raises(ValidationError, match="must be within"):
            validate_path_within(vault, vault / ".." / "escape.txt", name="user_path")

    def test_validate_path_within_refuses_symlink_escape(self, tmp_path: Path) -> None:
        vault = tmp_path / "vault"
        vault.mkdir()
        secret = tmp_path / "secret.txt"  # outside vault, inside tmp
        secret.write_text("s3cret")
        (vault / "escape").symlink_to(secret)
        with pytest.raises(ValidationError, match="must be within"):
            validate_path_within(vault, vault / "escape", name="user_path")

    def test_is_safe_path_blocks_traversal_and_metacharacters(self) -> None:
        assert is_safe_path("../etc/passwd") is False
        assert is_safe_path("/etc/passwd") is False
        assert is_safe_path("/root/.ssh/id_rsa") is False
        assert is_safe_path("out; rm -rf /") is False
        assert is_safe_path("a|b") is False
        assert is_safe_path("a&b") is False
        assert is_safe_path("a$b") is False
        assert is_safe_path("output/results.txt") is True
        assert is_safe_path("data/sample_001.vcf.gz") is True

    def test_sanitize_filename_neutralizes_dangerous_characters(self) -> None:
        assert sanitize_filename('report<>:"/\\|?*.txt') == "report_________.txt"
        assert sanitize_filename("\x00bad\x1fname") == "badname"
        assert sanitize_filename("  .hidden.  ") == "hidden"
        assert sanitize_filename("...") == "untitled"  # only dots -> must not be empty
        # Path separators become underscores, then leading/trailing dots are stripped.
        assert sanitize_filename("../..") == "_"


class TestValidationContracts:
    def test_validate_range_bounds_are_inclusive(self) -> None:
        validate_range(0.0, min_val=0.0, max_val=1.0, name="probability")
        validate_range(1.0, min_val=0.0, max_val=1.0, name="probability")
        with pytest.raises(ValidationError, match="probability must be >="):
            validate_range(-0.1, min_val=0.0, max_val=1.0, name="probability")
        with pytest.raises(ValidationError, match="probability must be <="):
            validate_range(1.5, min_val=0.0, max_val=1.0, name="probability")

    def test_validate_schema_refuses_missing_field_and_wrong_type(self) -> None:
        schema = {"sample_id": str, "depth": int}
        validate_schema({"sample_id": "S1", "depth": 30}, schema)
        with pytest.raises(ValidationError, match="missing required field: depth"):
            validate_schema({"sample_id": "S1"}, schema)
        with pytest.raises(ValidationError, match=r"data\.depth"):
            validate_schema({"sample_id": "S1", "depth": "30"}, schema)


# ---------------------------------------------------------------------------
# core.utils.optional_deps -- optional-dependency reporting contract
# ---------------------------------------------------------------------------


@pytest.fixture()
def clean_optional_dep_state():
    optional_deps.reset_warning_state()
    yield
    optional_deps.reset_warning_state()


class TestOptionalDependencyReporting:
    def test_message_format_and_once_only_issuance(self, caplog, clean_optional_dep_state) -> None:
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("seaborn", "enhanced plots", "basic plots used")
            optional_deps.warn_optional_dependency("seaborn", "enhanced plots", "basic plots used")

        records = [r for r in caplog.records if r.name == "metainformant.core.optional_deps"]
        assert len(records) == 1, "same (module, functionality) pair must warn exactly once"
        assert records[0].getMessage() == "seaborn not available, enhanced plots basic plots used"

    def test_distinct_functionality_warns_independently(self, caplog, clean_optional_dep_state) -> None:
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("anndata", "h5ad export", "export disabled")
            optional_deps.warn_optional_dependency("anndata", "cluster plots", "plots skipped")
        messages = [r.getMessage() for r in caplog.records if r.name == "metainformant.core.optional_deps"]
        assert len(messages) == 2

    def test_suppression_blocks_warnings_until_reenabled(self, caplog, clean_optional_dep_state) -> None:
        optional_deps.suppress_optional_warnings()
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("scanpy", "trajectory analysis", "analysis disabled")
        assert caplog.records == []
        assert optional_deps.get_warning_state()["suppress_warnings"] is True

        # reset_warning_state clears both suppression and the issued set.
        optional_deps.reset_warning_state()
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("scanpy", "trajectory analysis", "analysis disabled")
        assert any("scanpy not available" in r.getMessage() for r in caplog.records)


# ---------------------------------------------------------------------------
# core.execution.parallel -- worker limits and execution contracts
# ---------------------------------------------------------------------------


class TestParallelResourceLimits:
    def test_worker_recommendation_respects_cap_and_never_drops_below_one(self) -> None:
        capped = parallel.resource_aware_workers(task_type="io", max_cap=3)
        assert 1 <= capped <= 3
        assert parallel.resource_aware_workers(task_type="cpu", max_cap=1) == 1
        # Even a zero cap must not produce zero workers.
        assert parallel.resource_aware_workers(task_type="io", max_cap=0) == 1

    def test_thread_map_never_exceeds_max_workers(self) -> None:
        lock = threading.Lock()
        state = {"active": 0, "peak": 0}

        def work(x: int) -> int:
            with lock:
                state["active"] += 1
                state["peak"] = max(state["peak"], state["active"])
            time.sleep(0.005)
            with lock:
                state["active"] -= 1
            return x * 2

        results = parallel.thread_map(work, list(range(24)), max_workers=2)
        assert results == [x * 2 for x in range(24)]
        assert state["peak"] <= 2

    def test_thread_map_preserves_input_order_regardless_of_completion_order(self) -> None:
        items = list(range(8))

        def slow_early(x: int) -> int:
            # Later items finish first; ordering must still follow input index.
            time.sleep((len(items) - 1 - x) * 0.004)
            return x * 10

        assert parallel.thread_map(slow_early, items, max_workers=4) == [x * 10 for x in items]

    def test_thread_map_propagates_task_errors(self) -> None:
        def maybe_fail(x: int) -> int:
            if x == 3:
                raise ValueError(f"bad item {x}")
            return x

        with pytest.raises(ValueError, match="bad item 3"):
            parallel.thread_map(maybe_fail, list(range(6)), max_workers=2)

    def test_parallel_batch_flattens_results_in_input_order(self) -> None:
        items = list(range(10))
        out = parallel.parallel_batch(lambda batch: [x * x for x in batch], items, batch_size=3, max_workers=2)
        assert out == [x * x for x in items]

    def test_gather_results_separates_successes_from_errors(self) -> None:
        from concurrent.futures import ThreadPoolExecutor

        def ok() -> str:
            return "fine"

        def bad() -> str:
            raise RuntimeError("boom")

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(ok), pool.submit(bad)]
            successes, errors = parallel.gather_results(futures)

        assert successes == ["fine"]
        assert len(errors) == 1 and isinstance(errors[0], RuntimeError)

    def test_rate_limited_map_enforces_minimum_spacing_and_order(self) -> None:
        # 6 items at 50/s => nominal minimum ~0.10s; a lower bound is a stable assertion.
        start = time.monotonic()
        results = parallel.rate_limited_map(lambda x: x + 1, list(range(6)), max_per_second=50.0, max_workers=6)
        elapsed = time.monotonic() - start
        assert results == [1, 2, 3, 4, 5, 6]
        assert elapsed >= 0.08


# ---------------------------------------------------------------------------
# core.execution.workflow -- state-transition invariants
# ---------------------------------------------------------------------------


class TestWorkflowStateTransitions:
    def test_cycle_is_refused_and_no_step_executes(self) -> None:
        called: list[str] = []
        orch = BaseWorkflowOrchestrator({})
        orch.add_step("a", lambda: called.append("a") or "A", depends_on=["b"])
        orch.add_step("b", lambda: called.append("b") or "B", depends_on=["a"])

        result = orch.run_workflow()

        assert result["success"] is False
        assert result["execution_order"] == []
        assert any("Cycle" in error for error in result["errors"])
        assert called == [], "no step may execute when the DAG is cyclic"

    def test_unknown_dependency_is_refused(self) -> None:
        called: list[str] = []
        orch = BaseWorkflowOrchestrator({})
        orch.add_step("solo", lambda: called.append("solo") or "S", depends_on=["ghost"])

        result = orch.run_workflow()

        assert result["success"] is False
        assert any("unknown step 'ghost'" in error for error in result["errors"])
        assert called == []

    def test_failed_upstream_forces_downstream_skip_without_execution(self) -> None:
        child_called: list[str] = []

        def bad() -> str:
            raise RuntimeError("boom")

        def child() -> str:
            child_called.append("child")
            return "child-result"

        orch = BaseWorkflowOrchestrator({})
        orch.add_step("bad", bad)
        orch.add_step("child", child, depends_on=["bad"])
        orch.add_step("grandchild", lambda: "gc", depends_on=["child"])

        result = orch.run_workflow()

        assert result["success"] is False
        assert orch.get_step_status() == {"bad": "failed", "child": "skipped", "grandchild": "skipped"}
        assert child_called == [], "skipped steps must never run"
        assert any("skipped due to failed dependency" in error for error in result["errors"])
        assert any("boom" in error for error in result["errors"])
        # The failed step itself has no result recorded.
        assert "bad" not in result["results"]
        assert "child" not in result["results"]

    def test_step_lifecycle_transitions_and_reset(self) -> None:
        step = WorkflowStep("s", lambda: 42, {})
        assert step.status == "pending"

        assert step.execute() == 42
        assert step.status == "completed"
        assert step.result == 42
        assert step.duration() >= 0.0

        def bad() -> None:
            raise RuntimeError("exploded")

        failing = WorkflowStep("f", bad, {})
        with pytest.raises(RuntimeError, match="exploded"):
            failing.execute()
        assert failing.status == "failed"
        assert failing.error == "exploded"

        # reset() returns the step to a clean re-executable state.
        failing.reset()
        assert failing.status == "pending"
        assert failing.result is None
        assert failing.error is None

    def test_dependency_results_flow_into_downstream_kwargs(self) -> None:
        def consume(*, produce: dict) -> int:
            return produce["value"] * 2

        orch = BaseWorkflowOrchestrator({})
        orch.add_step("produce", lambda: {"value": 3})
        orch.add_step("consume", consume, depends_on=["produce"])

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"] == {"produce": {"value": 3}, "consume": 6}

    def test_successful_run_reports_all_steps_completed_in_topological_order(self) -> None:
        orch = BaseWorkflowOrchestrator({})
        # Dependent steps receive upstream results as keyword arguments named
        # after the dependency, so their callables must accept those kwargs.
        orch.add_step("download", lambda: "fastq")
        orch.add_step("qc", lambda download: f"clean:{download}", depends_on=["download"])
        orch.add_step("quantify", lambda qc: 7, depends_on=["qc"])

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["errors"] == []
        assert orch.get_step_status() == {"download": "completed", "qc": "completed", "quantify": "completed"}
        assert result["results"]["qc"] == "clean:fastq"
        # Topological constraint: dependencies precede dependents.
        order = result["execution_order"]
        assert order.index("download") < order.index("qc") < order.index("quantify")
