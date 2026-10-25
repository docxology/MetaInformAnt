from __future__ import annotations

import concurrent.futures
import gzip
import json
import threading
import time
from contextlib import ExitStack
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from metainformant.rna.retrieval.ena_staging import FirstPassDeferred, StagingDeadline, StagingError
from metainformant.rna.engine.streaming_orchestrator import _raw_validation_marker
from metainformant.rna.retrieval import ena_transport
from metainformant.rna.retrieval.ena_transport import ena_request_slot
from tests.rna.test_ena_staging import transfer_site
from tests.rna.test_ena_staging_boundaries import build_owner


@pytest.mark.parametrize("slots", [0, 1])
def test_acquisition_admission_prevents_file_queue_consuming_transfer_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, slots: int
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", str(slots))
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "5")
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    monkeypatch.setenv("AMALGKIT_ENA_API_CONNECTIONS", "2")
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    first_file_started = threading.Event()
    requests: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            if self.path.startswith("/api?"):
                accession = parse_qs(urlsplit(self.path).query)["accession"][0]
                body = (f"run_accession\tfastq_ftp\n{accession}\t"
                        f"http://127.0.0.1:{self.server.server_port}/{accession}.fastq.gz\n").encode()
            else:
                requests.append(self.path)
                first_file_started.set()
                time.sleep(3)
                body = payload
            try:
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                # The disabled-admission control deliberately expires native curl.
                return

        def log_message(self, format: str, *args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"http://127.0.0.1:{server.server_port}/api")
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=5)
    calls: dict[str, str] = {}

    def acquire(accession: str) -> bool:
        try:
            return owner.download_fastq(accession, tmp_path / "data/getfastq", False) is True
        except (FirstPassDeferred, StagingDeadline) as exc:
            calls[accession] = str(exc)
            return False

    try:
        started = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            first = pool.submit(acquire, "SRR101")
            assert first_file_started.wait(4)
            second = pool.submit(acquire, "SRR102")
            outcomes = [first.result(6), second.result(6)]
        assert outcomes == ([True, True] if slots else [True, False]), (outcomes, calls)
        assert time.monotonic() - started < 12
        assert requests == ["/SRR101.fastq.gz", "/SRR102.fastq.gz"]
        for accession, success in zip(("SRR101", "SRR102"), outcomes, strict=True):
            sample = tmp_path / "data/getfastq" / accession
            if success:
                assert (sample / f"{accession}.fastq.gz").read_bytes() == payload
                proof = json.loads(_raw_validation_marker(sample, accession).read_text())
                assert proof["layout_validation"].startswith("full_fastq_record_scan")
    finally:
        owner.db.close()
        server.shutdown()
        server.server_close()
        thread.join(5)
        assert not thread.is_alive()


@pytest.mark.parametrize("slots", [0, 1])
def test_one_second_service_budget_excludes_admission_but_not_file_queue(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, slots: int,
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", str(slots))
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "5")
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    monkeypatch.setenv("AMALGKIT_ENA_API_CONNECTIONS", "2")
    caplog.set_level("INFO", logger="metainformant.rna.engine.streaming_orchestrator")
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=1)
    sample = tmp_path / "data/getfastq/SRR123"
    sample.mkdir(parents=True)
    part = sample / "SRR123.fastq.gz.part"
    if not slots:
        part.write_bytes(b"retained")
    before = (part.read_bytes(), part.stat().st_mtime_ns) if part.exists() else None
    entered = threading.Event()
    record: dict = {}

    def acquire() -> bool:
        record["started"] = time.monotonic()
        entered.set()
        try:
            return owner.download_fastq("SRR123", sample.parent, False) is True
        except FirstPassDeferred as exc:
            record["reason"] = str(exc)
            return False
        finally:
            record["elapsed"] = time.monotonic() - record["started"]

    try:
        with transfer_site({"SRR123.fastq.gz": payload}) as site:
            monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", site.base + "/api")
            semaphore = owner._acquisition_semaphore
            with ExitStack() as ownership:
                if semaphore is not None:
                    assert semaphore.acquire(blocking=False)
                    ownership.callback(semaphore.release)
                with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                    ownership.enter_context(ena_request_slot(site.base + "/SRR123.fastq.gz", kind="file", deadline_at=time.monotonic() + 8))
                    try:
                        future = pool.submit(acquire)
                        assert entered.wait(2)
                        if not slots:
                            key = ("file", "127.0.0.1", urlsplit(site.base).port)
                            until = time.monotonic() + 0.9
                            with ena_transport._condition:
                                while not ena_transport._gates[key].queue and time.monotonic() < until:
                                    ena_transport._condition.wait(timeout=0.01)
                                assert ena_transport._gates[key].queue, site.requests
                        time.sleep(1.5)
                        assert not any(path == "/SRR123.fastq.gz" for _, path, _ in site.requests)
                        if slots:
                            assert not site.requests
                    finally:
                        ownership.close()
                    result = future.result(4)
            if slots:
                admitted = [entry for entry in caplog.records if entry.msg.startswith("Acquisition admitted for")]
                assert len(admitted) == 1
                queued = admitted[0].args[1]
                assert queued > 1.1, queued
                service_elapsed = record["elapsed"] - queued
                if result:
                    assert 0 < service_elapsed < 1, record
                    assert (sample / "SRR123.fastq.gz").read_bytes() == payload
                    proof = json.loads(_raw_validation_marker(sample, "SRR123").read_text())
                    assert proof["layout_validation"].startswith("full_fastq_record_scan")
                else:
                    # Native startup and HTTP processing consume the same one-second
                    # service clock. Expiry is correct only after a fresh full budget;
                    # an admission-clock regression would fail immediately instead.
                    assert service_elapsed >= 1, record
                    assert "timed out" in record["reason"], record
                    assert not (sample / "SRR123.fastq.gz").exists()
                    assert not _raw_validation_marker(sample, "SRR123").exists()
            else:
                assert result is False
                assert "file connection slot deadline" in record["reason"], record
                assert record["elapsed"] < 1.5, record
                assert (part.read_bytes(), part.stat().st_mtime_ns) == before
                assert not (sample / "SRR123.fastq.gz").exists()
                assert not _raw_validation_marker(sample, "SRR123").exists()
            # Prove release directly, without turning native subsecond latency into
            # a correctness requirement. Fresh-attempt success is exercised below.
            with ena_request_slot(site.base + "/SRR123.fastq.gz", kind="file", deadline_at=time.monotonic() + 1):
                pass
            if semaphore is not None:
                assert semaphore.acquire(blocking=False)
                semaphore.release()
    finally:
        (tmp_path / "one-second-control.json").write_text(json.dumps(record, indent=2))
        owner.db.close()


def test_retry_after_expiry_gets_the_same_full_service_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "1")
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=3)
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    sample = tmp_path / "data/getfastq/SRR123"
    try:
        with transfer_site({"SRR123.fastq.gz": payload}, slow=True, slow_bytes=4) as site:
            monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", site.base + "/api")
            started = time.monotonic()
            with pytest.raises(FirstPassDeferred, match="timed out"):
                owner.download_fastq("SRR123", sample.parent, False)
            assert time.monotonic() - started >= 3
            part = sample / "SRR123.fastq.gz.part"
            assert part.read_bytes() == payload[:len(part.read_bytes())]
            assert not _raw_validation_marker(sample, "SRR123").exists()
            site.release.set()
            # Neither attempt changes the configured allowance. A cached old
            # deadline would already be expired and this real retry would refuse.
            assert owner.download_first_pass_seconds == 3
            assert owner.download_fastq("SRR123", sample.parent, False) is True
            assert (sample / "SRR123.fastq.gz").read_bytes() == payload
            proof = json.loads(_raw_validation_marker(sample, "SRR123").read_text())
            assert proof["layout_validation"].startswith("full_fastq_record_scan")
    finally:
        owner.db.close()


def test_one_second_admitted_service_still_expires_and_preserves_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "1")
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=1)
    sample = tmp_path / "data/getfastq/SRR123"
    sample.mkdir(parents=True)
    part = sample / "SRR123.fastq.gz.part"
    part.write_bytes(payload[:4])
    try:
        with transfer_site({"SRR123.fastq.gz": payload}, slow=True, slow_bytes=4) as site:
            monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", site.base + "/api")
            started = time.monotonic()
            with pytest.raises(FirstPassDeferred, match="timed out"):
                owner.download_fastq("SRR123", tmp_path / "data/getfastq", False)
            assert time.monotonic() - started < 2.5
            retained = part.read_bytes()
            assert retained.startswith(payload[:4]) and retained == payload[:len(retained)]
            assert len(retained) < len(payload)
            assert not (sample / "SRR123.fastq.gz").exists()
            assert not _raw_validation_marker(sample, "SRR123").exists()
            assert owner._acquisition_semaphore.acquire(blocking=False)
            owner._acquisition_semaphore.release()
    finally:
        owner.db.close()


@pytest.mark.parametrize("variable,value", [
    ("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "-1"),
    ("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "257"),
    ("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "invalid"),
    ("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "0"),
    ("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "nan"),
    ("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "inf"),
])
def test_acquisition_admission_rejects_invalid_limits_before_creating_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, variable: str, value: str
) -> None:
    monkeypatch.setenv(variable, value)
    with pytest.raises(ValueError):
        build_owner(tmp_path, monkeypatch, staging=False, strict=True)
    assert not (tmp_path / "data/progress.db").exists()


def test_acquisition_queue_expiry_preserves_inputs_and_exception_releases_slot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "1")
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "0.05")
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True)
    sample = tmp_path / "data/getfastq/SRR101"
    sample.mkdir(parents=True)
    part = sample / "SRR101.fastq.gz.part"
    part.write_bytes(b"retained partial")
    before = part.read_bytes(), part.stat().st_mtime_ns
    try:
        assert owner._acquisition_semaphore.acquire(blocking=False)
        try:
            started = time.monotonic()
            with pytest.raises(StagingDeadline, match="admission queue"):
                owner.download_fastq("SRR101", sample.parent, False)
            assert time.monotonic() - started < 1
            assert (part.read_bytes(), part.stat().st_mtime_ns) == before
        finally:
            owner._acquisition_semaphore.release()
        with pytest.raises(StagingError, match="escapes"):
            owner.download_fastq("SRR101", tmp_path / "outside", False)
        assert owner._acquisition_semaphore.acquire(blocking=False)
        owner._acquisition_semaphore.release()
    finally:
        owner.db.close()
