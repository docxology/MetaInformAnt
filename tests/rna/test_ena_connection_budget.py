from __future__ import annotations

import concurrent.futures
import gzip
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from metainformant.rna.retrieval.ena_downloader import (
    ENADownloader,
    _remote_content_length,
)
from metainformant.rna.retrieval.ena_transport import ena_request_slot


@pytest.mark.parametrize("setting", ["-1", "257", "invalid"])
def test_invalid_connection_settings_fail_explicitly(monkeypatch: pytest.MonkeyPatch, setting: str) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", setting)
    with pytest.raises(ValueError):
        ENADownloader()
    with (
        pytest.raises(ValueError),
        ena_request_slot("https://example.invalid/file", kind="file", deadline_at=time.monotonic() + 1),
    ):
        pytest.fail("Invalid settings acquired a lease")


@pytest.mark.parametrize("deadline", [float("inf"), float("nan"), 0.0])
def test_request_deadlines_are_finite_and_not_expired(deadline: float) -> None:
    with (
        pytest.raises((ValueError, TimeoutError)),
        ena_request_slot("https://example.invalid/file", kind="file", deadline_at=deadline),
    ):
        pytest.fail("Invalid deadline acquired a lease")


def test_leases_release_on_exception_and_origins_and_metadata_are_independent(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    monkeypatch.setenv("AMALGKIT_ENA_API_CONNECTIONS", "1")
    url = "https://example.invalid/a"
    with (
        pytest.raises(RuntimeError, match="request failed"),
        ena_request_slot(url, kind="file", deadline_at=time.monotonic() + 2),
        ena_request_slot(url, kind="metadata", deadline_at=time.monotonic() + 2),
        ena_request_slot("https://other.invalid/a", kind="file", deadline_at=time.monotonic() + 2),
    ):
        raise RuntimeError("request failed")
    with ena_request_slot(url, kind="file", deadline_at=time.monotonic() + 2):
        monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "2")
        with (
            pytest.raises(ValueError, match="while requests own it"),
            ena_request_slot(url, kind="file", deadline_at=time.monotonic() + 2),
        ):
            pytest.fail("Live profile changed")
    with ena_request_slot(url, kind="file", deadline_at=time.monotonic() + 2):
        pass


def test_non_http_origins_are_rejected_when_limiting(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    with (
        pytest.raises(ValueError, match="HTTP origin"),
        ena_request_slot("file:///input", kind="file", deadline_at=time.monotonic() + 1),
    ):
        pytest.fail("Non-HTTP request acquired a lease")


def test_metadata_and_transfer_errors_retain_http_and_curl_status(
    tmp_path: Path,
) -> None:
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            self.send_response(503)
            self.send_header("Content-Length", "0")
            self.end_headers()

        def log_message(self, format: str, *args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    owner = ENADownloader(timeout=3, retries=0, api_retries=0, api_base_url=f"{base}/api")
    try:
        assert owner.get_fastq_urls("SRR123") == []
        assert owner.last_metadata_error is not None
        assert "curl_exit=22" in owner.last_metadata_error
        assert "http_status=503" in owner.last_metadata_error
        success, reason, files = owner._download_files([f"{base}/reads.fastq.gz"], tmp_path)
        assert not success and not files
        assert "curl_exit=22" in reason and "http_status=503" in reason
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)


@pytest.mark.parametrize("kind", ["file", "metadata", "head"])
@pytest.mark.parametrize("limit", [0, 2])
def test_connection_cap_is_shared_across_downloaders_with_unlimited_control(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str, limit: int
) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", str(limit))
    monkeypatch.setenv("AMALGKIT_ENA_API_CONNECTIONS", str(limit))
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    lock = threading.Lock()
    active = maximum = 0
    start = threading.Barrier(4)

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self) -> None:
            self.do_GET()

        def do_GET(self) -> None:
            nonlocal active, maximum
            with lock:
                active += 1
                maximum = max(maximum, active)
            try:
                time.sleep(0.3)
                if self.path.startswith("/api?"):
                    accession = parse_qs(urlsplit(self.path).query)["accession"][0]
                    body = f"run_accession\tfastq_ftp\n{accession}\thttp://127.0.0.1:{self.server.server_port}/{accession}.fastq.gz\n".encode()
                else:
                    body = payload
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if self.command != "HEAD":
                    self.wfile.write(body)
            finally:
                with lock:
                    active -= 1

        def log_message(self, format: str, *args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()

    def request(index: int) -> bool:
        start.wait(5)
        base = f"http://127.0.0.1:{server.server_port}"
        owner = ENADownloader(timeout=5, retries=0, api_retries=0, api_base_url=f"{base}/api")
        if kind == "metadata":
            return bool(owner.get_fastq_urls(f"SRR{index}"))
        if kind == "head":
            return _remote_content_length(f"{base}/reads.fastq.gz", timeout=5) == len(payload)
        return owner._download_files([f"{base}/reads.fastq.gz"], tmp_path / str(index))[0]

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            assert all(pool.map(request, range(4)))
        assert maximum == (limit or 4)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)
    assert not thread.is_alive()


def test_slot_wait_expires_without_http_or_partial_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "1")
    entered, release = threading.Event(), threading.Event()
    requests = []
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            requests.append(self.path)
            entered.set()
            if not release.wait(8):
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            try:
                self.wfile.write(payload)
            except BrokenPipeError:
                return

        def log_message(self, format: str, *args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    waiting = tmp_path / "waiting"
    waiting.mkdir()
    partial = waiting / "wait.fastq.gz.part"
    partial.write_bytes(b"retained-prefix")
    before = partial.read_bytes(), partial.stat().st_mtime_ns
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            first = pool.submit(ENADownloader(timeout=6, retries=0)._download_files,
                                [f"{base}/active.fastq.gz"], tmp_path / "active")
            try:
                assert entered.wait(3)
                assert _remote_content_length(f"{base}/waiting-head.fastq.gz", timeout=0.1) is None
                assert requests == ["/active.fastq.gz"]
                started = time.monotonic()
                outcome = ENADownloader(timeout=5, retries=0, run_timeout_seconds=1)._download_files(
                    [f"{base}/wait.fastq.gz"], waiting
                )
                assert not outcome[0]
                assert time.monotonic() - started < 2
                assert requests == ["/active.fastq.gz"]
                assert (partial.read_bytes(), partial.stat().st_mtime_ns) == before
            finally:
                release.set()
            assert first.result(timeout=3)[0]
    finally:
        release.set()
        server.shutdown()
        server.server_close()
        thread.join(5)
    assert not thread.is_alive()
