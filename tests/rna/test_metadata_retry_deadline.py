from __future__ import annotations

import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from metainformant.rna.retrieval.ena_downloader import ENADownloader


@contextmanager
def metadata_site(*, failures: int) -> Iterator[tuple[str, list[str]]]:
    calls: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            calls.append(self.path)
            failing = len(calls) <= failures
            body = b"" if failing else b"run_accession\tfastq_ftp\nSRR123\tftp.sra.ebi.ac.uk/reads.fastq.gz\n"
            self.send_response(500 if failing else 200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/api", calls
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)
        assert not thread.is_alive()


@pytest.mark.parametrize("budget", [0, 3])
def test_configured_metadata_retry_survives_transient_500(budget: int) -> None:
    with metadata_site(failures=1) as (base, calls):
        owner = ENADownloader(run_timeout_seconds=budget, api_retries=1,
                              api_retry_delay_seconds=1, api_base_url=base)
        started = time.monotonic()
        assert owner.get_fastq_urls("SRR123") == ["https://ftp.sra.ebi.ac.uk/reads.fastq.gz"]
        assert len(calls) == 2
        assert time.monotonic() - started < 3


@pytest.mark.parametrize("explicit", [False, True])
def test_metadata_retry_backoff_obeys_both_implicit_and_explicit_deadlines(explicit: bool) -> None:
    with metadata_site(failures=99) as (base, calls):
        owner = ENADownloader(run_timeout_seconds=1, api_retries=3,
                              api_retry_delay_seconds=10, api_base_url=base)
        started = time.monotonic()
        kwargs = {"deadline_at": started + 0.2} if explicit else {}
        assert owner.get_fastq_urls("SRR123", **kwargs) == []
        assert time.monotonic() - started < 1.5
        assert len(calls) == 1
        assert "http_status=500" in owner.last_metadata_error


def test_explicit_zero_retries_remains_one_request_in_bounded_mode() -> None:
    with metadata_site(failures=1) as (base, calls):
        owner = ENADownloader(run_timeout_seconds=3, api_retries=0, api_base_url=base)
        assert owner.get_fastq_urls("SRR123") == []
        assert len(calls) == 1
        assert "http_status=500" in owner.last_metadata_error


def test_expired_explicit_deadline_is_not_replaced_by_local_budget() -> None:
    with metadata_site(failures=0) as (base, calls):
        owner = ENADownloader(run_timeout_seconds=3, api_retries=1, api_base_url=base)
        assert owner.get_fastq_urls("SRR123", deadline_at=0.0) == []
        assert calls == []
