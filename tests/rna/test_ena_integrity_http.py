"""Real curl transfer controls with an unknown advertised remote size."""

from __future__ import annotations

import gzip
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from metainformant.rna.retrieval.ena_downloader import ENADownloader


def test_unknown_size_corruption_retries_fresh_and_preserves_witness(
    tmp_path: Path,
) -> None:
    requests: list[str | None] = []
    payload = gzip.compress(b"@read\nACGT\n+\n!!!!\n", mtime=0)

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            requests.append(self.headers.get("Range"))
            body = b"invalid gzip" if len(requests) == 1 else payload
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_port}/reads.fastq.gz"

    class LocalDownloader(ENADownloader):
        def get_fastq_urls(self, sample_id: str) -> list[str]:
            return [url]

    try:
        success, message, files = LocalDownloader(timeout=10, retries=0, integrity_retries=1).download_run(
            "local", tmp_path
        )
        assert success, message
        assert requests == [None, None]
        assert gzip.decompress(files[0].read_bytes()) == gzip.decompress(payload)
        assert (tmp_path / "reads.fastq.gz.part.invalid").read_bytes() == b"invalid gzip"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_paired_transfer_overlap_has_a_sequential_negative_control(tmp_path: Path) -> None:
    """Two real requests must arrive together; sequential curl cannot pass this gate."""
    barrier = threading.Barrier(2)
    payload = gzip.compress(b"@read\nACGT\n+\n!!!!\n", mtime=0)
    seen: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            seen.append(self.path)
            try:
                barrier.wait(timeout=1)
            except threading.BrokenBarrierError:
                self.send_error(503, "mate requests did not overlap")
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    class LocalDownloader(ENADownloader):
        def get_fastq_urls(self, sample_id: str) -> list[str]:
            return [f"http://127.0.0.1:{server.server_port}/reads_{mate}.fastq.gz" for mate in (1, 2)]

    try:
        sequential = LocalDownloader(timeout=3, retries=0, file_workers=1)
        assert not sequential.download_run("local", tmp_path / "sequential")[0]
        barrier.reset()
        seen.clear()
        concurrent = LocalDownloader(timeout=3, retries=0, file_workers=2)
        success, message, files = concurrent.download_run("local", tmp_path / "parallel")
        assert success, message
        assert sorted(seen) == ["/reads_1.fastq.gz", "/reads_2.fastq.gz"]
        assert [path.name for path in files] == ["reads_1.fastq.gz", "reads_2.fastq.gz"]
        assert all(gzip.decompress(path.read_bytes()) == gzip.decompress(payload) for path in files)
        before = [(path.read_bytes(), path.stat().st_mtime_ns) for path in files]
        seen.clear()
        assert concurrent.download_run("local", tmp_path / "parallel")[0]
        assert not seen
        assert before == [(path.read_bytes(), path.stat().st_mtime_ns) for path in files]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_parallel_failure_preserves_successful_mate_and_invalid_witness(tmp_path: Path) -> None:
    payload = gzip.compress(b"@read\nACGT\n+\n!!!!\n", mtime=0)
    requests: list[str] = []
    corrupt = True

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            requests.append(self.path)
            body = b"invalid gzip" if self.path.endswith("_2.fastq.gz") and corrupt else payload
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    class LocalDownloader(ENADownloader):
        def get_fastq_urls(self, sample_id: str) -> list[str]:
            return [f"http://127.0.0.1:{server.server_port}/reads_{mate}.fastq.gz" for mate in (1, 2)]

    try:
        downloader = LocalDownloader(timeout=3, retries=0, integrity_retries=0, file_workers=2)
        success, message, files = downloader.download_run("local", tmp_path)
        assert not success and "gzip integrity" in message
        assert [p.name for p in files] == ["reads_1.fastq.gz"]
        first_mtime = files[0].stat().st_mtime_ns
        assert (tmp_path / "reads_2.fastq.gz.part.invalid").read_bytes() == b"invalid gzip"
        requests.clear()
        corrupt = False
        success, message, files = downloader.download_run("local", tmp_path)
        assert success, message
        assert requests == ["/reads_2.fastq.gz"]
        assert files[0].stat().st_mtime_ns == first_mtime
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_file_worker_bounds_and_duplicate_targets_fail_before_transfer(tmp_path: Path) -> None:
    import pytest

    for invalid in (0, 3, True, 1.5):
        with pytest.raises(ValueError, match="file_workers"):
            ENADownloader(file_workers=invalid)

    class DuplicateDownloader(ENADownloader):
        def get_fastq_urls(self, sample_id: str) -> list[str]:
            return ["http://127.0.0.1:1/reads.gz"] * 2

    with pytest.raises(ValueError, match="distinct safe"):
        DuplicateDownloader(file_workers=2).download_run("local", tmp_path)
    assert not list(tmp_path.iterdir())


def test_interrupted_native_curl_resumes_from_retained_bytes(tmp_path: Path) -> None:
    payload = gzip.compress(b"@read\nACGT\n+\n!!!!\n" * 1000, mtime=0)
    requests: list[str | None] = []
    split = len(payload) // 2

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            byte_range = self.headers.get("Range")
            requests.append(byte_range)
            if byte_range is None:
                self.send_response(200)
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload[:split])
                self.close_connection = True
            else:
                assert byte_range == f"bytes={split}-"
                self.send_response(206)
                self.send_header("Content-Length", str(len(payload) - split))
                self.send_header("Content-Range", f"bytes {split}-{len(payload) - 1}/{len(payload)}")
                self.end_headers()
                self.wfile.write(payload[split:])

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    class LocalDownloader(ENADownloader):
        def get_fastq_urls(self, sample_id: str) -> list[str]:
            return [f"http://127.0.0.1:{server.server_port}/reads.fastq.gz"]

    try:
        success, message, files = LocalDownloader(timeout=10, retries=0, retry_delay_seconds=1).download_run(
            "local", tmp_path
        )
        assert success, message
        assert requests == [None, f"bytes={split}-"]
        assert files[0].read_bytes() == payload
        assert not (tmp_path / "reads.fastq.gz.part").exists()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)
