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
        success, message, files = LocalDownloader(
            timeout=10, retries=0, integrity_retries=1
        ).download_run("local", tmp_path)
        assert success, message
        assert requests == [None, None]
        assert gzip.decompress(files[0].read_bytes()) == gzip.decompress(payload)
        assert (
            tmp_path / "reads.fastq.gz.part.invalid"
        ).read_bytes() == b"invalid gzip"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)
