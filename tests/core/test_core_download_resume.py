"""Regression tests for ``core.io.download`` resume accounting, 416 handling and destination inference.

REAL IMPLEMENTATION POLICY: the HTTP tests run the real ``_download_http`` against a real HTTP server
on the loopback interface. The server answers HEAD, GET, byte ranges (206) and 416 for one resource
held in memory, and records each request so tests can check what was fetched. Destination tests use
``file://`` URLs and temporary directories.

Covered here:
- ``resume=False`` over an existing file reported the old file's bytes on top of the new ones
  (20480 for a 10240-byte file).
- A complete destination with ``resume=True`` sent a Range request past the end and failed on HTTP 416.
- A destination directory that does not exist yet was taken for a file name.
"""

from __future__ import annotations

import http.server
import threading
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pytest

from metainformant.core.io.download import (
    DownloadResult,
    _content_range_total,
    _download_http,
    _infer_dest_path,
    download_with_progress,
)


@dataclass
class Resource:
    """What the loopback server serves, and the requests it received."""

    body: bytes
    honour_range: bool = True
    head_length: bool = True
    content_range_on_416: bool = True
    calls: list[tuple[str, str | None]] = field(default_factory=list)  # (method, Range header)


class _ResourceHandler(http.server.BaseHTTPRequestHandler):
    resource: Resource  # set on a per-server subclass

    def log_message(self, format: str, *args: object) -> None:
        return

    def _send(self, status: int, headers: dict[str, str], body: bytes = b"") -> None:
        self.send_response(status)
        for name, value in headers.items():
            self.send_header(name, value)
        self.end_headers()
        if body and self.command != "HEAD":
            self.wfile.write(body)

    def do_HEAD(self) -> None:  # noqa: N802 - http.server dispatch name
        resource = self.resource
        resource.calls.append(("HEAD", self.headers.get("Range")))
        headers = {"Content-Length": str(len(resource.body))} if resource.head_length else {}
        self._send(200, headers)

    def do_GET(self) -> None:  # noqa: N802 - http.server dispatch name
        resource = self.resource
        range_header = self.headers.get("Range")
        resource.calls.append(("GET", range_header))
        size = len(resource.body)
        if range_header and resource.honour_range:
            start = int(range_header.removeprefix("bytes=").split("-")[0])
            if start >= size:
                headers = {"Content-Length": "0"}
                if resource.content_range_on_416:
                    headers["Content-Range"] = f"bytes */{size}"
                self._send(416, headers)
                return
            chunk = resource.body[start:]
            self._send(
                206, {"Content-Range": f"bytes {start}-{size - 1}/{size}", "Content-Length": str(len(chunk))}, chunk
            )
            return
        self._send(200, {"Content-Length": str(size)}, resource.body)


@pytest.fixture
def serve() -> Iterator[Callable[[Resource], str]]:
    """Start loopback servers for resources; return a function giving each one's URL."""
    servers: list[tuple[http.server.ThreadingHTTPServer, threading.Thread]] = []

    def start(resource: Resource) -> str:
        handler = type("Handler", (_ResourceHandler,), {"resource": resource})
        try:
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        except PermissionError as exc:  # sandboxes that forbid binding even loopback sockets
            pytest.skip(f"loopback sockets unavailable: {exc}")
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        servers.append((server, thread))
        return f"http://127.0.0.1:{server.server_port}/resource.bin"

    yield start
    for server, thread in servers:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


@pytest.fixture
def download(serve: Callable[[Resource], str]) -> Callable[..., DownloadResult]:
    def run(
        resource: Resource, dest: Path, *, resume: bool, progress_cb: Callable[..., None] | None = None
    ) -> DownloadResult:
        url = serve(resource)
        return _download_http(url, dest, chunk_size=1024, timeout=5, resume=resume, progress_cb=progress_cb)

    return run


def _resource(size: int = 10240, fill: int = ord("n")) -> bytes:
    return bytes([fill]) * size


class TestResumeAccounting:
    def test_fresh_download(self, tmp_path: Path, download: Callable[..., DownloadResult]) -> None:
        body = bytes(range(256)) * 40
        dest = tmp_path / "out.bin"
        result = download(Resource(body), dest, resume=True)
        assert result.success and not result.resumed
        assert result.bytes_downloaded == len(body) == result.total_bytes
        assert dest.read_bytes() == body

    def test_resume_false_over_an_existing_file_counts_only_the_new_bytes(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        """The old file was counted on top of the rewrite: 20480 reported for a 10240-byte file."""
        body = _resource(10240)
        dest = tmp_path / "out.bin"
        dest.write_bytes(_resource(10240, fill=ord("o")))

        result = download(Resource(body), dest, resume=False)

        assert result.success and not result.resumed
        assert result.bytes_downloaded == 10240
        assert dest.read_bytes() == body  # rewritten, not appended to

    def test_resume_false_over_a_longer_file_counts_only_the_new_bytes(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        body = _resource(2048)
        dest = tmp_path / "out.bin"
        dest.write_bytes(_resource(5000, fill=ord("o")))
        result = download(Resource(body), dest, resume=False)
        assert result.bytes_downloaded == 2048
        assert dest.read_bytes() == body

    def test_partial_file_is_resumed_with_a_range_request(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        body = bytes(range(256)) * 40
        dest = tmp_path / "out.bin"
        dest.write_bytes(body[:4000])
        resource = Resource(body)

        result = download(resource, dest, resume=True)

        assert result.success and result.resumed
        assert result.bytes_downloaded == len(body)
        assert dest.read_bytes() == body
        assert ("GET", "bytes=4000-") in resource.calls

    def test_server_that_ignores_range_restarts_from_zero(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        body = bytes(range(256)) * 40
        dest = tmp_path / "out.bin"
        dest.write_bytes(body[:4000])
        result = download(Resource(body, honour_range=False), dest, resume=True)
        assert result.success and not result.resumed
        assert result.bytes_downloaded == len(body)
        assert dest.read_bytes() == body

    def test_progress_reports_cumulative_bytes_when_resuming(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        body = bytes(range(256)) * 40
        dest = tmp_path / "out.bin"
        dest.write_bytes(body[:3000])
        seen: list[int] = []
        download(
            Resource(body),
            dest,
            resume=True,
            progress_cb=lambda *, bytes_written, total_bytes: seen.append(bytes_written),
        )
        assert seen and seen[0] > 3000 and seen[-1] == len(body)
        assert seen == sorted(seen)


class TestCompleteDestination:
    def test_complete_file_with_resume_is_already_complete(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        """A Range request starting at the end of the file only gets HTTP 416; the file is simply done."""
        body = _resource(10240)
        dest = tmp_path / "out.bin"
        dest.write_bytes(body)
        resource = Resource(body)

        result = download(resource, dest, resume=True)

        assert result.success
        assert result.error is None
        assert result.resumed  # the existing bytes were kept
        assert result.bytes_downloaded == len(body) == result.total_bytes
        assert dest.read_bytes() == body
        assert [method for method, _ in resource.calls] == ["HEAD"]  # nothing was fetched

    def test_complete_file_is_recognised_from_the_416_when_head_hides_the_length(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        body = _resource(10240)
        dest = tmp_path / "out.bin"
        dest.write_bytes(body)
        resource = Resource(body, head_length=False)

        result = download(resource, dest, resume=True)

        assert result.success and result.resumed
        assert result.bytes_downloaded == len(body)
        assert dest.read_bytes() == body
        assert ("GET", "bytes=10240-") in resource.calls  # the 416 itself proved the size

    def test_longer_local_file_is_replaced(self, tmp_path: Path, download: Callable[..., DownloadResult]) -> None:
        """A local file past the end of the resource is not a prefix of it: download the resource afresh."""
        body = _resource(4096)
        dest = tmp_path / "out.bin"
        dest.write_bytes(_resource(6000, fill=ord("o")))

        result = download(Resource(body), dest, resume=True)

        assert result.success and not result.resumed
        assert result.bytes_downloaded == 4096
        assert dest.read_bytes() == body

    def test_unverifiable_size_is_never_taken_for_complete(
        self, tmp_path: Path, download: Callable[..., DownloadResult]
    ) -> None:
        """No length from HEAD and no Content-Range on the 416: the size cannot be checked, so refetch."""
        body = _resource(4096)
        dest = tmp_path / "out.bin"
        dest.write_bytes(_resource(4096, fill=ord("o")))  # same length, different content

        result = download(Resource(body, head_length=False, content_range_on_416=False), dest, resume=True)

        assert result.success and not result.resumed
        assert dest.read_bytes() == body


class TestContentRangeTotal:
    @pytest.mark.parametrize(
        "header,expected",
        [
            ("bytes */100", 100),
            ("bytes 0-9/100", 100),
            ("bytes 4000-10239/10240", 10240),
            ("bytes 0-9/*", None),
            ("bytes 0-9", None),
            ("", None),
            (None, None),
            ("bytes 0-9/abc", None),
        ],
    )
    def test_total_is_read_from_the_header(self, header: str | None, expected: int | None) -> None:
        assert _content_range_total(header) == expected


class TestDestinationInference:
    SRC_URL = "https://example.org/files/reads.fastq.gz"

    def test_existing_directory_takes_the_name_from_the_url(self, tmp_path: Path) -> None:
        assert _infer_dest_path(self.SRC_URL, tmp_path) == tmp_path / "reads.fastq.gz"

    def test_missing_directory_with_a_trailing_separator_is_a_directory(self, tmp_path: Path) -> None:
        """``"out/"`` used to become a file called ``out`` because the directory did not exist yet."""
        missing = tmp_path / "not_yet_created"
        assert _infer_dest_path(self.SRC_URL, str(missing) + "/") == missing / "reads.fastq.gz"

    def test_url_without_a_file_name_gets_the_default_name(self, tmp_path: Path) -> None:
        assert _infer_dest_path("https://example.org/", str(tmp_path / "d") + "/") == tmp_path / "d" / "downloaded_file"

    def test_plain_missing_path_is_the_file_itself(self, tmp_path: Path) -> None:
        assert _infer_dest_path(self.SRC_URL, str(tmp_path / "renamed.fq")) == tmp_path / "renamed.fq"
        assert _infer_dest_path(self.SRC_URL, tmp_path / "renamed.fq") == tmp_path / "renamed.fq"

    def test_path_object_cannot_carry_a_trailing_separator(self, tmp_path: Path) -> None:
        # Documented limitation: ``Path('out/')`` is ``Path('out')``, so only strings can name a missing directory
        assert _infer_dest_path(self.SRC_URL, tmp_path / "out") == tmp_path / "out"

    def test_download_into_a_directory_that_does_not_exist_yet(self, tmp_path: Path) -> None:
        source = tmp_path / "source.bin"
        source.write_bytes(b"x" * 10240)
        target_dir = tmp_path / "new" / "nested"

        result = download_with_progress(
            f"file://{source}", str(target_dir) + "/", show_progress=False, max_retries=1, resume=False
        )

        assert result.success
        assert target_dir.is_dir()
        assert (target_dir / "source.bin").read_bytes() == source.read_bytes()
        assert result.dest_path == target_dir / "source.bin"

    def test_download_to_an_explicit_file_name_is_unchanged(self, tmp_path: Path) -> None:
        source = tmp_path / "source.bin"
        source.write_bytes(b"y" * 2048)
        dest = tmp_path / "chosen_name.bin"
        result = download_with_progress(f"file://{source}", dest, show_progress=False, max_retries=1, resume=False)
        assert result.success and dest.read_bytes() == source.read_bytes()
