"""Real HTTP download boundaries preserve existing genome packages on failure."""

from __future__ import annotations

import io
import threading
import zipfile
from collections.abc import Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from metainformant.dna.external.ncbi_datasets import DatasetsRequest, NCBIDatasetsClient, NCBIDatasetsError


@pytest.fixture
def package_server() -> Iterator[tuple[str, bytes]]:
    """Serve a real zip, a failed request, and a deliberately truncated HTTP body."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_STORED) as archive:
        archive.writestr("ncbi_dataset/data/genomic.fna", ">fixture\n" + "ACGT" * (512 * 1024))
    package = buffer.getvalue()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            if self.path == "/error":
                self.send_response(503)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_response(200)
            length = len(package) + (1024 if self.path == "/truncated" else 0)
            self.send_header("Content-Length", str(length))
            self.send_header("Content-Type", "application/zip")
            self.end_headers()
            self.wfile.write(package)
            self.wfile.flush()
            self.close_connection = True

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}", package
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
        assert not thread.is_alive()


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("failure", ["invalid_url", "error", "truncated"])
def test_failed_download_preserves_existing_package_and_cleans_only_owned_partial(
    tmp_path: Path, package_server: tuple[str, bytes], existing: bool, failure: str
) -> None:
    destination = tmp_path / "assembly.zip"
    previous = b"previous verified genome package"
    if existing:
        destination.write_bytes(previous)
    unrelated = tmp_path / ".assembly.zip.unrelated.part"
    unrelated.write_bytes(b"another transfer owns this file")
    url, _ = package_server

    with NCBIDatasetsClient(timeout=2) as client:
        client.base_url = "not-a-url" if failure == "invalid_url" else url
        with pytest.raises(NCBIDatasetsError):
            client.download_package(DatasetsRequest("/" + failure, {}), destination)

    assert unrelated.read_bytes() == b"another transfer owns this file"
    if existing:
        assert destination.read_bytes() == previous
    else:
        assert not destination.exists()
    assert set(tmp_path.iterdir()) == ({destination, unrelated} if existing else {unrelated})


@pytest.mark.parametrize("existing", [False, True])
def test_complete_download_atomically_publishes_the_package(
    tmp_path: Path, package_server: tuple[str, bytes], existing: bool
) -> None:
    destination = tmp_path / "assembly.zip"
    if existing:
        destination.write_bytes(b"old package")
    url, package = package_server

    with NCBIDatasetsClient(timeout=2) as client:
        client.base_url = url
        assert client.download_package(DatasetsRequest("/success", {}), destination) == destination

    assert destination.read_bytes() == package
    with zipfile.ZipFile(destination) as archive:
        assert archive.namelist() == ["ncbi_dataset/data/genomic.fna"]
    assert list(tmp_path.iterdir()) == [destination]


def test_failed_publication_preserves_destination_directory_and_removes_owned_partial(
    tmp_path: Path, package_server: tuple[str, bytes]
) -> None:
    destination = tmp_path / "assembly.zip"
    destination.mkdir()
    retained = destination / "retained.txt"
    retained.write_bytes(b"do not replace a directory")
    url, _ = package_server

    with NCBIDatasetsClient(timeout=2) as client:
        client.base_url = url
        with pytest.raises(OSError):
            client.download_package(DatasetsRequest("/success", {}), destination)

    assert retained.read_bytes() == b"do not replace a directory"
    assert list(tmp_path.iterdir()) == [destination]
