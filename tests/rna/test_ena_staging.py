from __future__ import annotations

import gzip
import hashlib
import os
import random
import threading
import time
from collections.abc import Callable, Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from uuid import uuid4

import pytest

from metainformant.rna.retrieval.ena_downloader import ENADownloader


@dataclass
class TransferSite:
    """Mutable request ledger and release event owned by one real HTTP server."""

    payloads: dict[str, bytes]
    slow: bool = False
    unknown_size: bool = False
    unknown_get_size: bool = False
    slow_metadata: bool = False
    on_metadata: Callable[[], None] | None = None
    head_extra: int = 0
    slow_bytes: int = 16384
    requests: list[tuple[str, str, str | None]] = field(default_factory=list)
    release: threading.Event = field(default_factory=threading.Event)
    base: str = ""


@contextmanager
def transfer_site(payloads: dict[str, bytes], **settings) -> Iterator[TransferSite]:
    site = TransferSite(payloads, **settings)

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self) -> None:
            site.requests.append(("HEAD", self.path, None))
            body = site.payloads[self.path.lstrip("/")]
            self.send_response(200)
            if not site.unknown_size:
                self.send_header("Content-Length", str(len(body) + site.head_extra))
            self.end_headers()

        def do_GET(self) -> None:
            byte_range = self.headers.get("Range")
            site.requests.append(("GET", self.path, byte_range))
            if self.path.startswith("/api?"):
                if site.on_metadata is not None:
                    site.on_metadata()
                urls = ";".join(f"{site.base}/{name}" for name in site.payloads)
                sizes = ";".join(str(len(body)) for body in site.payloads.values())
                md5 = ";".join(hashlib.md5(body, usedforsecurity=False).hexdigest() for body in site.payloads.values())
                body = f"run_accession\tfastq_ftp\tfastq_bytes\tfastq_md5\nSRR123\t{urls}\t{sizes}\t{md5}\n".encode()
                self.send_response(200)
            else:
                payload = site.payloads[self.path.lstrip("/")]
                offset = int(byte_range.removeprefix("bytes=").removesuffix("-")) if byte_range else 0
                body = payload[offset:]
                self.send_response(206 if byte_range else 200)
                if byte_range:
                    self.send_header(
                        "Content-Range",
                        f"bytes {offset}-{len(payload) - 1}/{len(payload)}",
                    )
            if not site.unknown_get_size:
                self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            try:
                if site.slow_metadata and self.path.startswith("/api?"):
                    for byte in body:
                        self.wfile.write(bytes([byte]))
                        self.wfile.flush()
                        if site.release.wait(0.1):
                            break
                elif site.slow and not self.path.startswith("/api?"):
                    self.wfile.write(body[: site.slow_bytes])
                    self.wfile.flush()
                    site.release.wait(10)
                    self.wfile.write(body[site.slow_bytes :])
                else:
                    self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                return

        def log_message(self, format: str, *args) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    site.base = f"http://127.0.0.1:{server.server_port}"
    thread.start()
    try:
        yield site
    finally:
        site.release.set()
        server.shutdown()
        server.server_close()
        thread.join(5)
        assert not thread.is_alive()


def test_first_pass_deadline_is_shared_by_paired_files_and_retains_bytes(
    tmp_path: Path,
) -> None:
    rng = random.Random(7)
    reads = ["".join(rng.choice("ACGT") for _ in range(100)) for _ in range(2000)]
    payload = gzip.compress(
        "".join(f"@r{i}\n{read}\n+\n{'I' * 100}\n" for i, read in enumerate(reads)).encode(),
        mtime=0,
    )
    with transfer_site({"SRR123_1.fastq.gz": payload, "SRR123_2.fastq.gz": payload}, slow=True) as site:
        downloader = ENADownloader(timeout=30, retries=0, run_timeout_seconds=1, file_workers=1)
        start = time.monotonic()
        success, _message, files = downloader._download_files(
            [f"{site.base}/{name}" for name in site.payloads], tmp_path
        )
        assert not success and not files
        assert time.monotonic() - start < 3
        retained = (tmp_path / "SRR123_1.fastq.gz.part").read_bytes()
        assert 0 < len(retained) <= 16384
        assert retained == payload[: len(retained)]
        assert not (tmp_path / "SRR123_2.fastq.gz.part").exists()


def test_transfer_size_bound_stops_oversized_response(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n" * 1000, mtime=0)
    with transfer_site({"SRR123.fastq.gz": payload}) as site:
        success, _message, files = ENADownloader(timeout=10, retries=0)._download_files(
            [f"{site.base}/SRR123.fastq.gz"],
            tmp_path,
            max_file_bytes={"SRR123.fastq.gz": 10},
        )
        assert not success and not files
        assert not (tmp_path / "SRR123.fastq.gz").exists()


@pytest.mark.parametrize("offset", [0, 32768])
def test_unknown_length_transfer_bound_includes_retained_prefix(tmp_path: Path, offset: int) -> None:
    payload = gzip.compress(random.Random(82).randbytes(128 * 1024), compresslevel=0, mtime=0)
    partial = tmp_path / "SRR123.fastq.gz.part"
    partial.write_bytes(payload[:offset])
    with transfer_site({"SRR123.fastq.gz": payload}, unknown_get_size=True) as site:
        result = ENADownloader(timeout=10, retries=0)._download_files(
            [f"{site.base}/SRR123.fastq.gz"],
            tmp_path,
            max_file_bytes={"SRR123.fastq.gz": 49152},
        )
    assert not result[0]
    retained = partial.read_bytes()
    assert offset <= len(retained) <= 49152
    assert retained == payload[: len(retained)]


@pytest.mark.parametrize("slots", [1, 2])
def test_publication_serialization_has_a_concurrent_negative_control(tmp_path: Path, slots: int) -> None:
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n" * 2_000_000, compresslevel=0, mtime=0)
    barrier = threading.Barrier(2)
    stopped = threading.Event()
    maximum = [0]

    def validate(directory: Path) -> list[Path]:
        barrier.wait(timeout=20)
        return list(directory.glob("*.fastq.gz"))

    def observe() -> None:
        while not stopped.wait(0.001):
            maximum[0] = max(maximum[0], len(list(tmp_path.glob("*/bulk/**/*.copying"))))

    watcher = threading.Thread(target=observe)
    watcher.start()
    lock = threading.Semaphore(slots)
    try:
        with (
            transfer_site({"SRR123.fastq.gz": payload}) as site,
            ThreadPoolExecutor(max_workers=2) as pool,
        ):
            jobs = [
                pool.submit(
                    stage_download,
                    site,
                    tmp_path / str(i),
                    validate=validate,
                    lock=lock,
                )
                for i in range(2)
            ]
            results = [job.result(timeout=30) for job in jobs]
        assert all(result is not None and result[0] for result in results)
    finally:
        stopped.set()
        watcher.join(5)
    assert maximum[0] == slots


@pytest.mark.external_tool
def test_real_cross_volume_publication_preserves_source_and_replay(
    tmp_path: Path,
) -> None:
    from metainformant.rna.retrieval.ena_staging import (
        StagingError,
        publish_verified_file,
    )

    bulk = os.environ.get("AMALGKIT_STAGING_TEST_BULK_ROOT")
    if not bulk:
        pytest.skip("Set AMALGKIT_STAGING_TEST_BULK_ROOT for an isolated cross-volume pilot")
    root = Path(bulk)
    assert root.is_dir() and root.stat().st_dev != tmp_path.stat().st_dev
    destination_dir = root / f"ssd-publication-pilot-{uuid4().hex}"
    destination_dir.mkdir()
    source = tmp_path / "source.fastq.gz"
    source.write_bytes(gzip.compress(random.Random(71).randbytes(8 * 1024**2), mtime=0))
    before = source.read_bytes(), source.stat().st_mtime_ns
    destination = destination_dir / source.name
    publish_verified_file(source, destination)
    assert (destination.read_bytes(), destination.stat().st_mtime_ns) == before
    assert (source.read_bytes(), source.stat().st_mtime_ns) == before
    publish_verified_file(source, destination)
    assert (destination.read_bytes(), destination.stat().st_mtime_ns) == before
    source.write_bytes(b"conflicting bytes")
    with pytest.raises(StagingError, match="conflicts"):
        publish_verified_file(source, destination)
    assert (destination.read_bytes(), destination.stat().st_mtime_ns) == before


def test_staging_publishes_only_validated_complete_files(tmp_path: Path) -> None:
    from metainformant.rna.retrieval.ena_staging import download_staged_ena

    payload = gzip.compress(b"@SRR123.1\nACGT\n+\n!!!!\n", mtime=0)
    canonical = tmp_path / "bulk/species/work/getfastq/SRR123"
    canonical.mkdir(parents=True)
    scratch = tmp_path / "ssd"
    reserved: list[int] = []
    released: list[int] = []
    validated: list[Path] = []

    def validate(directory: Path) -> list[Path]:
        assert directory.is_relative_to(scratch)
        assert not list(canonical.glob("*.fastq.gz"))
        paths = list(directory.glob("*.fastq.gz"))
        assert gzip.decompress(paths[0].read_bytes()).startswith(b"@SRR123")
        validated.extend(paths)
        return paths

    with transfer_site({"SRR123.fastq.gz": payload}) as site:
        outcome = download_staged_ena(
            ENADownloader(timeout=10, retries=0),
            [f"{site.base}/SRR123.fastq.gz"],
            canonical,
            scratch,
            reserve=lambda size: not reserved.append(size),
            release=released.append,
            promotion_lock=threading.Lock(),
            validate=validate,
            expected_paired=False,
        )
    assert outcome is not None and outcome[0], outcome
    assert outcome[2] == [canonical / "SRR123.fastq.gz"]
    assert outcome[2][0].read_bytes() == payload
    assert len(validated) == 1 and reserved == released


def stage_download(
    site: TransferSite,
    root: Path,
    *,
    downloader: ENADownloader | None = None,
    validate=None,
    reserve=None,
    lock=None,
):
    from metainformant.rna.engine.streaming_orchestrator import (
        _validated_local_fastq_inputs,
    )
    from metainformant.rna.retrieval.ena_staging import download_staged_ena

    canonical = root / "bulk/species/work/getfastq/SRR123"
    canonical.mkdir(parents=True, exist_ok=True)
    downloader = downloader or ENADownloader(timeout=10, retries=0, api_base_url=f"{site.base}/api")
    urls = downloader.get_fastq_urls("SRR123")
    return download_staged_ena(
        downloader,
        urls,
        canonical,
        root / "ssd",
        reserve=reserve or (lambda size: True),
        release=lambda size: None,
        promotion_lock=lock or threading.Lock(),
        expected_paired=False,
        validate=validate
        or (lambda path: _validated_local_fastq_inputs(path, "SRR123", False, assume_payload_valid=True)),
    )


def test_interrupted_ssd_transfer_resumes_with_bound_ena_checksum(
    tmp_path: Path,
) -> None:
    rng = random.Random(17)
    reads = ["".join(rng.choice("ACGT") for _ in range(100)) for _ in range(2000)]
    payload = gzip.compress(
        "".join(f"@r{i}\n{read}\n+\n{'I' * 100}\n" for i, read in enumerate(reads)).encode(),
        mtime=0,
    )
    with transfer_site({"SRR123.fastq.gz": payload}, slow=True) as site:
        first = ENADownloader(
            timeout=10,
            retries=0,
            run_timeout_seconds=1,
            api_base_url=f"{site.base}/api",
        )
        result = stage_download(site, tmp_path, downloader=first)
        assert result is not None and not result[0]
        partial = next((tmp_path / "ssd").rglob("*.part"))
        retained = partial.read_bytes()
        assert 0 < len(retained) <= 16384
        assert retained == payload[: len(retained)]
        assert not list((tmp_path / "bulk").rglob("*.fastq.gz"))
        site.slow = False
        site.release.set()
        result = stage_download(site, tmp_path)
        assert result is not None and result[0], result
        assert result[2][0].read_bytes() == payload
        requests = [row for row in site.requests if row[0] == "GET" and row[1].endswith(".fastq.gz")]
        assert [row[2] for row in requests] == [None, f"bytes={len(retained)}-"]
        assert not list((tmp_path / "ssd").rglob("*.part"))


@pytest.mark.parametrize("case", ["unknown", "capacity", "legacy_partial", "legacy_mate"])
def test_staging_declines_without_touching_bulk_inputs(tmp_path: Path, case: str) -> None:
    from metainformant.rna.retrieval.ena_staging import download_staged_ena

    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    canonical = tmp_path / "bulk/SRR123"
    canonical.mkdir(parents=True)
    if case == "legacy_partial":
        (canonical / "SRR123.fastq.gz.part").write_bytes(payload[:20])
    elif case == "legacy_mate":
        (canonical / "SRR123.fastq.gz").write_bytes(payload)
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in canonical.iterdir()}
    with transfer_site({"SRR123.fastq.gz": payload}, unknown_size=case == "unknown") as site:
        result = download_staged_ena(
            ENADownloader(timeout=10, retries=0),
            [f"{site.base}/SRR123.fastq.gz"],
            canonical,
            tmp_path / "ssd",
            reserve=lambda size: case != "capacity",
            release=lambda size: None,
            promotion_lock=threading.Lock(),
            validate=lambda path: list(path.glob("*.fastq.gz")),
            expected_paired=False,
        )
        assert result is None
        assert not any(row[0] == "GET" for row in site.requests)
    assert before == {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in canonical.iterdir()}


def test_publication_conflict_preserves_both_canonical_and_ssd_sources(
    tmp_path: Path,
) -> None:
    from metainformant.rna.retrieval.ena_staging import StagingError

    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    other = gzip.compress(b"@r\nTGCA\n+\n!!!!\n", mtime=0)
    destination = tmp_path / "bulk/species/work/getfastq/SRR123/SRR123.fastq.gz"

    def validate(directory: Path) -> list[Path]:
        destination.write_bytes(other)
        return list(directory.glob("*.fastq.gz"))

    with (
        transfer_site({"SRR123.fastq.gz": payload}) as site,
        pytest.raises(StagingError, match="conflicts"),
    ):
        stage_download(site, tmp_path, validate=validate)
    assert destination.read_bytes() == other
    assert next((tmp_path / "ssd").rglob("SRR123.fastq.gz")).read_bytes() == payload
    assert not list(destination.parent.glob("*.copying"))


def test_destination_directory_replacement_is_refused(tmp_path: Path) -> None:
    from metainformant.rna.retrieval.ena_staging import StagingError

    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    canonical = tmp_path / "bulk/species/work/getfastq/SRR123"

    def validate(directory: Path) -> list[Path]:
        canonical.rename(canonical.with_name("preserved-SRR123"))
        canonical.mkdir()
        return list(directory.glob("*.fastq.gz"))

    with (
        transfer_site({"SRR123.fastq.gz": payload}) as site,
        pytest.raises(StagingError, match="identity changed"),
    ):
        stage_download(site, tmp_path, validate=validate)
    assert not list(canonical.iterdir())
    assert next((tmp_path / "ssd").rglob("SRR123.fastq.gz")).read_bytes() == payload


def test_ssd_md5_mismatch_is_never_promoted(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    with transfer_site({"SRR123.fastq.gz": payload}) as site:
        downloader = ENADownloader(timeout=10, retries=0, api_base_url=f"{site.base}/api")
        urls = downloader.get_fastq_urls("SRR123")
        downloader.file_md5["SRR123.fastq.gz"] = "0" * 32
        from metainformant.rna.retrieval.ena_staging import download_staged_ena

        canonical = tmp_path / "bulk/SRR123"
        canonical.mkdir(parents=True)
        result = download_staged_ena(
            downloader,
            urls,
            canonical,
            tmp_path / "ssd",
            reserve=lambda size: True,
            release=lambda size: None,
            promotion_lock=threading.Lock(),
            validate=lambda directory: list(directory.glob("*.fastq.gz")),
            expected_paired=False,
        )
        assert result is not None and not result[0]
        assert not list(canonical.iterdir())
        assert list((tmp_path / "ssd").rglob("*.invalid"))


def test_declared_paired_unnumbered_input_keeps_bulk_path(tmp_path: Path) -> None:
    from metainformant.rna.retrieval.ena_staging import download_staged_ena

    canonical = tmp_path / "bulk/SRR123"
    canonical.mkdir(parents=True)
    result = download_staged_ena(
        ENADownloader(),
        ["https://example.invalid/SRR123.fastq.gz"],
        canonical,
        tmp_path / "ssd",
        reserve=lambda size: True,
        release=lambda size: None,
        promotion_lock=threading.Lock(),
        validate=lambda directory: [],
        expected_paired=True,
    )
    assert result is None and not (tmp_path / "ssd").exists()


def test_promotion_lock_wait_respects_first_pass_deadline(tmp_path: Path) -> None:
    from metainformant.rna.retrieval.ena_staging import StagingDeadline

    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    lock = threading.Lock()
    with lock, transfer_site({"SRR123.fastq.gz": payload}) as site:
        downloader = ENADownloader(
            timeout=10,
            retries=0,
            run_timeout_seconds=1,
            api_base_url=f"{site.base}/api",
        )
        started = time.monotonic()
        with pytest.raises(StagingDeadline):
            stage_download(site, tmp_path, downloader=downloader, lock=lock)
        assert time.monotonic() - started < 3
    assert not list((tmp_path / "bulk").rglob("*.fastq.gz"))
    assert next((tmp_path / "ssd").rglob("SRR123.fastq.gz")).read_bytes() == payload


@pytest.mark.external_tool
@pytest.mark.parametrize("staging", [False, True])
def test_native_first_pass_defers_then_ssd_retry_quantifies_and_reuses(
    native_acquisition_factory,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    staging: bool,
) -> None:
    from metainformant.rna.engine.acquisition_manifest import sha256_file
    from metainformant.rna.engine.acquisition_worker import run_manifest
    from metainformant.rna.engine.durable_quant import (
        DirectoryStore,
        lock_quantification,
        restore_quantification,
    )
    from metainformant.rna.engine.progress_db import ProgressDB
    from metainformant.rna.engine.provenance import classify_quantification

    fixture = native_acquisition_factory()
    source = fixture.work / "getfastq/SRR123/SRR123.fastq.gz"
    payload = source.read_bytes()
    source.rename(tmp_path / "preserved-input.fastq.gz")
    monkeypatch.setenv("AMALGKIT_LOCAL_ENA_SCRATCH_DIR", str(tmp_path / "ssd") if staging else "")
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "yes")
    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB", "0")
    monkeypatch.setenv("AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS", "1")
    monkeypatch.setenv("AMALGKIT_PIPELINE_DEDICATED_RECOVERY", "yes")
    kwargs = {
        "manifest_path": fixture.manifest,
        "data_root": fixture.data,
        "config_dir": fixture.config_dir,
        "workers": 2,
        "threads": 2,
        "fastq_threads": 1,
        "compression_threads": 1,
        "validation_slots": 1,
        "quant_slots": 1,
        "fasterq_slots": 1,
        "max_in_flight": 2,
    }
    restored = []

    def checkpoint(species: str, accession: str) -> None:
        store = DirectoryStore(tmp_path / "locked")
        config_hash, reference_hash = (
            sha256_file(fixture.config),
            sha256_file(fixture.index),
        )
        sample = fixture.work / "quant" / accession
        lock_quantification(
            store,
            "staged-test",
            species,
            accession,
            sample,
            expected_config_sha256=config_hash,
            expected_reference_index_sha256=reference_hash,
        )
        destination = tmp_path / f"restored-{len(restored)}" / accession
        restore_quantification(
            store,
            "staged-test",
            species,
            accession,
            destination,
            expected_config_sha256=config_hash,
            expected_reference_index_sha256=reference_hash,
            verified_config_path=fixture.config,
        )
        assert (
            classify_quantification(destination, accession, expected_config_path=fixture.config)["status"] == "current"
        )
        restored.append(accession)

    kwargs["on_quantified"] = checkpoint
    with transfer_site({"SRR123.fastq.gz": payload}, slow=True, slow_bytes=10) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        deferred = run_manifest(**kwargs)
        assert deferred["counts"] == {"deferred_download": 1}
        assert deferred["results"][0]["status"] == "deferred_download"
        assert not deferred["results"][0]["fallback_deferred"]
        assert not source.exists()
        assert not restored
        database = ProgressDB(fixture.data / "pipeline_progress.db")
        try:
            assert database.get_state("test_species", "SRR123") == "pending"
        finally:
            database.close()
        site.slow = False
        site.release.set()
        monkeypatch.setenv("AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS", "0")
        completed = run_manifest(**kwargs)
        assert completed["counts"] == {"quantified": 1, "newly_quantified": 1}
        quant = fixture.work / "quant/SRR123"
        before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in quant.iterdir() if p.is_file()}
        requests_before = len(site.requests)
        replay = run_manifest(**kwargs)
        assert replay["counts"] == {"quantified": 1, "reused": 1}
        assert len(site.requests) == requests_before
        assert before == {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in quant.iterdir() if p.is_file()}
        assert source.read_bytes() == payload
        assert not list((tmp_path / "ssd").rglob("*.part"))
        assert restored == ["SRR123", "SRR123"]


@pytest.mark.parametrize("case", ["empty", "symlink", "fifo", "missing_destination", "expired"])
def test_verified_publication_refuses_unsafe_sources(tmp_path: Path, case: str) -> None:
    import os

    from metainformant.rna.retrieval.ena_staging import (
        StagingError,
        publish_verified_file,
    )

    source = tmp_path / "source"
    source.write_bytes(b"valid input bytes")
    destination = tmp_path / "bulk/result.fastq.gz"
    destination.parent.mkdir()
    deadline = None
    if case == "empty":
        source.write_bytes(b"")
    elif case == "symlink":
        source.rename(tmp_path / "preserved")
        source.symlink_to(tmp_path / "preserved")
    elif case == "fifo":
        source.rename(tmp_path / "preserved")
        os.mkfifo(source)
    elif case == "missing_destination":
        destination = tmp_path / "unavailable/result.fastq.gz"
    else:
        deadline = time.monotonic() - 1
    with pytest.raises(StagingError):
        publish_verified_file(source, destination, deadline=deadline)
    assert not destination.exists()


def test_verified_publication_reuses_equal_bytes_without_touching_mtime(
    tmp_path: Path,
) -> None:
    from metainformant.rna.retrieval.ena_staging import publish_verified_file

    source, destination = tmp_path / "source", tmp_path / "destination"
    source.write_bytes(b"same current payload")
    destination.write_bytes(source.read_bytes())
    before = destination.stat().st_mtime_ns
    publish_verified_file(source, destination)
    assert destination.stat().st_mtime_ns == before
    assert source.exists()


def test_copy_source_mutation_refuses_publication_and_cleans_owned_temporary(
    tmp_path: Path,
) -> None:
    from metainformant.rna.retrieval.ena_staging import (
        StagingError,
        publish_verified_file,
    )

    source = tmp_path / "source"
    with source.open("wb") as handle:
        handle.truncate(64 * 1024**2)
    destination = tmp_path / "bulk/result.fastq.gz"
    destination.parent.mkdir()
    stop, changed = threading.Event(), threading.Event()

    def mutate() -> None:
        while not stop.wait(0.0005):
            if list(destination.parent.glob("*.copying")):
                with source.open("ab") as handle:
                    handle.write(b"changed")
                changed.set()
                return

    writer = threading.Thread(target=mutate)
    writer.start()
    try:
        with pytest.raises(StagingError, match="changed during copy"):
            publish_verified_file(source, destination)
        assert changed.wait(1)
    finally:
        stop.set()
        writer.join(5)
    assert not writer.is_alive() and not destination.exists()
    assert not list(destination.parent.glob("*.copying"))


@pytest.mark.parametrize(
    "case",
    ["identity", "lock_symlink", "lock_fifo", "identity_fifo", "oversize_partial"],
)
def test_saved_stage_corruption_is_refused(tmp_path: Path, case: str) -> None:
    from metainformant.rna.retrieval.ena_staging import StagingError

    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    with transfer_site({"SRR123.fastq.gz": payload}, slow=True, slow_bytes=10) as site:
        downloader = ENADownloader(
            timeout=10,
            retries=0,
            run_timeout_seconds=1,
            api_base_url=f"{site.base}/api",
        )
        result = stage_download(site, tmp_path, downloader=downloader)
        assert result is not None and not result[0]
        workspace = next((tmp_path / "ssd").iterdir())
        if case in {"lock_fifo", "identity_fifo"}:
            import os

            path = workspace / (".lock" if case == "lock_fifo" else "identity.json")
            path.rename(path.with_name(path.name + ".preserved"))
            os.mkfifo(path)
        elif case == "identity":
            (workspace / "identity.json").write_text("{}")
        elif case == "lock_symlink":
            lock = workspace / ".lock"
            lock.rename(workspace / ".preserved-lock")
            outside = tmp_path / "outside"
            outside.write_bytes(b"must remain unchanged")
            lock.symlink_to(outside)
        else:
            next(workspace.rglob("*.part")).write_bytes(b"x" * (len(payload) + 1))
        site.release.set()
        site.slow = False
        with pytest.raises(StagingError):
            stage_download(site, tmp_path)
    assert not list((tmp_path / "bulk").rglob("*.fastq.gz"))
    if case == "lock_symlink":
        assert outside.read_bytes() == b"must remain unchanged"


def test_real_scratch_ledger_refuses_stage_after_another_stage_reserves_capacity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import shutil

    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )
    from metainformant.rna.retrieval.ena_staging import download_staged_ena

    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB", "0")
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path))
    owner = StreamingPipelineOrchestrator(log_dir=tmp_path / "logs", db_path=tmp_path / "progress.db")
    scratch, canonical = tmp_path / "ssd", tmp_path / "bulk/SRR123"
    scratch.mkdir()
    canonical.mkdir(parents=True)
    free = shutil.disk_usage(scratch).free
    allocated = free // 2
    assert owner._reserve_local_scratch(scratch, allocated, purpose="quantification test reservation")
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    try:
        with transfer_site({"SRR123.fastq.gz": payload}) as site:
            downloader = ENADownloader(timeout=10, retries=0)
            downloader.file_sizes["SRR123.fastq.gz"] = free
            result = download_staged_ena(
                downloader,
                [f"{site.base}/SRR123.fastq.gz"],
                canonical,
                scratch,
                reserve=lambda size: owner._reserve_local_scratch(scratch, size, purpose="ENA test reservation"),
                release=lambda size: owner._release_local_scratch(size, purpose="ENA test reservation"),
                promotion_lock=threading.Lock(),
                validate=lambda path: [],
                expected_paired=False,
            )
            assert result is None and not any(row[0] == "GET" for row in site.requests)
            assert owner._local_scratch_reserved_bytes == allocated
    finally:
        owner._release_local_scratch(allocated, purpose="quantification test reservation")
        owner.db.close()
    assert owner._local_scratch_reserved_bytes == 0


@pytest.mark.external_tool
def test_interrupted_native_worker_releases_lock_retains_ssd_partial_and_resumes(
    native_acquisition_factory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    import os
    import signal
    import subprocess
    import sys

    from metainformant.rna.engine.acquisition_worker import run_manifest

    fixture = native_acquisition_factory()
    source = fixture.work / "getfastq/SRR123/SRR123.fastq.gz"
    rng = random.Random(23)
    lines = gzip.decompress(source.read_bytes()).decode().splitlines()
    for index in range(0, len(lines), 4):
        lines[index] += " " + "".join(rng.choice("ACGT") for _ in range(2000))
    payload = gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0)
    assert len(payload) > 16384
    source.rename(tmp_path / "preserved-input.fastq.gz")
    monkeypatch.setenv("AMALGKIT_LOCAL_ENA_SCRATCH_DIR", str(tmp_path / "ssd"))
    monkeypatch.setenv("AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB", "0")
    monkeypatch.setenv("AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS", "0")
    monkeypatch.setenv("AMALGKIT_PIPELINE_DEDICATED_RECOVERY", "yes")
    script = tmp_path / "worker.py"
    script.write_text(
        "import sys\nfrom pathlib import Path\nfrom metainformant.rna.engine.acquisition_worker import run_manifest\n"
        "run_manifest(manifest_path=Path(sys.argv[1]),data_root=Path(sys.argv[2]),config_dir=Path(sys.argv[3]),"
        "workers=1,threads=2,fastq_threads=1,compression_threads=1,validation_slots=1,quant_slots=1,"
        "fasterq_slots=1,max_in_flight=1)\n"
    )
    with transfer_site({"SRR123.fastq.gz": payload}, slow=True) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        with (tmp_path / "worker.log").open("w") as log:
            worker = subprocess.Popen(
                [
                    sys.executable,
                    str(script),
                    str(fixture.manifest),
                    str(fixture.data),
                    str(fixture.config_dir),
                ],
                stdout=log,
                stderr=log,
            )
            try:
                deadline = time.monotonic() + 15
                partial = None
                while time.monotonic() < deadline:
                    parts = list((tmp_path / "ssd").rglob("*.part"))
                    if parts and parts[0].stat().st_size > 0:
                        partial = parts[0]
                        break
                    assert worker.poll() is None, (tmp_path / "worker.log").read_text()
                    threading.Event().wait(0.01)
                assert partial is not None
                rows = [
                    line.split(None, 3)
                    for line in subprocess.run(
                        ["ps", "-axo", "pid=,ppid=,pgid=,command="],
                        capture_output=True,
                        text=True,
                        check=True,
                    ).stdout.splitlines()
                ]
                curl = [row for row in rows if int(row[1]) == worker.pid and "curl" in Path(row[3].split()[0]).name]
                assert len(curl) == 1 and int(curl[0][0]) == int(curl[0][2])
                os.kill(worker.pid, signal.SIGINT)
                os.killpg(int(curl[0][2]), signal.SIGTERM)
                assert worker.wait(timeout=15) != 0
                retained = partial.read_bytes()
                assert 0 < len(retained) <= site.slow_bytes
                assert retained == payload[: len(retained)]
                assert not source.exists()
            finally:
                site.release.set()
                if worker.poll() is None:
                    worker.kill()
                    worker.wait(timeout=5)
        site.slow = False
        recovered = run_manifest(
            manifest_path=fixture.manifest,
            data_root=fixture.data,
            config_dir=fixture.config_dir,
            workers=1,
            threads=2,
            fastq_threads=1,
            compression_threads=1,
            validation_slots=1,
            quant_slots=1,
            fasterq_slots=1,
            max_in_flight=1,
        )
        assert recovered["counts"] == {"quantified": 1, "newly_quantified": 1}
        assert any(row[2] == f"bytes={len(retained)}-" for row in site.requests)
        assert source.read_bytes() == payload
