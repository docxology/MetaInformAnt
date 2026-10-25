from __future__ import annotations

import errno
import gzip
import json
import os
import shutil
import threading
import time
from pathlib import Path

import pytest

from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    _promote_validated_fastq_inputs,
    _raw_validation_marker,
)
from metainformant.rna.retrieval.ena_downloader import ENADownloader
from metainformant.rna.retrieval.ena_staging import (
    FirstPassDeferred,
    StagingDeadline,
    StagingError,
    download_staged_ena,
    prepare_sample_directory,
    validate_fastq_records,
)
from tests.rna.test_ena_staging import transfer_site


def test_ssd_opt_in_alone_requires_records_for_retained_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    owner = build_owner(tmp_path, monkeypatch, staging=True, strict=False, budget=0)
    sample = tmp_path / "data/getfastq/SRR123"
    sample.mkdir(parents=True)
    payloads = {}
    for mate, key in ((1, "old"), (2, "other")):
        name = f"SRR123_{mate}.fastq.gz"
        payloads[name] = gzip.compress(
            f"@{key}/{mate}\nACGT\n+\n!!!!\n".encode(), mtime=0
        )
        (sample / name).write_bytes(payloads[name])
    before = {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in sample.iterdir()
    }
    with transfer_site(payloads) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        try:
            assert (
                owner.download_fastq(
                    "SRR123", sample.parent, True, defer_ncbi_fallback=True
                )
                is None
            )
            assert not _raw_validation_marker(sample, "SRR123").exists()
            for path, expected in before.items():
                assert (path.read_bytes(), path.stat().st_mtime_ns) == expected
        finally:
            owner.db.close()


@pytest.mark.external_tool
@pytest.mark.parametrize("conflict", [False, True])
@pytest.mark.parametrize("local_scratch", [False, True])
def test_native_sra_extraction_preserves_conflicting_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, conflict: bool, local_scratch: bool
) -> None:
    source_value = os.environ.get("AMALGKIT_NATIVE_SRA_TEST_SOURCE")
    if not source_value:
        pytest.skip(
            "Set AMALGKIT_NATIVE_SRA_TEST_SOURCE to a verified single-read SRA fixture"
        )
    source = Path(source_value).resolve(strict=True)
    monkeypatch.setenv(
        "AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR",
        str(tmp_path / "extract-scratch") if local_scratch else "",
    )
    monkeypatch.setenv("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", "0")
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=0)
    sample = tmp_path / "data/getfastq" / source.stem
    sample.mkdir(parents=True)
    existing = sample / f"{source.stem}_1.fastq.gz"
    if conflict:
        existing.write_bytes(gzip.compress(b"@old/1\nACGT\n+\n!!!!\n", mtime=0))
        before = existing.read_bytes(), existing.stat().st_mtime_ns
    try:
        assert owner._extract_sra_to_fastq_bounded(
            source, sample, sample.parent, source.stem, True
        ) is (not conflict)
        if conflict:
            assert (existing.read_bytes(), existing.stat().st_mtime_ns) == before
            assert not (sample / f"{source.stem}.fastq.gz").exists()
            retained = list(
                tmp_path.rglob(f"{source.stem}-conflict-*/output/*.fastq.gz")
            )
            assert retained and validate_fastq_records(retained)
            assert not _raw_validation_marker(sample, source.stem).exists()
        else:
            assert validate_fastq_records([sample / f"{source.stem}.fastq.gz"])
            assert json.loads(_raw_validation_marker(sample, source.stem).read_text())[
                "layout_validation"
            ].startswith("full_fastq_record_scan")
    finally:
        owner.db.close()


@pytest.mark.parametrize("budget,inline", [(0, False), (5, False), (0, True)])
def test_recovery_cannot_adopt_ena_rejected_mates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, budget: int, inline: bool
) -> None:
    owner = build_owner(
        tmp_path, monkeypatch, staging=False, strict=True, budget=budget
    )
    sample = tmp_path / "data/getfastq/SRR123"
    sample.mkdir(parents=True)
    first = gzip.compress(b"@left/1\nACGT\n+\n!!!!\n", mtime=0)
    second = gzip.compress(b"@different/2\nTGCA\n+\n!!!!\n", mtime=0)
    retained = sample / ("SRR123_1.fastq.gz" if inline else "SRR123_1.fastq.gz.part")
    retained.write_bytes(first if inline else first[:16])
    with transfer_site(
        {"SRR123_1.fastq.gz": first, "SRR123_2.fastq.gz": second}
    ) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        try:
            if inline:
                assert owner.download_fastq("SRR123", sample.parent, True) is False
            else:
                if budget:
                    with pytest.raises(FirstPassDeferred):
                        owner.download_fastq(
                            "SRR123", sample.parent, True, defer_ncbi_fallback=True
                        )
                else:
                    assert (
                        owner.download_fastq(
                            "SRR123", sample.parent, True, defer_ncbi_fallback=True
                        )
                        is None
                    )
                assert not _raw_validation_marker(sample, "SRR123").exists()
                assert (
                    owner._download_fastq_ncbi_only("SRR123", sample.parent, True)
                    is False
                )
                assert any(row[2] == "bytes=16-" for row in site.requests)
            assert not _raw_validation_marker(sample, "SRR123").exists()
            assert (sample / "SRR123_1.fastq.gz").read_bytes() == first
            assert (sample / "SRR123_2.fastq.gz").read_bytes() == second
        finally:
            owner.db.close()


@pytest.mark.parametrize("retained_key", ["old", "new"])
def test_extraction_promotion_validates_prospective_union_before_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, retained_key: str
) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "yes")
    sample, staged = tmp_path / "getfastq/SRR123", tmp_path / "staged"
    sample.mkdir(parents=True)
    staged.mkdir()
    existing = sample / "SRR123_1.fastq.gz"
    existing.write_bytes(
        gzip.compress(f"@{retained_key}/1\nACGT\n+\n!!!!\n".encode(), mtime=0)
    )
    inputs = []
    for mate in (1, 2):
        path = staged / f"SRR123_{mate}.fastq.gz"
        path.write_bytes(
            gzip.compress(f"@new/{mate}\nACGT\n+\n!!!!\n".encode(), mtime=0)
        )
        inputs.append(path)
    before = {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in [existing, *inputs]
    }
    outcome = _promote_validated_fastq_inputs(inputs, sample, "SRR123", True)
    assert outcome is (retained_key == "new")
    assert (existing.read_bytes(), existing.stat().st_mtime_ns) == before[existing]
    if retained_key == "old":
        assert not (sample / "SRR123_2.fastq.gz").exists()
        for path in inputs:
            assert (path.read_bytes(), path.stat().st_mtime_ns) == before[path]
        assert not _raw_validation_marker(sample, "SRR123").exists()


@pytest.mark.parametrize("state", ["present", "missing", "replaced"])
def test_directory_creation_cannot_recreate_a_disappeared_root(
    tmp_path: Path, state: str
) -> None:
    root = tmp_path / "selected"
    root.mkdir()
    identity = root.stat().st_dev, root.stat().st_ino
    sample = root / "species/work/getfastq/SRR123"
    if state == "present":
        prepare_sample_directory(sample, root, identity)
        prepare_sample_directory(sample, root, identity)
        assert sample.is_dir()
        return
    root.rename(tmp_path / "preserved-root")
    if state == "replaced":
        root.mkdir()
    with pytest.raises((StagingError, FileNotFoundError)):
        prepare_sample_directory(sample, root, identity)
    assert not sample.exists()
    assert not root.exists() or not list(root.iterdir())


def test_missing_transfer_directory_is_refused_without_recursive_creation(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "absent/root/SRR123"
    outcome = ENADownloader()._download_files(
        ["https://example.invalid/SRR123.fastq.gz"],
        destination,
        require_existing_output=True,
    )
    assert not outcome[0] and not (tmp_path / "absent").exists()


def test_root_replaced_during_real_metadata_request_refuses_bulk_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    owner = build_owner(tmp_path, monkeypatch, staging=False, strict=True, budget=5)
    root = tmp_path / "data"
    raw = root / "work/getfastq"
    raw.mkdir(parents=True)

    def replace_root() -> None:
        root.rename(tmp_path / "preserved-data")
        root.mkdir()

    with transfer_site(
        {"SRR123.fastq.gz": gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)},
        on_metadata=replace_root,
    ) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        try:
            with pytest.raises(StagingError, match="root identity"):
                owner.download_fastq("SRR123", raw, False)
            assert not list(root.iterdir())
            assert not any(row[1].endswith(".fastq.gz") for row in site.requests)
        finally:
            owner.db.close()


@pytest.mark.parametrize("retained", ["final", "partial"])
def test_bulk_resume_must_validate_paired_records_before_witness(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, retained: str
) -> None:
    owner = build_owner(tmp_path, monkeypatch)
    owner.download_first_pass_seconds = 5
    canonical = tmp_path / "data/getfastq/SRR123"
    canonical.mkdir(parents=True)
    first = gzip.compress(b"@left/1\nACGT\n+\n!!!!\n", mtime=0)
    second = gzip.compress(b"@different/2\nTGCA\n+\n!!!!\n", mtime=0)
    path = canonical / (
        "SRR123_1.fastq.gz" if retained == "final" else "SRR123_1.fastq.gz.part"
    )
    path.write_bytes(first if retained == "final" else first[:16])
    with transfer_site(
        {"SRR123_1.fastq.gz": first, "SRR123_2.fastq.gz": second}
    ) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        try:
            with pytest.raises(FirstPassDeferred):
                owner.download_fastq("SRR123", canonical.parent, True)
            assert not (canonical.parent.parent / "raw_validation/SRR123.json").exists()
            assert (canonical / "SRR123_1.fastq.gz").read_bytes() == first
            assert (canonical / "SRR123_2.fastq.gz").read_bytes() == second
            if retained == "partial":
                assert any(row[2] == "bytes=16-" for row in site.requests)
        finally:
            owner.db.close()


@pytest.mark.parametrize(
    "conflict", ["temporary_symlink", "marker_fifo", "oversized_marker"]
)
def test_unsafe_witness_objects_are_preserved_without_reading_or_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, conflict: str
) -> None:
    owner = build_owner(tmp_path, monkeypatch)
    canonical = tmp_path / "data/work/getfastq/SRR123"
    canonical.mkdir(parents=True)
    (canonical / "SRR123.fastq.gz").write_bytes(
        gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    )
    marker = canonical.parent.parent / "raw_validation/SRR123.json"
    marker.parent.mkdir()
    target = tmp_path / "protected"
    target.write_bytes(b"preserve these original bytes")
    temporary = marker.with_name(f".{marker.name}.tmp")
    if conflict == "temporary_symlink":
        temporary.symlink_to(target)
    elif conflict == "marker_fifo":
        os.mkfifo(marker)
    else:
        marker.write_bytes(b"x" * 65537)
    timer = None
    if conflict == "marker_fifo":

        def release_blocked_reader() -> None:
            try:
                descriptor = os.open(marker, os.O_WRONLY | os.O_NONBLOCK)
            except OSError as exc:
                if exc.errno != errno.ENXIO:
                    raise
            else:
                try:
                    os.write(descriptor, b"{}")
                finally:
                    os.close(descriptor)

        timer = threading.Timer(3, release_blocked_reader)
        timer.start()
    started = time.monotonic()
    try:
        with pytest.raises(StagingError):
            owner.download_fastq("SRR123", canonical.parent, False)
        assert time.monotonic() - started < 1
        assert target.read_bytes() == b"preserve these original bytes"
        if conflict == "temporary_symlink":
            assert temporary.is_symlink() and not marker.exists()
        elif conflict == "marker_fifo":
            import stat

            assert stat.S_ISFIFO(marker.stat().st_mode)
        else:
            assert marker.read_bytes() == b"x" * 65537
    finally:
        if timer is not None:
            timer.cancel()
            timer.join(5)
        owner.db.close()


def build_owner(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    staging: bool = True,
    strict: bool = False,
    budget: int = 1,
) -> StreamingPipelineOrchestrator:
    data = tmp_path / "data"
    data.mkdir()
    for key, value in {
        "AMALGKIT_DATA_ROOT": str(data),
        "AMALGKIT_LOCAL_ENA_SCRATCH_DIR": str(tmp_path / "ssd") if staging else "",
        "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS": "yes" if strict else "no",
        "AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB": "0",
        "AMALGKIT_MIN_EXTERNAL_FREE_GB": "0",
        "AMALGKIT_MIN_SYSTEM_FREE_GB": "0",
        "AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS": str(budget),
    }.items():
        monkeypatch.setenv(key, value)
    owner = StreamingPipelineOrchestrator(
        config_dir=tmp_path, log_dir=tmp_path / "logs", db_path=data / "progress.db"
    )
    owner._raw_validation_semaphore = threading.BoundedSemaphore(1)
    return owner


@pytest.mark.parametrize("budget", [0, 5])
def test_strict_bulk_profile_rejects_bad_mates_in_primary_and_full_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, budget: int
) -> None:
    owner = build_owner(
        tmp_path, monkeypatch, staging=False, strict=True, budget=budget
    )
    canonical = tmp_path / "data/getfastq/SRR123"
    canonical.mkdir(parents=True)
    first = gzip.compress(b"@left/1\nACGT\n+\n!!!!\n", mtime=0)
    second = gzip.compress(b"@different/2\nTGCA\n+\n!!!!\n", mtime=0)
    (canonical / "SRR123_1.fastq.gz.part").write_bytes(first[:16])
    with transfer_site(
        {"SRR123_1.fastq.gz": first, "SRR123_2.fastq.gz": second}
    ) as site:
        monkeypatch.setenv("AMALGKIT_ENA_API_BASE_URL", f"{site.base}/api")
        try:
            if budget:
                with pytest.raises(FirstPassDeferred):
                    owner.download_fastq(
                        "SRR123", canonical.parent, True, defer_ncbi_fallback=True
                    )
            else:
                assert (
                    owner.download_fastq(
                        "SRR123", canonical.parent, True, defer_ncbi_fallback=True
                    )
                    is None
                )
            assert owner.local_ena_scratch_dir is None
            assert not (canonical.parent.parent / "raw_validation/SRR123.json").exists()
            assert any(row[2] == "bytes=16-" for row in site.requests)
        finally:
            owner.db.close()


@pytest.mark.parametrize("component", ["sample", "parent"])
@pytest.mark.parametrize("broken", [False, True])
def test_opt_in_refuses_symlinks_before_any_canonical_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, component: str, broken: bool
) -> None:
    owner = build_owner(tmp_path, monkeypatch)
    parent = tmp_path / "data/getfastq"
    target = tmp_path / "data/target"
    if not broken:
        target.mkdir()
        (target / "preserved").write_bytes(b"original")
    link = parent if component == "parent" else parent / "SRR123"
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target, target_is_directory=True)
    try:
        with pytest.raises(StagingError, match="symlink"):
            owner.download_fastq("SRR123", parent, False)
        assert link.is_symlink() and link.readlink() == target
        assert target.exists() is (not broken)
        if not broken:
            assert (target / "preserved").read_bytes() == b"original"
    finally:
        owner.db.close()


@pytest.mark.parametrize(
    "wait", ["initial_validation", "bulk_capacity", "final_validation", "filesystem"]
)
def test_first_pass_bounds_real_waits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, wait: str
) -> None:
    owner = build_owner(tmp_path, monkeypatch)
    parent = tmp_path / "data/getfastq"
    parent.mkdir()
    timer = None
    if wait in {"initial_validation", "final_validation"}:
        owner._raw_validation_semaphore.acquire()
        timer = threading.Timer(3, owner._raw_validation_semaphore.release)
        timer.start()
    if wait == "filesystem":
        owner._fs_lock.acquire()
        timer = threading.Timer(3, owner._fs_lock.release)
        timer.start()
    if wait == "bulk_capacity":
        owner.min_external_free_gb = shutil.disk_usage(parent).free / 1024**3 + 1
    started = time.monotonic()
    try:
        with pytest.raises(StagingDeadline):
            if wait == "final_validation":
                owner._validate_acquisition_fastqs(
                    parent, "SRR123", False, started + 1, assume_payload_valid=True
                )
            else:
                owner.download_fastq("SRR123", parent, False)
        assert time.monotonic() - started < 2
        assert not list(parent.rglob("*.fastq.gz"))
    finally:
        if timer is not None:
            timer.join(5)
        owner.db.close()


def test_slow_metadata_body_obeys_absolute_budget(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    with transfer_site({"SRR123.fastq.gz": payload}, slow_metadata=True) as site:
        owner = ENADownloader(run_timeout_seconds=1, api_base_url=f"{site.base}/api")
        started = time.monotonic()
        assert owner.get_fastq_urls("SRR123") == []
        assert time.monotonic() - started < 2
        assert not owner.file_sizes and not owner.file_md5


@pytest.mark.parametrize(
    "failure",
    ["ids", "counts", "truncated", "quality", "mate_order", "casava_order", "empty"],
)
def test_checksum_correct_invalid_fastq_records_are_never_published(
    tmp_path: Path, failure: str
) -> None:
    first = b"@r/1\nACGT\n+\n!!!!\n"
    second = b"@r/2\nTGCA\n+\n!!!!\n"
    if failure == "ids":
        second = second.replace(b"@r/2", b"@other/2")
    elif failure == "counts":
        second += second
    elif failure == "truncated":
        second = b"@r/2\nTGCA\n+\n"
    elif failure == "quality":
        second = second.replace(b"!!!!", b"!!")
    elif failure == "mate_order":
        second = second.replace(b"/2", b"/1")
    elif failure == "casava_order":
        second = second.replace(b"@r/2", b"@r 1:N:0:1")
    else:
        second = b""
    canonical = tmp_path / "bulk/SRR123"
    canonical.mkdir(parents=True)
    with transfer_site(
        {
            "SRR123_1.fastq.gz": gzip.compress(first, mtime=0),
            "SRR123_2.fastq.gz": gzip.compress(second, mtime=0),
        }
    ) as site:
        downloader = ENADownloader(api_base_url=f"{site.base}/api")
        urls = downloader.get_fastq_urls("SRR123")
        result = download_staged_ena(
            downloader,
            urls,
            canonical,
            tmp_path / "ssd",
            reserve=lambda size: True,
            release=lambda size: None,
            promotion_lock=threading.Lock(),
            validate=lambda path: list(path.glob("*.fastq.gz")),
            expected_paired=True,
        )
    assert result is not None and not result[0]
    assert not list(canonical.iterdir())
    assert len(list((tmp_path / "ssd").rglob("*.fastq.gz"))) == 2
    assert not list((tmp_path / "bulk").rglob("raw_validation"))


def test_bulk_floor_is_rechecked_inside_publication_slot(tmp_path: Path) -> None:
    canonical = tmp_path / "bulk/SRR123"
    canonical.mkdir(parents=True)
    payload = gzip.compress(b"@r\nACGT\n+\n!!!!\n", mtime=0)
    lock = threading.Lock()
    entered = threading.Event()

    def validate(directory: Path) -> list[Path]:
        entered.set()
        return list(directory.glob("*.fastq.gz"))

    lock.acquire()
    with transfer_site({"SRR123.fastq.gz": payload}) as site:
        downloader = ENADownloader(api_base_url=f"{site.base}/api")
        urls = downloader.get_fastq_urls("SRR123")
        floor = shutil.disk_usage(canonical).total + 1

        def publish() -> None:
            with pytest.raises(StagingError, match="bulk capacity"):
                download_staged_ena(
                    downloader,
                    urls,
                    canonical,
                    tmp_path / "ssd",
                    reserve=lambda size: True,
                    release=lambda size: None,
                    promotion_lock=lock,
                    validate=validate,
                    expected_paired=False,
                    minimum_free_bytes=floor,
                )

        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(publish)
            assert entered.wait(5)
            assert not list(canonical.iterdir())
            lock.release()
            future.result(timeout=5)
    assert not list(canonical.iterdir())
    assert next((tmp_path / "ssd").rglob("*.fastq.gz")).read_bytes() == payload
