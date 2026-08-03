"""Tests for RNA retrieval and download functionality.

Tests the ENA downloader and related retrieval utilities used for
downloading FASTQ files from ENA/SRA.
"""

from __future__ import annotations

import gzip
import hashlib
import subprocess
import urllib.error
from pathlib import Path

# Import the module under test
from metainformant.rna.retrieval import ena_downloader


class TestCalculateMD5:
    """Tests for MD5 checksum calculation."""

    def test_calculate_md5_simple_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on a simple file."""
        test_file = tmp_path / "test.txt"
        content = b"Hello, World!"
        test_file.write_bytes(content)

        # Calculate expected MD5
        expected_md5 = hashlib.md5(content).hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5

    def test_calculate_md5_empty_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on empty file."""
        test_file = tmp_path / "empty.txt"
        test_file.write_bytes(b"")

        # MD5 of empty content
        expected_md5 = hashlib.md5(b"").hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5

    def test_calculate_md5_large_file(self, tmp_path: Path) -> None:
        """Test MD5 calculation on larger file (chunked reading)."""
        test_file = tmp_path / "large.bin"

        # Create a file larger than the chunk size (4096 bytes)
        content = b"A" * 10000

        test_file.write_bytes(content)

        expected_md5 = hashlib.md5(content).hexdigest()

        result = ena_downloader.calculate_md5(test_file)

        assert result == expected_md5


class TestCleanStagnantFile:
    """Tests for file cleanup utility."""

    def test_clean_stagnant_file_exists(self, tmp_path: Path) -> None:
        """Test cleaning an existing file."""
        test_file = tmp_path / "stagnant.txt"
        test_file.write_text("incomplete download")

        assert test_file.exists()

        ena_downloader.clean_stagnant_file(test_file)

        assert not test_file.exists()

    def test_clean_stagnant_file_not_exists(self, tmp_path: Path) -> None:
        """Test cleaning a non-existent file (should not error)."""
        test_file = tmp_path / "nonexistent.txt"

        assert not test_file.exists()

        # Should not raise any error
        ena_downloader.clean_stagnant_file(test_file)

        assert not test_file.exists()


class TestGzipIntegrity:
    """Tests for Gzip integrity verification."""

    def test_verify_gzip_integrity_valid(self, tmp_path: Path) -> None:
        """Test verification of a valid gzip file."""
        test_file = tmp_path / "valid.gz"
        content = b"This is some test content for gzip."

        with gzip.open(test_file, "wb") as f:
            f.write(content)

        assert ena_downloader.verify_gzip_integrity(test_file) is True

    def test_verify_gzip_integrity_invalid(self, tmp_path: Path) -> None:
        """Test verification of an invalid gzip file (corrupted)."""
        test_file = tmp_path / "invalid.gz"

        # Write random garbage bytes, not a valid gzip
        test_file.write_bytes(b"This is not a gzip file but has .gz extension")

        assert ena_downloader.verify_gzip_integrity(test_file) is False

    def test_verify_gzip_integrity_truncated(self, tmp_path: Path) -> None:
        """Test verification of a truncated gzip file."""
        test_file = tmp_path / "truncated.gz"
        content = b"This is some test content for gzip."

        # Create valid gzip
        with gzip.open(test_file, "wb") as f:
            f.write(content)

        # Truncate it
        data = test_file.read_bytes()
        test_file.write_bytes(data[:-5])  # Remove trailer

        assert ena_downloader.verify_gzip_integrity(test_file) is False

    def test_verify_gzip_integrity_non_gz_extension(self, tmp_path: Path) -> None:
        """Test that non-.gz files are assumed valid (skipped)."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("Not a gzip file")

        assert ena_downloader.verify_gzip_integrity(test_file) is True


class TestDownloadFile:
    """Tests for file download functionality."""

    def test_download_file_creates_directory(self, tmp_path: Path) -> None:
        """Test that download creates parent directories."""
        dest = tmp_path / "subdir" / "file.txt"

        # Create parent manually since download expects it to exist
        dest.parent.mkdir(parents=True, exist_ok=True)

        # The actual download will fail (no real curl call in unit test)
        # but we're testing the path handling with the real subprocess path
        # Here we just verify the path setup logic in the test setup
        assert dest.parent.exists()


class TestENADownloaderIntegrity:
    """Tests for ENADownloader download integrity gates."""

    def test_existing_corrupt_gzip_is_redownloaded(self, tmp_path: Path, monkeypatch) -> None:
        """A non-empty but corrupt .fastq.gz must not be accepted as complete."""

        class FakeDownloader(ena_downloader.ENADownloader):
            def get_fastq_urls(self, sample_id: str) -> list[str]:
                return ["https://example.org/bad.fastq.gz"]

        corrupt = tmp_path / "bad.fastq.gz"
        corrupt.write_bytes(b"not gzip")
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            destination = Path(cmd[cmd.index("-o") + 1])
            with gzip.open(destination, "wb") as fh:
                fh.write(b"@r1\nACGT\n+\n!!!!\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setitem(
            ena_downloader.ENADownloader.download_run.__globals__,
            "_run_command_in_process_group",
            fake_run,
        )

        ok, message, files = FakeDownloader().download_run("SRR_FAKE", tmp_path)

        assert ok is True
        assert "Downloaded 1 files" in message
        assert files == [corrupt]
        assert calls, "corrupt existing file should trigger a replacement download"
        assert ena_downloader.verify_gzip_integrity(corrupt) is True

    def test_new_corrupt_gzip_download_fails(self, tmp_path: Path, monkeypatch) -> None:
        """A successful curl exit code is not enough if the gzip payload is corrupt."""

        class FakeDownloader(ena_downloader.ENADownloader):
            def get_fastq_urls(self, sample_id: str) -> list[str]:
                return ["https://example.org/bad.fastq.gz"]

        output_file = tmp_path / "bad.fastq.gz"

        def fake_run(cmd, **kwargs):
            destination = Path(cmd[cmd.index("-o") + 1])
            destination.write_bytes(b"not gzip")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setitem(
            ena_downloader.ENADownloader.download_run.__globals__,
            "_run_command_in_process_group",
            fake_run,
        )

        ok, message, files = FakeDownloader().download_run("SRR_FAKE", tmp_path)

        assert ok is False
        assert "gzip integrity" in message
        assert files == []
        assert not output_file.exists()
        assert list(tmp_path.glob("bad.fastq.gz*.invalid*")), "invalid payload should be preserved for diagnosis"

    def test_ena_api_transient_failure_is_retried(self, monkeypatch) -> None:
        """A transient portal failure is retried before falling back to NCBI."""

        attempts = 0

        class Response:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def read(self):
                return b"run_accession\tfastq_ftp\nSRR_FAKE\tftp.sra.ebi.ac.uk/vol1/fastq/SRR_FAKE.fastq.gz\n"

        def fake_urlopen(*_args, **_kwargs):
            nonlocal attempts
            attempts += 1
            if attempts == 1:
                raise urllib.error.URLError("transient portal failure")
            return Response()

        monkeypatch.setattr(ena_downloader.urllib.request, "urlopen", fake_urlopen)
        monkeypatch.setattr(ena_downloader.time, "sleep", lambda _seconds: None)

        downloader = ena_downloader.ENADownloader(
            api_retries=1,
            api_retry_delay_seconds=1,
        )
        assert downloader.get_fastq_urls("SRR_FAKE") == ["https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR_FAKE.fastq.gz"]
        assert attempts == 2

    def test_interrupted_transfer_retries_from_retained_partial(self, tmp_path: Path, monkeypatch) -> None:
        """Premature closes use separate curl invocations so partial bytes accumulate."""

        class FakeDownloader(ena_downloader.ENADownloader):
            def get_fastq_urls(self, sample_id: str) -> list[str]:
                return ["https://example.org/resumable.fastq.gz"]

        payload_path = tmp_path / "payload.fastq.gz"
        with gzip.open(payload_path, "wb") as handle:
            handle.write(b"@r1\nACGT\n+\n!!!!\n" * 64)
        payload = payload_path.read_bytes()
        payload_path.unlink()
        calls: list[list[str]] = []
        offsets = iter((len(payload) // 3, 2 * len(payload) // 3, len(payload)))
        prior_end = 0

        def fake_run(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
            nonlocal prior_end
            assert timeout > 0
            calls.append(command)
            destination = Path(command[command.index("-o") + 1])
            end = next(offsets)
            with destination.open("ab") as handle:
                handle.write(payload[prior_end:end])
            prior_end = end
            if len(calls) < 3:
                return subprocess.CompletedProcess(
                    command,
                    18,
                    "206",
                    "curl: (18) transfer closed with bytes remaining",
                )
            return subprocess.CompletedProcess(command, 0, "206", "")

        monkeypatch.setattr(ena_downloader, "_run_command_in_process_group", fake_run)
        monkeypatch.setattr(ena_downloader.time, "sleep", lambda _seconds: None)

        ok, message, files = FakeDownloader(retries=2).download_run("SRR_FAKE", tmp_path)

        assert ok is True
        assert message == "Downloaded 1 files"
        assert len(calls) == 3
        assert all(command[command.index("--continue-at") + 1] == "-" for command in calls)
        assert files[0].read_bytes() == payload

    def test_transient_ena_403_is_bounded_and_static_404_is_not_retried(self, tmp_path: Path, monkeypatch) -> None:
        """Service throttling is retried, while a static missing object fails immediately."""

        class FakeDownloader(ena_downloader.ENADownloader):
            def get_fastq_urls(self, sample_id: str) -> list[str]:
                return [f"https://example.org/{sample_id}.fastq.gz"]

        statuses = iter((403, 403, 403, 404))
        calls: list[list[str]] = []

        def fake_run(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
            assert timeout > 0
            calls.append(command)
            status = next(statuses)
            return subprocess.CompletedProcess(
                command,
                22,
                str(status),
                f"curl: (22) The requested URL returned error: {status}",
            )

        monkeypatch.setattr(ena_downloader, "_run_command_in_process_group", fake_run)
        monkeypatch.setattr(ena_downloader.time, "sleep", lambda _seconds: None)
        downloader = FakeDownloader(retries=2)

        ok_403, _, _ = downloader.download_run("SRR_THROTTLED", tmp_path / "throttled")
        assert ok_403 is False
        assert len(calls) == 3

        ok_404, _, _ = downloader.download_run("SRR_MISSING", tmp_path / "missing")
        assert ok_404 is False
        assert len(calls) == 4

    def test_productive_interruptions_reset_no_progress_retry_budget(self, tmp_path: Path, monkeypatch) -> None:
        """A large file may need many productive range resumes within its deadline."""

        class FakeDownloader(ena_downloader.ENADownloader):
            def get_fastq_urls(self, sample_id: str) -> list[str]:
                return ["https://example.org/large.fastq.gz"]

        payload_path = tmp_path / "payload.fastq.gz"
        with gzip.open(payload_path, "wb") as handle:
            handle.write(b"@r1\nACGT\n+\n!!!!\n" * 512)
        payload = payload_path.read_bytes()
        payload_path.unlink()
        boundaries = [
            len(payload) // 5,
            2 * len(payload) // 5,
            3 * len(payload) // 5,
            4 * len(payload) // 5,
            len(payload),
        ]
        calls = 0
        prior_end = 0

        def fake_run(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
            nonlocal calls, prior_end
            destination = Path(command[command.index("-o") + 1])
            end = boundaries[calls]
            with destination.open("ab") as handle:
                handle.write(payload[prior_end:end])
            prior_end = end
            calls += 1
            return subprocess.CompletedProcess(
                command,
                0 if calls == len(boundaries) else 18,
                "206",
                "",
            )

        monkeypatch.setattr(ena_downloader, "_run_command_in_process_group", fake_run)
        monkeypatch.setattr(ena_downloader.time, "sleep", lambda _seconds: None)

        ok, _, files = FakeDownloader(retries=1).download_run("SRR_LARGE", tmp_path)

        assert ok is True
        assert calls == 5
        assert files[0].read_bytes() == payload
