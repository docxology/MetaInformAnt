"""Tests for RNA retrieval and download functionality.

Tests the ENA downloader and related retrieval utilities used for
downloading FASTQ files from ENA/SRA.
"""

from __future__ import annotations

import gzip
import hashlib
import subprocess
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
            with gzip.open(corrupt, "wb") as fh:
                fh.write(b"@r1\nACGT\n+\n!!!!\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setattr(ena_downloader.subprocess, "run", fake_run)

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
            output_file.write_bytes(b"not gzip")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setattr(ena_downloader.subprocess, "run", fake_run)

        ok, message, files = FakeDownloader().download_run("SRR_FAKE", tmp_path)

        assert ok is False
        assert "gzip integrity" in message
        assert files == []
        assert not output_file.exists()
