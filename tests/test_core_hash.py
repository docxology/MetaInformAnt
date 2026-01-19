"""Comprehensive tests for core.utils.hash module."""
from __future__ import annotations

from pathlib import Path
import tempfile

from metainformant.core.utils import hash as core_hash


class TestSha256Functions:
    """Tests for SHA256 hashing functions."""

    def test_sha256_bytes_basic(self) -> None:
        """Test sha256_bytes with basic input."""
        content = b"hello world"
        h = core_hash.sha256_bytes(content)
        assert len(h) == 64
        assert h.isalnum()
        # Verify determinism
        assert core_hash.sha256_bytes(content) == h

    def test_sha256_bytes_empty(self) -> None:
        """Test sha256_bytes with empty input."""
        h = core_hash.sha256_bytes(b"")
        assert len(h) == 64
        # SHA256 of empty string is well-known
        assert h == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

    def test_sha256_file_matches_bytes(self, tmp_path: Path) -> None:
        """Test sha256_file returns same hash as sha256_bytes for same content."""
        content = b"hello world\n"
        h_bytes = core_hash.sha256_bytes(content)
        
        path = tmp_path / "file.txt"
        path.write_bytes(content)
        h_file = core_hash.sha256_file(path)
        
        assert h_bytes == h_file

    def test_sha256_file_large_file(self, tmp_path: Path) -> None:
        """Test sha256_file handles large files correctly."""
        # Create a file larger than chunk size
        content = b"x" * (2 * 1024 * 1024)  # 2MB
        path = tmp_path / "large.bin"
        path.write_bytes(content)
        
        h = core_hash.sha256_file(path, chunk_size=512 * 1024)
        assert len(h) == 64

    def test_sha256_string_basic(self) -> None:
        """Test sha256_string function."""
        h = core_hash.sha256_string("hello")
        assert len(h) == 64
        # Should match bytes version
        assert h == core_hash.sha256_bytes(b"hello")

    def test_sha256_string_unicode(self) -> None:
        """Test sha256_string with unicode characters."""
        h = core_hash.sha256_string("héllo wörld 🌍")
        assert len(h) == 64


class TestDeterministicSeed:
    """Tests for deterministic_seed function."""

    def test_deterministic_seed_reproducible(self) -> None:
        """Test that same input produces same seed."""
        seed1 = core_hash.deterministic_seed("test_string")
        seed2 = core_hash.deterministic_seed("test_string")
        assert seed1 == seed2

    def test_deterministic_seed_different_inputs(self) -> None:
        """Test that different inputs produce different seeds."""
        seed1 = core_hash.deterministic_seed("input_a")
        seed2 = core_hash.deterministic_seed("input_b")
        assert seed1 != seed2

    def test_deterministic_seed_valid_range(self) -> None:
        """Test that seed is in valid range for random module."""
        for i in range(100):
            seed = core_hash.deterministic_seed(f"test_{i}")
            assert 0 <= seed < 2**31 - 1


class TestFileHashComparison:
    """Tests for file_hash_comparison function."""

    def test_file_hash_comparison_identical(self, tmp_path: Path) -> None:
        """Test comparison of identical files."""
        content = b"same content"
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_bytes(content)
        file2.write_bytes(content)
        
        assert core_hash.file_hash_comparison(file1, file2) is True

    def test_file_hash_comparison_different(self, tmp_path: Path) -> None:
        """Test comparison of different files."""
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_bytes(b"content a")
        file2.write_bytes(b"content b")
        
        assert core_hash.file_hash_comparison(file1, file2) is False


class TestHashDirectory:
    """Tests for hash_directory function."""

    def test_hash_directory_basic(self, tmp_path: Path) -> None:
        """Test hashing all files in a directory."""
        (tmp_path / "file1.txt").write_text("content1")
        (tmp_path / "file2.txt").write_text("content2")
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "file3.txt").write_text("content3")
        
        hashes = core_hash.hash_directory(tmp_path)
        
        assert "file1.txt" in hashes
        assert "file2.txt" in hashes
        assert "subdir/file3.txt" in hashes
        assert len(hashes) == 3

    def test_hash_directory_with_pattern(self, tmp_path: Path) -> None:
        """Test hashing with glob pattern."""
        (tmp_path / "file1.txt").write_text("content1")
        (tmp_path / "file2.py").write_text("content2")
        
        hashes = core_hash.hash_directory(tmp_path, pattern="*.txt")
        
        assert "file1.txt" in hashes
        assert "file2.py" not in hashes

    def test_hash_directory_empty(self, tmp_path: Path) -> None:
        """Test hashing empty directory."""
        hashes = core_hash.hash_directory(tmp_path)
        assert hashes == {}


class TestVerifyFileIntegrity:
    """Tests for verify_file_integrity function."""

    def test_verify_file_integrity_valid(self, tmp_path: Path) -> None:
        """Test verification with correct hash."""
        content = b"test content"
        path = tmp_path / "file.txt"
        path.write_bytes(content)
        expected_hash = core_hash.sha256_bytes(content)
        
        assert core_hash.verify_file_integrity(path, expected_hash) is True

    def test_verify_file_integrity_invalid(self, tmp_path: Path) -> None:
        """Test verification with incorrect hash."""
        path = tmp_path / "file.txt"
        path.write_bytes(b"actual content")
        wrong_hash = "a" * 64
        
        assert core_hash.verify_file_integrity(path, wrong_hash) is False

    def test_verify_file_integrity_missing_file(self, tmp_path: Path) -> None:
        """Test verification of non-existent file."""
        path = tmp_path / "nonexistent.txt"
        
        assert core_hash.verify_file_integrity(path, "a" * 64) is False
