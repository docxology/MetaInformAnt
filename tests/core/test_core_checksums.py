"""Comprehensive tests for metainformant.core.io.checksums module.

Tests cover:
- MD5 and SHA256 computation with known values
- Checksum verification (match and mismatch)
- Sidecar file creation and verification
- Batch processing with mixed existing/missing files
- Edge cases: empty files, large files, missing files, invalid algorithms
- Chunk size handling for large files

NO mocking - all tests use real files via tmp_path fixture.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from metainformant.core.io.checksums import (
    compute_checksums_batch,
    compute_md5,
    compute_sha256,
    verify_checksum,
    verify_checksum_file,
    write_checksum_file,
)


class TestBasicChecksum:
    """Test basic MD5 and SHA256 computation with known values."""

    def test_compute_md5_known_value(self, tmp_path: Path) -> None:
        """Test MD5 computation against known hash."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"hello world\n")

        result = compute_md5(test_file)

        # Known MD5 for "hello world\n"
        expected = "6f5902ac237024bdd0c176cb93063dc4"
        assert result == expected
        assert len(result) == 32

    def test_compute_sha256_known_value(self, tmp_path: Path) -> None:
        """Test SHA256 computation against known hash."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"hello world\n")

        result = compute_sha256(test_file)

        # Compute expected SHA256 for "hello world\n"
        expected = hashlib.sha256(b"hello world\n").hexdigest()
        assert result == expected
        assert len(result) == 64

    def test_compute_md5_empty_file(self, tmp_path: Path) -> None:
        """Test MD5 of empty file."""
        empty_file = tmp_path / "empty.txt"
        empty_file.touch()

        result = compute_md5(empty_file)

        # Known MD5 for empty string
        expected = "d41d8cd98f00b204e9800998ecf8427e"
        assert result == expected

    def test_compute_sha256_empty_file(self, tmp_path: Path) -> None:
        """Test SHA256 of empty file."""
        empty_file = tmp_path / "empty.txt"
        empty_file.touch()

        result = compute_sha256(empty_file)

        # Known SHA256 for empty string
        expected = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        assert result == expected

    def test_compute_md5_binary_file(self, tmp_path: Path) -> None:
        """Test MD5 computation on binary data."""
        binary_file = tmp_path / "binary.dat"
        binary_data = bytes(range(256))
        binary_file.write_bytes(binary_data)

        result = compute_md5(binary_file)

        # Verify against direct hashlib computation
        expected = hashlib.md5(binary_data).hexdigest()
        assert result == expected

    def test_compute_sha256_binary_file(self, tmp_path: Path) -> None:
        """Test SHA256 computation on binary data."""
        binary_file = tmp_path / "binary.dat"
        binary_data = bytes(range(256))
        binary_file.write_bytes(binary_data)

        result = compute_sha256(binary_file)

        # Verify against direct hashlib computation
        expected = hashlib.sha256(binary_data).hexdigest()
        assert result == expected

    def test_compute_md5_missing_file(self, tmp_path: Path) -> None:
        """Test MD5 computation raises FileNotFoundError for missing file."""
        missing_file = tmp_path / "missing.txt"

        with pytest.raises(FileNotFoundError, match="File not found"):
            compute_md5(missing_file)

    def test_compute_sha256_missing_file(self, tmp_path: Path) -> None:
        """Test SHA256 computation raises FileNotFoundError for missing file."""
        missing_file = tmp_path / "missing.txt"

        with pytest.raises(FileNotFoundError, match="File not found"):
            compute_sha256(missing_file)

    def test_compute_md5_directory(self, tmp_path: Path) -> None:
        """Test MD5 computation raises FileNotFoundError for directory."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()

        with pytest.raises(FileNotFoundError, match="File not found"):
            compute_md5(subdir)


class TestLargeFiles:
    """Test checksum computation on large files with different chunk sizes."""

    def test_compute_md5_large_file_default_chunk(self, tmp_path: Path) -> None:
        """Test MD5 computation on large file with default chunk size."""
        large_file = tmp_path / "large.bin"
        # Create 1MB file
        data = b"A" * (1024 * 1024)
        large_file.write_bytes(data)

        result = compute_md5(large_file)

        expected = hashlib.md5(data).hexdigest()
        assert result == expected

    def test_compute_sha256_large_file_default_chunk(self, tmp_path: Path) -> None:
        """Test SHA256 computation on large file with default chunk size."""
        large_file = tmp_path / "large.bin"
        # Create 1MB file
        data = b"B" * (1024 * 1024)
        large_file.write_bytes(data)

        result = compute_sha256(large_file)

        expected = hashlib.sha256(data).hexdigest()
        assert result == expected

    def test_compute_md5_large_file_small_chunk(self, tmp_path: Path) -> None:
        """Test MD5 computation with small chunk size."""
        large_file = tmp_path / "large.bin"
        data = b"C" * (100 * 1024)  # 100KB
        large_file.write_bytes(data)

        result = compute_md5(large_file, chunk_size=1024)

        expected = hashlib.md5(data).hexdigest()
        assert result == expected

    def test_compute_sha256_large_file_small_chunk(self, tmp_path: Path) -> None:
        """Test SHA256 computation with small chunk size."""
        large_file = tmp_path / "large.bin"
        data = b"D" * (100 * 1024)  # 100KB
        large_file.write_bytes(data)

        result = compute_sha256(large_file, chunk_size=1024)

        expected = hashlib.sha256(data).hexdigest()
        assert result == expected

    def test_compute_md5_large_file_large_chunk(self, tmp_path: Path) -> None:
        """Test MD5 computation with large chunk size."""
        large_file = tmp_path / "large.bin"
        data = b"E" * (100 * 1024)  # 100KB
        large_file.write_bytes(data)

        result = compute_md5(large_file, chunk_size=64 * 1024)

        expected = hashlib.md5(data).hexdigest()
        assert result == expected


class TestVerifyChecksum:
    """Test checksum verification against expected values."""

    def test_verify_checksum_md5_match(self, tmp_path: Path) -> None:
        """Test successful MD5 verification."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected = hashlib.md5(b"test data\n").hexdigest()
        result = verify_checksum(test_file, expected, algorithm="md5")

        assert result is True

    def test_verify_checksum_sha256_match(self, tmp_path: Path) -> None:
        """Test successful SHA256 verification."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected = hashlib.sha256(b"test data\n").hexdigest()
        result = verify_checksum(test_file, expected, algorithm="sha256")

        assert result is True

    def test_verify_checksum_md5_mismatch(self, tmp_path: Path) -> None:
        """Test MD5 verification returns False on mismatch."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        wrong_hash = "0" * 32
        result = verify_checksum(test_file, wrong_hash, algorithm="md5")

        assert result is False

    def test_verify_checksum_sha256_mismatch(self, tmp_path: Path) -> None:
        """Test SHA256 verification returns False on mismatch."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        wrong_hash = "0" * 64
        result = verify_checksum(test_file, wrong_hash, algorithm="sha256")

        assert result is False

    def test_verify_checksum_case_insensitive(self, tmp_path: Path) -> None:
        """Test checksum verification is case-insensitive."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected_lower = hashlib.sha256(b"test data\n").hexdigest()
        expected_upper = expected_lower.upper()

        assert verify_checksum(test_file, expected_lower, algorithm="sha256")
        assert verify_checksum(test_file, expected_upper, algorithm="sha256")

    def test_verify_checksum_with_whitespace(self, tmp_path: Path) -> None:
        """Test checksum verification strips whitespace from expected value."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected = hashlib.sha256(b"test data\n").hexdigest()
        expected_with_spaces = f"  {expected}  \n"

        result = verify_checksum(test_file, expected_with_spaces, algorithm="sha256")

        assert result is True

    def test_verify_checksum_sha1(self, tmp_path: Path) -> None:
        """Test verification with SHA1 algorithm."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected = hashlib.sha1(b"test data\n").hexdigest()
        result = verify_checksum(test_file, expected, algorithm="sha1")

        assert result is True

    def test_verify_checksum_sha512(self, tmp_path: Path) -> None:
        """Test verification with SHA512 algorithm."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        expected = hashlib.sha512(b"test data\n").hexdigest()
        result = verify_checksum(test_file, expected, algorithm="sha512")

        assert result is True

    def test_verify_checksum_invalid_algorithm(self, tmp_path: Path) -> None:
        """Test verification raises ValueError for unsupported algorithm."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        with pytest.raises(ValueError, match="Unsupported algorithm"):
            verify_checksum(test_file, "somehash", algorithm="invalid")

    def test_verify_checksum_missing_file(self, tmp_path: Path) -> None:
        """Test verification raises FileNotFoundError for missing file."""
        missing_file = tmp_path / "missing.txt"

        with pytest.raises(FileNotFoundError, match="File not found"):
            verify_checksum(missing_file, "somehash", algorithm="sha256")


class TestSidecarFiles:
    """Test sidecar checksum file creation and verification."""

    def test_write_checksum_file_sha256(self, tmp_path: Path) -> None:
        """Test writing SHA256 sidecar file."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        sidecar_path = write_checksum_file(test_file, algorithm="sha256")

        assert sidecar_path.exists()
        assert sidecar_path.name == "test.txt.sha256"
        assert sidecar_path.parent == tmp_path

        # Verify sidecar content format
        content = sidecar_path.read_text(encoding="utf-8")
        expected_hash = hashlib.sha256(b"test data\n").hexdigest()
        assert content == f"{expected_hash}  test.txt\n"

    def test_write_checksum_file_md5(self, tmp_path: Path) -> None:
        """Test writing MD5 sidecar file."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"test data\n")

        sidecar_path = write_checksum_file(test_file, algorithm="md5")

        assert sidecar_path.exists()
        assert sidecar_path.name == "test.txt.md5"

        # Verify sidecar content format
        content = sidecar_path.read_text(encoding="utf-8")
        expected_hash = hashlib.md5(b"test data\n").hexdigest()
        assert content == f"{expected_hash}  test.txt\n"

    def test_write_checksum_file_nested_path(self, tmp_path: Path) -> None:
        """Test writing sidecar file in nested directory."""
        subdir = tmp_path / "subdir" / "nested"
        subdir.mkdir(parents=True)
        test_file = subdir / "test.txt"
        test_file.write_bytes(b"nested data\n")

        sidecar_path = write_checksum_file(test_file, algorithm="sha256")

        assert sidecar_path.exists()
        assert sidecar_path.parent == subdir
        assert sidecar_path.name == "test.txt.sha256"

    def test_verify_checksum_file_sha256(self, tmp_path: Path) -> None:
        """Test verifying file against SHA256 sidecar."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"verify me\n")

        # Create sidecar file
        write_checksum_file(test_file, algorithm="sha256")

        # Verify
        result = verify_checksum_file(test_file)

        assert result is True

    def test_verify_checksum_file_md5(self, tmp_path: Path) -> None:
        """Test verifying file against MD5 sidecar."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"verify me\n")

        # Create sidecar file
        write_checksum_file(test_file, algorithm="md5")

        # Verify
        result = verify_checksum_file(test_file)

        assert result is True

    def test_verify_checksum_file_prefers_sha256(self, tmp_path: Path) -> None:
        """Test that verification prefers SHA256 when both sidecars exist."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"multi hash\n")

        # Create both sidecar files
        write_checksum_file(test_file, algorithm="sha256")
        write_checksum_file(test_file, algorithm="md5")

        # Modify MD5 sidecar to have wrong hash
        md5_sidecar = tmp_path / "test.txt.md5"
        md5_sidecar.write_text("0" * 32 + "  test.txt\n", encoding="utf-8")

        # Should still pass because SHA256 is checked first
        result = verify_checksum_file(test_file)

        assert result is True

    def test_verify_checksum_file_corrupted(self, tmp_path: Path) -> None:
        """Test verification fails when file is corrupted."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"original data\n")

        # Create sidecar
        write_checksum_file(test_file, algorithm="sha256")

        # Corrupt file
        test_file.write_bytes(b"corrupted data\n")

        # Verify should fail
        result = verify_checksum_file(test_file)

        assert result is False

    def test_verify_checksum_file_no_sidecar(self, tmp_path: Path) -> None:
        """Test verification raises FileNotFoundError when no sidecar exists."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"no sidecar\n")

        with pytest.raises(FileNotFoundError, match="No checksum sidecar file found"):
            verify_checksum_file(test_file)

    def test_verify_checksum_file_hash_only_format(self, tmp_path: Path) -> None:
        """Test verification works with hash-only sidecar format."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"simple hash\n")

        # Create sidecar with hash-only format (no filename)
        sidecar = tmp_path / "test.txt.sha256"
        expected_hash = hashlib.sha256(b"simple hash\n").hexdigest()
        sidecar.write_text(expected_hash, encoding="utf-8")

        result = verify_checksum_file(test_file)

        assert result is True

    def test_write_checksum_file_overwrites_existing(self, tmp_path: Path) -> None:
        """Test writing sidecar file overwrites existing sidecar."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"original\n")

        # Create first sidecar
        sidecar_path = write_checksum_file(test_file, algorithm="sha256")
        original_content = sidecar_path.read_text(encoding="utf-8")

        # Modify file
        test_file.write_bytes(b"modified\n")

        # Write new sidecar
        new_sidecar_path = write_checksum_file(test_file, algorithm="sha256")
        new_content = new_sidecar_path.read_text(encoding="utf-8")

        assert new_sidecar_path == sidecar_path
        assert new_content != original_content


class TestBatchChecksum:
    """Test batch checksum computation on multiple files."""

    def test_compute_checksums_batch_single_file(self, tmp_path: Path) -> None:
        """Test batch computation with single file."""
        test_file = tmp_path / "test.txt"
        test_file.write_bytes(b"batch test\n")

        result = compute_checksums_batch([test_file], algorithm="sha256")

        assert len(result) == 1
        assert str(test_file) in result
        expected = hashlib.sha256(b"batch test\n").hexdigest()
        assert result[str(test_file)] == expected

    def test_compute_checksums_batch_multiple_files(self, tmp_path: Path) -> None:
        """Test batch computation with multiple files."""
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file3 = tmp_path / "file3.txt"

        file1.write_bytes(b"data1\n")
        file2.write_bytes(b"data2\n")
        file3.write_bytes(b"data3\n")

        result = compute_checksums_batch([file1, file2, file3], algorithm="sha256")

        assert len(result) == 3
        assert str(file1) in result
        assert str(file2) in result
        assert str(file3) in result

        # Verify each hash
        assert result[str(file1)] == hashlib.sha256(b"data1\n").hexdigest()
        assert result[str(file2)] == hashlib.sha256(b"data2\n").hexdigest()
        assert result[str(file3)] == hashlib.sha256(b"data3\n").hexdigest()

    def test_compute_checksums_batch_md5(self, tmp_path: Path) -> None:
        """Test batch computation with MD5 algorithm."""
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"

        file1.write_bytes(b"md5 test 1\n")
        file2.write_bytes(b"md5 test 2\n")

        result = compute_checksums_batch([file1, file2], algorithm="md5")

        assert len(result) == 2
        assert result[str(file1)] == hashlib.md5(b"md5 test 1\n").hexdigest()
        assert result[str(file2)] == hashlib.md5(b"md5 test 2\n").hexdigest()

    def test_compute_checksums_batch_mixed_existing_missing(self, tmp_path: Path) -> None:
        """Test batch computation skips missing files."""
        file1 = tmp_path / "exists1.txt"
        file2 = tmp_path / "missing.txt"
        file3 = tmp_path / "exists2.txt"

        file1.write_bytes(b"exists 1\n")
        file3.write_bytes(b"exists 2\n")
        # file2 intentionally not created

        result = compute_checksums_batch([file1, file2, file3], algorithm="sha256")

        assert len(result) == 2
        assert str(file1) in result
        assert str(file2) not in result
        assert str(file3) in result

        assert result[str(file1)] == hashlib.sha256(b"exists 1\n").hexdigest()
        assert result[str(file3)] == hashlib.sha256(b"exists 2\n").hexdigest()

    def test_compute_checksums_batch_all_missing(self, tmp_path: Path) -> None:
        """Test batch computation returns empty dict when all files missing."""
        file1 = tmp_path / "missing1.txt"
        file2 = tmp_path / "missing2.txt"

        result = compute_checksums_batch([file1, file2], algorithm="sha256")

        assert result == {}

    def test_compute_checksums_batch_empty_list(self, tmp_path: Path) -> None:
        """Test batch computation with empty file list."""
        result = compute_checksums_batch([], algorithm="sha256")

        assert result == {}

    def test_compute_checksums_batch_nested_paths(self, tmp_path: Path) -> None:
        """Test batch computation with files in nested directories."""
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()

        file1 = dir1 / "file1.txt"
        file2 = dir2 / "file2.txt"

        file1.write_bytes(b"nested 1\n")
        file2.write_bytes(b"nested 2\n")

        result = compute_checksums_batch([file1, file2], algorithm="sha256")

        assert len(result) == 2
        assert str(file1) in result
        assert str(file2) in result


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_unicode_filename(self, tmp_path: Path) -> None:
        """Test checksum computation on file with unicode name."""
        test_file = tmp_path / "тест_файл_🐍.txt"
        test_file.write_bytes(b"unicode filename test\n")

        md5_result = compute_md5(test_file)
        sha256_result = compute_sha256(test_file)

        assert len(md5_result) == 32
        assert len(sha256_result) == 64

    def test_unicode_content(self, tmp_path: Path) -> None:
        """Test checksum computation on file with unicode content."""
        test_file = tmp_path / "unicode.txt"
        unicode_content = "Hello 世界 🌍\n"
        test_file.write_text(unicode_content, encoding="utf-8")

        result = compute_sha256(test_file)

        expected = hashlib.sha256(unicode_content.encode("utf-8")).hexdigest()
        assert result == expected

    def test_very_long_filename(self, tmp_path: Path) -> None:
        """Test checksum on file with very long name."""
        long_name = "a" * 200 + ".txt"
        test_file = tmp_path / long_name
        test_file.write_bytes(b"long filename\n")

        result = compute_md5(test_file)

        assert len(result) == 32

    def test_file_with_no_extension(self, tmp_path: Path) -> None:
        """Test sidecar file creation for file without extension."""
        test_file = tmp_path / "noext"
        test_file.write_bytes(b"no extension\n")

        sidecar_path = write_checksum_file(test_file, algorithm="sha256")

        assert sidecar_path.name == "noext.sha256"
        assert verify_checksum_file(test_file) is True

    def test_hidden_file(self, tmp_path: Path) -> None:
        """Test checksum on hidden file (starts with dot)."""
        hidden_file = tmp_path / ".hidden"
        hidden_file.write_bytes(b"hidden content\n")

        result = compute_sha256(hidden_file)

        expected = hashlib.sha256(b"hidden content\n").hexdigest()
        assert result == expected

    def test_symlink_to_file(self, tmp_path: Path) -> None:
        """Test checksum computation follows symlinks."""
        real_file = tmp_path / "real.txt"
        real_file.write_bytes(b"symlink test\n")

        link_file = tmp_path / "link.txt"
        link_file.symlink_to(real_file)

        result = compute_sha256(link_file)

        expected = hashlib.sha256(b"symlink test\n").hexdigest()
        assert result == expected

    def test_newline_variations(self, tmp_path: Path) -> None:
        """Test that different newline styles produce different hashes."""
        file_lf = tmp_path / "lf.txt"
        file_crlf = tmp_path / "crlf.txt"
        file_cr = tmp_path / "cr.txt"

        file_lf.write_bytes(b"line1\nline2\n")
        file_crlf.write_bytes(b"line1\r\nline2\r\n")
        file_cr.write_bytes(b"line1\rline2\r")

        hash_lf = compute_sha256(file_lf)
        hash_crlf = compute_sha256(file_crlf)
        hash_cr = compute_sha256(file_cr)

        # All should be different
        assert hash_lf != hash_crlf
        assert hash_lf != hash_cr
        assert hash_crlf != hash_cr

    def test_zero_byte_file(self, tmp_path: Path) -> None:
        """Test checksum on zero-byte file (same as empty file)."""
        zero_file = tmp_path / "zero.bin"
        zero_file.touch()

        md5_result = compute_md5(zero_file)
        sha256_result = compute_sha256(zero_file)

        assert md5_result == "d41d8cd98f00b204e9800998ecf8427e"
        assert sha256_result == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
