"""Comprehensive tests for metainformant.core.io.atomic module.

NO MOCKING - All tests use real file operations via tmp_path.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.core.io.atomic import (
    atomic_replace,
    atomic_write,
    safe_write_bytes,
    safe_write_text,
    temp_directory,
)


class TestAtomicWrite:
    """Tests for atomic_write context manager."""

    def test_successful_write_text_mode(self, tmp_path: Path) -> None:
        """Test successful atomic write in text mode."""
        target = tmp_path / "output.txt"
        content = "Hello, atomic world!"

        with atomic_write(target, mode="w") as f:
            f.write(content)

        assert target.exists()
        assert target.read_text() == content

    def test_successful_write_binary_mode(self, tmp_path: Path) -> None:
        """Test successful atomic write in binary mode."""
        target = tmp_path / "output.bin"
        content = b"\x00\x01\x02\x03\xff\xfe\xfd"

        with atomic_write(target, mode="wb") as f:
            f.write(content)

        assert target.exists()
        assert target.read_bytes() == content

    def test_nested_directory_creation(self, tmp_path: Path) -> None:
        """Test that atomic_write auto-creates nested directories."""
        target = tmp_path / "nested" / "deep" / "path" / "file.txt"
        content = "nested content"

        with atomic_write(target, mode="w") as f:
            f.write(content)

        assert target.exists()
        assert target.read_text() == content

    def test_exception_during_write_temp_cleaned_up(self, tmp_path: Path) -> None:
        """Test that temp file is cleaned up when exception occurs during write."""
        target = tmp_path / "target.txt"

        # Count files before
        files_before = list(tmp_path.iterdir())

        with pytest.raises(ValueError, match="intentional error"):
            with atomic_write(target, mode="w") as f:
                f.write("partial content")
                raise ValueError("intentional error")

        # Target should not exist
        assert not target.exists()

        # No temp files should remain
        files_after = list(tmp_path.iterdir())
        assert len(files_after) == len(files_before)

    def test_exception_does_not_modify_existing_file(self, tmp_path: Path) -> None:
        """Test that existing file is not modified if write fails."""
        target = tmp_path / "existing.txt"
        original_content = "original content"
        target.write_text(original_content)

        with pytest.raises(ValueError, match="intentional error"):
            with atomic_write(target, mode="w") as f:
                f.write("new content that should not be written")
                raise ValueError("intentional error")

        # Original content should be intact
        assert target.exists()
        assert target.read_text() == original_content

    def test_unicode_content(self, tmp_path: Path) -> None:
        """Test atomic write with unicode content."""
        target = tmp_path / "unicode.txt"
        content = "Hello 世界 🌍 Привет مرحبا"

        with atomic_write(target, mode="w", encoding="utf-8") as f:
            f.write(content)

        assert target.exists()
        assert target.read_text(encoding="utf-8") == content

    def test_custom_encoding(self, tmp_path: Path) -> None:
        """Test atomic write with custom encoding."""
        target = tmp_path / "latin1.txt"
        content = "Café résumé"

        with atomic_write(target, mode="w", encoding="latin-1") as f:
            f.write(content)

        assert target.exists()
        assert target.read_text(encoding="latin-1") == content

    def test_json_write_round_trip(self, tmp_path: Path) -> None:
        """Test atomic write with JSON data (realistic use case)."""
        target = tmp_path / "data.json"
        data = {
            "name": "METAINFORMANT",
            "values": [1, 2, 3, 4, 5],
            "nested": {"key": "value"},
        }

        with atomic_write(target, mode="w") as f:
            json.dump(data, f, indent=2)

        assert target.exists()
        with target.open() as f:
            loaded = json.load(f)
        assert loaded == data

    def test_large_content(self, tmp_path: Path) -> None:
        """Test atomic write with large content."""
        target = tmp_path / "large.txt"
        # 10MB of text
        content = "x" * (10 * 1024 * 1024)

        with atomic_write(target, mode="w") as f:
            f.write(content)

        assert target.exists()
        assert target.stat().st_size == len(content)
        assert target.read_text() == content

    def test_overwrite_existing_file(self, tmp_path: Path) -> None:
        """Test that atomic_write correctly overwrites existing files."""
        target = tmp_path / "overwrite.txt"
        target.write_text("old content")

        new_content = "new content"
        with atomic_write(target, mode="w") as f:
            f.write(new_content)

        assert target.read_text() == new_content


class TestAtomicReplace:
    """Tests for atomic_replace function."""

    def test_successful_replace(self, tmp_path: Path) -> None:
        """Test successful atomic replacement."""
        src = tmp_path / "source.txt"
        dst = tmp_path / "destination.txt"

        src.write_text("source content")
        dst.write_text("old destination content")

        atomic_replace(src, dst)

        assert not src.exists()
        assert dst.exists()
        assert dst.read_text() == "source content"

    def test_replace_creates_nested_directories(self, tmp_path: Path) -> None:
        """Test that atomic_replace creates destination directories."""
        src = tmp_path / "source.txt"
        dst = tmp_path / "nested" / "deep" / "destination.txt"

        src.write_text("content")

        atomic_replace(src, dst)

        assert not src.exists()
        assert dst.exists()
        assert dst.read_text() == "content"

    def test_replace_missing_source_raises_error(self, tmp_path: Path) -> None:
        """Test that replacing with missing source raises FileNotFoundError."""
        src = tmp_path / "nonexistent.txt"
        dst = tmp_path / "destination.txt"

        with pytest.raises(FileNotFoundError, match="Source file not found"):
            atomic_replace(src, dst)

    def test_replace_with_binary_content(self, tmp_path: Path) -> None:
        """Test atomic replace with binary content."""
        src = tmp_path / "source.bin"
        dst = tmp_path / "destination.bin"

        content = b"\x00\xff\xaa\x55" * 1000
        src.write_bytes(content)

        atomic_replace(src, dst)

        assert not src.exists()
        assert dst.exists()
        assert dst.read_bytes() == content

    def test_replace_nonexistent_destination(self, tmp_path: Path) -> None:
        """Test replacing when destination doesn't exist yet."""
        src = tmp_path / "source.txt"
        dst = tmp_path / "new_destination.txt"

        src.write_text("content")

        atomic_replace(src, dst)

        assert not src.exists()
        assert dst.exists()
        assert dst.read_text() == "content"


class TestSafeWriteText:
    """Tests for safe_write_text function."""

    def test_basic_text_write(self, tmp_path: Path) -> None:
        """Test basic text write."""
        target = tmp_path / "output.txt"
        content = "Safe text content"

        safe_write_text(target, content)

        assert target.exists()
        assert target.read_text() == content

    def test_unicode_text_write(self, tmp_path: Path) -> None:
        """Test writing unicode text."""
        target = tmp_path / "unicode.txt"
        content = "Unicode: 你好 мир 🎉 ñ é"

        safe_write_text(target, content)

        assert target.exists()
        assert target.read_text() == content

    def test_custom_encoding(self, tmp_path: Path) -> None:
        """Test safe_write_text with custom encoding."""
        target = tmp_path / "custom.txt"
        content = "Custom encoding: café"

        safe_write_text(target, content, encoding="latin-1")

        assert target.exists()
        assert target.read_text(encoding="latin-1") == content

    def test_nested_directory_creation(self, tmp_path: Path) -> None:
        """Test that safe_write_text creates nested directories."""
        target = tmp_path / "level1" / "level2" / "level3" / "file.txt"
        content = "nested"

        safe_write_text(target, content)

        assert target.exists()
        assert target.read_text() == content

    def test_multiline_content(self, tmp_path: Path) -> None:
        """Test writing multiline content."""
        target = tmp_path / "multiline.txt"
        content = "Line 1\nLine 2\nLine 3\n"

        safe_write_text(target, content)

        assert target.exists()
        assert target.read_text() == content

    def test_empty_string(self, tmp_path: Path) -> None:
        """Test writing empty string."""
        target = tmp_path / "empty.txt"
        content = ""

        safe_write_text(target, content)

        assert target.exists()
        assert target.read_text() == ""
        assert target.stat().st_size == 0


class TestSafeWriteBytes:
    """Tests for safe_write_bytes function."""

    def test_basic_bytes_write(self, tmp_path: Path) -> None:
        """Test basic binary write."""
        target = tmp_path / "output.bin"
        content = b"Binary content"

        safe_write_bytes(target, content)

        assert target.exists()
        assert target.read_bytes() == content

    def test_binary_data_with_null_bytes(self, tmp_path: Path) -> None:
        """Test writing binary data with null bytes."""
        target = tmp_path / "nulls.bin"
        content = b"\x00\x01\x02\x00\xff\xfe\x00"

        safe_write_bytes(target, content)

        assert target.exists()
        assert target.read_bytes() == content

    def test_large_binary_content(self, tmp_path: Path) -> None:
        """Test writing large binary content."""
        target = tmp_path / "large.bin"
        # 5MB of binary data
        content = bytes(range(256)) * (5 * 1024 * 1024 // 256)

        safe_write_bytes(target, content)

        assert target.exists()
        assert target.stat().st_size == len(content)
        assert target.read_bytes() == content

    def test_empty_bytes(self, tmp_path: Path) -> None:
        """Test writing empty bytes."""
        target = tmp_path / "empty.bin"
        content = b""

        safe_write_bytes(target, content)

        assert target.exists()
        assert target.read_bytes() == b""
        assert target.stat().st_size == 0

    def test_nested_directory_creation(self, tmp_path: Path) -> None:
        """Test that safe_write_bytes creates nested directories."""
        target = tmp_path / "dir1" / "dir2" / "file.bin"
        content = b"\xde\xad\xbe\xef"

        safe_write_bytes(target, content)

        assert target.exists()
        assert target.read_bytes() == content


class TestTempDirectory:
    """Tests for temp_directory context manager."""

    def test_temp_directory_exists_during_context(self, tmp_path: Path) -> None:
        """Test that temp directory exists during context."""
        temp_dir_path: Path | None = None

        with temp_directory() as tmp:
            temp_dir_path = tmp
            assert tmp.exists()
            assert tmp.is_dir()
            # Write a file to verify it's a real directory
            test_file = tmp / "test.txt"
            test_file.write_text("test content")
            assert test_file.exists()

        # Directory should be cleaned up after context
        assert temp_dir_path is not None
        assert not temp_dir_path.exists()

    def test_temp_directory_cleanup_default(self, tmp_path: Path) -> None:
        """Test that temp directory is cleaned up by default."""
        temp_dir_path: Path | None = None

        with temp_directory() as tmp:
            temp_dir_path = tmp
            # Create multiple files
            (tmp / "file1.txt").write_text("content1")
            (tmp / "file2.txt").write_text("content2")
            (tmp / "subdir").mkdir()
            (tmp / "subdir" / "file3.txt").write_text("content3")

        assert temp_dir_path is not None
        assert not temp_dir_path.exists()

    def test_temp_directory_no_cleanup(self, tmp_path: Path) -> None:
        """Test that temp directory is retained with cleanup=False."""
        temp_dir_path: Path | None = None

        with temp_directory(cleanup=False) as tmp:
            temp_dir_path = tmp
            test_file = tmp / "persistent.txt"
            test_file.write_text("persistent content")

        # Directory should still exist
        assert temp_dir_path is not None
        assert temp_dir_path.exists()
        assert temp_dir_path.is_dir()
        assert (temp_dir_path / "persistent.txt").read_text() == "persistent content"

        # Manual cleanup
        import shutil

        shutil.rmtree(temp_dir_path)

    def test_temp_directory_custom_prefix(self, tmp_path: Path) -> None:
        """Test temp directory with custom prefix."""
        with temp_directory(prefix="custom_test_") as tmp:
            assert tmp.exists()
            assert tmp.name.startswith("custom_test_")

    def test_temp_directory_nested_structure(self, tmp_path: Path) -> None:
        """Test creating nested structure in temp directory."""
        with temp_directory() as tmp:
            nested = tmp / "level1" / "level2" / "level3"
            nested.mkdir(parents=True)
            test_file = nested / "deep.txt"
            test_file.write_text("deep content")

            assert test_file.exists()
            assert test_file.read_text() == "deep content"

    def test_temp_directory_exception_still_cleans_up(self, tmp_path: Path) -> None:
        """Test that temp directory is cleaned up even when exception occurs."""
        temp_dir_path: Path | None = None

        with pytest.raises(ValueError, match="intentional"):
            with temp_directory() as tmp:
                temp_dir_path = tmp
                (tmp / "file.txt").write_text("content")
                raise ValueError("intentional")

        assert temp_dir_path is not None
        assert not temp_dir_path.exists()

    def test_temp_directory_multiple_files_and_dirs(self, tmp_path: Path) -> None:
        """Test temp directory with complex structure."""
        with temp_directory() as tmp:
            # Create files
            (tmp / "file1.txt").write_text("content1")
            (tmp / "file2.json").write_text('{"key": "value"}')

            # Create subdirectories with files
            subdir1 = tmp / "subdir1"
            subdir1.mkdir()
            (subdir1 / "nested1.txt").write_text("nested1")

            subdir2 = tmp / "subdir2"
            subdir2.mkdir()
            (subdir2 / "nested2.bin").write_bytes(b"\x00\xff")

            # Verify all exist
            assert (tmp / "file1.txt").exists()
            assert (tmp / "file2.json").exists()
            assert (subdir1 / "nested1.txt").exists()
            assert (subdir2 / "nested2.bin").exists()


class TestIntegrationScenarios:
    """Integration tests combining multiple atomic operations."""

    def test_atomic_write_then_replace(self, tmp_path: Path) -> None:
        """Test writing atomically then replacing with another file."""
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        final = tmp_path / "final.txt"

        # Write two files atomically
        with atomic_write(file1, mode="w") as f:
            f.write("content from file1")

        with atomic_write(file2, mode="w") as f:
            f.write("content from file2")

        # Replace final with file1
        atomic_replace(file1, final)
        assert final.read_text() == "content from file1"

        # Replace final with file2
        atomic_replace(file2, final)
        assert final.read_text() == "content from file2"

    def test_safe_write_in_temp_directory(self, tmp_path: Path) -> None:
        """Test using safe_write functions within temp directory."""
        with temp_directory() as tmp:
            text_file = tmp / "text.txt"
            binary_file = tmp / "binary.bin"

            safe_write_text(text_file, "temp text content")
            safe_write_bytes(binary_file, b"\xde\xad\xbe\xef")

            assert text_file.read_text() == "temp text content"
            assert binary_file.read_bytes() == b"\xde\xad\xbe\xef"

    def test_atomic_write_workflow(self, tmp_path: Path) -> None:
        """Test realistic workflow with atomic operations."""
        # Simulate processing pipeline
        input_file = tmp_path / "input.json"
        input_data = {"values": [1, 2, 3, 4, 5]}
        safe_write_text(input_file, json.dumps(input_data))

        # Process: read, transform, write atomically
        with input_file.open() as f:
            data = json.load(f)

        processed_data = {"values": [x * 2 for x in data["values"]]}

        output_file = tmp_path / "output.json"
        with atomic_write(output_file, mode="w") as f:
            json.dump(processed_data, f, indent=2)

        # Verify
        with output_file.open() as f:
            result = json.load(f)

        assert result == {"values": [2, 4, 6, 8, 10]}

    def test_error_recovery_preserves_original(self, tmp_path: Path) -> None:
        """Test that original file is preserved when update fails."""
        original_file = tmp_path / "important.json"
        original_data = {"version": 1, "data": "original"}
        safe_write_text(original_file, json.dumps(original_data))

        # Attempt to update but fail
        with pytest.raises(ValueError):
            with atomic_write(original_file, mode="w") as f:
                f.write('{"version": 2, "data": "corrupted"')
                raise ValueError("Simulated corruption")

        # Original should be intact
        with original_file.open() as f:
            recovered = json.load(f)

        assert recovered == original_data
