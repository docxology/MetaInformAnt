"""Comprehensive tests for core.io.paths module."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile

from metainformant.core.io import paths as core_paths


class TestExpandAndResolve:
    """Tests for expand_and_resolve function."""

    def test_expand_and_resolve_user(self, tmp_path: Path) -> None:
        """Test expansion of ~ to home directory."""
        home = tmp_path / "home"
        home.mkdir(parents=True, exist_ok=True)
        old_home = os.environ.get("HOME")
        os.environ["HOME"] = str(home)
        try:
            p = core_paths.expand_and_resolve("~/data")
            assert str(p).startswith(str(home))
        finally:
            if old_home is None:
                os.environ.pop("HOME", None)
            else:
                os.environ["HOME"] = old_home

    def test_expand_and_resolve_absolute(self) -> None:
        """Test that absolute path is returned unchanged but resolved."""
        p = core_paths.expand_and_resolve("/tmp/test")
        assert p.is_absolute()
        assert str(p).startswith("/tmp") or str(p).startswith("/private/tmp")

    def test_expand_and_resolve_relative(self) -> None:
        """Test that relative path becomes absolute."""
        p = core_paths.expand_and_resolve("relative/path")
        assert p.is_absolute()


class TestIsWithin:
    """Tests for is_within function."""

    def test_is_within_subdirectory(self, tmp_path: Path) -> None:
        """Test that subdirectory is within parent."""
        root = tmp_path / "root"
        sub = root / "a" / "b"
        sub.mkdir(parents=True, exist_ok=True)
        assert core_paths.is_within(sub, root)

    def test_is_within_not_contained(self, tmp_path: Path) -> None:
        """Test that parent is not within child."""
        root = tmp_path / "root"
        sub = root / "a"
        sub.mkdir(parents=True, exist_ok=True)
        assert not core_paths.is_within(root, sub)

    def test_is_within_sibling(self, tmp_path: Path) -> None:
        """Test that sibling directories are not within each other."""
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()
        assert not core_paths.is_within(dir1, dir2)


class TestEnsureDirectory:
    """Tests for ensure_directory function."""

    def test_ensure_directory_creates(self, tmp_path: Path) -> None:
        """Test creating new directory."""
        new_dir = tmp_path / "new" / "nested" / "dir"
        core_paths.ensure_directory(new_dir)
        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_ensure_directory_existing(self, tmp_path: Path) -> None:
        """Test with existing directory (should not error)."""
        existing = tmp_path / "existing"
        existing.mkdir()
        core_paths.ensure_directory(existing)  # Should not raise
        assert existing.exists()


class TestPrepareFilePath:
    """Tests for prepare_file_path function."""

    def test_prepare_file_path_creates_parent(self, tmp_path: Path) -> None:
        """Test that parent directories are created."""
        file_path = tmp_path / "a" / "b" / "file.txt"
        core_paths.prepare_file_path(file_path)
        assert file_path.parent.exists()


class TestIsSafePath:
    """Tests for is_safe_path function."""

    def test_is_safe_path_normal(self) -> None:
        """Test normal paths are safe."""
        assert core_paths.is_safe_path("data/file.txt")
        assert core_paths.is_safe_path("output/results.json")

    def test_is_safe_path_traversal(self) -> None:
        """Test path traversal attempts are blocked."""
        assert not core_paths.is_safe_path("../etc/passwd")
        assert not core_paths.is_safe_path("data/../../secret")


class TestGetFileExtension:
    """Tests for get_file_extension function."""

    def test_get_file_extension_simple(self) -> None:
        """Test simple extension."""
        assert core_paths.get_file_extension("file.txt") == ".txt"

    def test_get_file_extension_double(self) -> None:
        """Test double extension."""
        ext = core_paths.get_file_extension("file.tar.gz")
        assert ext == ".gz"

    def test_get_file_extension_none(self) -> None:
        """Test file without extension."""
        assert core_paths.get_file_extension("README") == ""


class TestChangeExtension:
    """Tests for change_extension function."""

    def test_change_extension_basic(self) -> None:
        """Test basic extension change."""
        result = core_paths.change_extension("file.txt", ".json")
        assert str(result) == "file.json"  # Returns Path

    def test_change_extension_without_dot(self) -> None:
        """Test extension without leading dot."""
        result = core_paths.change_extension("file.txt", "json")
        assert str(result) == "file.json"  # Returns Path


class TestFindFilesByExtension:
    """Tests for find_files_by_extension function."""

    def test_find_files_by_extension(self, tmp_path: Path) -> None:
        """Test finding files by extension."""
        (tmp_path / "a.txt").touch()
        (tmp_path / "b.txt").touch()
        (tmp_path / "c.py").touch()
        
        results = core_paths.find_files_by_extension(tmp_path, ".txt")
        
        assert len(results) == 2
        assert all(str(r).endswith(".txt") for r in results)


class TestGetFileSize:
    """Tests for get_file_size function."""

    def test_get_file_size_basic(self, tmp_path: Path) -> None:
        """Test file size calculation."""
        file = tmp_path / "test.txt"
        content = b"hello world"
        file.write_bytes(content)
        
        size = core_paths.get_file_size(file)
        assert size == len(content)

    def test_get_file_size_missing(self, tmp_path: Path) -> None:
        """Test missing file returns 0."""
        missing = tmp_path / "nonexistent.txt"
        assert core_paths.get_file_size(missing) == 0


class TestGetDirectorySize:
    """Tests for get_directory_size function."""

    def test_get_directory_size(self, tmp_path: Path) -> None:
        """Test directory size calculation."""
        (tmp_path / "a.txt").write_bytes(b"a" * 100)
        (tmp_path / "b.txt").write_bytes(b"b" * 200)
        
        size = core_paths.get_directory_size(tmp_path)
        assert size == 300


class TestSanitizeFilename:
    """Tests for sanitize_filename function."""

    def test_sanitize_filename_basic(self) -> None:
        """Test basic sanitization."""
        result = core_paths.sanitize_filename("normal_file.txt")
        assert result == "normal_file.txt"

    def test_sanitize_filename_special_chars(self) -> None:
        """Test removal of special characters."""
        result = core_paths.sanitize_filename("file<>:\"/\\|?*.txt")
        assert "<" not in result
        assert ">" not in result


class TestCreateTempFile:
    """Tests for create_temp_file function."""

    def test_create_temp_file_basic(self) -> None:
        """Test temp file path creation."""
        path = core_paths.create_temp_file(suffix=".txt", prefix="test_")
        assert str(path).endswith(".txt")
        assert "test_" in str(path)
        # Path should not exist yet
        assert not path.exists()


class TestProjectDirectories:
    """Tests for project directory functions."""

    def test_get_project_root_returns_path(self) -> None:
        """Test get_project_root returns a Path."""
        root = core_paths.get_project_root()
        assert isinstance(root, Path)
        assert root.is_absolute()

    def test_get_data_dir_returns_path(self) -> None:
        """Test get_data_dir returns a Path."""
        data_dir = core_paths.get_data_dir()
        assert isinstance(data_dir, Path)

    def test_get_cache_dir_returns_path(self) -> None:
        """Test get_cache_dir returns a Path."""
        cache_dir = core_paths.get_cache_dir()
        assert isinstance(cache_dir, Path)

    def test_get_logs_dir_returns_path(self) -> None:
        """Test get_logs_dir returns a Path."""
        logs_dir = core_paths.get_logs_dir()
        assert isinstance(logs_dir, Path)

    def test_get_temp_dir_returns_path(self) -> None:
        """Test get_temp_dir returns a Path."""
        temp_dir = core_paths.get_temp_dir()
        assert isinstance(temp_dir, Path)


class TestResolvePath:
    """Tests for resolve_path function."""

    def test_resolve_path_basic(self) -> None:
        """Test path resolution."""
        resolved = core_paths.resolve_path(".")
        assert isinstance(resolved, Path)
        assert resolved.is_absolute()
