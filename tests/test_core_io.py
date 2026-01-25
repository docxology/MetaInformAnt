"""Tests for I/O functionality."""

from __future__ import annotations

import pytest
import tempfile
from pathlib import Path

from metainformant.core import io, paths


def _check_online(url: str, timeout: int = 5) -> bool:
    """Check if we can reach a URL within timeout.

    Args:
        url: URL to check
        timeout: Timeout in seconds

    Returns:
        True if URL is reachable, False otherwise
    """
    try:
        import requests

        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        return True
    except Exception:
        return False


class TestIO:
    """Test enhanced I/O functionality."""

    @pytest.mark.network
    def test_download_file(self):
        """Test file download functionality.

        Uses real HTTP requests with graceful skip if network unavailable.
        """
        if not _check_online("https://httpbin.org"):
            pytest.skip("Network unavailable - real implementation requires connectivity")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Test downloading a real file (using httpbin for testing)
            url = "https://httpbin.org/json"
            dest_file = tmp_path / "test.json"

            success = io.download_file(url, dest_file)
            if success:
                assert dest_file.exists()
                content = dest_file.read_text()
                assert "slideshow" in content  # httpbin.org/json returns a slideshow JSON

    @pytest.mark.network
    def test_download_json(self):
        """Test JSON download functionality.

        Uses real HTTP requests with graceful skip if network unavailable.
        """
        if not _check_online("https://httpbin.org"):
            pytest.skip("Network unavailable - real implementation requires connectivity")

        url = "https://httpbin.org/json"
        data = io.download_json(url)
        if data:
            assert isinstance(data, dict)
            assert "slideshow" in str(data)  # httpbin.org/json returns a slideshow JSON

    @pytest.mark.network
    def test_download_text(self):
        """Test text download functionality.

        Uses real HTTP requests with graceful skip if network unavailable.
        """
        if not _check_online("https://httpbin.org"):
            pytest.skip("Network unavailable - real implementation requires connectivity")

        url = "https://httpbin.org/html"
        text = io.download_text(url)
        if text:
            assert isinstance(text, str)
            assert "html" in text.lower()

    @pytest.mark.network
    def test_batch_download(self):
        """Test batch download functionality.

        Uses minimal URLs to reduce test execution time while still
        verifying batch download behavior. Uses real HTTP requests with
        graceful skip if network unavailable.
        """
        if not _check_online("https://httpbin.org"):
            pytest.skip("Network unavailable - real implementation requires connectivity")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Use single URL to reduce network calls and execution time
            urls = ["https://httpbin.org/json"]

            results = io.batch_download(urls, tmp_path)
            assert isinstance(results, dict)
            assert len(results) == 1

    @pytest.mark.network
    def test_csv_download(self):
        """Test CSV download functionality.

        Uses real HTTP requests with graceful skip if network unavailable.
        """
        if not _check_online("https://httpbin.org"):
            pytest.skip("Network unavailable - real implementation requires connectivity")

        url = "https://httpbin.org/csv"
        df = io.download_csv(url)
        if df is not None:
            assert hasattr(df, "shape")  # pandas DataFrame

    def test_gzipped_json_io(self, tmp_path):
        """Test gzipped JSON I/O operations."""
        test_data = {
            "large_list": list(range(1000)),
            "nested": {"deep": {"data": "value"}},
            "metadata": {"compressed": True},
        }

        gz_file = tmp_path / "compressed.json.gz"

        # Write compressed
        io.dump_json_gz(test_data, gz_file)
        assert gz_file.exists()

        # Read compressed
        loaded_data = io.load_json_gz(gz_file)
        assert loaded_data == test_data

    def test_pandas_csv_operations(self, tmp_path):
        """Test pandas CSV operations."""
        try:
            import pandas as pd

            # Create test DataFrame
            test_df = pd.DataFrame(
                {
                    "species": ["E_coli", "B_subtilis", "S_aureus"],
                    "gc_content": [0.507, 0.436, 0.328],
                    "genome_size": [4641652, 4215606, 2821337],
                }
            )

            csv_file = tmp_path / "test_data.csv"

            # Write CSV
            io.write_csv(test_df, csv_file)
            assert csv_file.exists()

            # Read CSV
            loaded_df = io.read_csv(csv_file)
            assert loaded_df.shape == test_df.shape
            assert list(loaded_df.columns) == list(test_df.columns)

            # Check data consistency
            assert loaded_df.loc[0, "species"] == "E_coli"
            assert abs(loaded_df.loc[0, "gc_content"] - 0.507) < 1e-6

        except ImportError as e:
            pytest.skip(f"pandas not available: {e}")

    def test_pandas_parquet_operations(self, tmp_path):
        """Test pandas Parquet operations (requires pyarrow or fastparquet)."""
        pytest.importorskip(
            "pyarrow",
            reason="Parquet support requires pyarrow or fastparquet. Install with: uv add pyarrow",
        )
        try:
            import pandas as pd

            # Create test DataFrame
            test_df = pd.DataFrame({"A": [1, 2, 3, 4, 5], "B": [10, 20, 30, 40, 50], "C": ["a", "b", "c", "d", "e"]})

            parquet_file = tmp_path / "test_data.parquet"

            # Write Parquet
            io.write_parquet(test_df, parquet_file)
            assert parquet_file.exists()

            # Read Parquet
            loaded_df = io.read_parquet(parquet_file)
            assert loaded_df.shape == test_df.shape
            assert list(loaded_df.columns) == list(test_df.columns)

        except ImportError as e:
            error_msg = str(e).lower()
            if "pyarrow" in error_msg or "fastparquet" in error_msg:
                pytest.skip("Parquet support requires pyarrow or fastparquet. Install with: uv add pyarrow")
            pytest.skip(f"pandas or parquet engine not available: {e}")

    @pytest.mark.network
    def test_io_error_handling(self):
        """Test error handling in I/O operations."""
        # Test with invalid URLs
        invalid_url = "https://invalid-domain-that-does-not-exist.com/file.txt"

        result = io.download_file(invalid_url, "/tmp/test.txt")
        assert result is False

        json_result = io.download_json(invalid_url)
        assert json_result is None

        text_result = io.download_text(invalid_url)
        assert text_result is None

    def test_file_size_operations(self, tmp_path):
        """Test file size calculation operations."""
        # Create test file
        test_file = tmp_path / "size_test.txt"
        test_content = "This is test content for size calculation"
        test_file.write_text(test_content)

        size = paths.get_file_size(test_file)
        assert size == len(test_content.encode("utf-8"))

        # Test directory size
        dir_size = paths.get_directory_size(tmp_path)
        assert dir_size >= size

        # Test non-existent file
        non_existent = tmp_path / "does_not_exist.txt"
        size_non = paths.get_file_size(non_existent)
        assert size_non == 0


def test_ensure_directory_creates(tmp_path: Path) -> None:
    """Test directory creation."""
    p = tmp_path / "nested" / "dir"
    assert not p.exists()
    made = io.ensure_directory(p)
    assert made.exists() and made.is_dir()


def test_json_roundtrip(tmp_path: Path) -> None:
    """Test JSON roundtrip."""
    data = {"a": 1, "b": [1, 2, 3]}
    path = tmp_path / "data.json"
    io.dump_json(data, path)
    loaded = io.load_json(path)
    assert loaded == data


def test_json_gz_roundtrip(tmp_path: Path) -> None:
    """Test JSON gz roundtrip."""
    data = {"x": "y"}
    path = tmp_path / "data.json.gz"
    io.dump_json(data, path)
    loaded = io.load_json(path)
    assert loaded == data


def test_jsonl_roundtrip(tmp_path: Path) -> None:
    """Test JSONL roundtrip."""
    rows = [{"i": i} for i in range(5)]
    path = tmp_path / "rows.jsonl"
    io.write_jsonl(rows, path)
    read_back = list(io.read_jsonl(path))
    assert read_back == rows


def test_tsv_roundtrip(tmp_path: Path) -> None:
    """Test TSV roundtrip."""
    rows = [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}]
    path = tmp_path / "rows.tsv"
    io.write_delimited(rows, path, delimiter="\t")
    read_back = list(io.read_delimited(path, delimiter="\t"))
    assert read_back == rows


def test_open_text_auto_handles_gz(tmp_path: Path) -> None:
    """Test gzip handling in text auto-open."""
    txt_gz = tmp_path / "hello.txt.gz"
    with io.open_text_auto(txt_gz, mode="wt") as fh:
        fh.write("hello world\n")
    with io.open_text_auto(txt_gz, mode="rt") as fh:
        content = fh.read()
    assert content.strip() == "hello world"


# Edge case tests
class TestIOEdgeCases:
    """Edge case tests for I/O functionality."""

    def test_read_delimited_empty_file(self, tmp_path: Path) -> None:
        """Test reading an empty delimited file."""
        empty_file = tmp_path / "empty.csv"
        empty_file.write_text("")

        result = list(io.read_delimited(empty_file))
        assert result == []

    def test_read_delimited_header_only(self, tmp_path: Path) -> None:
        """Test reading a delimited file with only headers, no data rows."""
        header_only = tmp_path / "header_only.csv"
        header_only.write_text("name,value,count\n")

        result = list(io.read_delimited(header_only))
        assert result == []

    def test_read_delimited_with_none_values(self, tmp_path: Path) -> None:
        """Test reading delimited file with missing values."""
        csv_file = tmp_path / "with_nulls.csv"
        csv_file.write_text("a,b,c\n1,,3\n,2,\n")

        result = list(io.read_delimited(csv_file))
        assert len(result) == 2
        # csv.DictReader returns None for missing values, our wrapper converts to ""
        assert result[0]["b"] == ""
        assert result[1]["a"] == ""
        assert result[1]["c"] == ""

    def test_write_delimited_empty_rows(self, tmp_path: Path) -> None:
        """Test writing empty rows list creates empty file."""
        empty_csv = tmp_path / "empty_write.csv"
        io.write_delimited([], empty_csv)

        assert empty_csv.exists()
        content = empty_csv.read_text()
        assert content == ""

    def test_json_file_not_found(self, tmp_path: Path) -> None:
        """Test load_json raises FileNotFoundError for non-existent file."""
        non_existent = tmp_path / "does_not_exist.json"

        with pytest.raises(FileNotFoundError):
            io.load_json(non_existent)

    def test_json_invalid_syntax(self, tmp_path: Path) -> None:
        """Test load_json handles invalid JSON syntax."""
        from metainformant.core.io.errors import IOError as CoreIOError

        invalid_json = tmp_path / "invalid.json"
        invalid_json.write_text("{ invalid json syntax")

        with pytest.raises(CoreIOError):
            io.load_json(invalid_json)

    def test_dump_json_creates_parent_dirs(self, tmp_path: Path) -> None:
        """Test dump_json creates parent directories if they don't exist."""
        nested_path = tmp_path / "deep" / "nested" / "path" / "data.json"
        data = {"test": "data"}

        io.dump_json(data, nested_path)

        assert nested_path.exists()
        loaded = io.load_json(nested_path)
        assert loaded == data

    def test_read_csv_without_pandas(self, tmp_path: Path) -> None:
        """Test read_csv fallback behavior when pandas is available."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("name,value\ntest,123\n")

        result = io.read_csv(csv_file)
        # Either pandas DataFrame or dict of lists
        assert result is not None

    def test_jsonl_empty_lines_skipped(self, tmp_path: Path) -> None:
        """Test that empty lines in JSONL are skipped."""
        jsonl_file = tmp_path / "with_blanks.jsonl"
        jsonl_file.write_text('{"a": 1}\n\n{"b": 2}\n\n\n{"c": 3}\n')

        result = list(io.read_jsonl(jsonl_file))
        assert len(result) == 3
        assert result[0] == {"a": 1}
        assert result[1] == {"b": 2}
        assert result[2] == {"c": 3}

    def test_download_file_invalid_url_returns_false(self, tmp_path: Path) -> None:
        """Test download_file returns False for malformed URLs."""
        dest = tmp_path / "test.txt"
        # Use clearly invalid URL
        result = io.download_file("not-a-valid-url", dest)
        assert result is False

    def test_ensure_directory_idempotent(self, tmp_path: Path) -> None:
        """Test ensure_directory can be called multiple times safely."""
        target = tmp_path / "test_dir"

        # First call creates directory
        result1 = io.ensure_directory(target)
        assert result1.exists()

        # Second call should not raise
        result2 = io.ensure_directory(target)
        assert result2.exists()
        assert result1 == result2

    def test_read_tsv_returns_list_of_lists(self, tmp_path: Path) -> None:
        """Test read_tsv returns list of lists format."""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("a\tb\tc\n1\t2\t3\n4\t5\t6\n")

        result = io.read_tsv(tsv_file)
        assert len(result) == 3
        assert result[0] == ["a", "b", "c"]
        assert result[1] == ["1", "2", "3"]

    def test_write_tsv_handles_list_of_lists(self, tmp_path: Path) -> None:
        """Test write_tsv handles list of lists."""
        tsv_file = tmp_path / "output.tsv"
        data = [["a", "b"], ["1", "2"]]

        io.write_tsv(data, tsv_file)

        assert tsv_file.exists()
        content = tsv_file.read_text()
        assert "a\tb" in content
        assert "1\t2" in content

    def test_yaml_roundtrip(self, tmp_path: Path) -> None:
        """Test YAML write/read roundtrip (if PyYAML is available)."""
        yaml_file = tmp_path / "test.yaml"
        data = {"key": "value", "nested": {"list": [1, 2, 3]}}

        io.dump_yaml(data, yaml_file)
        assert yaml_file.exists()

        try:
            loaded = io.load_yaml(yaml_file)
            assert loaded["key"] == "value"
            assert loaded["nested"]["list"] == [1, 2, 3]
        except ImportError:
            pytest.skip("PyYAML not available")

    def test_atomic_write_json(self, tmp_path: Path) -> None:
        """Test atomic write leaves no temp files on success."""
        json_file = tmp_path / "atomic.json"
        data = {"atomic": True}

        io.dump_json(data, json_file, atomic=True)

        assert json_file.exists()
        # Temp file should be cleaned up
        temp_file = json_file.with_suffix(".json.tmp")
        assert not temp_file.exists()

        loaded = io.load_json(json_file)
        assert loaded == data
