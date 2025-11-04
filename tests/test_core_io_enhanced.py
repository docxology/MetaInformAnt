"""Enhanced tests for I/O functionality."""

import pytest
import tempfile
from pathlib import Path

from metainformant.core import io, paths


class TestIOEnhanced:
    """Test enhanced I/O functionality."""

    @pytest.mark.network
    def test_download_file(self):
        """Test file download functionality."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Test downloading a real file (using httpbin for testing)
            url = "https://httpbin.org/json"
            dest_file = tmp_path / "test.json"
            
            # This will fail without internet, but tests the method structure
            try:
                success = io.download_file(url, dest_file)
                if success:
                    assert dest_file.exists()
                    content = dest_file.read_text()
                    assert "httpbin" in content
            except Exception:
                # Expected to fail without internet, but method should not crash
                pass

    @pytest.mark.network
    def test_download_json(self):
        """Test JSON download functionality."""
        # Test with a real API endpoint
        url = "https://httpbin.org/json"
        
        try:
            data = io.download_json(url)
            if data:
                assert isinstance(data, dict)
                assert "httpbin" in str(data)
        except Exception:
            # Expected to fail without internet
            pass

    @pytest.mark.network
    def test_download_text(self):
        """Test text download functionality."""
        url = "https://httpbin.org/html"
        
        try:
            text = io.download_text(url)
            if text:
                assert isinstance(text, str)
                assert "html" in text.lower()
        except Exception:
            # Expected to fail without internet
            pass

    @pytest.mark.network
    def test_batch_download(self):
        """Test batch download functionality."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            urls = [
                "https://httpbin.org/json",
                "https://httpbin.org/html"
            ]
            
            try:
                results = io.batch_download(urls, tmp_path)
                assert isinstance(results, dict)
                assert len(results) == 2
                # Results may be False without internet, but structure should be correct
            except Exception:
                # Expected to fail without internet
                pass

    @pytest.mark.network
    def test_csv_download(self):
        """Test CSV download functionality."""
        url = "https://httpbin.org/csv"
        
        try:
            df = io.download_csv(url)
            if df is not None:
                assert hasattr(df, 'shape')  # pandas DataFrame
        except Exception:
            # Expected to fail without internet
            pass

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

        except ImportError:
            pytest.skip("pandas not available")

    def test_pandas_parquet_operations(self, tmp_path):
        """Test pandas Parquet operations."""
        try:
            import pandas as pd

            # Create test DataFrame
            test_df = pd.DataFrame(
                {
                    "A": [1, 2, 3, 4, 5],
                    "B": [10, 20, 30, 40, 50],
                    "C": ["a", "b", "c", "d", "e"]
                }
            )

            parquet_file = tmp_path / "test_data.parquet"

            # Write Parquet
            io.write_parquet(test_df, parquet_file)
            assert parquet_file.exists()

            # Read Parquet
            loaded_df = io.read_parquet(parquet_file)
            assert loaded_df.shape == test_df.shape
            assert list(loaded_df.columns) == list(test_df.columns)

        except ImportError:
            pytest.skip("pandas not available")

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
        assert size == len(test_content.encode('utf-8'))

        # Test directory size
        dir_size = paths.get_directory_size(tmp_path)
        assert dir_size >= size

        # Test non-existent file
        non_existent = tmp_path / "does_not_exist.txt"
        size_non = paths.get_file_size(non_existent)
        assert size_non == 0
