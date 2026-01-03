"""Comprehensive tests for core functionality modules.

These tests expand coverage for core utilities that form the foundation
of the METAINFORMANT toolkit. All tests follow the NO_MOCKING policy
and use real implementations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core import cache, config, hash, io, logging, parallel, paths, text

try:
    from metainformant.core import db

    DB_AVAILABLE = True
except ImportError:
    db = None
    DB_AVAILABLE = False


class TestCoreCache:
    """Test cache functionality with real file operations."""

    def test_simple_json_cache(self, tmp_path):
        """Test basic JSON caching operations."""
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        # Create a cache instance
        json_cache = cache.JsonCache(cache_dir)

        # Test cache miss (key doesn't exist)
        result = json_cache.get("test_key")
        assert result is None

        # Test cache set
        test_data = {"key1": "value1", "numbers": [1, 2, 3]}
        json_cache.set("test_key", test_data)

        # Test cache hit
        cached_result = json_cache.get("test_key")
        assert cached_result == test_data

        # Test cache update
        test_data["key2"] = "value2"
        json_cache.set("test_key", test_data)

        updated_result = json_cache.get("test_key")
        assert updated_result == test_data
        assert "key2" in updated_result

    def test_cache_with_ttl(self, tmp_path):
        """Test cache with time-to-live functionality."""
        import time

        cache_dir = tmp_path / "ttl_cache"
        cache_dir.mkdir()

        # Create a cache instance with 1 second TTL
        json_cache = cache.JsonCache(cache_dir, ttl_seconds=1)
        test_data = {"timestamp": time.time()}

        # Set cache
        json_cache.set("test_key", test_data)

        # Should be fresh immediately
        result = json_cache.get("test_key")
        assert result == test_data

        # Verify cache size
        assert json_cache.size() == 1


class TestCoreConfig:
    """Test configuration management."""

    def test_config_loading(self, tmp_path):
        """Test loading configuration from files."""
        config_file = tmp_path / "test_config.yaml"
        config_data = """
        database:
          host: localhost
          port: 5432
          name: testdb
        analysis:
          threads: 4
          memory_gb: 8
        """
        config_file.write_text(config_data)

        # Test loading YAML config
        loaded_config = config.load_config_file(config_file)
        assert "database" in loaded_config
        assert loaded_config["database"]["host"] == "localhost"
        assert loaded_config["analysis"]["threads"] == 4

    def test_env_override(self):
        """Test environment variable configuration override."""
        import os

        # Set test environment variable
        test_key = "METAINFORMANT_TEST_VALUE"
        test_value = "override_value"
        os.environ[test_key] = test_value

        try:
            # Test env var retrieval
            result = config.get_env_or_default("METAINFORMANT_TEST_VALUE", "default")
            assert result == test_value

            # Test default when env var not set
            result = config.get_env_or_default("NONEXISTENT_VAR", "default")
            assert result == "default"
        finally:
            # Clean up
            os.environ.pop(test_key, None)


class TestCoreDatabase:
    """Test database utilities."""

    @pytest.mark.skipif(not DB_AVAILABLE, reason="Database module not available")
    def test_connection_string_building(self):
        """Test database connection string construction."""
        params = {"host": "localhost", "port": 5432, "database": "testdb", "user": "testuser", "password": "testpass"}

        conn_str = db.build_postgres_url(**params)
        assert "postgresql://" in conn_str
        assert "testuser:testpass@localhost:5432/testdb" in conn_str

    @pytest.mark.skipif(not DB_AVAILABLE, reason="Database module not available")
    def test_safe_connection_params(self):
        """Test parameter sanitization."""
        unsafe_params = {"host": "localhost; DROP TABLE users;", "database": "test'db", "user": "test user"}

        safe_params = db.sanitize_connection_params(unsafe_params)

        # Should remove or escape dangerous characters
        assert "DROP TABLE" not in safe_params["host"]
        assert "'" not in safe_params["database"]


class TestCoreHashing:
    """Test cryptographic hashing utilities."""

    def test_sha256_consistency(self):
        """Test SHA256 hashing produces consistent results."""
        test_data = b"Hello, METAINFORMANT!"

        hash1 = hash.sha256_bytes(test_data)
        hash2 = hash.sha256_bytes(test_data)

        assert hash1 == hash2
        assert len(hash1) == 64  # SHA256 hex string length
        assert isinstance(hash1, str)

        # Different data should produce different hashes
        different_data = b"Different data"
        hash3 = hash.sha256_bytes(different_data)
        assert hash1 != hash3

    def test_file_hashing(self, tmp_path):
        """Test hashing file contents."""
        test_file = tmp_path / "test_file.txt"
        test_content = "Test file content for hashing"
        test_file.write_text(test_content)

        file_hash = hash.sha256_file(test_file)
        content_hash = hash.sha256_bytes(test_content.encode("utf-8"))

        assert file_hash == content_hash

    def test_deterministic_seeding(self):
        """Test deterministic hash-based seeding."""
        seed_data = "reproducible_seed"

        seed1 = hash.deterministic_seed(seed_data)
        seed2 = hash.deterministic_seed(seed_data)

        assert seed1 == seed2
        assert isinstance(seed1, int)
        assert 0 <= seed1 <= 2**31 - 1  # Valid random seed range


class TestCoreIO:
    """Expanded I/O functionality tests."""

    def test_json_with_compression(self, tmp_path):
        """Test JSON I/O with compression using gzip extension."""
        test_data = {
            "large_list": list(range(1000)),
            "nested": {"deep": {"data": "value"}},
            "metadata": {"compressed": True},
        }

        json_file = tmp_path / "compressed.json.gz"

        # Write compressed (automatic based on .gz extension)
        io.dump_json(test_data, json_file)
        assert json_file.exists()

        # Read compressed (automatic based on .gz extension)
        loaded_data = io.load_json(json_file)
        assert loaded_data == test_data

    def test_csv_operations(self, tmp_path):
        """Test CSV reading and writing."""
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

        # Read CSV (no index column expected since we write without index)
        loaded_df = io.read_csv(csv_file)
        assert loaded_df.shape == test_df.shape
        assert list(loaded_df.columns) == list(test_df.columns)

        # Check data consistency
        assert loaded_df.loc[0, "species"] == "E_coli"
        assert abs(loaded_df.loc[0, "gc_content"] - 0.507) < 1e-6

    def test_tsv_operations(self, tmp_path):
        """Test TSV (tab-separated values) operations."""
        test_data = [
            ["gene_id", "expression", "condition"],
            ["GENE001", "150.5", "control"],
            ["GENE002", "200.3", "treatment"],
            ["GENE003", "75.8", "control"],
        ]

        tsv_file = tmp_path / "expression.tsv"

        # Write TSV
        io.write_tsv(test_data, tsv_file)

        # Read TSV
        loaded_data = io.read_tsv(tsv_file)
        assert len(loaded_data) == 4  # Including header
        assert loaded_data[0] == test_data[0]
        assert loaded_data[1][1] == "150.5"


class TestCorePaths:
    """Test path manipulation utilities."""

    def test_path_expansion(self, tmp_path):
        """Test path expansion and resolution."""
        # Test relative path resolution
        rel_path = "data/test.txt"
        expanded = paths.expand_and_resolve(rel_path)
        assert expanded.is_absolute()

        # Test home directory expansion
        home_path = "~/test.txt"
        expanded_home = paths.expand_and_resolve(home_path)
        assert expanded_home.is_absolute()
        assert str(expanded_home).startswith(str(Path.home()))

    def test_safe_path_creation(self, tmp_path):
        """Test safe directory and file creation."""
        test_dir = tmp_path / "nested" / "directory"

        # Create nested directory safely
        paths.ensure_directory(test_dir)
        assert test_dir.exists()
        assert test_dir.is_dir()

        # Test file path preparation
        test_file = test_dir / "test_file.txt"
        paths.prepare_file_path(test_file)
        assert test_file.parent.exists()

    def test_path_validation(self):
        """Test path validation and security."""
        # Test valid paths
        assert paths.is_safe_path("/tmp/valid_file.txt")
        assert paths.is_safe_path("relative/path.txt")

        # Test potentially unsafe paths (path traversal)
        assert not paths.is_safe_path("../../../etc/passwd")
        assert not paths.is_safe_path("/etc/shadow")

    def test_file_extension_handling(self):
        """Test file extension utilities."""
        # Test extension detection
        assert paths.get_file_extension("data.json") == ".json"
        assert paths.get_file_extension("archive.tar.gz") == ".gz"
        assert paths.get_file_extension("no_extension") == ""

        # Test extension changing
        new_path = paths.change_extension("data.csv", ".tsv")
        assert str(new_path).endswith(".tsv")


class TestCoreText:
    """Test text processing utilities."""

    def test_advanced_slugify(self):
        """Test advanced text slugification."""
        # Basic slugification
        assert text.slugify("Hello World") == "hello-world"

        # Special characters (Greek letters removed as expected)
        assert text.slugify("Protein α-helix β-sheet") == "protein-helix-sheet"

        # Scientific notation
        assert text.slugify("E. coli K-12") == "e-coli-k-12"

        # Multiple spaces and punctuation
        assert text.slugify("Multiple   spaces!!!") == "multiple-spaces"

    def test_text_cleaning(self):
        """Test text cleaning utilities."""
        dirty_text = "  \t\nExtra   whitespace\n\n  "
        clean_text = text.clean_whitespace(dirty_text)
        assert clean_text == "Extra whitespace"

        # Test removing control characters
        text_with_controls = "Normal text\x00\x01\x02"
        cleaned = text.remove_control_chars(text_with_controls)
        assert cleaned == "Normal text"

    def test_biological_text_processing(self):
        """Test biology-specific text processing."""
        # Gene name standardization
        gene_variants = ["BRCA1", "brca1", "Brca1", "BRCA-1"]
        standardized = [text.standardize_gene_name(g) for g in gene_variants]
        assert len(set(standardized)) == 1  # Should all be the same

        # Species name formatting
        species_raw = "escherichia coli"
        formatted = text.format_species_name(species_raw)
        assert formatted == "Escherichia coli"

        # Sequence ID cleaning
        messy_id = ">gi|123456|ref|NM_001234.1| some gene [Homo sapiens]"
        clean_id = text.clean_sequence_id(messy_id)
        assert clean_id.startswith("NM_001234.1")
        assert "[Homo sapiens]" not in clean_id


class TestCoreLogging:
    """Test logging configuration and utilities."""

    def test_logger_setup(self, tmp_path):
        """Test logger configuration."""
        log_file = tmp_path / "test.log"

        # Set up logger
        logger = logging.setup_logger("test_logger", log_file, level="INFO")

        # Test logging
        logger.info("Test info message")
        logger.warning("Test warning message")

        # Verify log file creation
        assert log_file.exists()

        # Read log contents
        log_content = log_file.read_text()
        assert "Test info message" in log_content
        assert "Test warning message" in log_content

    def test_structured_logging(self):
        """Test structured logging with metadata."""
        # Create in-memory log handler for testing
        import io
        import logging as stdlib_logging

        log_stream = io.StringIO()
        handler = stdlib_logging.StreamHandler(log_stream)

        logger = stdlib_logging.getLogger("test_structured")
        logger.addHandler(handler)
        logger.setLevel(stdlib_logging.INFO)

        # Log with structured data
        logging.log_with_metadata(
            logger, "Analysis completed", {"samples": 100, "genes": 20000, "runtime_seconds": 45.2}
        )

        log_output = log_stream.getvalue()
        assert "Analysis completed" in log_output
        assert "samples" in log_output


class TestCoreMethods:
    """Test all the new methods added to core modules."""

    def test_enhanced_hash_utilities(self, tmp_path):
        """Test enhanced hash utilities."""
        # Test string hashing
        hash1 = hash.sha256_string("test1")
        hash2 = hash.sha256_string("test2")
        assert len(hash1) == 64
        assert hash1 != hash2

        # Test file comparison
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("identical content")
        file2.write_text("identical content")

        assert hash.file_hash_comparison(file1, file2)

        # Different content
        file2.write_text("different content")
        assert not hash.file_hash_comparison(file1, file2)

        # Test directory hashing
        dir_hashes = hash.hash_directory(tmp_path)
        assert "file1.txt" in dir_hashes
        assert "file2.txt" in dir_hashes

        # Test file integrity verification
        expected_hash = hash.sha256_file(file1)
        assert hash.verify_file_integrity(file1, expected_hash)
        assert not hash.verify_file_integrity(file1, "wrong_hash")

    def test_enhanced_path_utilities(self, tmp_path):
        """Test enhanced path utilities."""
        # Test file finding by extension
        txt_files = paths.find_files_by_extension(tmp_path, "txt")
        assert isinstance(txt_files, list)

        # Test file size
        test_file = tmp_path / "size_test.txt"
        test_file.write_text("test content")
        size = paths.get_file_size(test_file)
        assert size > 0

        # Test directory size
        dir_size = paths.get_directory_size(tmp_path)
        assert dir_size > 0

        # Test filename sanitization
        unsafe = "file<name>test.txt"
        safe = paths.sanitize_filename(unsafe)
        assert "<" not in safe
        assert ">" not in safe
        assert safe == "file_name_test.txt"

        # Test temp file creation (function returns path, we create the file)
        temp_file = paths.create_temp_file(suffix=".txt", directory=tmp_path)
        temp_file.write_text("test content")
        assert temp_file.exists()
        assert temp_file.parent == tmp_path

    def test_enhanced_text_utilities(self):
        """Test enhanced text utilities."""
        # Test number extraction
        text_with_numbers = "Hello 123 world 456.78 more 999"
        numbers = text.extract_numbers(text_with_numbers)
        assert numbers == [123.0, 456.78, 999.0]

        # Test text truncation
        long_text = "This is a very long text that should be truncated"
        truncated = text.truncate_text(long_text, 30)
        assert len(truncated) <= 30
        assert truncated.endswith("...")

        # Test word counting
        word_count = text.count_words("Hello world this is a test")
        assert word_count == 6

        # Test email extraction
        text_with_emails = "Contact us at test@example.com or admin@company.org"
        emails = text.extract_email_addresses(text_with_emails)
        assert "test@example.com" in emails
        assert "admin@company.org" in emails

    @pytest.mark.skip("API functions not implemented in cache module")
    def test_enhanced_cache_utilities(self, tmp_path):
        """Test enhanced cache utilities."""
        cache_dir = tmp_path / "enhanced_cache"

        # Test cache info
        info = cache.get_cache_info(cache_dir)
        assert not info["exists"]  # Directory doesn't exist yet

        # Create cache directory and add data
        cache.cache_json(cache_dir, "test1", {"data": "value1"})
        cache.cache_json(cache_dir, "test2", {"data": "value2"})

        # Check cache info
        info = cache.get_cache_info(cache_dir)
        assert info["exists"]
        assert info["total_files"] == 2

        # Test cache clearing
        cache.clear_cache_dir(cache_dir)
        info = cache.get_cache_info(cache_dir)
        assert info["total_files"] == 0

    def test_enhanced_io_utilities(self, tmp_path):
        """Test enhanced I/O utilities."""
        # Test gzipped JSON
        test_data = {"test": "data", "numbers": [1, 2, 3]}
        gz_path = tmp_path / "test.json.gz"

        io.dump_json_gz(test_data, gz_path)
        loaded = io.load_json_gz(gz_path)
        assert loaded == test_data

        # Test CSV operations (if pandas available)
        try:
            import pandas as pd
            df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
            csv_path = tmp_path / "test.csv"

            io.write_csv(df, csv_path)
            loaded_df = io.read_csv(csv_path)
            assert isinstance(loaded_df, pd.DataFrame)

            # Test gzipped CSV
            gz_csv_path = tmp_path / "test.csv.gz"
            io.write_csv(df, gz_csv_path)
            loaded_gz_df = io.read_csv(gz_csv_path)
            assert isinstance(loaded_gz_df, pd.DataFrame)

        except ImportError:
            # Pandas not available, skip CSV tests
            pass

        # Test Parquet operations (if pandas available)
        try:
            import pandas as pd
            df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
            parquet_path = tmp_path / "test.parquet"

            io.write_parquet(df, parquet_path)
            loaded_df = io.read_parquet(parquet_path)
            assert isinstance(loaded_df, pd.DataFrame)

        except ImportError:
            # Pandas not available, skip Parquet tests
            pass

    def test_enhanced_parallel_utilities(self):
        """Test enhanced parallel utilities."""
        # Test CPU count
        cpu_count = parallel.cpu_count()
        assert cpu_count > 0

        # Test unordered mapping
        def square(x):
            return x * x

        results = parallel.thread_map_unordered(square, [1, 2, 3, 4, 5])
        assert len(results) == 5
        assert 1 in results  # 1^2 = 1
        assert 4 in results  # 2^2 = 4
        assert 9 in results  # 3^2 = 9

        # Test batch processing
        def process_batch(batch):
            return [x * 2 for x in batch]

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        results = parallel.parallel_batch(process_batch, data, batch_size=3)
        expected = [2, 4, 6, 8, 10, 12, 14, 16]
        assert sorted(results) == sorted(expected)
