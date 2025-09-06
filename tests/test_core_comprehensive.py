"""Comprehensive tests for core functionality modules.

These tests expand coverage for core utilities that form the foundation
of the METAINFORMANT toolkit. All tests follow the NO_MOCKING policy
and use real implementations.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from metainformant.core import cache, config, hash, io, logging, paths, text

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
        cache_file = tmp_path / "test_cache.json"

        # Test cache miss (file doesn't exist)
        result = cache.get_json_cache(cache_file, default={})
        assert result == {}

        # Test cache set
        test_data = {"key1": "value1", "numbers": [1, 2, 3]}
        cache.set_json_cache(cache_file, test_data)

        # Test cache hit
        cached_result = cache.get_json_cache(cache_file)
        assert cached_result == test_data

        # Test cache update
        test_data["key2"] = "value2"
        cache.set_json_cache(cache_file, test_data)

        updated_result = cache.get_json_cache(cache_file)
        assert updated_result == test_data
        assert "key2" in updated_result

    def test_cache_with_ttl(self, tmp_path):
        """Test cache with time-to-live functionality."""
        import time

        cache_file = tmp_path / "ttl_cache.json"
        test_data = {"timestamp": time.time()}

        # Set cache with short TTL
        cache.set_json_cache(cache_file, test_data)

        # Should be fresh immediately
        result = cache.get_json_cache(cache_file, max_age_seconds=1.0)
        assert result == test_data

        # Sleep briefly and test expiration (if TTL is implemented)
        time.sleep(0.1)
        # Cache should still be valid for 1 second TTL
        result = cache.get_json_cache(cache_file, max_age_seconds=1.0)
        assert result is not None


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

        # Read CSV
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
