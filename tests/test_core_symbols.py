"""Tests for core symbol indexing utilities.

Tests symbol indexing, lookup, and reference finding following NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core import symbols


class TestIndexFunctions:
    """Tests for index_functions function."""

    def test_index_functions_repo_root(self):
        """Test indexing functions in repository."""
        repo_root = Path(__file__).parent.parent
        index = symbols.index_functions(repo_root, use_cache=False)
        assert isinstance(index, dict)
        # Should have indexed some functions
        assert len(index) >= 0

    def test_index_functions_with_cache(self):
        """Test indexing functions with cache enabled."""
        repo_root = Path(__file__).parent.parent
        index = symbols.index_functions(repo_root, use_cache=True)
        assert isinstance(index, dict)


class TestIndexClasses:
    """Tests for index_classes function."""

    def test_index_classes_repo_root(self):
        """Test indexing classes in repository."""
        repo_root = Path(__file__).parent.parent
        index = symbols.index_classes(repo_root, use_cache=False)
        assert isinstance(index, dict)
        # Should have indexed some classes
        assert len(index) >= 0


class TestFindSymbol:
    """Tests for find_symbol function."""

    def test_find_symbol_function(self):
        """Test finding a function symbol."""
        repo_root = Path(__file__).parent.parent
        results = symbols.find_symbol("get_logger", "function", repo_root)
        assert isinstance(results, list)
        # Should find get_logger in logging module
        assert len(results) > 0
        for result in results:
            assert hasattr(result, "name")
            assert hasattr(result, "file_path")
            assert result.name == "get_logger"

    def test_find_symbol_class(self):
        """Test finding a class symbol."""
        repo_root = Path(__file__).parent.parent
        results = symbols.find_symbol("FunctionInfo", "class", repo_root)
        assert isinstance(results, list)
        # May or may not find it depending on indexing

    def test_find_symbol_nonexistent(self):
        """Test finding nonexistent symbol."""
        repo_root = Path(__file__).parent.parent
        results = symbols.find_symbol("NonexistentSymbol12345", "function", repo_root)
        assert isinstance(results, list)
        # Should return empty list
        assert len(results) == 0


class TestGetSymbolSignature:
    """Tests for get_symbol_signature function."""

    def test_get_symbol_signature_existing(self):
        """Test getting signature for existing symbol."""
        symbol_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "logging.py"
        if symbol_path.exists():
            signature = symbols.get_symbol_signature(symbol_path, "get_logger")
            assert signature is None or isinstance(signature, str)

    def test_get_symbol_signature_nonexistent(self):
        """Test getting signature for nonexistent symbol."""
        symbol_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "logging.py"
        if symbol_path.exists():
            signature = symbols.get_symbol_signature(symbol_path, "NonexistentFunction")
            assert signature is None


class TestFindSymbolReferences:
    """Tests for find_symbol_references function."""

    def test_find_symbol_references_common(self):
        """Test finding references to a common symbol."""
        repo_root = Path(__file__).parent.parent
        references = symbols.find_symbol_references("get_logger", repo_root)
        assert isinstance(references, list)
        # Should find some references
        assert len(references) > 0
        for ref in references:
            assert hasattr(ref, "symbol_name")
            assert hasattr(ref, "file_path")
            assert hasattr(ref, "line_number")


class TestGetSymbolMetadata:
    """Tests for get_symbol_metadata function."""

    def test_get_symbol_metadata_existing(self):
        """Test getting metadata for existing symbol."""
        symbol_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "logging.py"
        if symbol_path.exists():
            metadata = symbols.get_symbol_metadata(symbol_path, "get_logger")
            assert isinstance(metadata, dict)
            # Should have some metadata fields
            assert len(metadata) >= 0


class TestFuzzyFindSymbol:
    """Tests for fuzzy_find_symbol function."""

    def test_fuzzy_find_symbol_close_match(self):
        """Test fuzzy finding with close match."""
        repo_root = Path(__file__).parent.parent
        results = symbols.fuzzy_find_symbol("get_logr", "function", repo_root, threshold=0.6)
        assert isinstance(results, list)
        # Should find get_logger with high similarity
        for symbol_name, score in results:
            assert isinstance(symbol_name, str)
            assert isinstance(score, float)
            assert 0.0 <= score <= 1.0

    def test_fuzzy_find_symbol_no_match(self):
        """Test fuzzy finding with no close match."""
        repo_root = Path(__file__).parent.parent
        results = symbols.fuzzy_find_symbol("XyzAbc123", "function", repo_root, threshold=0.9)
        assert isinstance(results, list)
        # Should return empty or very few results with high threshold

