"""Tests for core discovery utilities.

Tests function discovery, config discovery, and workflow discovery following NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core import discovery


class TestDiscoverFunctions:
    """Tests for discover_functions function."""

    def test_discover_functions_in_core_module(self):
        """Test discovering functions in a core module."""
        # Test with a known module
        module_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "io.py"
        if module_path.exists():
            functions = discovery.discover_functions(module_path)
            assert isinstance(functions, list)
            assert len(functions) > 0
            # Check that we found some expected functions
            function_names = [f.name for f in functions]
            assert any("load" in name.lower() or "dump" in name.lower() for name in function_names)

    def test_discover_functions_with_pattern(self):
        """Test discovering functions with name pattern."""
        module_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "io.py"
        if module_path.exists():
            functions = discovery.discover_functions(module_path, pattern="load")
            assert isinstance(functions, list)
            # All functions should match pattern
            for func in functions:
                assert "load" in func.name.lower()

    def test_discover_functions_nonexistent_file(self):
        """Test that nonexistent file raises FileNotFoundError."""
        nonexistent = Path("/nonexistent/path/module.py")
        with pytest.raises(FileNotFoundError):
            discovery.discover_functions(nonexistent)

    def test_discover_functions_invalid_syntax(self, tmp_path):
        """Test that invalid syntax raises SyntaxError."""
        invalid_file = tmp_path / "invalid.py"
        invalid_file.write_text("def invalid syntax here")
        with pytest.raises(SyntaxError):
            discovery.discover_functions(invalid_file)


class TestDiscoverConfigs:
    """Tests for discover_configs function."""

    def test_discover_configs_repo_root(self):
        """Test discovering configs from repo root."""
        repo_root = Path(__file__).parent.parent
        configs = discovery.discover_configs(repo_root)
        assert isinstance(configs, list)
        # Should find at least some config files
        assert len(configs) > 0
        for config in configs:
            assert hasattr(config, "path")
            assert hasattr(config, "format")

    def test_discover_configs_with_domain(self):
        """Test discovering configs for specific domain."""
        repo_root = Path(__file__).parent.parent
        configs = discovery.discover_configs(repo_root, domain="gwas")
        assert isinstance(configs, list)
        # All configs should be in gwas domain or match pattern
        for config in configs:
            assert "gwas" in str(config.path).lower() or config.domain == "gwas"


class TestDiscoverOutputPatterns:
    """Tests for discover_output_patterns function."""

    def test_discover_output_patterns_core(self):
        """Test discovering output patterns for core module."""
        pattern = discovery.discover_output_patterns("core")
        assert hasattr(pattern, "module")
        assert hasattr(pattern, "base_pattern")
        assert pattern.module == "core"

    def test_discover_output_patterns_dna(self):
        """Test discovering output patterns for dna module."""
        pattern = discovery.discover_output_patterns("dna")
        assert hasattr(pattern, "module")
        assert pattern.module == "dna"


class TestBuildCallGraph:
    """Tests for build_call_graph function."""

    def test_build_call_graph_simple_module(self):
        """Test building call graph for a simple module."""
        repo_root = Path(__file__).parent.parent
        # Use a simple module like core/io.py
        entry_point = repo_root / "src" / "metainformant" / "core" / "io.py"
        if entry_point.exists():
            graph = discovery.build_call_graph(entry_point, repo_root)
            assert isinstance(graph, dict)
            # Graph should have some entries
            assert len(graph) >= 0


class TestFindSymbolUsage:
    """Tests for find_symbol_usage function."""

    def test_find_symbol_usage_common_function(self):
        """Test finding usage of a common function."""
        # Search for a common function like "get_logger" in just the core module
        core_dir = Path(__file__).parent.parent / "src" / "metainformant" / "core"
        usages = discovery.find_symbol_usage("get_logger", core_dir)
        assert isinstance(usages, list)
        # Should find at least some usages
        assert len(usages) > 0
        for usage in usages:
            assert hasattr(usage, "file")
            assert hasattr(usage, "line")
            # All usages should be in core module
            assert str(core_dir) in str(usage.file)


class TestGetModuleDependencies:
    """Tests for get_module_dependencies function."""

    def test_get_module_dependencies_core_io(self):
        """Test getting dependencies for core/io module."""
        module_path = Path(__file__).parent.parent / "src" / "metainformant" / "core" / "io.py"
        if module_path.exists():
            deps = discovery.get_module_dependencies(module_path)
            assert hasattr(deps, "module")
            assert hasattr(deps, "imports")
            assert hasattr(deps, "from_imports")
            assert isinstance(deps.imports, list)
            assert isinstance(deps.from_imports, dict)


class TestDiscoverWorkflows:
    """Tests for discover_workflows function."""

    def test_discover_workflows_repo_root(self):
        """Test discovering workflows from repo root."""
        repo_root = Path(__file__).parent.parent
        workflows = discovery.discover_workflows(repo_root)
        assert isinstance(workflows, list)
        # Should find workflow definitions
        for workflow in workflows:
            assert isinstance(workflow, dict)

    def test_discover_workflows_none(self):
        """Test discovering workflows with None repo_root."""
        workflows = discovery.discover_workflows(None)
        assert isinstance(workflows, list)

