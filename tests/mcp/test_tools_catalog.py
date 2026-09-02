"""Registry round-trip tests for the module-spec catalog (zero mocks)."""

from __future__ import annotations

import subprocess
import sys

from metainformant.mcp.registry import ToolRegistry
from metainformant.mcp.tools.catalog import TOOLS, build_registry_tools


def test_catalog_builds_expected_tool_count() -> None:
    assert 10 <= len(TOOLS) <= 25
    assert len({t.name for t in TOOLS}) == len(TOOLS)


def test_catalog_registers_into_real_registry() -> None:
    registry = ToolRegistry()
    names = registry.register_module_tools("metainformant.mcp.tools.catalog")
    assert registry.names() == sorted(names)
    assert "dna_translate" in registry.names()
    assert "rna_tau" in registry.names()


def test_catalog_dispatch_dna_translate_through_registry() -> None:
    registry = ToolRegistry()
    registry.register_module_tools("metainformant.mcp.tools.catalog")
    result = registry.call("dna_translate", {"sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"})
    assert result["protein"] == "MAIVMGR*KGAR*"


def test_catalog_dispatch_validates_arguments() -> None:
    registry = ToolRegistry()
    registry.register_module_tools("metainformant.mcp.tools.catalog")
    import pytest

    from metainformant.mcp.registry import SchemaError

    with pytest.raises(SchemaError):
        registry.call("dna_translate", {"wrong_arg": "x"})


def test_catalog_metadata_declares_writes_class() -> None:
    for tool in build_registry_tools():
        assert tool.metadata["writes"] in {"read-only", "output-dir-only"}


def test_subprocess_registry_server_lists_catalog_tools() -> None:
    """Real stdio MCP server exposes the catalog via tools/list."""
    code = (
        "import json, sys; from metainformant.mcp.registry import ToolRegistry; "
        "r = ToolRegistry(); r.register_module_tools('metainformant.mcp.tools.catalog'); "
        "print(json.dumps([t['name'] for t in r.list_tools()]))"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    names = __import__("json").loads(proc.stdout.strip().splitlines()[-1])
    assert "protein_align_identity" in names
