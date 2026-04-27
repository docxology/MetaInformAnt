# MCP Module Technical Specification

## Module: `metainformant.mcp`

**Status:** Minimal (40% complete)
**Version:** 0.3.0

## Public API

### `server.py` — Entry Point

```python
def main() -> None:
    """Start MCP server (stdio or SSE transport)."""
```

**CLI:**
```bash
python -m metainformant.mcp.server --transport stdio
python -m metainformant.mcp.server --transport sse --port 8080
```

### Tool Registration

```python
from metainformant.mcp import register_tool

@register_tool(
    name="amalgkit_status",
    description="Get RNA pipeline status",
    input_schema={"type": "object", "properties": {...}}
)
def amalgkit_status(species: str) -> dict:
    ...
```

## Transport Modes

1. **STDIO** (default): Pipe-based for desktop LLMs
2. **SSE** (planned): HTTP endpoint for web-based LLMs

## See Also

- [MCP Specification](https://spec.modelcontextprotocol.io)
- [RNA Module](../rna/amalgkit/)

---

**Changelog**

v0.3.0 (2025-04-15): Initial scaffold, STDIO transport, amalgkit_status tool
