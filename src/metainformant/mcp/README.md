# MCP Package

This package contains MCP-adjacent helper code for METAINFORMANT. The current
public interface is a standalone monitor script; a full MCP server is not
implemented in this checkout.

## Public Surface

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor
```

```python
from metainformant.mcp.tools import amalgkit_monitor
```

## Available Modules

| Module | Purpose | Status |
|--------|---------|--------|
| `tools.amalgkit_monitor` | Inspect local Amalgkit/RNA workflow progress | Implemented |

## Deferred Interfaces

Do not document or depend on these names until implementation and tests land:

- `metainformant.mcp.server`
- `run_workflow` as an MCP tool
- `list_outputs` as an MCP tool

## Related

- [docs/mcp/](../../../docs/mcp/)
- [RNA package](../rna/)
