# MCP

The MCP area is currently a small helper package for monitoring local
RNA/Amalgkit runs. It is not yet a complete Model Context Protocol server.

## Current Command

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor
```

Use this command to inspect local Amalgkit process/log status from a terminal or
from wrapper scripts.

## Current Python API

```python
from metainformant.mcp.tools import amalgkit_monitor

progress = amalgkit_monitor.parse_log_progress("output/amalgkit/species/logs")
```

## Planned, Not Available

The following names appear in older task notes but are not implemented:

- `metainformant.mcp.server`
- `run_workflow`
- `list_outputs`
- stdio/SSE MCP transports

## Development Checklist

Before documenting a true MCP server as available, add source implementation,
tool registration tests, transport tests, and a CLI smoke test.

## Related

- [README](README.md)
- [SPEC](SPEC.md)
- [RNA workflow docs](../rna/)
