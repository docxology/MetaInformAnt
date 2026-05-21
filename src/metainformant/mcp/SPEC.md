# MCP Package Specification

## Module

`metainformant.mcp`

## Status

Experimental helper package. There is no MCP server entry point in the current
checkout.

## Implemented Interface

### `metainformant.mcp.tools.amalgkit_monitor`

Command:

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor
```

Responsibilities:

- inspect local Amalgkit/RNA pipeline processes;
- parse available log progress;
- return a small status payload suitable for terminal use and future adapters.

## Package Exports

```python
from metainformant.mcp import tools
from metainformant.mcp.tools import amalgkit_monitor
```

## Non-Goals For This Pass

- no `metainformant.mcp.server` module;
- no stdio/SSE MCP transport;
- no MCP tool registry;
- no `run_workflow` or `list_outputs` MCP tools.

Those interfaces must stay undocumented as available until they have source
implementation and tests.
