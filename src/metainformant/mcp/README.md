# MCP Module (Model Context Protocol)

METAINFORMANT tools exposed via Model Context Protocol for LLM-based development environments.

## Overview

METAINFORMANT tools exposed via Model Context Protocol for LLM-based development environments.


## Table of Contents

- [Quick Start](#quick-start)
  - [Claude Desktop](#claude-desktop)
  - [Cursor](#cursor)
- [Available Tools](#available-tools)
- [See Also](#see-also)

## Quick Start

### Claude Desktop

1. Install: `uv pip install -e .`
2. Create `~/.config/Claude/mcp_config.json`:
```json
{
  "mcpServers": {
    "metainformant": {
      "command": "uv",
      "args": ["run", "python", "-m", "metainformant.mcp.server"]
    }
  }
}
```
3. Restart Claude Desktop

### Cursor

Settings → AI Agents → MCP Servers → Add:
- Command: `uv`
- Args: `["run", "python", "-m", "metainformant.mcp.server"]`

## Available Tools

| Tool | Description | Status |
|------|-------------|--------|
| `amalgkit_status` | Check RNA pipeline status for species | [DONE] Implemented |
| `run_workflow` | Trigger METAINFORMANT workflow | [PARTIAL] |
| `list_outputs` | List output files for module | [PLANNED] |

## See Also

- [RNA Module](../rna/) — Amalgkit pipeline
- [MCP SPEC](SPEC.md) — Technical API reference
