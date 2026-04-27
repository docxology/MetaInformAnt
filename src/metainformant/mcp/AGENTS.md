# Agent Directives: mcp

**Context**: Model Context Protocol (MCP) tool implementations for METAINFORMANT.

## Capabilities

MCP-compliant tools that expose METAINFORMANT functionality to LLM-based development environments: pipeline status monitoring, workflow orchestration queries, and data retrieval.

## Subpackages

### tools/

| File | Purpose |
|------|---------|
| `amalgkit_monitor.py` | Amalgkit pipeline status monitoring tool |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
- Tools must be stateless and idempotent
