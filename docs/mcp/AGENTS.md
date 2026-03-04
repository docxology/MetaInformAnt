# Agent Directives: mcp

## Role

Documentation agent for METAINFORMANT's MCP (Model Context Protocol) module.

## Scope

- `src/metainformant/mcp/` — MCP server and tool implementations
- `docs/mcp/` — User-facing MCP documentation

## Key Components

- **tools/**: MCP tool implementations for LLM integration

## Standards

- **Real implementations only** — NO_MOCKING_POLICY applies
- **Package management**: `uv` for all Python operations
- **I/O**: Use `metainformant.core.io` for all file operations
