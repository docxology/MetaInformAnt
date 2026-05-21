# Agent Directives: mcp

**Context**: MCP-adjacent helper package for METAINFORMANT.

## Capabilities

The current checkout provides a standalone Amalgkit monitor. It does not provide
an MCP server, workflow orchestration tool, or output-listing tool yet.

## Subpackages

### tools/

| File | Purpose |
|------|---------|
| `amalgkit_monitor.py` | Amalgkit pipeline status monitoring tool |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
- Tools must be stateless and idempotent
- Do not document `metainformant.mcp.server`, `run_workflow`, or `list_outputs` as available until implemented and tested
