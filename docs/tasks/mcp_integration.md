# MCP Integration Task Notes

Status: draft. The current checkout does not implement a METAINFORMANT MCP
server. Use this page as planning context, not as runnable setup instructions.

## Current Runnable Interface

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor
```

This standalone monitor inspects local RNA/Amalgkit workflow progress. It is the
only implemented MCP-adjacent command today.

## Unavailable Until Implemented

- `metainformant.mcp.server`
- Claude Desktop/Cursor MCP server config for METAINFORMANT
- `run_workflow` MCP tool
- `list_outputs` MCP tool
- stdio or SSE transport commands

## Future Acceptance Criteria

1. Add a real server module and transport tests.
2. Add a tool registry with JSON-serializable request/response contracts.
3. Implement and test each tool before documenting it as available.
4. Update [docs/mcp/](../mcp/) and this task guide with runnable commands.

Until then, workflow execution remains available through the Python APIs,
`uv run metainformant gwas run --help`, and the domain scripts under `scripts/`.
