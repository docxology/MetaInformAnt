# MCP Notes

METAINFORMANT currently ships a lightweight MCP-adjacent helper module, not a
full Model Context Protocol server.

## Current Interface

The implemented command is the standalone Amalgkit monitor:

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor
```

It inspects explicitly selected local RNA/Amalgkit workflow state and reports
process/log diagnostics alongside database- and receipt-backed readiness.
The helper can be imported from `metainformant.mcp.tools.amalgkit_monitor`.

## Not Yet Implemented

The checkout does not provide:

- `metainformant.mcp.server`
- MCP stdio or SSE transports
- registered MCP tools named `run_workflow` or `list_outputs`

Keep examples and integrations on this page limited to the standalone monitor
until a real server module and tests are added.

The monitor is not a completion oracle: it separates executable readiness,
cohort readiness, descriptive analysis, and biological inference, and always
withholds the last field. Use `--no-process-scan` when a lock/receipt/database
snapshot is preferred over process inspection.

## See Also

- [MCP SPEC](SPEC.md)
- [RNA docs](../rna/)
- [Source README](../../src/metainformant/mcp/README.md)
