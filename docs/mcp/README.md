# MCP Notes

METAINFORMANT ships a Model Context Protocol server: `metainformant.mcp.server`
implements a stdio JSON-RPC 2.0 server with a declarative, schema-validated
tool registry (standard library only).

## Server

```bash
uv run python -m metainformant.mcp.server
```

Implemented methods: `initialize`, `tools/list`, `tools/call`,
`resources/list`, and `resources/read`, plus graceful shutdown via EOF or the
`exit` notification. Read-only resources: `metainformant://capabilities`
(registered tools and schemas) and `metainformant://methods` (connection
documentation).

## Registry

`metainformant.mcp.registry` provides `Tool`, `ToolRegistry`, and `SchemaError`
with JSON-schema validation at registration and argument-validation at call
time. `build_default_registry()` bundles the `amalgkit_monitor` adapter plus
the 20 `TOOL_SPEC` tools from `metainformant.mcp.tools` (core, dna, gwas,
math, protein, rna, visualization).

## Standalone monitor

The Amalgkit monitor still runs standalone and feeds the default registry:

```bash
uv run python -m metainformant.mcp.tools.amalgkit_monitor \
  --data-root "$AMALGKIT_DATA_ROOT" --no-process-scan
```

Flags: `--data-root` (default `AMALGKIT_DATA_ROOT` or `output/amalgkit`),
`--log-file` (explicit log file path), and `--no-process-scan` for a
lock/receipt/database-only snapshot.

The monitor is not a completion oracle: it separates executable readiness,
cohort readiness, descriptive analysis, and biological inference, and always
withholds the last field.

## See Also

- [MCP SPEC](SPEC.md)
- [RNA docs](../rna/)
- [Source README](../../src/metainformant/mcp/README.md)
