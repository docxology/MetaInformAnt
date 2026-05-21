# Agent Directives: mcp (Model Context Protocol)

## Role

Maintain truthful documentation for METAINFORMANT's MCP-adjacent helper package.
The current checkout has a standalone monitor script, not a protocol server.

## Module Scope

| File/Dir | Purpose |
|-----------|---------|
| `tools/` | Standalone helper modules |
| `__init__.py` | Package exports |

## Key Source Files

- **Tools**: `src/metainformant/mcp/tools/amalgkit_monitor.py` — local Amalgkit/RNA progress monitor

## Related Documentation

- **Module guide**: [index.md](index.md) — current status and future checklist
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **API spec**: [SPEC.md](SPEC.md) — current helper surface and unavailable interfaces
- **Amalgkit**: [../rna/amalgkit/](../rna/amalgkit/) — The pipeline that MCP tools monitor

## Rules & Constraints

1. Do not document `metainformant.mcp.server` until the module exists.
2. Do not document `run_workflow` or `list_outputs` as MCP tools until they are implemented and tested.
3. Helper output must be JSON-serializable if it is intended for future adapters.
4. Keep tools stateless and idempotent.

## Cross-Module Dependencies

- **rna**: RNA/Amalgkit workflow logs and process state are the current monitoring target.
- **core.io**: use canonical helpers for domain data file I/O.

## Development Notes

Current implementation: standalone `amalgkit_monitor` helper only.
