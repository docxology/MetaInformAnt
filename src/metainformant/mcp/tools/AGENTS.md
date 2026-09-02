# AGENTS.md — `src/metainformant/mcp/tools`

Standalone MCP-adjacent tool modules for the `mcp` domain module. Each
`*_tools.py` module exposes a module-level `TOOL_SPEC` (or `ALL_SPECS` list)
following the shared contract in `_spec.py`: `name`, `description`,
`input_schema` (JSON-schema object), `handler` (`**kwargs -> dict`), and
`writes` (`"read-only"` or `"output-dir-only"`).

`catalog.py` adapts every spec into `metainformant.mcp.registry.Tool` records
for `ToolRegistry.register_module_tools("metainformant.mcp.tools.catalog")`;
the stdio server currently registers via the registry lane's
`metainformant.mcp.tool_adapters` (amalgkit_monitor only), so exposing the
catalog through the server is a one-entry addition to `TOOLS_MODULES` in
`metainformant/mcp/__init__.py`.

Invariants:
- Handlers are deterministic; identical inputs -> identical outputs.
- Read-only tools touch only caller-supplied paths; writing tools create
  files only under a caller-supplied `output_dir`
  (`_spec.validate_output_dir`).
- RNA campaign/analysis outputs are DESCRIPTIVE ONLY until the evidence
  manifest freezes (no inferential statistics in tool outputs).
- `amalgkit_monitor.py` predates the spec contract and is wired by
  `metainformant/mcp/tool_adapters.py`; do not duplicate it in the catalog.

Tests: `tests/mcp/test_tools_specs.py` (contract),
`test_tools_handlers.py` (behavior, real synthetic data),
`test_tools_subprocess.py` (fresh-interpreter round trips),
`test_tools_catalog.py` (registry round trip + stdio listing).
Run one test file per pytest invocation.

Repo-wide policy: see the repository-root `AGENTS.md`.
