# tools

Standalone MCP-adjacent tool modules for `metainformant.mcp`.

- Contract: `_spec.py` (`TOOL_SPEC` shape; `validate_output_dir`, `read_table`, `dump_json` helpers).
- Modules: `core_tools`, `dna_tools`, `gwas_tools`, `math_tools`, `protein_tools`, `rna_tools`, `visualization_tools` (one `TOOL_SPEC`/`ALL_SPECS` each).
- Bridge: `catalog.py` -> `metainformant.mcp.registry.Tool` for `register_module_tools`.
- Legacy: `amalgkit_monitor.py` (wired separately by `metainformant/mcp/tool_adapters.py`).

Every tool declares `writes: "read-only" | "output-dir-only"`; deterministic
outputs; descriptive-only statistics on RNA campaign surfaces.
