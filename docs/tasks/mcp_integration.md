# MCP Integration Quick Reference

Use METAINFORMANT from Claude Desktop, Cursor, Windsurf, or any MCP-compatible LLM. Exposes bioinformatics workflows as JSON-RPC tools.

## When to Use

Use `mcp_integration` when you want LLM-driven workflow orchestration (natural language → pipeline execution) — not for direct CLI use or programmatic Python API calls.

## Table of Contents

- [Setup (3 minutes)](#setup-3-minutes)
- [Example Prompts](#example-prompts)
- [Available Tools (v0.3.0)](#available-tools-v030)
- [Debugging](#debugging)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

1. **Install metainformant with MCP extras:**

```bash
uv pip install -e ".[mcp]"
```

2. **Create MCP config:**

**Claude Desktop** (`~/.config/Claude/mcp_config.json`):
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

**Cursor** (Settings → AI → MCP Servers):
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

3. **Restart your LLM app**.

## Example Prompts

**Status query:**
> "Is the Amel RNA-seq pipeline done?"

→ Calls `amalgkit_status(species="Apis mellifera")` → Returns progress + output location

**Run analysis:**
> "Run GWAS on the filtered VCF with phenotype file `phenotypes.tsv`"

→ Calls `run_workflow(module="gwas", config={"vcf": "...", "pheno": "..."})`

**Get results:**
> "List all output files from the last run"

→ Calls `list_outputs(workflow_id="...")` → Returns file list with sizes

**Create plot:**
> "Make a Manhattan plot of the GWAS results"

→ Future tool → Returns plot file path

## Available Tools (v0.3.0)

```bash
# List tools in terminal
uv run python -m metainformant.mcp.server --list

# Manually invoke tool
echo '{"tool":"amalgkit_status","arguments":{"species":"Apis mellifera"}}' \
  | uv run python -m metainformant.mcp.server --transport stdio
```

2. **Create MCP config:**

**Claude Desktop** (`~/.config/Claude/mcp_config.json`):
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

**Cursor** (Settings → AI → MCP Servers):
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

3. **Restart your LLM app.**

## Example Prompts

**Status query:**
> "Is the Amel RNA-seq pipeline done?"

→ Calls `amalgkit_status(species="Apis mellifera")` → Returns progress + output location

**Run analysis:**
> "Run GWAS on the filtered VCF with phenotype file `phenotypes.tsv`"

→ Calls `run_workflow(module="gwas", config={"vcf": "...", "pheno": "..."})`

**Get results:**
> "List all output files from the last run"

→ Calls `list_outputs(workflow_id="...")` → Returns file list with sizes

**Create plot:**
> "Make a Manhattan plot of the GWAS results"

→ Future tool → Returns plot file path

## Available Tools (v0.3.0)

```bash
# List tools in terminal
uv run python -m metainformant.mcp.server --list

# Manually invoke tool
echo '{"tool":"amalgkit_status","arguments":{"species":"Apis mellifera"}}' \
  | uv run python -m metainformant.mcp.server --transport stdio
```

## Advanced Examples

### Custom tool registration (extend MCP)
```python
# In your custom metainformant extension module
from metainformant.mcp import register_tool, ToolSpec

@register_tool(
    name="run_custom_analysis",
    description="Run custom QC on RNA-seq samples",
    parameters={
        "type": "object",
        "properties": {
            "samples_dir": {"type": "string", "description": "Path to FASTQ directory"},
            "min_quality": {"type": "number", "default": 20}
        }
    }
)
def run_custom_analysis(samples_dir: str, min_quality: float = 20.0) -> dict:
    """Custom analysis exposed to MCP clients."""
    from metainformant.quality import qc
    results = qc.fastq_quality_check(samples_dir, min_q=min_quality)
    return {"passed": results.passed, "failed": results.failed, "report": results.report_path}
```

### Multi-tool workflow via LLM orchestration
```bash
# User tells Claude: "Compare RNA-seq expression vs GWAS hits for bee immunity genes"
# Claude auto-generates and executes:

# Tool 1: get_gwas_hits(threshold=1e-5)
# Tool 2: load_expression_matrix("output/amellifera/counts.tsv")
# Tool 3: subset_genes(gwas_hits)
# Tool 4: plot_volcano(expression_subset)
# Tool 5: save_report("immunity_comparison.html")
```

### MCP over HTTP (SSE) for web integration
```bash
# Start MCP server with SSE transport (for web clients)
uv run python -m metainformant.mcp.server --transport sse --port 8080

# In JavaScript frontend:
const client = new EventSource("http://localhost:8080/events");
client.onmessage = (event) => {
  const data = JSON.parse(event.data);
  console.log("Tool call result:", data);
};
```

## Expected Output

### Tool list JSON
```json
{
  "tools": [
    {
      "name": "amalgkit_status",
      "description": "Get status of-running RNA-seq pipeline",
      "inputSchema": {
        "type": "object",
        "properties": {
          "species": {"type": "string", "description": "Species name"}
        },
        "required": ["species"]
      }
    },
    {
      "name": "run_workflow",
      "description": "Execute METAINFORMANT workflow",
      "inputSchema": {...}
    }
  ]
}
```

### amalgkit_status tool response
```json
{
  "species": "Apis mellifera",
  "status": "COMPLETE",
  "output_dir": "/data/output/amellifera/",
  "elapsed_hours": 4.3,
  "samples_processed": 8,
  "files": {
    "counts.tsv": "11.2 MB",
    "deseq2_results.tsv": "0.8 MB",
    "pca.png": "0.4 MB"
  }
}
```

### run_workflow tool response (streaming)
```json
{
  "workflow_id": "gwas-20260426-001",
  "status": "RUNNING",
  "progress": {"completed": 3, "total": 16, "pct": 18.75},
  "current_step": "Chromosome 3 association testing",
  "logs_url": "https://console.cloud.google.com/logs/..."
}
```

## Debugging

```bash
# Run server in foreground with debug logging
METAINFORMANT_LOG_LEVEL=DEBUG \
uv run python -m metainformant.mcp.server --transport stdio

# Test with mcp-cli (install via npm)
mcp list-tools --transport stdio --command "uv run python -m metainformant.mcp.server"

# Inspect JSON-RPC messages directly
echo '{"jsonrpc":"2.0","id":1,"method":"tools/list","params":{}}' \
  | uv run python -m metainformant.mcp.server --transport stdio
```

### Debug session example
```bash
$ METAINFORMANT_LOG_LEVEL=DEBUG uv run python -m metainformant.mcp.server --transport stdio
[DEBUG] 2026-04-26 14:02:01 MCP server starting (transport=stdio)
[DEBUG] 2026-04-26 14:02:01 Loading tool registry from metainformant.mcp.tools
[INFO ] 2026-04-26 14:02:01 Registered 3 tools: amalgkit_status, run_workflow, list_outputs
[DEBUG] 2026-04-26 14:02:01 Waiting for JSON-RPC requests on stdin...
<-- {"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"amalgkit_status",...}}
--> {"jsonrpc":"2.0","id":1,"result":{"status":"RUNNING","progress":0.42,...}}
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `"tool not found"` error in Claude Desktop | MCP server not running or `mcp_config.json` syntax error | Verify `uv run python -m metainformant.mcp.server --list` works; validate JSON with `jq . mcp_config.json` |
| `Connection refused` on MCP startup | Port 8080 already in use (SSE mode) | Kill existing server (`pkill -f metainformant.mcp.server`); use `--port 8081` |
| No tools appear in Cursor | Cursor's MCP integration stale or sandboxed | Restart Cursor completely; check `~/.config/Cursor/mcp_config.json`; ensure `uv` in PATH |
| Timeout after 300s (5min) | Long-running pipeline exceeds default MCP timeout | Increase `METAINFORMANT_MCP_TIMEOUT=3600` env var; or use async tool that returns job_id immediately |
| `stdin/stdout pipes closed` | Claude Desktop stdout buffer filled (server blocking on print) | Flush logs; use non-blocking I/O; reduce log verbosity |
| Tool returns `None` or empty | Underlying function raised exception not caught | Check server logs (DEBUG level); run same function directly in Python REPL to reproduce |

---

**Related:** [MCP module docs](../mcp/index.md) | [MCP SPEC](../mcp/SPEC.md) | [Claude Desktop MCP setup](https://docs.anthropic.com/claude/docs/claude-desktop-mcp) | [Cursor MCP docs](https://docs.cursor.com/context/model-context-protocol)
