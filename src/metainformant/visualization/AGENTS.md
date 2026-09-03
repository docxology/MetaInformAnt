# Agent Directives: visualization

**Context**: Visualization and plotting utilities module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `dimred`, `information`, `math_plots`, `quality`, `quality_assessment`, `quality_omics`, `simulation_plots`
- `config/` — exports: `palettes`, `themes`
- `dashboards/` — exports: `composite`, `interactive`
- `genomics/` — exports: `expression`, `genomics`, `networks`, `trees`
- `interactive_dashboards/` — exports: `dashboards`
- `plots/` — exports: `animations`, `basic`, `general`, `multidim`, `specialized`

## Shared chart conventions

- `config/conventions.py` is the canonical home for the Okabe-Ito palette,
  category hatches, log-axis zero marks, and deterministic figure export
  (`save_figure_deterministic`). All module plotting code routes figure saves
  through it; other palette names (`WONG`, theme cyclers) are derived from it.
  The single-source invariant is pinned by `tests/visualization/test_palette_single_source.py`.

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/visualization/](../../../docs/visualization/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
- **Dna module**: [../dna/AGENTS.md](../dna/AGENTS.md)
- **Rna module**: [../rna/AGENTS.md](../rna/AGENTS.md)
- **Gwas module**: [../gwas/AGENTS.md](../gwas/AGENTS.md)
- **Quality module**: [../quality/AGENTS.md](../quality/AGENTS.md)
