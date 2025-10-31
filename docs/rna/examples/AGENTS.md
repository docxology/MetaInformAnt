# AI Agents in RNA Example Documentation

This document records AI assistance for the RNA workflow example guides.

## AI Contributions

### Scenario Design
**Documentation Agent** assembled the tutorial sequence covering:
- Representative amalgkit-driven runs for *Pogonomyrmex barbatus*
- Step-by-step operator checklists for metadata, download, quantification, and merge phases
- Cross-references to reusable configuration patterns

### Content Drafting
**Documentation Agent** prepared:
- Narrative walkthroughs with CLI and Python snippets
- Quick reference tables summarizing required inputs and expected outputs
- Troubleshooting notes for common failure modes in high-throughput RNA processing

### Technical Validation
**Code Assistant Agent** reviewed:
- Command sequences for compatibility with `uv run` workflows
- Integration points with `scripts/rna/batch_ena.py`
- Consistency with repository I/O conventions (`output/` placement, `metainformant.core.io` usage)

## Maintenance Approach
- Regenerate command captures using real amalgkit runs before updating screenshots or logs
- Keep example datasets aligned with test fixtures in `tests/data/rna`
- Flag any AI-assisted changes for human verification prior to publishing


