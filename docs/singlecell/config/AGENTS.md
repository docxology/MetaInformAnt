# AI Agents in Single-Cell Configuration Documentation

This document captures AI participation in drafting and maintaining single-cell configuration references.

## AI Contributions

### Configuration Cataloging
**Documentation Agent** compiled:
- Template sections for QC, normalization, clustering, and trajectory setups
- Checklists for environment variables and dataset metadata requirements
- Cross-links between documentation, YAML configs, and CLI usage

### Technical Review
**Code Assistant Agent** verified:
- Alignment between documented parameters and `metainformant.singlecell` function signatures
- Compatibility of example commands with the `uv` execution model
- Consistency with repository path and I/O rules (`output/` destinations, cache usage)

## Maintenance Notes
- Reconfirm defaults whenever single-cell modules gain new options or change behavior
- Annotate AI-authored updates for human review prior to merging
- Sync configuration guidance with examples in `docs/singlecell/` and workflows in `src/metainformant/singlecell`

