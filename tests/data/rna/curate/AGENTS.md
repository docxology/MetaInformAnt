# AI Agents in RNA Curation Test Data

This document summarizes AI support for the RNA curation fixtures used in testing.

## AI Contributions

### Data Assembly
**Code Assistant Agent** produced:
- Synthetic expression matrices and metadata JSON aligned with curation pipelines
- Edge-case datasets targeting quality control branches and batch-effect handling
- Expected-results fixtures for regression testing of the curation stack

### Documentation
**Documentation Agent** drafted:
- The README outlining directory structure and usage patterns
- Notes on data provenance and maintenance needs
- Cross-links to RNA workflow documentation

### Validation
**Code Assistant Agent** validated:
- Schema compatibility with `metainformant.rna.curate`
- Deterministic checksums for regression comparisons
- Coverage of error conditions invoked by integration tests

## Maintenance Checklist
- Rebuild synthetic datasets when pipeline parameters change significantly
- Keep fixtures lightweight while still exercising failure paths
- Mark AI-authored updates for human review before merging

