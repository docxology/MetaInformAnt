# AI Agents in Apis mellifera Curation Tables

This document captures AI assistance specific to the tabular fixtures for *Apis mellifera* curation tests.

## AI Contributions

### Fixture Generation
**Code Assistant Agent** generated:
- `Apis_mellifera.metadata.tsv` with representative sample annotations and edge-case values
- `Apis_mellifera.uncorrected.tc.tsv` expressing raw transcript counts for normalization workflows
- Data blends designed to trigger QC thresholds inside `metainformant.rna.curate`

### Documentation
**Documentation Agent** authored the accompanying README to:
- Explain row/column semantics and expected formatting
- Describe relationships between metadata and count matrices
- Provide maintenance reminders tied to RNA workflow evolution

### Validation
**Code Assistant Agent** ensured:
- TSV schemas stay in sync with parser expectations
- Deterministic ordering for regression comparisons
- Compatibility with downstream analytics used in integration tests

## Maintenance Notes
- Refresh these tables when reference genomes or feature annotations change materially
- Keep the datasets compact to preserve fast test execution
- Mark AI-driven updates for human verification prior to merging

