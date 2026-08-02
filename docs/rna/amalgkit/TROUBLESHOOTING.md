# Amalgkit Pipeline Troubleshooting

This guide covers common issues, performance bottlenecks, and environment configurations encountered during large-scale RNA-seq runs (8k+ samples).

## Database & IO Issues

### IO Contention (Hanging Queries)
**Problem**: The `output/amalgkit` directory becomes IO-saturated during quantification, causing `sqlite3` queries on `pipeline_progress.db` to hang or deadlock.

**Solution**: 
1. **Read-Only Mode**: Always use URI-based read-only connections for monitoring.
   ```python
   conn = sqlite3.connect('file:output/amalgkit/pipeline_progress.db?mode=ro', uri=True)
   ```
2. **Copy to TMP**: For heavy analysis, copy the DB to a local fast drive (SSD or `/tmp`).
   ```bash
   cp output/amalgkit/pipeline_progress.db /tmp/pipeline_status.db
   ```
3. **Avoid Pandas `read_sql` on live DBs**: Standard `pandas.read_sql` can trigger locking issues. Prefer raw `sqlite3` cursors with short timeouts.

## Tool & Environment Setup

### SRA Toolkit: `fasterq-dump` Not Found
**Problem**: Amalgkit fails to find `fasterq-dump` despite the SRA Toolkit being installed.

**Solution**: Ensure the binary path is explicitly in the environment within the Dockerfile or orchestration script.
```bash
export PATH=$PATH:/usr/local/ncbi/sra-tools/bin
```

### Native Amalgkit taxonomy-cache preparation
**Problem**: Metadata or acquisition startup spends excessive time rebuilding the
shared NCBI taxonomy cache, or multiple species contend for the same cache.

**Solution**: Use the native Amalgkit taxonomy cache under the selected data
root and let the current launcher share it across species.
```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
export AMALGKIT_SHARED_DOWNLOAD_DIR="$AMALGKIT_DATA_ROOT/shared/ncbi_taxonomy"
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```
The launcher records the cache provenance beside the native database. Do not
seed a second taxonomy implementation or copy an incomplete database into a
species work directory.

## Data Integrity

### Tissue Patch Inconsistency
**Problem**: Inconsistent BioProject-to-Tissue mappings in `config/amalgkit/tissue_patches.yaml`.

**Solution**: 
- Use `uv run pytest tests/rna/test_rna_tissue_normalization.py` to validate
  patches before bulk processing.
- Ensure no duplicate BioProject entries exist in the `NCBI Batch Patches` section.

### Ortholog Generation Failures
**Problem**: Missing OrthoDB mappings for new species.

**Solution**: Refer to [Ortholog Generation Guide](ortholog_generation.md). Ensure the automated extraction script handles species name synonyms correctly (e.g., `Apis mellifera` vs `apis_mellifera`).
