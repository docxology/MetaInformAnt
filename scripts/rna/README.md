# RNA-seq Processing Scripts

Scripts for RNA-seq data processing, including SRA download and quantification.

## Available Scripts

### `batch_ena.py`
**ENA-based parallel downloader and quantifier**

Fast RNA-seq processing using ENA (European Nucleotide Archive) for downloads:
- Direct FASTQ downloads (no SRA conversion)
- 100-500X faster than NCBI
- Parallel downloads and quantifications
- Automatic cleanup

**Configuration:**
```python
MAX_CONCURRENT_DOWNLOADS = 5  # Parallel sample downloads
MAX_CONCURRENT_QUANTS = 3     # Parallel quantifications
KALLISTO_THREADS = 3          # Threads per kallisto process
```

**Usage:**
```bash
# Run from project root
cd /path/to/metainformant
python3 scripts/rna/batch_ena.py
```

**Requirements:**
- wget (for downloads)
- kallisto (for quantification)
- Sample list in `output/amalgkit/{species}/remaining_samples.txt`
- Kallisto index in `output/amalgkit/{species}/work/index/`

### `restart_batch.sh`
**Convenience script to restart batch processing**

Automatically handles:
- Background execution
- PID tracking
- Log file management

**Usage:**
```bash
bash scripts/rna/restart_batch.sh
```

## Output Structure

All outputs go to `output/amalgkit/{species}/`:
- `work/quant/` - Quantification results (abundance.tsv)
- `work/metadata/` - Filtered metadata
- `logs/` - Processing logs
- `*.json`, `*.log` - Status files

See `docs/rna/examples/` for complete documentation.

## Examples

See `docs/rna/examples/pbarbatus_analysis.md` for a complete workflow example.

