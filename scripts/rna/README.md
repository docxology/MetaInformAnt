# RNA-seq Processing Scripts

Scripts for RNA-seq data processing, including SRA download, quantification, and multi-species workflows.

## Directory Structure

```
scripts/rna/
├── amalgkit/                      # Amalgkit-specific workflow scripts
│   ├── run_amalgkit.sh           # Comprehensive pipeline orchestrator
│   ├── verify_workflow.sh        # Workflow validation
│   └── README.md                 # Detailed amalgkit documentation
├── batch_ena.py                  # Fast ENA parallel downloader
├── cleanup_quantified_sra.sh     # Safe deletion of quantified SRA files
├── force_fasterq.sh              # SRA FASTQ processing
├── force_fasterq_parallel.sh     # Parallel FASTQ processing
├── list_unquantified.sh          # Report unquantified samples
├── monitor_amalgkit_progress.sh  # Multi-species progress monitor
├── monitor_workflow.py           # Real-time monitoring dashboard
├── process_one_srr.sh            # Single SRR processing
├── quant_downloaded_samples.py   # Quantify downloaded samples
├── run_multi_species.py          # Multi-species with cross-species analysis (RECOMMENDED)
├── test_pbarbatus_workflow.py    # P. barbatus workflow testing
├── test_single_species.py        # Single species testing
├── test_skip_logic.py            # Skip logic verification
├── verify_skip_logic.sh          # Comprehensive skip verification
├── README.md                     # This file
└── AGENTS.md                     # AI agent documentation
```

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

### Multi-Species Workflow Script

#### `run_multi_species.py` ⭐ **RECOMMENDED**
Complete workflow for all species with batched processing and cross-species analysis:

**Features:**
- **Auto-discovery**: Finds all `config/amalgkit/amalgkit_*.yaml` files (excludes template)
- **Batched processing**: Downloads 8 samples → Quantifies → Deletes FASTQs → Repeats
- **Disk-friendly**: Only ~16-40 GB peak disk usage (8 samples worth of FASTQs)
- **Fast**: Parallel processing within each batch (8 threads by default)
- **Resume**: Automatically skips already-quantified samples
- **Cross-species analysis**: Runs CSTMM and CSCA for multi-species comparisons
- **Complete pipeline**: metadata → select → getfastq → quant → merge → curate → sanity

**Usage:**
```bash
# Process all discovered species (pbarbatus, cfloridanus, mpharaonis, sinvicta)
python3 scripts/rna/run_multi_species.py
```

**Current Species Configured:**
- Pogonomyrmex barbatus (Red harvester ant) - 83 samples
- Camponotus floridanus (Florida carpenter ant) - 307 samples
- Monomorium pharaonis (Pharaoh ant) - 100 samples
- Solenopsis invicta (Red fire ant) - 354 samples

**Processing Mode:**
```
Batch 1: [Download 8] → [Quantify 8] → [Delete FASTQs] → Disk freed
Batch 2: [Download 8] → [Quantify 8] → [Delete FASTQs] → Disk freed
...
```

**Performance:**
- ~8x faster than strict sequential processing
- Disk usage bounded regardless of cohort size
- Tested with 300+ samples per species

### Monitoring Scripts

#### `monitor_workflow.py`
Real-time monitoring dashboard with:
- Species progress tracking
- Disk usage monitoring
- Running process detection
- Batch activity logs

**Usage:**
```bash
python3 scripts/rna/monitor_workflow.py
```

#### `monitor_amalgkit_progress.sh`
Multi-species progress monitor (simple version):

**Usage:**
```bash
bash scripts/rna/monitor_amalgkit_progress.sh
```

### `list_unquantified.sh`
**Generate reports of samples needing quantification**

Identifies samples with downloaded SRA files but no quantification output:
- Scans all species directories
- Reports size and sample count
- Creates sample lists in output/amalgkit/

**Usage:**
```bash
bash scripts/rna/list_unquantified.sh
```

**Output:**
- `output/amalgkit/{species}_unquantified.txt` - Sample lists
- Console report with sizes and counts

### `cleanup_quantified_sra.sh`
**Safe deletion of quantified SRA files**

Reclaims disk space by removing SRA files after successful quantification:
- Verifies quantification completion before deletion
- Detailed logging of operations
- Safe: skips unquantified samples

**Usage:**
```bash
# Preview what will be deleted (dry run)
bash scripts/rna/cleanup_quantified_sra.sh

# Execute cleanup
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

### Testing Scripts

#### `test_pbarbatus_workflow.py`
Test script for P. barbatus workflow:
- Single sample end-to-end testing
- Disk usage verification
- Workflow execution with sample limits

**Usage:**
```bash
python3 scripts/rna/test_pbarbatus_workflow.py --check-only
python3 scripts/rna/test_pbarbatus_workflow.py --max-samples 3
```

#### `test_single_species.py`
Single species workflow testing:
```bash
python3 scripts/rna/test_single_species.py
```

#### `test_skip_logic.py` and `verify_skip_logic.sh`
Verify skip logic for already-quantified samples:
```bash
python3 scripts/rna/test_skip_logic.py
bash scripts/rna/verify_skip_logic.sh
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

