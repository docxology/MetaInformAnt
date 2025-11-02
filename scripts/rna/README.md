# RNA-seq Processing Scripts

Scripts for RNA-seq data processing, including SRA download, quantification, and multi-species workflows.

## Directory Structure

```
scripts/rna/
├── amalgkit/                       # Amalgkit-specific workflow scripts
│   ├── run_amalgkit.sh            # Comprehensive pipeline orchestrator
│   ├── verify_workflow.sh         # Workflow validation
│   └── README.md                  # Detailed amalgkit documentation
├── workflow_ena_integrated.py     # Integrated ENA download + quantification (PRODUCTION) ⭐
├── download_ena_robust.py         # Robust ENA downloader with retry logic ⭐
├── run_multi_species.py           # Multi-species with cross-species analysis (legacy SRA-based)
├── check_environment.py           # Environment validation
├── cleanup_quantified_sra.sh      # Safe deletion of FASTQ files after quantification
├── monitor_comprehensive.py       # Comprehensive real-time monitoring
├── monitor_workflow.py            # Real-time monitoring dashboard
├── monitor_amalgkit_progress.sh   # Simple progress monitor
├── list_unquantified.sh           # Report unquantified samples
├── quant_downloaded_samples.py    # Quantify already-downloaded samples
├── README.md                      # This file
└── AGENTS.md                      # AI agent documentation
```

## Available Scripts

### `workflow_ena_integrated.py` ⭐ **PRODUCTION**
**Robust integrated download + quantification workflow**

Production-ready workflow bypassing SRA Toolkit issues:

**Features:**
- **Direct ENA downloads**: Fetches FASTQs directly from European Nucleotide Archive API
- **Robust retry logic**: Uses wget with --continue for resume capability and automatic retries
- **Batched processing**: Download N samples → Quantify → Delete FASTQs → Repeat
- **Disk-friendly**: Only one batch of FASTQs on disk at a time
- **Auto-detection**: Handles both single-end and paired-end data
- **Resume support**: Skips already-quantified samples automatically

**Why ENA over SRA Toolkit:**
- SRA Toolkit downloads fail frequently (~100% failure rate on large samples)
- ENA provides direct FASTQ files (no SRA→FASTQ conversion needed)
- wget --continue allows proper resumption after network interruptions
- Much more reliable for large-scale downloads

**Usage:**
```bash
# Full workflow with default settings (12 samples/batch, 12 threads)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12

# Test with 3 samples
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 3 \
  --threads 8 \
  --max-samples 3

# Resume (skip download, only quantify existing FASTQs)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --skip-download
```

**Performance:**
- Download: ~6 minutes per 3 samples (varies by sample size)
- Quantification: ~36 seconds per sample (single-end, 8 threads)
- Cleanup: Instant FASTQ deletion after quantification
- Peak disk: ~1.5 GB per sample (batch size × 1.5 GB)

**Requirements:**
- wget (for robust downloads)
- kallisto (for quantification)
- Kallisto index must exist (built once per species)

### `download_ena_robust.py` ⭐
**Standalone robust ENA downloader**

Used internally by `workflow_ena_integrated.py`, but can be used standalone:

**Usage:**
```bash
python3 scripts/rna/download_ena_robust.py \
  --metadata output/amalgkit/cfloridanus/work/metadata/metadata.tsv \
  --out-dir output/amalgkit/cfloridanus/fastq \
  --threads 12 \
  --max-retries 3
```

### `run_multi_species.py`
**Legacy multi-species workflow (SRA Toolkit-based)**

Alternative workflow using SRA Toolkit instead of ENA direct downloads:
- Auto-activation and environment management
- Batched processing (10 samples at a time)
- Cross-species analysis (CSTMM, CSCA)
- Complete pipeline automation

**Note:** Consider using `workflow_ena_integrated.py` for better reliability.

**Usage:**
```bash
# Process all discovered species
python3 scripts/rna/run_multi_species.py
```

### Monitoring Scripts

#### `monitor_comprehensive.py` ⭐
**Comprehensive real-time workflow monitor**

Tracks all 4 species simultaneously with detailed progress:
- Sample counts (quantified/total)
- Current batch numbers
- Downloading sample counts
- FASTQ directory sizes
- Time estimates

**Usage:**
```bash
python3 scripts/rna/monitor_comprehensive.py
```

#### `monitor_workflow.py`
**Alternative monitoring dashboard**

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
**Simple progress monitor**

Lightweight bash-based monitoring (no Python dependencies):

**Usage:**
```bash
bash scripts/rna/monitor_amalgkit_progress.sh

# Or watch mode
watch -n 60 bash scripts/rna/monitor_amalgkit_progress.sh
```

### Utility Scripts

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
**Safe deletion of FASTQ files after quantification**

Reclaims disk space by removing FASTQ files after successful quantification.
Works with both SRA Toolkit downloads (.sra files) and ENA downloads (.fastq.gz files):
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

All test scripts have been moved to `tests/test_rna_ena_workflow.py`.
See that file for comprehensive integration tests of the ENA workflow.

## Disk Space Management

The ENA-based workflow automatically manages disk space through batched processing:

**Batched Processing:**
- Downloads N samples from ENA (parallel, robust)
- Quantifies all downloaded samples with kallisto
- Deletes FASTQ files immediately after quantification
- Repeats with next batch
- Peak usage: ~1.5 GB per sample × batch size

**FASTQ Cleanup:**
- Automatic cleanup in `workflow_ena_integrated.py`
- Manual cleanup available via `cleanup_quantified_sra.sh`
- Quantification files retained permanently (~2 MB per sample)

## Output Structure

All outputs go to `output/amalgkit/{species}/`:
- `quant/` - Quantification results (abundance.tsv, ~2 MB per sample)
- `work/metadata/` - Filtered metadata
- `fastq/` - FASTQ files (automatically cleaned after quantification)
- `work/index/` - Kallisto index files
- `logs/` - Processing logs from workflow_ena_integrated.py

See `docs/rna/examples/` for complete documentation.

## Troubleshooting

**ENA Downloads:**
- The workflow uses wget with automatic retry and resume
- If downloads fail, re-run the workflow - it will resume automatically
- Check network connectivity if many downloads fail

**Virtual Environment:**
- `run_multi_species.py` auto-activates venv if available
- If venv missing, script provides setup instructions
- Manual activation: `source .venv/bin/activate`

## Examples

See `docs/rna/examples/pbarbatus_analysis.md` for a complete workflow example.

