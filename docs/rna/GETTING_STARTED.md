# Getting Started with RNA-seq Analysis

Complete guide for setting up and running RNA-seq workflows with METAINFORMANT, from single-species to multi-species production workflows.

## Quick Start

### 1. Installation

```bash
cd /path/to/MetaInformAnt

# Create virtual environment
uv venv
source .venv/bin/activate

# Install metainformant
uv pip install -e .

# Install amalgkit
uv pip install git+https://github.com/kfuku52/amalgkit

# Install system dependencies
sudo apt-get install -y sra-toolkit kallisto fastp seqkit
```

### 2. Verify Installation

```bash
# Check Python packages
python3 -c "import metainformant; print('✓ metainformant')"
amalgkit --version

# Check command-line tools
command -v amalgkit && echo "✓ amalgkit CLI"
command -v fasterq-dump && echo "✓ fasterq-dump"
command -v kallisto && echo "✓ kallisto"

# Set NCBI email (required)
export NCBI_EMAIL="your.email@example.com"
```

### 3. Run Your First Workflow

**Recommended**: Use `run_workflow.py` for complete end-to-end execution. The workflow automatically handles genome setup: if a genome configuration is provided in the YAML file, it will automatically download the genome, prepare the transcriptome, and build the kallisto index before starting sample processing.

```bash
# Full end-to-end workflow (all steps)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Check status at any time
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status
```

**End-to-End Workflow** (executed automatically):
1. **Automatic genome setup** (if genome config exists): Download → Prepare transcriptome → Build kallisto index
2. **Metadata retrieval**: Download sample metadata from NCBI SRA
3. **Per-sample processing**: For each sample: download → quantify → delete FASTQ (as configured with `keep_fastq: no`)
4. **Post-processing**: Integrate, merge, cstmm, curate, csca, sanity

See [scripts/rna/README.md](../../scripts/rna/README.md) for complete usage documentation.

## System Requirements

- **Python 3.11+** with virtual environment (`.venv/` or `/tmp/metainformant_venv` on ext6 filesystems)
- **uv** (Python package manager) - **REQUIRED** for all setup and package management
- **amalgkit** (installed in venv via `uv pip install`)
- **wget** (for ENA downloads) or **SRA Toolkit** (for SRA downloads)
- **kallisto** or **salmon** (for quantification)
- **fastp** and **seqkit** (for quality control)

**Note**: On external drives or filesystems with limitations (e.g., ext6 without symlink support), the system automatically handles setup. See [External Drive Setup Guide](EXTERNAL_DRIVE_SETUP.md) for details.

## Detailed Setup

### Python Environment Setup

The project uses `uv` for fast, reliable Python package management:

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
cd /path/to/MetaInformAnt
uv venv

# Activate virtual environment
source .venv/bin/activate

# Install metainformant in editable mode
uv pip install -e .

# Install amalgkit
uv pip install git+https://github.com/kfuku52/amalgkit
```

**Verification:**
```bash
python3 -c "import metainformant; print('metainformant OK')"
amalgkit --version  # Should show: amalgkit version 0.12.20 or higher
```

### System Dependencies

#### SRA Toolkit (for SRA downloads)
```bash
# Debian/Ubuntu
sudo apt-get update
sudo apt-get install -y sra-toolkit

# Verify installation
fasterq-dump --version
prefetch --version
```

#### RNA-seq Quantification Tools

**kallisto (Recommended):**
```bash
sudo apt-get install -y kallisto
kallisto version  # Should show: kallisto, version 0.48.0 or higher
```

**salmon (Alternative):**
```bash
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar -xzf salmon-1.10.0_linux_x86_64.tar.gz
export PATH=$PATH:$(pwd)/salmon-latest_linux_x86_64/bin
salmon --version
```

#### Quality Control Tools
```bash
sudo apt-get install -y fastp seqkit

# Verify installations
fastp --version  # Should show 0.23.2 or higher
seqkit version   # Should show 2.3.0 or higher
```

### Environment Configuration

#### Required Environment Variables

```bash
# Add to ~/.bashrc or ~/.zshrc
export NCBI_EMAIL="your.email@example.com"  # Required for NCBI API

# Optional: Configure data directories
export METAINFORMANT_OUTPUT="/path/to/MetaInformAnt/output"
export METAINFORMANT_DATA="/path/to/MetaInformAnt/data"
```

## Running Workflows

### Single Species End-to-End Workflow (Recommended)

**Use `run_workflow.py` for complete end-to-end execution**:
```bash
# Full workflow (all 11 steps)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge
```

**Features:**
- Complete end-to-end execution via `execute_workflow()`
- Automatic genome setup (if genome config exists)
- Per-sample processing: download → quantify → delete FASTQ
- All 11 amalgkit steps in correct order
- Status checking and cleanup operations
- Resume support (automatically skips completed steps)

### Multi-Species Workflows

**For multiple species**, run `run_workflow.py` separately for each species config:

```bash
# Process each species separately with parallel downloads
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species2.yaml
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species3.yaml
```

**Or run in parallel (background):**
```bash
nohup python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species1.yaml > logs/species1.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species2.yaml > logs/species2.log 2>&1 &
```

**Parallel downloads** are controlled via `num_download_workers` in each config file:

```yaml
steps:
  getfastq:
    num_download_workers: 16  # Number of parallel download processes
    threads: 24
```

**Features:**
- ✅ **Immediate processing**: Each sample: download → immediately quantify → immediately delete FASTQs
- ✅ **Parallel downloads**: Configurable via `num_download_workers` in config file
- ✅ **Maximum disk efficiency**: Only one sample's FASTQs exist at a time
- ✅ **100% reliability**: ENA-based downloads (vs ~0% for SRA Toolkit)
- ✅ **Virtual environment**: Auto-activates if available
- ✅ **Cloud acceleration**: AWS, GCP, NCBI enabled by default
- ✅ **Progress tracking**: Status checking via `--status` flag

**Configuration Options:**
- `num_download_workers`: Number of parallel download processes per species (default: 16)
- `threads`: Threads for quantification (default: 24)

**Note**: For multiple species, run `run_workflow.py` in parallel (using `nohup` or `screen`) or sequentially. Each workflow processes one species at a time with configurable parallelism.

## Monitoring Progress

### Status Checking

```bash
# Check workflow status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Detailed status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

### Check Running Processes

```bash
# Count workflow processes
ps aux | grep "run_workflow\|amalgkit" | grep -v grep

# Count active downloads
ps aux | grep wget | grep -v grep | wc -l

# Show all RNA workflow processes
ps aux | grep -E "(run_workflow|amalgkit|wget)" | grep -v grep
```

### Check Individual Species Logs

```bash
# Find latest logs
ls -lt output/workflow_*.log | head -4

# Tail specific species
tail -f output/workflow_cfloridanus_*.log
tail -f output/workflow_pbarbatus_*.log

# Check for errors
grep -i error output/workflow_*.log
grep -i failed output/workflow_*.log
```

### Check Disk Usage

```bash
# Overall amalgkit directory
du -sh output/amalgkit/

# Per species
du -sh output/amalgkit/*/

# FASTQ directories (should stay small with immediate processing)
du -sh output/amalgkit/*/fastq/

# Quantification results
du -sh output/amalgkit/*/quant/
```

## Verification & Status Checking

### Quick Verification

Check workflow outputs manually:

```bash
cd output/amalgkit/<species_name>

# Verify quantification results
ls -lh quant/*/abundance.tsv

# Verify merged expression matrix
ls -lh merged/*.tsv

# Check sanity check results
cat sanity/*.txt

# View curate QC reports
ls -lh curate/*.tsv
```

### Sanity Check

**Purpose**: Validate integrity of all workflow outputs

```bash
cd output/amalgkit/<species_name>
amalgkit sanity --out_dir work --all
```

**Expected Output**:
```
Quant outputs found for all SRA IDs in --metadata (inferred)
amalgkit sanity: end
```

**Exit Code**: 0 (success)

**Files Created**:
- `work/sanity/SRA_IDs_without_fastq.txt` - Lists samples without FASTQ files (expected if FASTQs deleted after quant)
- NO `SRA_IDs_without_quant.txt` file = all samples validated ✅

**Verification**:
```bash
# Should NOT exist (or be empty) if all samples are good
ls work/sanity/SRA_IDs_without_quant.txt 2>/dev/null || echo "✅ All samples validated!"
```

### Curate Step

**Purpose**: Quality control, outlier removal, and visualization generation

```bash
cd output/amalgkit/<species_name>
amalgkit curate --out_dir work --batch_effect_alg no
```

**Why `--batch_effect_alg no`?**  
- The `sva` R package has system-level compilation issues
- Using `no` skips SVA batch correction but still generates all other outputs
- This is perfectly acceptable for most analyses

**Expected Output**:
```
Removing samples with mapping rate of 0.2
Mapping rate cutoff: 20%
No entry removed due to low mapping rate.
Iteratively checking within-sample_group correlation
No batch effect correction was performed.
Round: 2 : # before = N : # after = N
Writing summary files for <species_name>
transcriptome_curation.r: Completed.
```

**Files Created** (17 total):

#### Visualizations (6 PDFs)
```
work/curate/<species_name>/plots/
├── <species_name>.0.original.pdf
├── <species_name>.0.original.no.pdf
├── <species_name>.1.mapping_cutoff.pdf
├── <species_name>.1.mapping_cutoff.no.pdf
├── <species_name>.2.correlation_cutoff.pdf
└── <species_name>.2.correlation_cutoff.no.pdf
```

Each PDF contains:
- Hierarchical clustering dendrogram
- Correlation heatmap
- PCA plot
- QC metrics

#### Data Tables (7 TSV files)
```
work/curate/<species_name>/tables/
├── <species_name>.uncorrected.tc.tsv                    # Original expression matrix
├── <species_name>.uncorrected.sample_group.mean.tsv     # Sample group means
├── <species_name>.no.tc.tsv                             # Final curated expression (⭐ PRIMARY OUTPUT)
├── <species_name>.no.sample_group.mean.tsv              # Final sample means
├── <species_name>.no.tau.tsv                            # Tissue specificity scores
├── <species_name>.no.correlation_statistics.tsv         # QC correlation statistics
└── <species_name>.metadata.tsv                          # Sample metadata with QC metrics
```

## Understanding the Output

### Directory Structure

```
output/amalgkit/
├── {species}/
│   ├── quant/              # Quantification results (~2 MB per sample)
│   │   ├── SRR1234567/
│   │   │   ├── abundance.tsv
│   │   │   └── run_info.json
│   │   └── ...
│   ├── fastq/              # Temporary FASTQs (auto-deleted)
│   ├── work/
│   │   ├── metadata/
│   │   └── index/          # Kallisto index
│   └── logs/               # Workflow logs
```

### Key Files

- **Per-sample quantification**: `output/amalgkit/{species}/quant/{SRR_ID}/abundance.tsv`
- **Metadata**: `output/amalgkit/{species}/work/metadata/metadata.tsv`
- **Kallisto index**: `output/amalgkit/{species}/work/index/transcriptome.idx`
- **Logs**: `output/amalgkit/{species}/logs/`

### Progress Indicators

**Completed sample:**
- Quantification directory exists: `output/amalgkit/{species}/quant/{SRR_ID}/`
- Contains: `abundance.tsv` and `run_info.json`
- FASTQ files deleted from `fastq/{SRR_ID}/`

**Downloading sample:**
- FASTQ directory exists: `output/amalgkit/{species}/fastq/{SRR_ID}/`
- Contains `.fastq.gz` files (growing)
- wget process visible in `ps aux | grep wget`

**Failed sample:**
- Check logs: `grep {SRR_ID} output/workflow_{species}_*.log`
- Workflow will retry automatically on restart

## Performance Characteristics

### Expected Timeline

**Per sample:**
- Download: 2-8 minutes (varies by sample size, 1-4 GB)
- Quantification: 30-90 seconds (with kallisto, 12 threads)
- Cleanup: <1 second
- **Total**: ~7.5 minutes average

**Full workflow:**
- Small species (50-100 samples): ~6-12 hours
- Medium species (100-300 samples): ~12-24 hours
- Large species (300+ samples): ~24-48 hours

### Resource Usage

**CPU:**
- Parallel downloads: Configurable via `num_download_workers` (default: 16)
- Quantification: Configurable via `threads` (default: 24)
- Recommend: 64+ cores for smooth parallel processing of multiple species

**RAM:**
- ~4 GB per workflow
- ~16 GB total for 4 parallel workflows

**Disk:**
- Peak: ~2-10 GB per sample (only one sample's FASTQs at a time with immediate processing)
- Total peak: ~50-100 GB (temporary, auto-cleaned per sample)
- Final: ~2-3 GB per species (quantification files only)

**Network:**
- Bandwidth: ~40 GB downloading at any time
- ENA servers: highly reliable, automatic resume

## Troubleshooting

### "amalgkit: command not found"
```bash
# Ensure virtual environment is activated
source .venv/bin/activate

# Reinstall if needed with uv
uv pip install --force-reinstall git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3
```

### "fasterq-dump: command not found"
```bash
# Install SRA Toolkit
sudo apt-get install -y sra-toolkit

# Or add to PATH if manually installed
export PATH=$PATH:/path/to/sratoolkit/bin
```

### "Permission denied" errors
```bash
# Ensure output directory is writable
chmod -R u+w output/

# Run script from repository root
cd /path/to/MetaInformAnt
python3 scripts/rna/run_workflow.py ...
```

### External Drive and Filesystem Issues

**See [External Drive Setup Guide](EXTERNAL_DRIVE_SETUP.md) for comprehensive documentation.**

Common issues on external drives (ext6 filesystems):
- Symlink errors when installing packages
- Virtual environment creation failures
- UV cache directory issues

**All issues are automatically handled** by the setup utilities, but see the guide for manual troubleshooting.

### PEP 668 "externally-managed-environment" error
This is expected on modern Linux systems. Always use virtual environments with uv:
```bash
# Create venv with uv
uv venv

# Activate venv
source .venv/bin/activate

# Install packages with uv (no activation needed if using --python flag)
uv pip install -e . --python .venv/bin/python3
```

### Downloads Failing
- Check network connectivity: `ping -c 4 8.8.8.8`
- Verify NCBI_EMAIL is set: `echo $NCBI_EMAIL`
- Try reducing threads: `--threads 8`
- Check disk space: `df -h /`
- Test ENA connectivity: `wget --spider https://www.ebi.ac.uk/ena/`

### Disk Space Issues
- Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
- Reduce `num_download_workers` in config file for fewer parallel downloads
- Use immediate per-sample processing (default): only one sample's FASTQs exist at a time

### Workflow Not Starting

**Check virtual environment:**
```bash
# Verify venv exists (check both locations)
ls -la .venv/ || ls -la /tmp/metainformant_venv/

# Verify amalgkit installed
source .venv/bin/activate || source /tmp/metainformant_venv/bin/activate
amalgkit --version
```

**Check dependencies:**
```bash
# Verify wget
which wget
wget --version

# Verify kallisto
which kallisto
kallisto version
```

### Processes Stuck or Hung

**Kill and restart:**
```bash
# Find process IDs
ps aux | grep "run_workflow\|amalgkit" | grep -v grep

# Kill specific workflow
kill <PID>

# Or kill all workflows (use with caution)
pkill -f run_workflow

# Restart (will resume from last completed sample)
nohup python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml > logs/species.log 2>&1 &
```

### Quantification Failing

**Check kallisto index:**
```bash
# Verify index exists
ls -lh output/amalgkit/*/work/index/*.idx

# Rebuild if needed (automatic on first run)
# Just restart the workflow
```

**Check FASTQ quality:**
```bash
# List downloaded FASTQs
ls -lh output/amalgkit/*/fastq/*/

# Check if files are complete (not 0 bytes)
find output/amalgkit/*/fastq -name "*.fastq.gz" -size 0
```

### Curate completes in 0 seconds

**Cause**: Existing outputs detected, no reprocessing needed

**Solution**: This is normal! To regenerate:
```bash
rm -rf work/curate/*
amalgkit curate --out_dir work --batch_effect_alg no
```

### R package errors (Rtsne, amap, vegan, sva, RUVSeq)

**Cause**: System-level gcc compilation issues (documented)

**Solution**: Already handled via patches to `amalgkit/curate.r`:
- Changed `library()` to `require()` for graceful degradation
- Added fallbacks for missing functions
- Core visualizations still work perfectly ✅

**Impact**: Minimal - only advanced optional features affected

## Next Steps

1. **Learn about workflows**: See [workflow.md](workflow.md) for planning and execution
2. **Understand steps**: See [Step Documentation](amalgkit/steps/README.md) for individual workflow steps
3. **Configure workflows**: See [CONFIGURATION.md](CONFIGURATION.md) for configuration options
4. **Discover species**: See [DISCOVERY.md](DISCOVERY.md) for species discovery
5. **View examples**: See [EXAMPLES.md](EXAMPLES.md) for real-world analysis examples

## Getting Help

- **Documentation**: See [README.md](README.md) for navigation
- **Workflow Issues**: See [workflow.md](workflow.md#common-issues-and-solutions)
- **Orchestration**: See [ORCHESTRATION.md](ORCHESTRATION.md)

## Related Documentation

- **[workflow.md](workflow.md)**: Workflow planning and execution
- **[Step Documentation](amalgkit/steps/README.md)**: Individual workflow steps
- **[CONFIGURATION.md](CONFIGURATION.md)**: Configuration management
- **[ORCHESTRATION.md](ORCHESTRATION.md)**: Orchestrator overview
- **[DISCOVERY.md](DISCOVERY.md)**: Species discovery system
- **[EXAMPLES.md](EXAMPLES.md)**: Real-world analysis examples

## See Also

- **Source Code**: `src/metainformant/rna/` - Implementation details
- **Module Documentation**: `src/metainformant/rna/README.md` - API reference
- **Tests**: `tests/test_rna_*.py` - Comprehensive test coverage
