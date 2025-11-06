# Amalgkit Pipeline Commands Reference

**Species-Specific Example** - Complete command reference for running the amalgkit RNA-seq workflow for Pogonomyrmex barbatus.

> **For general step documentation**: See [amalgkit/steps/](amalgkit/steps/)  
> **For workflow guide**: See [workflow.md](workflow.md)  
> **For API reference**: See [API.md](API.md)

This document provides a complete walkthrough of all 11 amalgkit steps using *Pogonomyrmex barbatus* as an example. For general step documentation applicable to all species, see the [step documentation](amalgkit/steps/README.md).

**Configuration File**: `config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml`

**Base Paths**:
- Work Directory: `output/amalgkit/pogonomyrmex_barbatus/work`
- Log Directory: `output/amalgkit/pogonomyrmex_barbatus/logs`
- FastQ Directory: `output/amalgkit/pogonomyrmex_barbatus/fastq`
- Quant Directory: `output/amalgkit/pogonomyrmex_barbatus/quant`

---

## Step 0: Environment Setup and Verification

### 0.1: Check Environment
```bash
python3 scripts/rna/check_environment.py
```

### 0.2: Verify Genome Setup
```bash
python3 scripts/rna/setup_genome.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --verify-only
```

### 0.3: Check Current Status
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --status \
  --detailed
```

---

## Step 1: Metadata (Discover Samples from SRA)

**Purpose**: Query NCBI SRA for RNA-seq samples matching the search criteria.

**Command**:
```bash
amalgkit metadata \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work \
  --search_string '"Pogonomyrmex barbatus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]' \
  --entrez_email "${NCBI_EMAIL:-DanielAriFriedman@gmail.com}" \
  --redo yes
```

**Output**: 
- `output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv`

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps metadata
```

---

## Step 2: Config (Generate Configuration Files)

**Purpose**: Generate tool configuration files (kallisto, salmon, etc.) for downstream steps.

**Command**:
```bash
amalgkit config \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work \
  --config base
```

**Output**:
- `output/amalgkit/pogonomyrmex_barbatus/work/config_base/` (directory with config files)

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps config
```

---

## Step 3: Select (Filter Samples)

**Purpose**: Filter samples based on quality criteria (read count, base count, read length) and tissue mappings.

**Command**:
```bash
amalgkit select \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --config_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/config_base
```

**Note**: The `config_dir` path is automatically resolved to an absolute path by the workflow.

**Output**:
- `output/amalgkit/pogonomyrmex_barbatus/work/metadata/pivot_qualified.tsv`
- `output/amalgkit/pogonomyrmex_barbatus/work/metadata/pivot_selected.tsv`

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps select
```

---

## Step 4: Getfastq (Download FASTQ Files)

**Purpose**: Download SRA files and convert them to FASTQ format. **THIS IS THE LONGEST STEP** (can take hours or days).

**Command**:
```bash
amalgkit getfastq \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/fastq \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/pivot_qualified.tsv \
  --threads 24 \
  --aws yes \
  --gcp yes \
  --ncbi yes \
  --pfd no \
  --fastp yes
```

**Note**: 
- `pfd: no` disables parallel-fastq-dump (prevents dependency check error)
- `fastp: yes` enables quality filtering
- Cloud sources (AWS, GCP, NCBI) are all enabled for faster downloads

**Output**:
- `output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq/SRR*/` (one directory per sample)
  - `SRR*_1.fastq.gz` (paired-end read 1)
  - `SRR*_2.fastq.gz` (paired-end read 2)
  - `fastp.json` and `fastp.html` (quality reports)

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps getfastq
```

**WARNING**: This step can take hours or days depending on the number of samples and network speed.

---

## Step 5: Integrate (Integrate FASTQ Files with Metadata)

**Purpose**: Add FASTQ file paths to metadata after download.

**Command**:
```bash
amalgkit integrate \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --fastq_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/fastq
```

**Output**:
- Updated `metadata.tsv` with FASTQ file paths

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps integrate
```

**Note**: This step can be run before getfastq (will skip if no FASTQ files), but should be run again after getfastq to integrate downloaded files.

---

## Step 6: Quant (Quantify with Kallisto)

**Purpose**: Quantify transcript abundances using kallisto.

**Command**:
```bash
amalgkit quant \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/quant \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --threads 24 \
  --redo no \
  --keep_fastq no \
  --build_index yes
```

**Note**:
- `build_index: yes` builds kallisto index if not present
- `keep_fastq: no` deletes FASTQ files after quantification to save space
- `redo: no` skips already quantified samples

**Output**:
- `output/amalgkit/pogonomyrmex_barbatus/quant/SRR*/abundance.tsv` (per-sample quantification)

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps quant
```

---

## Step 7: Merge (Merge Quantification Results)

**Purpose**: Combine per-sample quantification results into a single matrix.

**Command**:
```bash
amalgkit merge \
  --out /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/merged/merged_abundance.tsv \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/merged \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv
```

**Output**:
- `output/amalgkit/pogonomyrmex_barbatus/merged/merged_abundance.tsv` (gene × sample matrix)

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps merge
```

---

## Step 8: CSTMM (Cross-Species TMM Normalization)

**Purpose**: Perform cross-species TMM (Trimmed Mean of M-values) normalization.

**Command**:
```bash
amalgkit cstmm \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/cstmm \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --threads 12
```

**Output**:
- Normalized expression matrices in `output/amalgkit/pogonomyrmex_barbatus/cstmm/`

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps cstmm
```

---

## Step 9: Curate (Quality Control)

**Purpose**: Detect outliers and unwanted biases in the dataset.

**Command**:
```bash
amalgkit curate \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/curate \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --threads 12
```

**Output**:
- Quality control reports and outlier detection results in `output/amalgkit/pogonomyrmex_barbatus/curate/`

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps curate
```

---

## Step 10: CSCA (Cross-Species Correlation Analysis)

**Purpose**: Perform cross-species correlation analysis and generate plots.

**Command**:
```bash
amalgkit csca \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/csca \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --threads 12
```

**Output**:
- Correlation analysis results and plots in `output/amalgkit/pogonomyrmex_barbatus/csca/`

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps csca
```

---

## Step 11: Sanity (Final Validation)

**Purpose**: Perform final integrity checks on all outputs.

**Command**:
```bash
amalgkit sanity \
  --out_dir /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work \
  --metadata /media/q/ext6/github/MetaInformAnt/output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv
```

**Output**:
- Validation report

**Via Workflow Script**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps sanity
```

---

## Step 12: Final Status Check

**Purpose**: Verify all steps completed successfully.

**Command**:
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --status \
  --detailed
```

---

## Running the Complete Pipeline

### Option 1: Run All Steps Sequentially
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

### Option 2: Run Specific Steps
```bash
python3 scripts/rna/run_workflow.py \
  --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml \
  --steps metadata,config,select,getfastq,integrate,quant,merge,cstmm,curate,csca,sanity
```

### Option 3: Run Steps Individually (Recommended for Debugging)
```bash
# Step 1: Metadata
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps metadata

# Step 2: Config
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps config

# Step 3: Select
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps select

# Step 4: Getfastq (LONGEST STEP - can take hours/days)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq

# Step 5: Integrate
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps integrate

# Step 6: Quant
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps quant

# Step 7: Merge
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps merge

# Step 8: CSTMM
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps cstmm

# Step 9: Curate
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps curate

# Step 10: CSCA
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps csca

# Step 11: Sanity
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps sanity
```

---

## Important Notes

1. **Path Resolution**: All paths in the commands are automatically resolved to absolute paths by the workflow system. The commands shown use absolute paths for clarity.

2. **Environment Variables**: 
   - `NCBI_EMAIL` is used for Entrez API queries (set via environment or config)
   - Default: `DanielAriFriedman@gmail.com`

3. **Thread Configuration**:
   - Default: 12 threads (from config)
   - Getfastq and Quant steps: 24 threads (overridden in config)

4. **Dependencies**:
   - `parallel-fastq-dump` is disabled (`pfd: no`) to avoid dependency check errors
   - `fastp` is enabled for quality filtering
   - Cloud sources (AWS, GCP, NCBI) are enabled for faster downloads

5. **Workflow Order**:
   - Metadata → Config → Select → Getfastq → Integrate → Quant → Merge → CSTMM → Curate → CSCA → Sanity
   - Integrate can be run before getfastq (will skip if no FASTQ files), but should be run again after getfastq

6. **Log Files**: All steps write logs to `output/amalgkit/pogonomyrmex_barbatus/logs/` with timestamps.

---

## Troubleshooting

### Select Step Fails with "No such file or directory: './config'"
- **Fixed**: The workflow now automatically resolves `config_dir` to an absolute path
- Ensure the `config` step completed successfully before running `select`

### Getfastq Step Fails with "No such file or directory: 'parallel-fastq-dump'"
- **Fixed**: When `pfd: no` is set, `pfd_exe` is removed from parameters to prevent dependency checks
- Verify `pfd: no` is set in your config file

### Path Issues
- All paths are automatically resolved relative to the repository root
- If you see path errors, check that you're running from the repository root directory

## See Also

### Step Documentation
- **[amalgkit/steps/README.md](amalgkit/steps/README.md)** - All step documentation
- **[amalgkit/steps/metadata.md](amalgkit/steps/metadata.md)** - Metadata step details
- **[amalgkit/steps/quant.md](amalgkit/steps/quant.md)** - Quantification step details

### Documentation
- **[API Reference](API.md)** - Complete function documentation
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Configuration Guide](CONFIGURATION.md)** - Configuration management
- **[Main Index](README.md)** - RNA domain master index

