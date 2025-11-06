# Complete Workflow Sequence for Pogonomyrmex barbatus

This document provides the sequential commands to run the complete RNA-seq workflow from start to finish.

## Quick Start: Full Workflow

```bash
# 1. Check environment
python3 scripts/rna/check_environment.py

# 2. Verify genome setup
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only

# 3. Run full workflow (all steps)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# 4. Check final status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

## Detailed Sequential Steps

### Step 1: Environment Check
Verify all dependencies are available:
```bash
python3 scripts/rna/check_environment.py
```

**Expected output:** All checks should pass (virtual environment, metainformant, amalgkit, SRA Toolkit, kallisto)

### Step 2: Genome Setup
Check if genome, transcriptome, and kallisto index are ready:
```bash
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only
```

**If genome is not ready, run full setup:**
```bash
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

This will:
- Download genome from NCBI
- Prepare transcriptome FASTA
- Build kallisto index

### Step 3: Check Current Status
See what's already done and what needs to be done:
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

### Step 4: Run Workflow Steps

#### Option A: Run All Steps at Once
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

This runs all steps sequentially:
1. `metadata` - Discover samples from SRA
2. `integrate` - Prepare metadata
3. `config` - Generate config files
4. `select` - Filter samples
5. `getfastq` - Download FASTQ files
6. `quant` - Quantify with kallisto
7. `merge` - Merge quantification results
8. `cstmm` - Cross-species transcriptome mapping (optional)
9. `curate` - Quality control
10. `csca` - Cross-species comparative analysis (optional)
11. `sanity` - Validate outputs

#### Option B: Run Steps Sequentially (Recommended for Large Workflows)

**Step 4a: Metadata**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps metadata
```
Discovers RNA-seq samples from NCBI SRA.

**Step 4b: Integrate**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps integrate
```
Prepares and integrates metadata.

**Step 4c: Config**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps config
```
Generates configuration files for downstream steps.

**Step 4d: Select**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps select
```
Filters samples based on quality criteria.

**Step 4e: Getfastq**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq
```
Downloads FASTQ files from ENA/SRA. This is the longest step (can take hours/days for many samples).

**Step 4f: Quant**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps quant
```
Quantifies expression with kallisto. Requires kallisto index (from genome setup).

**Step 4g: Merge**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps merge
```
Merges quantification results from all samples into a single matrix.

**Step 4h: Curate**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps curate
```
Performs quality control and curation of merged results.

**Step 4i: Sanity**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps sanity
```
Validates workflow outputs and checks for errors.

### Step 5: Cleanup Operations (As Needed)

**Clean up partial/failed downloads:**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-partial
```
Removes samples with partial downloads that aren't quantified.

**Quantify downloaded samples and cleanup FASTQs:**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
```
Finds downloaded samples that haven't been quantified, quantifies them, and deletes FASTQ files to free space.

**Fix abundance file naming:**
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --fix-abundance-naming
```
Creates symlinks from `abundance.tsv` to `{SRR}_abundance.tsv` for merge compatibility.

### Step 6: Final Status Check
Verify everything completed successfully:
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

## Workflow Execution Tips

1. **Monitor Progress**: Use `--status` frequently to check progress
2. **Resume After Interruption**: The workflow automatically skips completed steps
3. **Disk Space**: Use `--cleanup-unquantified` periodically to free space
4. **Parallel Downloads**: The `getfastq` step uses multiple threads (configured in YAML)
5. **Logs**: Check `output/amalgkit/pogonomyrmex_barbatus/logs/` for detailed logs

## Expected Output Locations

- **Metadata**: `output/amalgkit/pogonomyrmex_barbatus/work/metadata/`
- **FASTQ files**: `output/amalgkit/pogonomyrmex_barbatus/fastq/`
- **Quantification**: `output/amalgkit/pogonomyrmex_barbatus/quant/`
- **Merged results**: `output/amalgkit/pogonomyrmex_barbatus/merged/`
- **Curated results**: `output/amalgkit/pogonomyrmex_barbatus/curate/`
- **Logs**: `output/amalgkit/pogonomyrmex_barbatus/logs/`

## Troubleshooting

- **Environment issues**: Run `check_environment.py` first
- **Genome not found**: Run `setup_genome.py` to download and prepare
- **Download failures**: Use `--cleanup-partial` then retry `getfastq` step
- **Disk space**: Use `--cleanup-unquantified` to free space
- **Status unclear**: Use `--status --detailed` for comprehensive information

