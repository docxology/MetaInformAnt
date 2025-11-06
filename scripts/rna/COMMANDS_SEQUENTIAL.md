# Sequential Commands for Pogonomyrmex barbatus Workflow

Complete list of commands to run the workflow from start to finish.

## Quick Start (All Steps at Once)

```bash
# Run all steps sequentially
bash scripts/rna/run_pogonomyrmex_barbatus_workflow.sh
```

## Individual Commands (Run One at a Time)

### Step 1: Check Environment
```bash
python3 scripts/rna/check_environment.py
```

### Step 2: Verify Genome Setup
```bash
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only
```

### Step 3: Check Current Status
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

### Step 4: Metadata (Discover Samples from SRA)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps metadata
```

### Step 5: Integrate (Prepare Metadata)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps integrate
```
*Note: Will skip gracefully if no FASTQ files exist yet (expected before getfastq)*

### Step 6: Config (Generate Config Files)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps config
```

### Step 7: Select (Filter Samples)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps select
```

### Step 8: Getfastq (Download FASTQ Files)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq
```
**⚠️ WARNING: This is the longest step - can take hours or days depending on number of samples**

### Step 9: Integrate (Now with FASTQ Files)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps integrate
```
*Now that FASTQ files exist, this will actually integrate metadata with FASTQ stats*

### Step 10: Quant (Quantify with Kallisto)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps quant
```

### Step 11: Merge (Merge Quantification Results)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps merge
```

### Step 12: Curate (Quality Control)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps curate
```

### Step 13: Sanity (Validate Outputs)
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps sanity
```

### Step 14: Final Status Check
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

## Optional Cleanup Operations (Run as Needed)

### Clean Up Partial/Failed Downloads
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-partial
```

### Quantify Downloaded Samples and Cleanup FASTQs
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --cleanup-unquantified
```

### Fix Abundance File Naming
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --fix-abundance-naming
```

## Alternative: Run All Steps at Once

If you want to run all steps in one command (not recommended for large workflows):

```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

This will run: metadata → integrate → config → select → getfastq → quant → merge → cstmm → curate → csca → sanity

## Notes

- **Resume Support**: The workflow automatically skips completed steps, so you can resume after interruptions
- **Progress Monitoring**: Use `--status` frequently to check progress
- **Disk Space**: Use `--cleanup-unquantified` periodically to free space
- **Logs**: Check `output/amalgkit/pogonomyrmex_barbatus/logs/` for detailed logs
- **Getfastq Step**: This is the longest step - monitor progress and be patient

## Expected Output Locations

- **Metadata**: `output/amalgkit/pogonomyrmex_barbatus/work/metadata/`
- **FASTQ files**: `output/amalgkit/pogonomyrmex_barbatus/fastq/`
- **Quantification**: `output/amalgkit/pogonomyrmex_barbatus/quant/`
- **Merged results**: `output/amalgkit/pogonomyrmex_barbatus/merged/`
- **Curated results**: `output/amalgkit/pogonomyrmex_barbatus/curate/`
- **Logs**: `output/amalgkit/pogonomyrmex_barbatus/logs/`

