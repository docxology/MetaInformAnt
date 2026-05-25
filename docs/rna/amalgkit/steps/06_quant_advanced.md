# amalgkit quant: Advanced Usage & Best Practices

Advanced workflows, cleanup strategies, and real-world examples for the `amalgkit quant` step.

**See**: [06_quant.md](06_quant.md) for core usage and parameters.

---

## Automatic FASTQ Cleanup After Quantification

FASTQ files are deleted immediately after successful quantification. The `abundance.tsv` file is the canonical proof that a sample completed — no separate marker files are used.

> **Note (v0.2.7):** `.safely_removed` marker files were deleted in v0.2.7 and added to `.gitignore`. They are no longer created.

### Automatic Cleanup in Workflow

When using the METAINFORMANT workflow (`execute_workflow()` or `run_workflow.py`), FASTQ files are automatically deleted after quantification:

1. **Download**: Sample FASTQ files are downloaded via `getfastq`
2. **Quantify**: Sample is quantified using `quantify_sample()` from `metainformant.rna.engine.workflow_steps`
3. **Delete**: FASTQ files are automatically deleted using `delete_sample_fastqs()` from `metainformant.rna.engine.sra_extraction`

This per-sample workflow ensures maximum disk efficiency — only one sample's FASTQ files exist at any time.

### Manual Cleanup

For manual processing or recovery of individual samples:

```python-snippet
from metainformant.rna.engine.workflow_steps import quantify_sample
from metainformant.rna.engine.sra_extraction import delete_sample_fastqs
from pathlib import Path

# Quantify sample
success, message, abundance_path = quantify_sample(
    sample_id="SRR14740514",
    metadata_rows=sample_rows,
    quant_params=quant_params,
    log_dir=log_dir,
)

# Delete FASTQ files after successful quantification
if success and abundance_path and abundance_path.exists():
    fastq_dir = Path("output/amalgkit/pogonomyrmex_barbatus/fastq")
    delete_sample_fastqs("SRR14740514", fastq_dir)
```

### Batch Cleanup

The `cleanup_unquantified_samples()` function processes all downloaded but unquantified samples:

```python-snippet
from metainformant.rna.orchestration import cleanup_unquantified_samples
from pathlib import Path

config_path = Path("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
quantified, failed = cleanup_unquantified_samples(config_path)
```

---

## Best Practices

### 1. Always Clean FASTQs After Quantification

```bash
# Good: Save massive disk space
--clean_fastq yes

# Bad: Wastes 100s of GBs
--clean_fastq no
```

**Reasoning**: FASTQ files are 10-50GB per sample. After quantification, you only need the ~1MB abundance.tsv files.

### 2. Verify Reference Transcriptome

```bash
# Before quantifying 100s of samples, test with one
amalgkit quant --batch 1 --threads 8

# Check alignment rate
cat output/work/quant/SRR*/run_info.json | grep "p_pseudoaligned"

# Good: >60%
# Acceptable: 40-60%
# Concerning: <40% (check reference)
```

### 3. Use Appropriate Threading

```bash
# Single sample: use all cores
--threads 16

# Multiple samples in parallel: divide cores
# 4 samples × 4 threads = 16 cores total
```

### 4. Monitor Quantification Progress

```bash
# Check how many samples completed
find output/work/quant -name "abundance.tsv" | wc -l

# Check for failed samples
find output/work/quant -type d | while read dir; do
    if [ ! -f "$dir/abundance.tsv" ]; then
        echo "Failed: $dir"
    fi
done
```

---

## Real-World Examples

### Example 1: Apis mellifera (83 Samples)

```bash
amalgkit quant \
  --out_dir output/amalgkit/apis_mellifera/work \
  --metadata output/amalgkit/apis_mellifera/work/metadata/pivot_qualified.tsv \
  --index_dir output/amalgkit/apis_mellifera/work/index \
  --threads 16 \
  --clean_fastq yes
```

**Result**: 83 samples quantified in ~8 hours (serial), 64.5% average alignment rate, ~350GB disk saved.

### Example 2: Pogonomyrmex barbatus (120 Samples, HPC)

```bash
sbatch --array=1-120 --cpus-per-task=8 --mem=8G quant.sh

# In quant.sh:
amalgkit quant \
  --out_dir output/amalgkit/pogonomyrmex_barbatus/work \
  --batch ${SLURM_ARRAY_TASK_ID} \
  --threads 8 \
  --clean_fastq yes
```

**Result**: All 120 samples completed in 30 minutes on HPC cluster.

### Example 3: METAINFORMANT Workflow Integration

```python
from metainformant.rna.engine.workflow import execute_workflow, load_workflow_config

cfg = load_workflow_config("config/amalgkit/amalgkit_apis_mellifera.yaml")
execute_workflow(cfg)  # quant runs automatically after getfastq/integrate
```

---

## References

- **kallisto**: https://pachterlab.github.io/kallisto/
- **kallisto paper**: https://www.nature.com/articles/nbt.3519

**See Also**: [06_quant.md](06_quant.md) | [06_quant_troubleshooting.md](06_quant_troubleshooting.md) | [Steps Index](README.md)
