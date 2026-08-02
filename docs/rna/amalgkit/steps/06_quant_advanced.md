# amalgkit quant: Advanced Usage & Best Practices

Advanced workflows, cleanup strategies, and real-world examples for the `amalgkit quant` step.

**See**: [06_quant.md](06_quant.md) for core usage and parameters.

---

## Automatic FASTQ Cleanup After Quantification

The project retains FASTQ files during quantification. After a non-empty
abundance table and exact current-method provenance sidecar are written,
`scripts/rna/reclaim_quantified_raw.py` removes the sample's raw inputs and
writes an audit manifest.

Amalgkit may leave `.safely_removed` markers from older runs; they are
redundant and are removed by the same guarded utility after current evidence
is present.

### Automatic Cleanup in Workflow

When using the current streaming workflow, FASTQ files are reclaimed only
after current quantification evidence is written:

1. **Download**: Sample FASTQ files are downloaded via `getfastq`
2. **Quantify**: Sample is quantified using `quantify_sample()` from `metainformant.rna.engine.workflow_steps`
3. **Reclaim**: `reclaim_quantified_raw.py` removes only validated raw inputs
   for samples with current provenance

The per-sample workflow bounds concurrent FASTQ trees by the selected sample
worker profile. When workers exceed one, disk planning must use the selected
worker count and sample-size bound because multiple raw-read trees may coexist.

### Manual Cleanup

For manual processing or recovery of individual samples, run the same
provenance-aware lane and then use its audited reclamation command:

```bash
uv run python scripts/rna/reclaim_quantified_raw.py \
  --data-root "$AMALGKIT_DATA_ROOT" --species apis_mellifera --execute
```

### Batch Cleanup

The `cleanup_unquantified_samples()` function processes all downloaded but unquantified samples:

```python-snippet
from metainformant.rna.engine.workflow_cleanup import cleanup_unquantified_samples
from pathlib import Path

config_path = Path("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
quantified, failed = cleanup_unquantified_samples(config_path)
```

---

## Best Practices

### 1. Use the guarded reclamation contract

```bash
# Keep inputs until the current sidecar exists
--clean_fastq no

# Then run an audited reclamation pass
uv run python scripts/rna/reclaim_quantified_raw.py \
  --data-root "$AMALGKIT_DATA_ROOT" --execute
```

**Reasoning**: FASTQ files are large, but deleting them before current
provenance is written can make a failed quantification contract ambiguous.

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
  --clean_fastq no
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
  --clean_fastq no
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
