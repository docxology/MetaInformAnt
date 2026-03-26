# METAINFORMANT Troubleshooting Guide

Advanced troubleshooting guide for complex issues in METAINFORMANT workflows.

## Memory Issues with Large Datasets

### Problem: Out of Memory (OOM) errors with large VCF/FASTQ files

**Symptoms**: `MemoryError`, `OSError: cannot allocate memory`, process killed

**Solutions**:

1. **Process in chunks**:
```python
from metainformant.core.io import read_jsonl

# Read large JSONL
for record in read_jsonl("large_file.jsonl"):
    process_record(record)
    # Optionally track memory and clear periodically
```

2. **Use streaming for VCF**:
```python
from metainformant.gwas.quality import parse_vcf_full

# Parse VCF - processes in memory
# For very large files, consider using external tools:
# bcftools view -H large.vcf.gz | chunk processing
vcf_data = parse_vcf_full("large.vcf.gz")
```

3. **Reduce working set**:
```python
# Filter to essential fields only
vcf_data = parse_vcf_full(vcf_path, fields=["CHROM", "POS", "ID", "GT"])
```

### Problem: GPU memory exhausted during deep learning

**Solutions**:

```python
# Enable gradient checkpointing
model.enable_gradient_checkpointing()

# Use mixed precision
from torch.cuda.amp import autocast
with autocast():
    outputs = model(inputs)

# Reduce batch size
trainer = Trainer(batch_size=16)  # Reduce from 32
```

## Parallel Processing Failures

### Problem: Thread/Process Pool Crashes

**Symptoms**: `RuntimeError: pool not started`, deadlocks, zombie processes

**Solutions**:

1. **Check worker count**:
```python
from metainformant.core.parallel import ParallelProcessor

# Don't exceed CPU count
import os
max_workers = min(8, os.cpu_count() or 4)
processor = ParallelProcessor(max_workers=max_workers)
```

2. **Handle exceptions in workers**:
```python
def safe_process(item):
    try:
        return process_item(item)
    except Exception as e:
        logger.error(f"Failed to process {item}: {e}")
        return None

results = parallel.run_parallel(safe_process, items, max_workers=4)
results = [r for r in results if r is not None]
```

3. **Use spawn start method (Linux)**:
```python
import multiprocessing as mp
if __name__ == '__main__':
    mp.set_start_method('spawn', force=True)
```

### Problem: Deadlock in SQLite Progress Database

**Symptoms**: Pipeline hangs, database locked errors

**Solutions**:

The RNA pipeline uses SQLite for progress tracking. High-throughput runs can cause deadlocks:

1. **Use read-only connections for status queries**:
```python
# Instead of Python (can cause WAL deadlock), use shell sqlite3 binary
import subprocess
result = subprocess.run(
    ["sqlite3", "pipeline_progress.db", 
     "SELECT COUNT(*) FROM samples WHERE status='quantified'"],
    capture_output=True, text=True
)
```

2. **Use the monitoring script**:
```bash
python3 scripts/rna/check_pipeline_status.py
```

3. **Copy DB to /tmp for intensive queries**:
```bash
cp pipeline_progress.db /tmp/status_db.sqlite
sqlite3 /tmp/status_db.sqlite "SELECT * FROM samples LIMIT 10;"
```

## External Tool Integration Issues

### Problem: amalgkit CLI not found

**Symptoms**: `FileNotFoundError: amalgkit not found`, `CommandNotFoundError`

**Solutions**:

1. **Check installation**:
```bash
which amalgkit
amalgkit --version
```

2. **Add to PATH**:
```bash
# Add to ~/.bashrc
export PATH="$HOME/.local/bin:$PATH"

# Or use full path in config
amalgkit:
  cli_path: "/home/user/.local/bin/amalgkit"
```

3. **Install R dependencies**:
```bash
# Check R packages
Rscript -e "library(tidyr); library(dplyr)"

# Install missing
Rscript -e "install.packages(c('tidyr', 'dplyr'))"
```

### Problem: bcftools/vcftools errors

**Symptoms**: `bcftools: error while loading shared libraries`, malformed VCF

**Solutions**:

1. **Verify tool installation**:
```bash
bcftools --version
vcftools --version
```

2. **Check for compatible versions**:
```bash
# bcftools >= 1.10 required for some features
bcftools --version | head -1
```

3. **Validate input files**:
```python
from metainformant.gwas.quality import validate_vcf

errors = validate_vcf("input.vcf.gz")
if errors:
    print(f"VCF validation errors: {errors}")
```

### Problem: Kallisto index errors

**Symptoms**: `Error: index file not found`, `Error: invalid index`

**Solutions**:

1. **Generate index**:
```bash
kallisto index -i transcripts.idx transcripts.fasta
```

2. **Check index compatibility**:
```bash
kallisto version
# Ensure index was created with same version
```

3. **Use correct index path in config**:
```yaml
quant:
  index: "path/to/transcripts.idx"
```

## Network and Download Issues

### Problem: ENA/NCBI download failures

**Symptoms**: `ConnectionError`, `HTTP 404`, timeout during download

**Solutions**:

1. **Use SRA fallback**:
```python
# Configure SRA fallback in workflow
config = {
    "getfastq": {
        "prefer_ena": True,
        "fallback_to_sra": True,
        "max_retries": 3
    }
}
```

2. **Check network connectivity**:
```bash
curl -s https://www.ebi.ac.uk/ena/portal/api/ | head
curl -s https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ | head
```

3. **Use fasterq-dump explicitly**:
```python
from metainformant.rna.sra_extraction import download_sra_run

# Direct SRA download
download_sra_run(
    accession="SRR14740514",
    output_dir="output/fastq",
    tool="fasterq-dump"
)
```

### Problem: NCBI rate limiting

**Symptoms**: `429 Too Many Requests`, `Rate limit exceeded`

**Solutions**:

1. **Configure rate limiting**:
```yaml
ncbi:
  rate_limit_delay: 1.0  # seconds between requests
  max_retries: 5
```

2. **Set email for polite pool**:
```bash
export NCBI_EMAIL="your.email@example.com"
```

3. **Use caching**:
```python
from metainformant.core.cache import JsonCache

cache = JsonCache("output/ncbi_cache", ttl_seconds=86400)  # 24h
```

## Data Format Issues

### Problem: Parsing errors in VCF/BCF files

**Symptoms**: `ValueError: could not convert string to float`, malformed records

**Solutions**:

1. **Validate VCF format**:
```python
from metainformant.gwas.quality import validate_vcf

errors = validate_vcf(vcf_path)
```

2. **Filter problematic variants**:
```python
vcf_data = parse_vcf_full(vcf_path)
vcf_data = filter_variants_by_quality(vcf_data, min_qual=30)
```

3. **Use bgzip + tabix index**:
```bash
bgzip -c input.vcf > input.vcf.gz
tabix -p vcf input.vcf.gz
```

### Problem: FASTQ quality encoding issues

**Symptoms**: `Quality score out of range`, incorrect GC content

**Solutions**:

1. **Detect quality encoding**:
```python
from metainformant.quality.fastq import detect_quality_encoding

encoding = detect_quality_encoding(fastq_path)
print(f"Detected quality encoding: {encoding}")
```

2. **Convert encoding**:
```bash
# Convert Sanger to Illumina
fastq_quality_convert -i input.fastq -o output.fastq -s sanger-to-illumina
```

## Environment and Configuration Issues

### Problem: Virtual environment issues

**Symptoms**: `ModuleNotFoundError`, incorrect package versions

**Solutions**:

1. **Recreate virtual environment**:
```bash
rm -rf .venv
uv venv
source .venv/bin/activate
uv pip install -e .
```

2. **Check Python version**:
```bash
python --version  # Should be 3.11+
```

3. **Verify package installation**:
```bash
uv pip list | grep metainformant
```

### Problem: Configuration file errors

**Symptoms**: `YAML error`, missing required fields

**Solutions**:

1. **Validate YAML syntax**:
```bash
python -c "import yaml; yaml.safe_load(open('config.yaml'))"
```

2. **Use schema validation**:
```python
from metainformant.core.validation import validate_config_schema

errors = validate_config_schema(config, "gwas_schema")
```

3. **Check required fields**:
```python
required = ["work_dir", "threads", "genome"]
for field in required:
    if field not in config:
        print(f"Missing required field: {field}")
```

## Performance Optimization

### Problem: Slow analysis on large datasets

**Symptoms**: Analysis taking hours/days longer than expected

**Solutions**:

1. **Enable parallel processing**:
```python
results = parallel.run_parallel(
    analyze_function, 
    data_items, 
    max_workers=8  # Adjust based on CPU
)
```

2. **Use appropriate data structures**:
```python
import numpy as np

# Fast vectorized operations
data = np.array(data_list)
result = np.mean(data)  # Instead of sum(data)/len(data)
```

3. **Profile to find bottlenecks**:
```python
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()
# Your analysis code
profiler.disable()

stats = pstats.Stats(profiler)
stats.sort_stats('cumulative').print_stats(20)
```

### Problem: Disk space issues

**Symptoms**: `No space left on device`, cannot write files

**Solutions**:

1. **Use cleanup options**:
```bash
# RNA workflow cleanup
python3 scripts/rna/run_workflow.py config.yaml --cleanup-unquantified
```

2. **Set custom temp directory**:
```bash
export TMPDIR="/path/to/large/disk/tmp"
```

3. **Monitor disk usage**:
```python
from metainformant.core.paths import get_directory_size

size = get_directory_size("output/")
print(f"Directory size: {size / 1e9:.2f} GB")
```

## Getting Help

### Collecting diagnostic information

When reporting issues, include:

1. **Version information**:
```bash
uv run metainformant --version
python --version
```

2. **Error traceback**:
```bash
uv run python -c "import metainformant; ..." 2>&1
```

3. **Configuration (sanitized)**:
```bash
cat config.yaml | grep -v password  # Remove sensitive info
```

4. **System info**:
```bash
uname -a
free -h
nproc
```

### Additional Resources

- [FAQ](docs/FAQ.md) - Common questions and answers
- [Documentation Guide](DOCUMENTATION_GUIDE.md) - How to use documentation
- [GitHub Issues](https://github.com/username/metainformant/issues) - Report bugs
- [Discussions](https://github.com/username/metainformant/discussions) - Ask questions