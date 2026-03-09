# amalgkit quant: Troubleshooting Guide

Common issues and solutions for the `amalgkit quant` step.

**See**: [06_quant.md](06_quant.md) for core usage and parameters.

---

## Issue: "kallisto: command not found"

**Solutions**:
```bash
# Install kallisto (system tool)
conda install -c bioconda kallisto
# Or using system package manager:
sudo apt-get install -y kallisto

# Verify installation
which kallisto
kallisto version

# Note: For Python packages, use uv pip install (primary method)
```

## Issue: Index file not found

```
Error: Could not find index file
```

**Diagnosis**:
```bash
# Check index directory
ls -lh output/amalgkit/work/index/

# Check species name in metadata
cut -f2 output/amalgkit/work/metadata/pivot_qualified.tsv | sort -u
# Should match index filename (with spaces → underscores)
```

**Solutions**:
1. Build index:
   ```bash
   amalgkit quant --build_index yes --fasta_dir output/work/fasta
   ```

2. Download pre-built index and place in correct location:
   ```bash
   cp Apis_mellifera_transcripts.idx output/work/index/
   ```

3. Specify index directory explicitly:
   ```bash
   --index_dir /path/to/indices
   ```

## Issue: Low alignment rate

```json
{
    "p_pseudoaligned": 15.2  // <30% is concerning
}
```

**Causes**:
1. Wrong reference transcriptome (different species/version)
2. Low-quality sequencing data
3. Contamination from other organisms
4. Incorrect library type

**Solutions**:
1. Verify reference matches species:
   ```bash
   # Check what's in your index
   grep "^>" output/work/fasta/Apis_mellifera_rna.fasta | head
   ```

2. Check FASTQ quality:
   ```bash
   # Inspect fastp reports
   cat output/work/getfastq/SRR*/fastp.json | grep "total_reads"
   ```

3. Try different reference version:
   ```bash
   # Download different annotation release
   datasets download genome taxon 7460 --include rna
   ```

## Issue: Out of memory

```
kallisto: std::bad_alloc
```

**Solutions**:
1. Increase memory allocation:
   ```bash
   #SBATCH --mem=16G  # Instead of 8G
   ```

2. Reduce number of concurrent jobs

3. Use smaller index (subset of transcriptome)

## Issue: FASTQ files missing

```
Error: Could not find FASTQ files for SRR12345678
```

**Solutions**:
1. Verify getfastq completed:
   ```bash
   ls output/work/getfastq/SRR12345678/
   ```

2. Check if already cleaned (from previous quant run):
   ```bash
   # Check if abundance files exist
   ls output/work/quant/SRR12345678/abundance.tsv
   # FASTQs deleted after successful quant with --clean_fastq yes
   ```

3. Re-run getfastq:
   ```bash
   amalgkit getfastq --id SRR12345678 --out_dir output/work
   ```

---

**See Also**: [06_quant.md](06_quant.md) | [06_quant_advanced.md](06_quant_advanced.md)
