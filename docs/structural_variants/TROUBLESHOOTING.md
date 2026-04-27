# Troubleshooting: Structural Variants

## Common failure modes

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
`OSError: Manta binary not found` | Manta not on `PATH` | Install Manta or set `SV_MANTA_PATH` env var |
`Zero SVs detected from 30× WGS` | BAM not indexed | `samtools index sample.bam` before calling |
`ValueError: reference sequence not found` | FASTA lacks .fai index | `samtools faidx ref.fa` |
`PermissionError` on temp dir | /tmp full or no exec perms | Set `SV_TMPDIR=/large_tmp` with enough space |
`ImportError: No module named 'cyvcf2'` | Missing optional extra | `uv pip install metainformant[sv]` |
`bgzip: command not found` | tabix/bgzip not installed | `conda install -c bioconda htslib` |

## Diagnostic checklist

1. Validate BAM
   ```bash
   samtools quickcheck -q all sample.bam
   samtools idxstats sample.bam | head
   ```

2. Test Manta config generation
   ```bash
   configManta.py --referenceFasta hg38.fa --bam sample.bam --runDir manta_test/
   ls manta_test/  # expect manta_workflow.py inside
   ```

3. Restrict to chr1 as sanity
   ```bash
   metainformant sv detect --bam sample.bam --region chr1:1-1000000 --out test_sv
   ```

4. Inspect VCF header
   ```bash
   bcftools view -h merged.vcf.gz | grep -E '^##|#CHROM'
   ```

## Known caller-specific quirks
- Manta: END for INS is POS+1; our normalizer adjusts back to POS-based length.
- DELLY: writes sample name from BAM header; if BAM renamed after alignment run, pass `--sample NEWNAME`.
- LUMPY: requires `sambamba` or `samtools` to extract split reads; error 'no split reads found' may mean BAM not marked for duplicates.
- CNVkit: target BED must match exactly the capture design used; off-target bins need `cnvkit.py export` first.

## Getting help
Provide: OS, Python version, metainformant version, caller versions (`manta --version`), full command line, and 10 kbp mini-BAM (chr subset) with which the bug can be reproduced.