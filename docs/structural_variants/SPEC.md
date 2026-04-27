# Technical Specification: Structural Variants

## Public API

```python
def detect_svs(bam: str, reference: str, callers: list[str], output_dir: str) -> pd.DataFrame
def annotate_svs(sv_calls: pd.DataFrame, gencode_gff: str, repeatmasker: str, gnomad_sv: str) -> pd.DataFrame
def merge_caller_outputs(calls_dict: dict[str,str], strategy: str, min_callers: int) -> list[SVRecord]
def filter_svs(annotated: pd.DataFrame, *, min_size: int, quality_threshold: float, exclude_common: bool) -> pd.DataFrame
def aggregate_cohort(vcf_dir: str, metadata: str) -> CohortMatrix
def burden_test(cohort: CohortMatrix, gene_intervals: str, phenotype: str, weight_col: str) -> pd.DataFrame
```

## Data model

`StructuralVariant` fields: chrom, start (1-based), end, sv_type, caller, quality, filters, info, sample_id.
`SVRecord` adds: genes, gene_types, is_coding, is_known, repeat_mask, impact_score.

## Merging logic

Reciprocal overlap threshold default 0.8:
  merged.start = min(sv.start for sv in cluster)
  merged.end   = max(sv.end for sv in cluster)
  merged.quality = Σ(weight_i × quality_i) / Σ(weights)
  merged.sv_type = majority({sv.sv_type})
  merged.info['CALLERS'] = ','.join(sorted(set(sv.caller for sv in cluster)))

## Annotation pipeline steps

1. Gene overlap — pybedtools.intersect against GENCODE GFF; collects all overlapping
   transcripts and the most severe functional class among them.
2. RepeatMasker — UCSC `rmsk.txt` converted to BED; classification: LINE, SINE, LTR, DNA.
3. gnomAD-SV — frequencies for population-specific rarity filter.
4. Impact scoring — rule-based:

   | Feature | Points |
   |---------|--------|
   | Overlaps protein-coding gene | +3 |
   | Frameshift / exon-loss DEL/INS | +2 |
   | Present in ClinGen Dosage Sensitivity | +2 |
   | Found in gnomAD-SV at AF >0.1 % | -1 |
   | Overlaps segmental duplication | -1 |
   | In telomere/centromere | -2 |

   Total ≥4 → High impact; 2–3 Medium; 0–1 Low.

## Configuration

YAML schema fully described in CONFIGURATION.md.

## Input formats

Accepted: BAM/CRAM (coordinate-sorted+indexed), VCF (caller output, bgzip + tabix).

## Output formats

| Format | Tool | Purpose |
|--------|------|---------|
| VCF.gz | all callers | downstream compatibility with GATK, bcftools |
| TSV | annotation | easy pandas loading, Excel viewing |
| BED | integration | overlaps with other interval sets (peaks, regulatory) |
| MAF | cBioPortal | mutation mapper format |

## Performance numbers (hg38 WGS, 30×)

- Manta: 2.5 h (32 GB), DELLY: 45 m (16 GB), LUMPY: 20 m (8 GB)
- Annotation: 45 s for 10 k SV calls.
- Merge: 1–2 min for 3-callers ensemble.

## Limitations

- Cannot call SVs <50 bp (use indel caller).
- Translocation detection requires at least one split-read supporting both breakends.
- Very high-coverage WGS (>80×) may cause Manta OOM; use `--maxInputSize` to limit.

## Future work

- Sniffles2 (long-read) integration, graph-based genotyping (SViply), ONT-specific pattern mining.