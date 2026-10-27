# Technical Specification: Structural Variants

## Public API

Four of the five subpackages (`detection`, `annotation`, `filtering`, `visualization`) only expose their modules, so import functions from the module path, for example `from metainformant.structural_variants.detection.sv_calling import call_structural_variants`. `population` also re-exports its six functions (`from metainformant.structural_variants.population import genotype_sv_population`).

```python-snippet
# detection.sv_calling
def call_structural_variants(alignments: list[dict], min_support: int = 3, insert_size_stats: InsertSizeStats | None = None, min_mapq: int = 20, min_sv_size: int = 50) -> list[StructuralVariant]
# detection.cnv
def detect_cnv_from_depth(depth_data: dict[str, Any], window_size: int = 1000, method: str = "segmentation", significance: float = 0.01, ploidy: int = 2, min_segment_bins: int = 3, merge_distance: int = 1000) -> dict[str, CNVResult]
# detection.breakpoints
def refine_breakpoints(variants: list[dict], reads: list[dict], window: int = 200, min_clip: int = 5) -> list[BreakpointPair]
# annotation.overlap
def annotate_gene_overlap(variants: list[dict], gene_db: list[dict] | list[GenomicInterval]) -> list[dict]
# annotation.functional_impact
def predict_functional_impact(variant: dict, gene_annotations: list[dict], haploinsufficiency_db: dict[str, float] | None = None, tad_boundaries: list[dict] | None = None) -> FunctionalImpact
# filtering.quality_filter
def filter_by_quality(variants: list[dict], min_qual: float = 20.0, min_support: int = 3, min_mapq: float = 0.0, min_gq: float = 0.0) -> tuple[list[dict], FilterStats]
# filtering.merge
def merge_callsets(callsets: dict[str, list[dict]], min_overlap: float = 0.5, type_match: bool = True) -> tuple[list[MergedVariant], MergeStats]
# population.sv_population
def genotype_sv_population(sv_calls: list[dict], samples: list[str], method: str = "depth") -> dict
def sv_association_test(genotypes: list[int], phenotypes: list[float], covariates: list[list[float]] | None = None) -> dict
```

Further functions (`segment_coverage`, `call_cnv_states`, `calculate_log2_ratio`, `detect_split_reads`, `detect_discordant_pairs`, `classify_sv_type`, `genotype_sv`, `annotate_regulatory_overlap`, `find_nearest_gene`, `filter_by_size`, `filter_by_frequency`, `apply_blacklist`, `survivor_merge`, `deduplicate_variants`, `sv_allele_frequency`, `sv_population_structure`, `sv_ld_analysis`, `merge_sv_callsets` and the `plot_*` functions) are listed in [ARCHITECTURE.md](ARCHITECTURE.md).

## Data model

`StructuralVariant` (`detection.sv_calling`) fields: chrom, start (0-based), end (0-based), sv_type (`SVType`: DEL, DUP, INV, TRA, INS, BND, UNKNOWN), size, quality, evidence (`SVEvidence`: split_reads, discordant_pairs, depth_support, evidence_reads, breakpoint1, breakpoint2, chrom1, chrom2, strand1, strand2), genotype, filter_status, chrom2, info.
`MergedVariant` (`filtering.merge`) fields: chrom, start, end, sv_type, size, n_callers, caller_names, support_variants, consensus_quality, genotype.
The filters, annotation and merge functions take plain dicts with the keys `chrom`, `start`, `end`, `sv_type` and, depending on the function, `quality` (or `qual`), `support` (or `n_support`), `size`, `mapq`, `gq`, `af`. Use `dataclasses.asdict()` to turn `StructuralVariant` objects into such dicts.

## Merging logic

`merge_callsets` default reciprocal overlap threshold 0.5 (`deduplicate_variants` uses 0.8 within a single callset). Two variants from the same caller are never compared directly, but clusters are formed transitively (union-find), so a call from another caller can bridge two same-caller calls into one cluster:

```
merged.start = median(sv.start for sv in cluster)
merged.end   = median(sv.end for sv in cluster)
merged.consensus_quality = mean(sv.quality) + 10 * log10(n_callers + 1) * 10
merged.sv_type = majority({sv.sv_type})
merged.genotype = majority vote of the non-missing genotypes
merged.caller_names = unique callers in the cluster
```

`survivor_merge` matches on breakpoint distance (`max_distance`, default 1000 bp) instead of overlap and keeps clusters seen by at least `min_callers` callers.

## Annotation pipeline steps

1. Gene overlap: `annotate_gene_overlap` builds an `IntervalIndex` from the gene entries you pass (dicts with `chrom`, `start`, `end`, `name` or `GenomicInterval` objects) and adds `overlapping_genes`, `gene_overlaps` (`OverlapResult` per gene: overlap in bp, overlap fractions of variant and gene, relationship) and `n_genes_affected`. `annotate_regulatory_overlap` does the same for regulatory elements.
2. Dosage sensitivity: `assess_dosage_sensitivity` classifies a gene as haploinsufficient (HI ≥ 0.5, pLI ≥ 0.9 or LOEUF < 0.35) or triplosensitive (TS ≥ 0.5) from score dictionaries you supply. `predict_functional_impact` passes only the haploinsufficiency dictionary.
3. TAD disruption: `predict_tad_disruption` flags boundaries (midpoint of each `chrom`/`start`/`end` entry) with a rule that depends on the SV type: a DEL spanning the midpoint or with a breakpoint within `boundary_window` (10 kb), an INV/TRA/BND with a breakpoint within the window, a DUP spanning the midpoint ([ARCHITECTURE.md](ARCHITECTURE.md) section 8.4).
4. Impact scoring, rule-based: `predict_functional_impact` assigns an impact level and type:

   | Level | Type | Condition |
   |-------|------|-----------|
   | HIGH | `dosage_loss` | DEL/HOMODEL overlapping a haploinsufficient gene |
   | HIGH | `gene_disruption` | DEL/HOMODEL overlapping a coding gene |
   | HIGH | `gene_fusion_candidate` | TRA overlapping a gene |
   | MODERATE | `gene_duplication`, `gene_inversion`, `gene_insertion` | DUP, INV or INS overlapping a gene |
   | MODERATE | `tad_disruption` | A TAD boundary is disrupted |
   | LOW | `large_intergenic` | Size > 100 kb and no rule above matched |
   | LOW | `intergenic` | No overlapping genes and size > 10 kb |
   | MODIFIER | `benign` | Anything else |

   `score_pathogenicity` then combines six weighted lines of evidence into a [0, 1] score: gene impact 30%, size 15%, SV type 15%, impact level 15%, population frequency 15% and TAD disruption 10%.

## Configuration

Thresholds are function arguments. Five `SV_` environment variables (`SV_MIN_SUPPORT`, `SV_MIN_QUAL`, `SV_CNV_SIGNIFICANCE`, `SV_CNV_DEL_THRESHOLD`, `SV_CNV_DUP_THRESHOLD`) override the matching argument at call time. There is no YAML schema; see [CONFIGURATION.md](CONFIGURATION.md).

## Input formats

- Alignment records: dicts with `name`, `chrom`, `pos` (0-based), `mapq`, `cigar` (a CIGAR string, or a list of `(op_code, length)` tuples with BAM integer codes such as 0 = M and 4 = S; `refine_breakpoints` accepts strings only) and optionally `read_length`, `mate_chrom`, `mate_pos`, `insert_size`, `is_reverse`, `mate_is_reverse`, `sa_tag`, `seq`. You build them from BAM/CRAM yourself (pysam is not a dependency of the detection code).
- Depth arrays: `{chrom: per-bin depth}` for `detect_cnv_from_depth`.
- Variant, gene, TAD, blacklist and population-frequency data: Python dicts and lists (see CONFIGURATION.md for the layouts).
- VCF: only `survivor_merge` reads VCF files (through pysam when installed, otherwise a basic text parser); VCF writing uses `dna.variation.variants.write_vcf`.

## Output formats

| Format | Produced by | Purpose |
|--------|-------------|---------|
| Python dataclasses and dicts | all functions | `StructuralVariant`, `CNVResult`, `FunctionalImpact`, `MergedVariant`, `FilterStats`, population dicts |
| JSON | `core.io.dump_json` on `dataclasses.asdict()` output | Machine-readable results |
| VCF | `dna.variation.variants.write_vcf` (not part of this module; see EXAMPLES.md Example 14) | Downstream compatibility with bcftools and similar tools |
| PNG/PDF | `visualization.plots.plot_*` | Figures; the format follows the output file extension |

## Performance

The module ships no benchmark suite. [PERFORMANCE.md](PERFORMANCE.md) describes what drives the cost of each stage and gives indicative timings measured on synthetic data.

## Limitations

- Input is in-memory dict records and arrays: the module has no BAM/CRAM reader and no streaming mode.
- `call_structural_variants` drops calls shorter than `min_sv_size` (default 50 bp), except translocations. Insertion calls always have a breakpoint offset below 50, so they only appear when `min_sv_size` is lowered, and no inserted sequence is reported.
- Translocations come from split reads whose SA tag points to another chromosome, or from inter-chromosomal discordant pairs.
- CBS segmentation is pure Python and slows quickly with the number of bins per chromosome; use a window size that leaves a few thousand bins.
- `predict_functional_impact` evaluates haploinsufficiency only (no triplosensitivity); the pathogenicity score is a triage heuristic, not a clinical classification.

## Not implemented

Long-read-specific SV calling in this module, graph-based genotyping and integrated annotation against RepeatMasker, GENCODE or gnomAD-SV files (you supply any annotation data as Python objects). Long-read SV calling from split alignments is in `metainformant.longread.analysis.structural`.
