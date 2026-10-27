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
def predict_functional_impact(variant: dict, gene_annotations: list[dict], haploinsufficiency_db: dict[str, float] | None = None, tad_boundaries: list[dict] | None = None, triplosensitivity_db: dict[str, float] | None = None, pli_db: dict[str, float] | None = None, loeuf_db: dict[str, float] | None = None) -> FunctionalImpact
# filtering.quality_filter
def filter_by_quality(variants: list[dict], min_qual: float = 20.0, min_support: int = 3, min_mapq: float = 0.0, min_gq: float = 0.0) -> tuple[list[dict], FilterStats]
# filtering.merge
def merge_callsets(callsets: dict[str, list[dict]], min_overlap: float = 0.5, type_match: bool = True) -> tuple[list[MergedVariant], MergeStats]
# population.sv_population
def genotype_sv_population(sv_calls: list[dict], samples: list[str], method: str = "depth") -> dict
def sv_association_test(genotypes: list[int], phenotypes: list[float], covariates: list[list[float]] | None = None) -> dict
```

Further functions (`segment_coverage`, `call_cnv_states`, `calculate_log2_ratio`, `detect_split_reads`, `detect_discordant_pairs`, `classify_sv_type`, `genotype_sv`, `annotate_regulatory_overlap`, `find_nearest_gene`, `filter_by_size`, `filter_by_frequency`, `apply_blacklist`, `survivor_merge`, `deduplicate_variants`, `sv_allele_frequency`, `sv_population_structure`, `sv_ld_analysis`, `merge_sv_callsets` and the `plot_*` functions) and the `ReadIndex` and `IntervalIndex` classes are listed in [ARCHITECTURE.md](ARCHITECTURE.md).

## Data model

`StructuralVariant` (`detection.sv_calling`) fields: chrom, start (0-based), end (0-based), sv_type (`SVType`: DEL, DUP, INV, TRA, INS, BND, UNKNOWN), size, quality, evidence (`SVEvidence`: split_reads, discordant_pairs, depth_support, evidence_reads, breakpoint1, breakpoint2, chrom1, chrom2, strand1, strand2), genotype, filter_status, chrom2, info, pos2.

- Calls other than translocations (DEL, DUP, INV, INS, UNKNOWN) occupy `chrom:start-end`; their `chrom2` is `""` and `pos2` is `None`.
- A translocation (TRA) joins two loci on different chromosomes and has no interval: `chrom` and `start` give the first breakpoint, `chrom2` and `pos2` the second, `end` repeats `start`, and `size` is 0.
- `size == 0` means "not measured". It is 0 for TRA and for INS, because split-read evidence does not measure an inserted length.

`MergedVariant` (`filtering.merge`) fields: chrom, start, end, sv_type, size, n_callers, caller_names (sorted), support_variants, consensus_quality, genotype.
The filters, annotation and merge functions take plain dicts with the keys `chrom`, `start`, `end`, `sv_type` and, depending on the function, `quality` (or `qual`), `support` (or `n_support`), `size`, `mapq`, `gq`, `af`. Use `dataclasses.asdict()` to turn `StructuralVariant` objects into such dicts.

## Merging logic

`merge_callsets` default reciprocal overlap threshold 0.5 (`deduplicate_variants` uses 0.8 within a single callset). Two variants from the same caller are never compared directly, but clusters are formed transitively (union-find), so a call from another caller can bridge two same-caller calls into one cluster:

```
merged.start = median(sv.start for sv in cluster)
merged.end   = median(sv.end for sv in cluster)
merged.consensus_quality = mean(sv.quality) + 10 * log10(n_callers + 1) * 10
merged.sv_type = majority({sv.sv_type})
merged.genotype = majority vote of the non-missing genotypes (a tie goes to the genotype seen first)
merged.caller_names = sorted unique callers in the cluster
```

`survivor_merge` matches on breakpoint distance (`max_distance`, default 1000 bp) instead of overlap and keeps clusters seen by at least `min_callers` callers. It takes the median start and end and the majority type as above, and combines the other fields as follows:

```
merged.genotype = majority vote of the non-missing genotypes (allele order ignored);
                  a tie goes to the genotype with the larger summed call quality, then to the one that sorts first
merged.consensus_quality = mean(sv.quality) + 10 * n_callers * log10(n_callers + 1)
```

Genotypes come from each dict's `genotype` key, or from the first sample of a VCF file (a phased `0|1` becomes `0/1`); `./.` means that no call has one. Two callers with calls of quality 80 (genotype `0/1`) and 40 (genotype `1/1`) give genotype `0/1` and `consensus_quality` 69.54 (a mean of 60 plus a bonus of 20 · log10(3)).

## Annotation pipeline steps

1. Gene overlap: `annotate_gene_overlap` builds an `IntervalIndex` from the gene entries you pass (dicts with `chrom`, `start`, `end`, `name` or `GenomicInterval` objects) and adds `overlapping_genes`, `gene_overlaps` (`OverlapResult` per gene: overlap in bp, overlap fractions of variant and gene, relationship) and `n_genes_affected`. `annotate_regulatory_overlap` does the same for regulatory elements. The index is a max-end segment tree, so it finds every overlapping interval however long, and `query_nearest` considers every interval of the chromosome.
2. Dosage sensitivity: `assess_dosage_sensitivity` classifies a gene as haploinsufficient (HI ≥ 0.5, pLI ≥ 0.9 or LOEUF < 0.35) or triplosensitive (TS ≥ 0.5) from score dictionaries you supply. `predict_functional_impact` takes all four dictionaries (`haploinsufficiency_db`, `triplosensitivity_db`, `pli_db`, `loeuf_db`; each optional) and marks an overlapped gene as dosage-sensitive when it is haploinsufficient and the variant is a DEL or HOMODEL, or triplosensitive and the variant is a DUP. Without any dictionary no gene is marked.
3. TAD disruption: `predict_tad_disruption` flags boundaries (midpoint of each `chrom`/`start`/`end` entry) with a rule that depends on the SV type: a DEL spanning the midpoint or with a breakpoint within `boundary_window` (10 kb), an INV/TRA/BND with a breakpoint within the window, a DUP spanning the midpoint ([ARCHITECTURE.md](ARCHITECTURE.md) section 8.4). It reads the variant's `chrom`, `start` and `end` only, so a translocation (`end == start`) is tested at its first breakpoint.
4. Impact scoring, rule-based: `predict_functional_impact` assigns an impact level and type:

   | Level | Type | Condition |
   |-------|------|-----------|
   | HIGH | `dosage_loss` | DEL/HOMODEL overlapping a haploinsufficient gene |
   | HIGH | `dosage_gain` | DUP overlapping a triplosensitive gene |
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

- Alignment records: dicts with `name`, `chrom`, `pos` (0-based), `mapq`, `cigar` (a CIGAR string, or a list of `(op_code, length)` tuples with BAM integer codes such as 0 = M and 4 = S; `call_structural_variants` and `refine_breakpoints` accept both) and optionally `read_length`, `mate_chrom`, `mate_pos`, `insert_size`, `is_reverse`, `mate_is_reverse`, `sa_tag`, `seq`. `call_structural_variants` ignores `seq`; `refine_breakpoints` reads `is_reverse` (breakpoint strand), `seq` (microhomology and inserted sequence) and `sa_tag` (inserted sequence). You build the records from BAM/CRAM yourself (pysam is not a dependency of the detection code).
- Depth arrays: `{chrom: per-bin depth}` for `detect_cnv_from_depth`, finite and non-negative (zero is an uncovered bin). Log2 ratios are not depths: a negative or non-finite value raises `ValueError` naming the chromosome.
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
- `call_structural_variants` applies `min_sv_size` (default 50 bp) to DEL, DUP, INV and UNKNOWN calls, measured as the breakpoint distance. Translocations have no length and the inserted length of an insertion is not measured, so TRA and INS calls are exempt from it and have `size` 0. It reports no inserted sequence; `refine_breakpoints` derives one from reads that carry `seq` and `sa_tag`.
- Translocations come from split reads whose SA tag points to another chromosome, or from inter-chromosomal discordant pairs. A TRA call carries its second breakpoint in `chrom2` and `pos2`, but the functions that take variant dicts read `chrom`, `start` and `end` only: `survivor_merge` joins translocations whose first breakpoints are close whatever their partner, `merge_callsets` and `deduplicate_variants` give a zero-length call a reciprocal overlap of 0 (so with a threshold above 0 they never join or collapse translocations), `filter_by_frequency` matches no database entry to a zero-length call unless `min_reciprocal_overlap` is 0, and `predict_tad_disruption` tests the first breakpoint only. Only `plot_circos` and `refine_breakpoints` use `pos2`.
- `call_structural_variants` genotypes with `genotype_sv`'s default `min_mapq=20`, whatever `min_mapq` you pass; the reference reads it counts are the reads that span the 200 bp window around the first breakpoint.
- CBS segmentation is vectorised with NumPy, but its cost still grows with the square of the number of bins per chromosome; use a window size that leaves tens of thousands of bins at most (see [PERFORMANCE.md](PERFORMANCE.md)).
- `refine_breakpoints` derives `microhomology` with a mirror-match heuristic on the read sequence around each clip junction, not from the reference, so it is not true microhomology.
- The pathogenicity score is a triage heuristic, not a clinical classification.

## Not implemented

Long-read-specific SV calling in this module, graph-based genotyping and integrated annotation against RepeatMasker, GENCODE or gnomAD-SV files (you supply any annotation data as Python objects). Long-read SV calling from split alignments is in `metainformant.longread.analysis.structural`.
