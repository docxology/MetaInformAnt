# Structural Variants Module

Structural variant analysis: CNV detection via circular binary segmentation, SV calling from split/discordant reads, annotation, filtering, population genotyping, and visualization.

## Architecture

```mermaid
graph TD
    subgraph "Structural Variants Module"
        DT[detection/] --> |cnv.py| CN[CNV Detection (CBS)]
        DT --> |sv_calling.py| SC[SV Calling (Split/Discordant)]
        DT --> |breakpoints.py| BP[Breakpoint Refinement]

        AN[annotation/] --> |overlap.py| OV[Gene & Regulatory Overlap]
        AN --> |functional_impact.py| FI[Dosage Sensitivity / TAD]

        FL[filtering/] --> |quality_filter.py| QF[Quality & Size Filtering]
        FL --> |merge.py| MG[Multi-Caller Consensus]

        PO[population/] --> |sv_population.py| SP[Population Genotyping]

        VZ[visualization/] --> |plots.py| PL[Circos, Coverage, CNV Plots]
    end

    DT --> AN
    AN --> FL
    FL --> PO
```

## Submodules

| Module | Purpose |
|--------|---------|
| [`detection/`](detection/) | CNV detection (CBS), SV calling (split/discordant reads), breakpoint refinement |
| [`annotation/`](annotation/) | Gene/regulatory overlap annotation, functional impact, TAD disruption |
| [`filtering/`](filtering/) | Quality filtering, size filtering, blacklist regions, multi-caller merging |
| [`population/`](population/) | Population genotyping, allele frequency, association testing, PCA, LD |
| [`visualization/`](visualization/) | Circos plots, coverage tracks, size distributions, CNV profiles |

## Key Capabilities

### CNV Detection

| Class/Function | Description |
|----------------|-------------|
| `detect_cnv_from_depth` | Detect CNVs from read depth via circular binary segmentation |
| `CNVSegment` | Dataclass for a segment with log2 ratio, state, and copy number |
| `call_cnv_states` | Assign DEL/DUP/NEUTRAL/AMP/HOMODEL from log2 ratios |
| `merge_adjacent_segments` | Merge neighboring segments with similar copy number |
| `calculate_log2_ratio` | Compute GC-corrected log2 ratios from depth profiles |

### SV Calling

| Class/Function | Description |
|----------------|-------------|
| `call_structural_variants` | Full SV calling pipeline from aligned reads |
| `SVType` | Enum: DEL, DUP, INV, TRA, INS, BND |
| `StructuralVariant` | Dataclass with position, type, evidence, and quality |
| `detect_split_reads` | Extract split-read SV evidence from alignments |
| `detect_discordant_pairs` | Extract discordant read-pair SV evidence |
| `genotype_sv` | Genotype a single SV in one sample |

### Annotation

| Class/Function | Description |
|----------------|-------------|
| `annotate_gene_overlap` | Find genes overlapping structural variants |
| `annotate_regulatory_overlap` | Find regulatory elements overlapping SVs |
| `find_nearest_gene` | Find closest gene to an SV breakpoint |
| `IntervalIndex` | Efficient interval overlap queries via sorted index |

### Filtering

| Function | Description |
|----------|-------------|
| `filter_by_quality` | Remove low-quality SV calls by score threshold |
| `filter_by_size` | Filter SVs by minimum/maximum size |
| `filter_by_frequency` | Remove common population variants |
| `apply_blacklist` | Exclude calls in known problematic regions |

### Population Analysis

| Function | Description |
|----------|-------------|
| `genotype_sv_population` | Population-scale SV genotyping across samples |
| `sv_allele_frequency` | Compute allele frequencies per SV |
| `sv_association_test` | Test SV-phenotype association (linear/logistic) |
| `sv_population_structure` | PCA-based population structure from SV genotypes |
| `sv_ld_analysis` | Linkage disequilibrium between SVs and SNPs |
| `merge_sv_callsets` | Merge SV calls from multiple samples/callers |

## Quick Start

```python
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality
from metainformant.structural_variants.population.sv_population import sv_allele_frequency

# Detect CNVs from read depth
cnv_result = detect_cnv_from_depth(
    depths=[120, 115, 60, 55, 58, 110, 125],
    chrom="chr1",
    bin_size=10000,
)

# Call SVs from aligned reads
sv_calls = call_structural_variants(reads=aligned_reads, chrom="chr1")

# Filter and compute population frequencies
filtered = filter_by_quality(sv_calls, min_quality=20.0)
freq = sv_allele_frequency(sv_calls=filtered, samples=sample_ids)
```

## Integration

SV calls integrate with GWAS and expression analysis for functional interpretation:

```python
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap
from metainformant.structural_variants.population.sv_population import sv_association_test

# Annotate SVs with gene overlaps
annotated = annotate_gene_overlap(sv_calls=filtered, gene_intervals=gene_db)

# Test association with phenotype
assoc = sv_association_test(genotypes=gt_matrix, phenotypes=pheno_values)
```

## Related

- [docs/structural_variants/](../../../docs/structural_variants/) - Detection, annotation, filtering docs
- [metainformant.dna](../dna/) - Upstream alignment and variant calling
- [metainformant.gwas](../gwas/) - GWAS integration for SV-phenotype associations
- Config prefix: `SV_`
