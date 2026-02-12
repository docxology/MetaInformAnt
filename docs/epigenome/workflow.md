# Epigenome Workflow Orchestration

End-to-end workflow pipelines for methylation, ChIP-seq, and ATAC-seq analysis. Provides configuration management via the `EpigenomeConfig` dataclass and multi-assay integration.

## Key Concepts

**Workflow orchestration** chains individual analysis steps (loading, QC, peak calling, annotation) into reproducible pipelines with consistent configuration and output paths.

**EpigenomeConfig** centralizes all parameters for an epigenome analysis session: input paths, output directory, quality thresholds, and assay-specific settings.

**Multi-assay integration** combines results from methylation, ChIP-seq, and ATAC-seq into a unified view of the epigenomic landscape, identifying concordant and discordant regulatory signals.

## Configuration

### `EpigenomeConfig`

Dataclass for workflow configuration.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | `str` | `"output/epigenome"` | Output directory |
| `genome` | `str` | `"hg38"` | Reference genome |
| `bin_size` | `int` | `200` | Genomic bin size in bp |
| `n_states` | `int` | `10` | Number of chromatin states |
| `min_peak_score` | `float` | `5.0` | Minimum peak score |
| `max_q_value` | `float` | `0.05` | Maximum q-value for peaks |
| `min_fold_enrichment` | `float` | `2.0` | Minimum fold enrichment |
| `dmr_min_diff` | `float` | `0.2` | Minimum beta difference for DMRs |
| `dmr_window_size` | `int` | `1000` | Sliding window for DMR detection |

Environment variable prefix: `EPI_` (e.g., `EPI_BIN_SIZE=500`).

## Function Reference

### `run_methylation_workflow(config, methylation_data) -> Dict`

Execute a complete methylation analysis pipeline:

1. Load and validate CpG methylation data.
2. Compute beta values with regularization.
3. Generate per-chromosome summaries.
4. (If two groups) Find differentially methylated regions.
5. Identify CpG islands.
6. Calculate methylation entropy.
7. Write results and summary report.

**Returns** dict with `beta_values`, `chromosome_summary`, `dmrs`, `cpg_islands`, `entropy`, and `report_path`.

### `run_chipseq_workflow(config, peak_data) -> Dict`

Execute a ChIP-seq analysis pipeline:

1. Load peaks from narrowPeak/broadPeak files.
2. Filter by quality thresholds (score, q-value, fold enrichment).
3. Compute peak statistics.
4. (If two conditions) Find overlapping peaks and differential regions.
5. Calculate enrichment in genomic features.
6. Run motif discovery in peak sequences.
7. Write results.

**Returns** dict with `filtered_peaks`, `statistics`, `overlaps`, `enrichment`, `motifs`, and `report_path`.

### `run_atacseq_workflow(config, atac_data) -> Dict`

Execute an ATAC-seq analysis pipeline:

1. Load and filter accessibility peaks.
2. Analyze fragment size distribution and nucleosome positioning.
3. Calculate TSS enrichment as QC metric.
4. Detect transcription factor binding sites.
5. (If two conditions) Compare accessibility between conditions.
6. Write results.

**Returns** dict with `peaks`, `nucleosome_fractions`, `tss_enrichment`, `tf_sites`, `comparison`, and `report_path`.

### `integrate_epigenome_results(methylation_results, chipseq_results, atacseq_results) -> Dict`

Combine results from all three assay types into an integrated epigenomic annotation:

- Correlate methylation levels with chromatin accessibility.
- Identify regions with concordant signals (e.g., low methylation + open chromatin + active histone marks).
- Flag discordant regions for further investigation.
- Generate a unified regulatory element catalog.

**Returns** dict with `concordant_regions`, `discordant_regions`, `regulatory_catalog`, and `integration_summary`.

## Usage Examples

```python
from metainformant.epigenome import (
    EpigenomeConfig,
    run_methylation_workflow,
    run_chipseq_workflow,
    run_atacseq_workflow,
    integrate_epigenome_results,
)

# Configure the analysis
config = EpigenomeConfig(
    output_dir="output/epigenome/sample1",
    genome="hg38",
    n_states=8,
    min_fold_enrichment=2.0,
    dmr_min_diff=0.25,
)

# Run individual assay workflows
meth_results = run_methylation_workflow(config, methylation_data)
chip_results = run_chipseq_workflow(config, chipseq_peaks)
atac_results = run_atacseq_workflow(config, atac_peaks)

# Integrate across assays
integrated = integrate_epigenome_results(
    meth_results, chip_results, atac_results
)
print(f"Concordant regions: {len(integrated['concordant_regions'])}")
print(f"Regulatory catalog size: {len(integrated['regulatory_catalog'])}")
```

## Configuration

Environment variable prefix: `EPI_`

All workflow parameters can be overridden via environment variables following the pattern `EPI_PARAMETER_NAME`.

## Related Modules

- `metainformant.epigenome.methylation` -- DNA methylation analysis
- `metainformant.epigenome.chipseq` -- ChIP-seq peak analysis
- `metainformant.epigenome.atacseq` -- ATAC-seq accessibility analysis
- `metainformant.epigenome.chromatin_state` -- chromatin state learning
- `metainformant.epigenome.peak_calling` -- de novo peak calling
