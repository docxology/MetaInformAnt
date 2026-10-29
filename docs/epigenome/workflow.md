# Epigenome Workflow Orchestration

Directory-driven workflows for methylation, ChIP-seq, and ATAC-seq data, implemented in `metainformant.epigenome.workflow.workflow`. Provides configuration management via the `EpigenomeConfig` dataclass and a multi-assay integration step.

## Key Concepts

**Workflow orchestration** chains loading, filtering, summary statistics and report writing for every matching file in an input directory, with consistent configuration and output paths. The workflows summarize each sample; they do not call peaks or annotate them, and they do not run the differential, enrichment, motif or TSS analyses of the assay modules.

**EpigenomeConfig** collects the thresholds and options for a session. Only some of its fields are read by the workflow functions (see the field table below).

**Multi-assay integration** relates the three assays pairwise: where methylation sites fall relative to ChIP-seq and ATAC-seq peaks, how methylation differs between accessible and closed chromatin, and how much the ChIP-seq and ATAC-seq peak sets overlap.

## Configuration

### `EpigenomeConfig`

Dataclass for workflow configuration. Every field is an optional keyword argument; an unknown field name raises `TypeError`. `__post_init__` validates the numeric ranges noted below.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_dir` | `str \| Path` | `"data/epigenome"` | Not read by the workflow functions (they take `input_dir` as an argument) |
| `output_dir` | `str \| Path` | `"output/epigenome"` | Not read by the workflow functions (they take `output_dir` as an argument) |
| `temp_dir` | `str \| Path` | `"temp/epigenome"` | Not used |
| `methylation_threshold` | `float` | `0.8` | Validated to lie in 0--1; not otherwise used |
| `chipseq_qvalue_threshold` | `float` | `0.05` | ChIP-seq workflow: peaks with a q-value above it are dropped, and peaks without a q-value are kept. It is an ordinary FDR limit on the q-value as a probability (a q column of 2 in the file is 0.01); the peak score plays no part. Validated to lie in 0--1 |
| `atacseq_qvalue_threshold` | `float` | `0.05` | ATAC-seq workflow: the same q-value filter for the ATAC-seq peaks; the accessibility score (`signal_value`) plays no part. Validated to lie in 0--1 |
| `min_peak_length` | `int` | `50` | Minimum peak length in bp (ChIP-seq and ATAC-seq workflows) |
| `max_peak_length` | `int` | `10000` | Maximum peak length in bp (ChIP-seq and ATAC-seq workflows) |
| `min_methylation_coverage` | `int` | `5` | Minimum read coverage for a site to be loaded by the methylation workflow (BEDgraph files carry no read counts, so it removes none of their sites) |
| `min_chip_coverage` | `int` | `10` | Not used |
| `min_atac_coverage` | `int` | `5` | Not used |
| `max_distance_for_integration` | `int` | `1000` | Distance in bp within which a methylation site counts as proximal or near a peak in `integrate_epigenome_results` |
| `correlation_threshold` | `float` | `0.3` | Validated to lie in -1 to 1 and copied into the integration results; no correlation is computed |
| `run_methylation`, `run_chipseq`, `run_atacseq`, `integrate_results` | `bool` | `True` | Not read; call the workflow functions you need directly |
| `generate_reports` | `bool` | `True` | Write the text report of each workflow and of the integration step |
| `n_threads` | `int` | `4` | Validated to lie in 1--64; the workflows do not use threads |
| `memory_limit_gb` | `float` | `8.0` | Validated to lie in 0.1--1024; not otherwise used |
| `chunk_size` | `int` | `1000000` | Not used |

Environment-variable overrides (an `EPI_` prefix) are not implemented: `EPI_*` variables have no effect. Set fields in code or load them from a file with `load_epigenome_config`.

### `load_epigenome_config(config_path=None) -> EpigenomeConfig`

Build an `EpigenomeConfig` from a YAML, TOML or JSON file whose top-level keys are field names. Without `config_path` it uses the first of `config/epigenome.yaml`, `config/epigenome.yml` and `config/epigenome.json` that exists relative to the working directory. If no file is found the defaults are returned; if a file is found but cannot be loaded (for example because it contains an unknown key) a warning is logged and the defaults are returned.

## Function Reference

### `run_methylation_workflow(input_dir, output_dir, config=None) -> Dict`

Summarize every methylation file found under `input_dir` (searched recursively for `*.bedgraph`, `*.bedGraph`, `*.BEDGRAPH`, `*.bg` and `*.cov`). `config` defaults to `EpigenomeConfig()`; `output_dir` is created if needed. `input_dir` must exist, and `PipelineError` is raised if it holds no matching file.

1. Load each file with `load_methylation_bedgraph` (`.bedgraph`, `.bg`) or `load_methylation_cov` (`.cov`), keeping sites with at least `min_methylation_coverage` reads. BEDgraph files are read with `level_scale="auto"`, so a file with a level above 1 is read as percentages (with a warning). `.cov` files must have the layout `load_methylation_cov` reads (see the methylation page).
2. Compute `calculate_methylation_statistics` for each sample. A file that fails to load is recorded in `errors` and the remaining files are still processed.
3. Write `<sample>_methylation.json` for each sample and `methylation_summary.json`.
4. If `generate_reports` is true, write `methylation_report.txt` for the sites of all samples combined.

Beta values, DMRs, CpG islands and entropy are not computed here; call the functions of `metainformant.epigenome.assays.methylation` on the loaded sites.

**Returns** dict with `workflow_type` (`"methylation"`), `input_dir`, `output_dir`, `config` (the configuration as a dict), `processed_files`, `statistics` (`total_samples`, `total_sites` and `sample_statistics`, one entry with `sample`, `file` and `statistics` per sample), `sites` (one dict per loaded site of every processed sample, with `sample`, `chrom`, `position`, `methylation_level` and `coverage`), `errors`, and `report_path` when a report was written. `integrate_epigenome_results` reads `sites`.

### `run_chipseq_workflow(input_dir, output_dir, config=None) -> Dict`

Summarize every ChIP-seq peak file found under `input_dir` (searched recursively for `*.narrowPeak`, `*.broadPeak`, `*.bed` and the lower-case `.narrowpeak` / `.broadpeak` spellings). The peak layout of each file is detected from its name and its first data line: a name containing "narrow" means narrowPeak, one containing "broad" broadPeak, a `.bed` extension BED, and any other name narrowPeak. The layout is never richer than the first data line supports (a narrowPeak-named file whose first data line has nine columns is read as broadPeak, one with fewer as BED, and a broadPeak-named file with fewer than nine columns as BED), so a mislabelled file is read for what it holds instead of silently yielding no peaks. `PipelineError` is raised if no file is found.

1. Load each file with `load_chip_peaks`.
2. Filter: drop peaks whose q-value is above `chipseq_qvalue_threshold`, then keep peaks between `min_peak_length` and `max_peak_length` bp. The q-value filter is an ordinary FDR threshold on the q-value as a probability (a q column of 2 in the file is 0.01 and passes a threshold of 0.01). A peak exactly on the threshold is kept: the comparison allows a relative 1e-9 for the round-off of the `-log10` column. A peak without a q-value (a BED file, or `-1` in the column) cannot be judged and is kept, and a warning counts such peaks. A q-value above 1 is not a probability and draws a warning, but peaks read from files never have one. Peaks are never dropped for their score, and file order is preserved.
3. Compute `calculate_peak_statistics` for each sample. A file that fails is recorded in `errors`.
4. Write `<sample>_chipseq.json`, the filtered peaks as `<sample>_peaks.narrowPeak` (written by `save_chip_peaks` in the narrowPeak layout, with `-1` for a missing p-value, q-value or summit), and `chipseq_summary.json`.
5. If `generate_reports` is true, write `chipseq_report.txt` for the filtered peaks of all samples.

Overlap, enrichment and motif analyses are not run here; call the functions of `metainformant.epigenome.assays.chipseq`.

**Returns** dict with `workflow_type` (`"chipseq"`), `input_dir`, `output_dir`, `config`, `processed_files`, `statistics` (`total_samples`, `total_peaks` and `sample_statistics`, one entry with `sample`, `file`, `statistics` and `peaks_count` per sample), `peaks` (one dict per kept peak of every processed sample, with `sample`, `chrom`, `start`, `end`, `score`, `signal_value` and `q_value`), `errors`, and `report_path` when a report was written. `integrate_epigenome_results` reads `peaks`.

### `run_atacseq_workflow(input_dir, output_dir, config=None) -> Dict`

Summarize every ATAC-seq peak file under `input_dir`. File discovery and format detection are the same as for `run_chipseq_workflow`, so keep ChIP-seq and ATAC-seq files in separate directories. `PipelineError` is raised if no file is found.

1. Load each file with `load_atac_peaks`.
2. Filter: drop peaks whose q-value is above `atacseq_qvalue_threshold`, then keep peaks between `min_peak_length` and `max_peak_length` bp. The q-value filter works as in `run_chipseq_workflow`: an ordinary FDR threshold on the q-value as a probability, with a relative 1e-9 tolerance for the `-log10` round-off. Peaks without a q-value are kept with a warning, the accessibility score plays no part, and file order is preserved.
3. Compute `calculate_atac_statistics` for each sample, which includes the peak-length-based nucleosome fractions (`nfr_peak_fraction` and similar), not fragment sizes. A file that fails is recorded in `errors`.
4. Write `<sample>_atacseq.json`, the filtered peaks as `<sample>_peaks.narrowPeak` (written by `save_atac_peaks`, with the same layout), and `atacseq_summary.json`.
5. If `generate_reports` is true, write `atacseq_report.txt`.

TSS enrichment and condition comparison are not run here; call `identify_tss_enrichment` and `compare_atac_conditions` from `metainformant.epigenome.assays.atacseq`.

**Returns** dict with `workflow_type` (`"atacseq"`), `input_dir`, `output_dir`, `config`, `processed_files`, `statistics` (`total_samples`, `total_peaks` and `sample_statistics`), `peaks` (one dict per kept peak of every processed sample, with the same keys as the ChIP-seq `peaks`), `errors`, and `report_path` when a report was written. `integrate_epigenome_results` reads `peaks`.

### `integrate_epigenome_results(methylation_results, chipseq_results, atacseq_results, output_dir, config=None) -> Dict`

Combine the result dicts of the three workflows with pairwise association analyses. An analysis runs when both of its inputs report `statistics["total_samples"] > 0`:

- `methylation_chipseq`: classifies each methylation site as overlapping a ChIP-seq peak, proximal (within `max_distance_for_integration` bp) or distal; a site whose chromosome has no peak is counted in none of the three classes. A per-chromosome interval index supports these queries; missing inputs produce warning findings rather than inferred associations.
- `methylation_atacseq`: counts methylation sites inside, near and outside accessible regions and compares their mean methylation (`mean_methylation_accessible`, `mean_methylation_near_accessible`, `mean_methylation_closed`).
- `chipseq_atacseq`: how much the ChIP-seq peaks and the ATAC-seq accessible regions coincide, with the peaks of all samples of an assay pooled. Its `statistics` are `total_chip_peaks` and `total_atac_peaks`; `chip_peaks_overlapping_atac` and `atac_peaks_overlapping_chip` (peaks sharing at least one base with the other assay) and their rates `chip_overlap_rate` and `atac_overlap_rate`; and the base-pair overlap of the two peak sets: `chip_bases` and `atac_bases` (each assay's peaks merged first, so no base counts twice), `shared_bases`, `union_bases` and `jaccard_index` (`shared_bases / union_bases`). A single finding of type `"overlap"` summarizes them. Entries without a chromosome, `start` or `end`, or with `end <= start`, are ignored.

Methylation-to-peak association uses closed intervals: both endpoints count as overlapping, and `start == end` represents one genomic position. The ChIP-seq/ATAC-seq base-overlap summary instead requires positive half-open extents as described above.

The analyses read lists from the input dicts, which the `run_*_workflow` functions return: `sites` in the methylation results (dicts with `chrom` or `chromosome`, `position` or `pos`, and optionally `methylation_level` or `beta`), and `peaks` in the ChIP-seq and ATAC-seq results (dicts with `chrom` or `chromosome`, `start` and `end`; `accessible_regions` is accepted as an alternative key for ATAC-seq). Pass the workflow output as it is, or supply lists of that shape. A missing or empty list gives only a warning finding for the analyses that need it.

**Returns** dict with `integration_type`, `output_dir`, `config`, `methylation_summary`, `chipseq_summary` and `atacseq_summary` (the `statistics` of each input), `integrated_analyses` (analysis name mapped to a dict with `analysis_type`, `max_distance`, `correlation_threshold`, `findings` and `statistics`, which is empty when the analysis only reports a warning), `errors`, and `report_path` when `generate_reports` is true. It writes `integrated_epigenome_results.json` and, with reports enabled, `epigenome_integration_report.txt` into `output_dir`. There are no concordant/discordant region calls or regulatory element catalog.

## Usage Examples

```python
from metainformant.epigenome.workflow.workflow import (
    EpigenomeConfig,
    integrate_epigenome_results,
    run_atacseq_workflow,
    run_chipseq_workflow,
    run_methylation_workflow,
)

# Configure the analysis
config = EpigenomeConfig(
    min_methylation_coverage=10,
    chipseq_qvalue_threshold=0.01,
    min_peak_length=100,
    max_distance_for_integration=2000,
)

# Run the assay workflows on directories of input files
out = "output/epigenome/sample1"
meth_results = run_methylation_workflow("data/epigenome/methylation", f"{out}/methylation", config)
chip_results = run_chipseq_workflow("data/epigenome/chipseq", f"{out}/chipseq", config)
atac_results = run_atacseq_workflow("data/epigenome/atacseq", f"{out}/atacseq", config)
print(f"Methylation sites: {meth_results['statistics']['total_sites']}")
print(f"ChIP-seq peaks kept: {chip_results['statistics']['total_peaks']}")
print(f"ATAC-seq peaks kept: {atac_results['statistics']['total_peaks']}")

# Integrate across assays: the workflow results carry the `sites` and `peaks` lists the integration reads
integrated = integrate_epigenome_results(
    meth_results, chip_results, atac_results, f"{out}/integrated", config
)
for name, analysis in integrated["integrated_analyses"].items():
    print(name, [finding["type"] for finding in analysis["findings"]])
```

## Configuration

Environment-variable overrides are not implemented; see the `EpigenomeConfig` field table and `load_epigenome_config` above.

## Related Modules

- `metainformant.epigenome.assays.methylation` -- DNA methylation analysis
- `metainformant.epigenome.assays.chipseq` -- ChIP-seq peak analysis
- `metainformant.epigenome.assays.atacseq` -- ATAC-seq accessibility analysis
- `metainformant.epigenome.chromatin_state` -- chromatin state learning
- `metainformant.epigenome.peak_calling` -- de novo peak calling
