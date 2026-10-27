# Epigenome Workflow Orchestration

Directory-driven workflows for methylation, ChIP-seq, and ATAC-seq data, implemented in `metainformant.epigenome.workflow.workflow`. Provides configuration management via the `EpigenomeConfig` dataclass and a multi-assay integration step.

## Key Concepts

**Workflow orchestration** chains loading, filtering, summary statistics and report writing for every matching file in an input directory, with consistent configuration and output paths. The workflows summarize each sample; they do not call peaks or annotate them, and they do not run the differential, enrichment, motif or TSS analyses of the assay modules.

**EpigenomeConfig** collects the thresholds and options for a session. Only some of its fields are read by the workflow functions (see the field table below).

**Multi-assay integration** relates the three assays pairwise: where methylation sites fall relative to ChIP-seq and ATAC-seq peaks, and how methylation differs between accessible and closed chromatin.

## Configuration

### `EpigenomeConfig`

Dataclass for workflow configuration. Every field is an optional keyword argument; an unknown field name raises `TypeError`. `__post_init__` validates the numeric ranges noted below.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_dir` | `str \| Path` | `"data/epigenome"` | Not read by the workflow functions (they take `input_dir` as an argument) |
| `output_dir` | `str \| Path` | `"output/epigenome"` | Not read by the workflow functions (they take `output_dir` as an argument) |
| `temp_dir` | `str \| Path` | `"temp/epigenome"` | Not used |
| `methylation_threshold` | `float` | `0.8` | Validated to lie in 0--1; not otherwise used |
| `chipseq_qvalue_threshold` | `float` | `0.05` | ChIP-seq workflow: peaks with a q-value above it are dropped (peaks without a q-value are kept), and the same number is passed to `filter_peaks_by_score` as the minimum peak score. Validated to lie in 0--1 |
| `atacseq_qvalue_threshold` | `float` | `0.05` | ATAC-seq workflow: peaks with a q-value above it are dropped; a peak without a q-value is kept only if its accessibility score (`signal_value`) is at least this value. Validated to lie in 0--1 |
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

1. Load each file with `load_methylation_bedgraph` (`.bedgraph`, `.bg`) or `load_methylation_cov` (`.cov`), keeping sites with at least `min_methylation_coverage` reads. `.cov` files must have the layout `load_methylation_cov` reads (see the methylation page).
2. Compute `calculate_methylation_statistics` for each sample. A file that fails to load is recorded in `errors` and the remaining files are still processed.
3. Write `<sample>_methylation.json` for each sample and `methylation_summary.json`.
4. If `generate_reports` is true, write `methylation_report.txt` for the sites of all samples combined.

Beta values, DMRs, CpG islands and entropy are not computed here; call the functions of `metainformant.epigenome.assays.methylation` on the loaded sites.

**Returns** dict with `workflow_type` (`"methylation"`), `input_dir`, `output_dir`, `config` (the configuration as a dict), `processed_files`, `statistics` (`total_samples`, `total_sites` and `sample_statistics`, one entry with `sample`, `file` and `statistics` per sample), `errors`, and `report_path` when a report was written.

### `run_chipseq_workflow(input_dir, output_dir, config=None) -> Dict`

Summarize every ChIP-seq peak file found under `input_dir` (searched recursively for `*.narrowPeak`, `*.broadPeak`, `*.bed` and the lower-case `.narrowpeak` / `.broadpeak` spellings). The peak format is chosen from the file name: a name containing "narrow" is read as narrowPeak, one containing "broad" as broadPeak, anything else (including plain `.bed` files, which therefore yield no peaks unless they have ten columns) as narrowPeak. `PipelineError` is raised if no file is found.

1. Load each file with `load_chip_peaks`.
2. Filter: `filter_peaks_by_score` with `min_score=chipseq_qvalue_threshold`, then drop peaks whose q-value exceeds `chipseq_qvalue_threshold`, then keep peaks between `min_peak_length` and `max_peak_length` bp.
3. Compute `calculate_peak_statistics` for each sample. A file that fails is recorded in `errors`.
4. Write `<sample>_chipseq.json`, the filtered peaks as `<sample>_peaks.narrowPeak` (written by `save_chip_peaks`, so the columns after the strand are not in standard narrowPeak order), and `chipseq_summary.json`.
5. If `generate_reports` is true, write `chipseq_report.txt` for the filtered peaks of all samples.

Overlap, enrichment and motif analyses are not run here; call the functions of `metainformant.epigenome.assays.chipseq`.

**Returns** dict with `workflow_type` (`"chipseq"`), `input_dir`, `output_dir`, `config`, `processed_files`, `statistics` (`total_samples`, `total_peaks` and `sample_statistics`, one entry with `sample`, `file`, `statistics` and `peaks_count` per sample), `errors`, and `report_path` when a report was written.

### `run_atacseq_workflow(input_dir, output_dir, config=None) -> Dict`

Summarize every ATAC-seq peak file under `input_dir`. File discovery and format detection are the same as for `run_chipseq_workflow`, so keep ChIP-seq and ATAC-seq files in separate directories. `PipelineError` is raised if no file is found.

1. Load each file with `load_atac_peaks`.
2. Filter: keep peaks whose q-value is at most `atacseq_qvalue_threshold` (a peak without a q-value needs an accessibility score of at least that value), then peaks between `min_peak_length` and `max_peak_length` bp.
3. Compute `calculate_atac_statistics` for each sample, which includes the peak-length-based nucleosome fractions (`nfr_peak_fraction` and similar), not fragment sizes. A file that fails is recorded in `errors`.
4. Write `<sample>_atacseq.json`, the filtered peaks as `<sample>_peaks.narrowPeak` (written by `save_atac_peaks`, with the same column layout), and `atacseq_summary.json`.
5. If `generate_reports` is true, write `atacseq_report.txt`.

TSS enrichment and condition comparison are not run here; call `identify_tss_enrichment` and `compare_atac_conditions` from `metainformant.epigenome.assays.atacseq`.

**Returns** dict with `workflow_type` (`"atacseq"`), `input_dir`, `output_dir`, `config`, `processed_files`, `statistics` (`total_samples`, `total_peaks` and `sample_statistics`), `errors`, and `report_path` when a report was written.

### `integrate_epigenome_results(methylation_results, chipseq_results, atacseq_results, output_dir, config=None) -> Dict`

Combine the result dicts of the three workflows with pairwise association analyses. An analysis runs when both of its inputs report `statistics["total_samples"] > 0`:

- `methylation_chipseq`: classifies each methylation site as overlapping a ChIP-seq peak, proximal (within `max_distance_for_integration` bp) or distal; a site whose chromosome has no peak is counted in none of the three classes.
- `methylation_atacseq`: counts methylation sites inside, near and outside accessible regions and compares their mean methylation (`mean_methylation_accessible`, `mean_methylation_near_accessible`, `mean_methylation_closed`).
- `chipseq_atacseq`: a placeholder that returns a single finding of type `"placeholder"`.

The first two read optional lists from the input dicts: `sites` in the methylation results (dicts with `chrom` or `chromosome`, `position` or `pos`, and optionally `methylation_level` or `beta`), and `peaks` in the ChIP-seq and ATAC-seq results (dicts with `chrom` or `chromosome`, `start` and `end`; `accessible_regions` is accepted as an alternative key for ATAC-seq). The `run_*_workflow` functions do not return these lists, so with their output as-is both analyses return only a warning finding; add the lists yourself to get statistics.

**Returns** dict with `integration_type`, `output_dir`, `config`, `methylation_summary`, `chipseq_summary` and `atacseq_summary` (the `statistics` of each input), `integrated_analyses` (analysis name mapped to a dict with `analysis_type`, `max_distance`, `correlation_threshold`, `findings` and, for the first two analyses, `statistics`), `errors`, and `report_path` when `generate_reports` is true. It writes `integrated_epigenome_results.json` and, with reports enabled, `epigenome_integration_report.txt` into `output_dir`. There are no concordant/discordant region calls or regulatory element catalog.

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

# The workflow results carry no per-site or per-peak lists; add them to get association statistics
meth_results["sites"] = [{"chrom": "chr1", "position": 1100, "methylation_level": 0.1}]
chip_results["peaks"] = [{"chrom": "chr1", "start": 1000, "end": 1400}]
atac_results["peaks"] = [{"chrom": "chr1", "start": 1050, "end": 1250}]

# Integrate across assays
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
