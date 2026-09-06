# Specification: longread

## Scope

PacBio and Oxford Nanopore long-read analysis: signal I/O, quality assessment, assembly, methylation calling, haplotype phasing, and structural variant detection.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: io, quality, analysis, assembly, methylation, phasing, workflow, visualization, utils
- **Key Concepts**: `LongReadOrchestrator`, `PipelineStep`, `PipelineResult`

## API Definition

### Exports — `workflow/`

- `LongReadOrchestrator` — DAG execution engine (`__init__(config, output_dir)`,
  `run_qc_pipeline(reads)`, `run_assembly_pipeline(reads)`,
  `run_methylation_pipeline(reads)`, `run_sv_pipeline(alignments)`,
  `run_full_pipeline(reads, alignments)`, `run_pipeline(name, input_data)`)
- `PipelineStep` — Named step with `function`, `params`, and `depends_on`
- `PipelineResult` — Aggregated step results, success flag, and timings
- `pipelines.get_qc_pipeline_config` / `get_assembly_pipeline_config` /
  `get_methylation_pipeline_config` / `get_sv_pipeline_config` —
  Configuration factories with validated parameters
- `pipelines.load_pipeline_config` — YAML/JSON config loading
- `pipelines.validate_pipeline_config` — Config validation (steps, dependencies, bounds)
- `reporting.generate_qc_report`, `reporting.export_report` — Report generation (JSON/text/HTML)

### Exports — `io/`

- `fast5.read_fast5` — HDF5 FAST5 signal/basecall reading
- `bam.read_long_read_bam`, `bam.extract_methylation_tags` — BAM reading and MM/ML tag decoding
- `formats.fast5_to_fastq`, `formats.convert_pod5_to_fast5`, `formats.write_paf` — Format conversion

### Exports — `quality/`

- `metrics.read_length_stats`, `metrics.quality_score_distribution`, `metrics.calculate_n50`,
  `metrics.estimate_accuracy` — N50/Nx, length and Phred statistics
- `filtering.filter_by_length`, `filter_by_quality`, `trim_adapters`, `split_chimeric_reads`

### Exports — `analysis/`

- `modified_bases.call_5mc`, `call_6ma`, `aggregate_methylation`, `differential_methylation`
- `structural.detect_sv_from_long_reads`, `phase_structural_variants` — Split-read SV calling
- `phasing.phase_reads`, `build_haplotype_blocks`, `tag_reads_by_haplotype` — Dataclass-based phasing

### Exports — `assembly/`

- `overlap.minimizer_sketch`, `find_overlaps`, `filter_contained_reads`
- `consensus.generate_consensus`, `calculate_consensus_quality`
- `hybrid.hybrid_assemble`, `correct_with_short_reads`

### Exports — `methylation/` and `phasing/`

- `methylation.calling.call_methylation_from_signal`, `aggregate_methylation`,
  `detect_dmrs`, `compute_methylation_stats` — Signal-level calling and aggregation
- `phasing.haplotyping.phase_reads`, `build_phase_blocks`, `haplotag_reads` — Dict-based phasing API

### Exports — `utils/`

- `batch.process_batch`, `batch_filter_reads`, `batch_compute_metrics` — Batched/parallel processing helpers
- `summary.generate_qc_summary`, `generate_assembly_summary`, `generate_methylation_summary`,
  `generate_sv_summary`, `build_run_summary`, `export_run_summary`, `compare_run_summaries`, `RunSummary`

### Exports — `visualization/`

- `plots.plot_read_length_histogram`, `plot_quality_vs_length`, `plot_dotplot`,
  `plot_alignment_view`, `plot_methylation_track`, `plot_phasing_blocks` — Matplotlib/seaborn plot generation
