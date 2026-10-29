# AGENTS.md — `MetaInformAnt/tests/longread/`

tests for the longread module.
Files (verified 2026-10-28): __init__.py, test_assembly_overlaps.py, test_longread.py, test_longread_adapter_matching.py, test_longread_analysis.py, test_longread_assembly.py, test_longread_bam_tags.py, test_longread_bayesian_methylation.py, test_longread_determinism.py, test_longread_fast5_mux.py, test_longread_io.py, test_longread_methylation.py, test_longread_methylation_phasing_quality.py, test_longread_mm_ml_spec.py, test_longread_numpy_lengths.py, test_longread_phasing_json_blocks.py, test_longread_phasing_pipeline.py, test_longread_placeholder_values.py, test_longread_quality.py, test_longread_structural_merge.py, test_longread_true_median.py, test_longread_visualization.py.

## Coverage by file
- `test_assembly_overlaps.py`: reverse-complement (antiparallel) overlaps in `find_overlaps` and the frames of minimizer chaining.
- `test_longread.py`: the public API of every subpackage on synthetic data: quality metrics and filtering, methylation, structural variants, phasing, assembly (overlap, consensus, hybrid), I/O and visualization.
- `test_longread_adapter_matching.py`: the default adapter table (no entry shorter than `minimum_adapter_length()`; the renamed `Illumina_Nextera_read1` and the deprecated `PacBio_barcoded` alias), nested adapters (a longer adapter cannot take an exact shorter one's site by spending edits), adapter-free reads left untouched, adapters at any offset and in both orientations, the substitution, insertion and deletion edit budget, `adapter_false_positive_bound` against an exhaustive count and against random sequence, exact-only search with a warning for short adapters, and argument validation.
- `test_longread_analysis.py`: `BatchResult`, `process_batch`, `batch_filter_reads` and `batch_compute_metrics`.
- `test_longread_assembly.py`: the summaries in `utils.summary`: `generate_assembly_summary`, `generate_sv_summary`, `RunSummary`, `build_run_summary`, `export_run_summary` and `compare_run_summaries`.
- `test_longread_bam_tags.py`: MM/ML decoding for every `ML` form a BAM reader returns (`array('B')`, numpy arrays, lists, tuples, bytes), the strand in the `extract_methylation_tags` key (`C+m` and `C-m` stay separate), and the workflow's methylation features built from array tags.
- `test_longread_mm_ml_spec.py`: MM/ML decoding against the SAM tags specification and the htslib reference: read orientation (FLAG 0x10), the `-` strand, multi-code entries, the `.`/`?`/no flag, malformed and truncated tags, positions without a reference coordinate, a randomized check against an independent decoder, and the workflow's methylation features.
- `test_longread_bayesian_methylation.py`: the `bayesian` model of `call_methylation_from_signal`: the single-read closed form and the prior shift, the posterior against exhaustive enumeration, pooling by `(chrom, position, strand)`, order independence, NaN and extreme LLRs, and prior validation.
- `test_longread_determinism.py`: an AST guard against unseeded random number generation in the package, and input-order keys in parallel `process_batch` results.
- `test_longread_fast5_mux.py`: `mux` read from the `start_mux` attribute, for multi-read and legacy single-read FAST5 files written with h5py.
- `test_longread_io.py`: error paths of `read_long_read_bam`, `convert_pod5_to_fast5` and `fast5_to_fastq`, and the mean Phred quality helper in `formats`.
- `test_longread_methylation.py`: `generate_methylation_summary`.
- `test_longread_methylation_phasing_quality.py`: the public functions of `methylation.calling`, `phasing.haplotyping` and `quality.filtering` (length and quality filters, `trim_adapters`, `detect_adapters`), and `workflow.pipelines.load_pipeline_config`.
- `test_longread_numpy_lengths.py`: numpy integer and whole-number float lengths (fractional, NaN and bool lengths raise) in `read_length_stats`, `calculate_throughput`, `calculate_n50`, `batch_compute_metrics` and the read-length histogram, and `generate_qc_summary` / `generate_assembly_summary` on numpy arrays.
- `test_longread_phasing_json_blocks.py`: phase blocks reloaded from JSON (text `variant_haplotypes` keys) tagging the same reads, key and allele validation, and read alleles that are not 0 or 1.
- `test_longread_phasing_pipeline.py`: `phase_reads` blocks passed straight to `haplotag_reads` (`variant_haplotypes`, best-supported block, inclusive block end, the warning for blocks without the map), the `build_phase_blocks` map from `haplotype` and `allele`, `analysis.phasing.phase_reads` from `aligned_pairs` and with `min_variant_quality`, and `tag_reads_by_haplotype` block selection.
- `test_longread_placeholder_values.py`: placeholder values and silent defaults removed from methylation calling, phase blocks, consensus, SV genotypes, the QC report, `process_batch` and the differential step of the methylation pipeline.
- `test_longread_quality.py`: `generate_qc_summary`, `QCReport`, `generate_qc_report`, `generate_run_summary`, short reads passing through `split_chimeric_reads`, and the text export of `RunSummary`.
- `test_longread_structural_merge.py`: merged structural variant calls: breakend mates (`mate_chromosome`, median `mate_position`, chained mate positions), `inserted_sequence` (soft and hard clips, the representative of a merged insertion), and output that does not depend on the hash seed (`read_names`, `info["sources"]`, `compare_run_summaries` platforms).
- `test_longread_true_median.py`: true medians in merged SV calls, `generate_sv_summary`, the SV pipeline summary and the median line of the read-length histogram.
- `test_longread_visualization.py`: the workflow layer (`PipelineStep`, `PipelineResult`, pipeline configurations and their validation, `LongReadOrchestrator`) and plots written for empty input; the plot functions themselves are tested in `test_longread.py`.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/longread`
Repo-wide policy: see the repository-root `AGENTS.md`.
