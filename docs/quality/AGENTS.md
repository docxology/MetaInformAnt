# Agent Directives: docs/quality

## Role

Documentation for the data quality assessment module covering FASTQ quality analysis (Phred scores, adapter contamination, overrepresented sequences), contamination detection (k-mer based, cross-species), assembly validation (QUAST-style metrics, misassembly detection), and comprehensive quality reporting (MultiQC integration, HTML reports).

## Module Scope

Quality control workflows:
- **Raw read QC**: FastQC-style per-base quality, per-sequence quality, adapter content, duplication levels, overrepresented k-mers
- **Contamination detection**: Taxonomy-based k-mer screening (Kraken2-like), cross-species contamination, vector/adapter masking
- **Assembly QC**: N50, L50, misassembly detection (relocation, inversion), BUSCO completeness, k-mer spectra
- **BAM/CRAM QC**: Mapping rate, insert size distribution, coverage depth, duplicate rate, index metrics
- **Reporting**: MultiQC aggregation, HTML report generation, JSON metrics, threshold-based pass/fail

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/quality/analysis/` | Quality metric calculation (Phred averaging, adapter detection, k-mer counting, contamination screening) |
| `src/metainformant/quality/io/` | FASTQ/BAM parsing, sequence streaming, low-memory k-mer sketches |
| `src/metainformant/quality/reporting/` | HTML report generation, MultiQC.yaml emitter, JSON serialization |
| `src/metainformant/quality/contamination/` | Cross-species contamination detection (using k-mer indexing against reference panel) |
| `scripts/quality/` | CLI tools: `fastq_qc.py`, `bam_qc.py`, `assembly_qc.py`, `contam_detect.py` |

## Cross-Module Dependencies

### Upstream
- **Core I/O**: `metainformant.core.io` — FASTQ/BAM reading, file handling
- **Core config**: `metainformant.core.utils.config` — QC thresholds (min_quality, min_length, adapter patterns)

### Downstream (consumers of QC metrics)
- **RNA**: `metainformant.rna.engine` — Reads pass QC before alignment; low-quality samples flagged
- **DNA**: `metainformant.dna.variants` — BAM QC required before variant calling
- **GWAS**: `metainformant.gwas.qc` — Sample/genotype missingness filters based on QC
- **Cloud**: `metainformant.cloud` — QC checkpoints in deployment status

## Rule Constraints
1. **Streaming required** — All QC operations must support `stream=True` for >100GB FASTQs.
2. **Thresholds configurable** — QC pass/fail cutoffs in config, not hard-coded.
3. **Report localization** — HTML reports should support i18n strings; default English.
4. **Interoperability** — Output MultiQC-compatible YAML so external MultiQC run can aggregate across tools.

## Related Tasks & Guides
- **Task reference**: N/A (QC integrated into other task pages)
- **Module guide**: `docs/quality/index.md` — Overview and quick start
- **FASTQ QC**: `docs/quality/fastq.md` — Per-base quality, adapter trimming, k-mer analysis
- **Contamination**: `docs/quality/contamination.md` — Cross-species screening methods
- **Assembly**: `docs/quality/assembly.md` — QUAST metrics, misassembly detection
- **BAM QC**: `docs/quality/bam.md` — Mapping rate, coverage, duplicate stats
- **Plot integration**: `docs/visualization/plots.md` — Quality control figures (box, per-base quality, GC)

## Maintenance Notes
- **FastQC compatibility**: QC metrics broadly compatible with FastQC 0.11+ outputs for MultiQC integration.
- **Reference genomes**: Contamination detection requires k-mer index of common contaminants (phiX, adapter, host genomes). Indices stored in `data/contamination/`.
- **Parallelization**: QC is embarrassingly parallel across samples; use `parallel_map()` from core.
- **Benchmarks**: Run `tests/benchmark/test_qc_speed.py` on 30GB FASTQ to ensure <10 min wall time.

## AI Assistant Guidance
When adding QC metric:
1. Implement in `src/metainformant/quality/analysis/<metric>.py`
2. Add to `src/metainformant/quality/reporting/report.py` aggregation
3. Write test in `tests/test_quality_<metric>.py` with small FASTQ fixture
4. Add MultiQC section in `src/metainformant/quality/reporting/multiqc.py`
5. Document in `docs/quality/<metric>.md` with before/after example
