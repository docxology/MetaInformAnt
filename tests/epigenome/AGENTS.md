# AGENTS.md — `MetaInformAnt/tests/epigenome/`

tests for the epigenome module.

## Files

Module-level tests, one file per sub-package:

- `test_epigenome.py` — BEDgraph track reading and beta-value summaries from the `tests/data/epigenome` fixtures.
- `test_epigenome_analysis.py` — `GenomicTrack`, track loading, saving and merging, statistics, region extraction and comparison.
- `test_epigenome_assays.py` — methylation sites and beta tables, `ChIPPeak` and `ATACPeak` statistics, overlap, merging and condition comparison.
- `test_epigenome_chromatin.py` — chromatin state learning, assignment, interpretation, enrichment and segmentation.
- `test_epigenome_peak_calling.py` — peak calling, merging, filtering, FRiP, differential peaks, local lambda and summit refinement.
- `test_epigenome_visualization.py` — the plotting functions, rendered with the Agg backend.
- `test_epigenome_workflow.py` — `EpigenomeConfig`, config loading, integration and the methylation workflow.

Behavior-specific tests:

- `test_epigenome_peak_io.py` — narrowPeak, broadPeak and BED round trips for ChIP-seq and ATAC-seq: column layout, summit offset, `-1` placeholders, `-log10` p- and q-values, invalid values and all-or-nothing saves.
- `test_chipseq_overlap_enrichment.py` — `find_overlapping_peaks` against a brute-force comparison, the log-space Poisson enrichment, overlapping motif occurrences and zero p- and q-values in merging.
- `test_chipseq_motif_scan.py` — IUPAC motif scanning of peak windows on both strands.
- `test_atacseq_tf_scan.py` — `find_tf_binding_sites` on planted motifs: both strands, IUPAC codes, expected sites, enrichment and argument errors.
- `test_epigenome_methylation_deep.py` — BEDgraph and Bismark `.cov` loaders and exports, DMRs, CpG islands, entropy and the report.
- `test_epigenome_methylation_fixes.py` — entropy bins, DMR p-values against an exact integer computation, and BEDgraph level scales.
- `test_peak_calling_fixes.py` — `merge_peaks` across chromosomes and merged q-values, and the `chrom` argument of the peak callers.
- `test_tracks_merge.py` — `merge_tracks` operations, including a real intersection checked against a position-by-position model.
- `test_workflow_fixes.py` — BED input, the q-value filter, the `sites` and `peaks` lists, and the ChIP-seq/ATAC-seq integration.
- `test_epigenome_determinism.py` — no unseeded randomness and no dependence on the hash seed in result order.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/epigenome` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
