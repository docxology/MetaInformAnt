# RNA Splicing

Alternative splicing detection, quantification, and isoform analysis for RNA-seq data. Includes splice junction detection from CIGAR strings, event classification, PSI computation, differential splicing testing, splice site scoring, isoform quantification via EM, and splice graph construction.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `detection`, `isoforms`, `splice_analysis`, `splice_sites` |
| `detection.py` | Backward-compatible shim re-exporting splice_sites and splice_analysis |
| `splice_sites.py` | Junction detection from aligned reads and splice site strength scoring |
| `splice_analysis.py` | Event classification, PSI computation, differential splicing |
| `isoforms.py` | EM-based isoform quantification, splice graphs, diversity analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `splice_sites.detect_splice_junctions()` | Identify junctions from aligned reads via CIGAR parsing |
| `splice_sites.compute_splice_site_strength()` | Score splice sites using position weight matrices |
| `splice_analysis.classify_splicing_events()` | Classify junctions by splicing event type |
| `splice_analysis.compute_psi()` | Compute Percent Spliced In with confidence interval |
| `splice_analysis.differential_splicing()` | Test for differential splicing between conditions |
| `splice_analysis.find_novel_junctions()` | Discover unannotated splice junctions |
| `isoforms.quantify_isoforms()` | EM algorithm for transcript abundance estimation |
| `isoforms.build_isoform_graph()` | Construct splice graph from exons and junctions |
| `isoforms.enumerate_isoforms()` | Enumerate possible isoforms by graph traversal |
| `isoforms.compute_isoform_diversity()` | Shannon entropy and effective isoform count |
| `isoforms.compare_isoform_usage()` | Compare isoform usage between conditions |

## Usage

```python
from metainformant.rna.splicing import splice_sites, splice_analysis, isoforms

junctions = splice_sites.detect_splice_junctions(alignments, min_reads=3)
psi = splice_analysis.compute_psi(inclusion_reads=30, exclusion_reads=10)
abundance = isoforms.quantify_isoforms(read_assignments, isoform_models)
```
