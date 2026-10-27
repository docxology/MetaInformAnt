# Detection

Structural variant detection from sequencing data, providing SV calling from split/discordant reads, copy number variation detection via circular binary segmentation, and breakpoint refinement to base-pair resolution.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports cnv, sv_calling, breakpoints submodules |
| `sv_calling.py` | SV calling from split-read and discordant read-pair evidence |
| `cnv.py` | CNV detection via circular binary segmentation (CBS) of read depth |
| `breakpoints.py` | Breakpoint refinement, microhomology detection, clustering |

## Key Functions

| Function | Description |
|----------|-------------|
| `sv_calling.call_structural_variants()` | Call SVs (DEL, DUP, INV, TRA, INS) from aligned reads; a TRA call carries its second breakpoint in `chrom2`/`pos2`, an INS call has `size` 0 |
| `sv_calling.detect_split_reads()` | Detect split-read evidence for SVs |
| `sv_calling.detect_discordant_pairs()` | Detect discordant read-pair evidence |
| `sv_calling.classify_sv_type()` | Classify SV type from evidence signatures |
| `sv_calling.genotype_sv()` | Genotype an SV using read support; accepts a read list or a `ReadIndex` |
| `sv_calling.ReadIndex` | Per-chromosome index of reads by start position, built once and shared by all genotyping calls |
| `cnv.detect_cnv_from_depth()` | Detect CNVs from read depth profiles (finite, non-negative depths) |
| `cnv.segment_coverage()` | Segment coverage data using CBS algorithm |
| `cnv.call_cnv_states()` | Assign copy number states to segments |
| `cnv.calculate_log2_ratio()` | Compute GC-corrected log2 ratios |
| `breakpoints.refine_breakpoints()` | Refine breakpoint positions to base-pair resolution; reports strand, microhomology and inserted sequence per breakpoint |
| `breakpoints.detect_microhomology()` | Detect microhomology at breakpoint junctions |
| `breakpoints.cluster_breakpoints()` | Cluster nearby breakpoints from multiple samples |

## Usage

```python
from dataclasses import asdict

from metainformant.structural_variants.detection import sv_calling, cnv, breakpoints

svs = sv_calling.call_structural_variants(alignments, min_support=3)
cnv_results = cnv.detect_cnv_from_depth(depth_data, window_size=1000)  # {chrom: CNVResult}
refined = breakpoints.refine_breakpoints([asdict(sv) for sv in svs], alignments)
```
