# Specification: detection

## 🎯 Scope
Structural variant detection subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `breakpoints.py`
- `cnv.py`
- `sv_calling.py`

### Function contracts (behavior pinned by tests)

- `sv_calling.classify_sv_type(evidence)` makes a call only when strand and
  direction evidence support it: TRA for inter-chromosomal pairs; INS for
  split reads with nearly coincident breakpoints; INV for same-strand
  pairs with actual strand characters; DUP for everted (upstream '-',
  downstream '+') pairs after normalising to genomic order; DEL for
  head-to-head ('+', '-') pairs with split-read support or a span of at
  least `_DEL_MIN_BREAKPOINT_DISTANCE`. Weak evidence (short-span FR pairs
  without split support, unknown strand characters) returns
  `SVType.UNKNOWN`; the classifier never manufactures a deletion.
- `sv_calling.call_structural_variants(...)` keeps each breakpoint on its own chromosome. A
  translocation (`SVType.TRA`) carries them as `chrom`/`start` and `chrom2`/`pos2`, with
  `end == start` and `size == 0`; every other call occupies `chrom:start-end` and has
  `pos2 is None` and `chrom2 == ""`. The two breakpoints of an insertion nearly coincide (less than
  50 bp apart) and nothing in split-read evidence measures the inserted length, so an `INS` call
  has `size == 0` ("not measured") and `min_sv_size` is not applied to it; `min_sv_size` filters
  DEL, DUP, INV and UNKNOWN calls on their breakpoint distance. Genotyping uses
  `genotype_sv`'s default `min_mapq=20` whatever `min_mapq` the caller passes.
- `sv_calling.genotype_sv(variant, reads)` counts the reads that span `start - 100 .. start + 100`
  on `variant.chrom` and are not among the call's evidence reads. `reads` is a list of read
  dictionaries or a `sv_calling.ReadIndex` built from them once: the answers are identical, but
  the index replaces a pass over every read per call with a binary search.
- `sv_calling.StructuralVariant.evidence.evidence_reads`, like every list built from a set in this
  package, is sorted, so a run does not depend on the hash seed.
- `breakpoints.refine_breakpoints(variants, reads)` reads the second breakpoint from
  `breakpoint2`, else `pos2`, else `end`, on `chrom2` (absent, `None` or `""` means `chrom`). Reads are
  found by the reference coordinate of their soft clip, so a read that reaches the window from far
  away counts. The supporting reads of a breakpoint (those clipping at its consensus position)
  give, for both breakpoints: `strand` (the majority `is_reverse`; `'+'` without orientation or on
  a tie), `microhomology` (the longest `detect_microhomology` result over the window of each
  read's `seq` centred on its own junction) and `inserted_sequence` (the read bases between the
  primary and the single supplementary alignment of `sa_tag`, in the orientation of the
  alignment at that breakpoint; empty unless `seq`, `is_reverse`, `cigar` and a one-entry
  `sa_tag` located at the other breakpoint are all present and consistent).
- `cnv.detect_cnv_from_depth(depth_data, ...)` raises `ValueError` naming the chromosome, the
  number of negative and of non-finite values, and the likely cause (log2 ratios passed instead
  of depths) for any analysed chromosome whose depths are negative or not finite; zero depth is
  valid. Chromosomes shorter than `min_segment_bins` are skipped before this check.
- `cnv.segment_coverage(...)` scores every arc of at most half the segment, vectorised over the
  arcs of one length with the same floating-point operations in the same order as a scalar scan
  of the arcs (the first arc in `(i, j)` order wins a tie). The work still grows with the square of
  the segment length. `min_width` must be at least 1.
