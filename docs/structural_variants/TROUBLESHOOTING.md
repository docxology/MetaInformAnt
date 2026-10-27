# Troubleshooting: Structural Variants

The module works on in-memory Python records. It starts no external callers and reads no BAM
files (only `survivor_merge` reads VCF paths), so failures come from the records you pass in,
from the thresholds, or from a missing optional package (numpy, matplotlib, pysam).

## Errors

| Message | Cause | Fix |
|---------|-------|-----|
| `ValueError: depth_data['chr1'] has N negative and M non-finite depth value(s) ...` | `detect_cnv_from_depth` received log2 ratios, or NaN or infinite values | Pass read depths (finite, `>= 0`; 0 is an uncovered bin). For a tumor/normal comparison pass `2**ratio`; segment log2 ratios with `segment_coverage` and `call_cnv_states` instead |
| `ValueError: population_db layout is not supported: none of its N entries could be used ...` | `filter_by_frequency` could read no entry of `population_db` (for example a dict keyed by chromosome) | Use a dict `{"chrom:start-end[-svtype]": af}` or a list of dicts with `chrom`, `start`, `end`, `af` (or `allele_frequency`) and optionally `sv_type`, with `af` in [0, 1] |
| `Frequency filter: skipped N of M population_db entries ...` (warning) | Some entries are in an unsupported layout or have an `af` outside [0, 1] | Fix the entries named in the message; the readable ones are used |
| `ValueError: min_width must be at least 1` | `segment_coverage(min_width=0)` or less | Pass `min_width >= 1` |
| `ValueError: Unsupported method ...` | `detect_cnv_from_depth(method=...)` is not `"segmentation"` | The only method is `"segmentation"` |
| `ValueError: Unknown genotyping method` | `genotype_sv_population(method=...)` is neither `"depth"` nor `"split"` | Use one of the two |
| `ValueError: Length mismatch: N genotypes vs M phenotypes` | `sv_association_test` got lists of different lengths | Pass one genotype and one phenotype per sample, in the same order. The genotypes are those of a single SV: `gt_matrix[i]` |
| `ImportError: numpy is required ...` / `matplotlib is required ...` | Optional package missing | `uv pip install numpy` or `uv pip install matplotlib` |
| `ValueError: invalid literal for int()` or `could not convert string to float` at the start of a call | An `SV_` environment variable does not parse as its type (an int for `SV_MIN_SUPPORT`, a float for the others) | Fix or unset it (see [CONFIGURATION.md](CONFIGURATION.md)) |

## Unexpected results

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `call_structural_variants` returns no calls | Records without a `mapq` key are dropped (a missing `mapq` counts as 0, below `min_mapq`); reads carry no soft clip of at least 20 bp and no discordant mate information; support is below `min_support`; `SV_MIN_SUPPORT` is set in the environment and overrides the argument | Check the keys of your records (`name`, `chrom`, `pos`, `mapq`, `cigar`, and `mate_chrom`, `mate_pos`, `insert_size`, `is_reverse`, `mate_is_reverse` for discordant pairs, `sa_tag` for split reads); run `detect_split_reads` and `detect_discordant_pairs` on the records to count the evidence; `env \| grep '^SV_'` |
| `WARNING: Could not estimate insert size, using defaults (mean=400, std=100)` | No record has an `insert_size` between 1 and 9,999 bp on the same chromosome as its mate | Add `insert_size` to the records, or pass `insert_size_stats=InsertSizeStats(mean=..., std=...)` |
| A call is missing although it has support | A DEL, DUP, INV or UNKNOWN call whose breakpoint distance is below `min_sv_size` (default 50) is dropped | Lower `min_sv_size`; translocations and insertions are not subject to it |
| An `INS` call has `size` 0 | Split-read evidence does not measure an inserted length, so `size` 0 means "not measured" | Use `Breakpoint.inserted_sequence` from `refine_breakpoints` when the reads allow it; `filter_by_size` keeps these calls |
| A `TRA` call has `start == end` | A translocation joins two loci on different chromosomes: `chrom`/`start` are the first breakpoint and `chrom2`/`pos2` the second | Read `chrom2` and `pos2` for the partner; `refine_breakpoints` and `plot_circos` use them |
| Annotation, merging or frequency filtering ignores a translocation's partner locus | These functions read `chrom`, `start` and `end` only. `merge_callsets` and `deduplicate_variants` give a zero-length call a reciprocal overlap of 0, so at their default thresholds they never join or collapse translocations; `survivor_merge` joins them by their first breakpoint alone | Compare translocations yourself on `(chrom, start, chrom2, pos2)` |
| With `min_mapq` below 20, reads with MAPQ from `min_mapq` to 19 give evidence but are not counted as reference reads | `call_structural_variants` genotypes with `genotype_sv`'s default `min_mapq=20` | Use `min_mapq` of 20 or more, or call `genotype_sv(variant, reads, min_mapq=...)` yourself |
| `refine_breakpoints` returns `support` 0 and the input position | No soft clip of at least `min_clip` bp lies within `window` of that breakpoint | Check that the reads are on the same chromosome (`chrom`, and `chrom2` for the second breakpoint), that their `cigar` contains `S`, and that `window` and `min_clip` fit; for a translocation give `chrom2` and `pos2` |
| `inserted_sequence` is empty | It is filled only when a supporting read has `seq`, `is_reverse`, `cigar` and an `sa_tag` with exactly one entry located at the other breakpoint, and both CIGARs add up to the length of `seq` | Add the missing keys; hard-clipped records (their `seq` is shorter than the CIGAR) give no sequence |
| `microhomology` is empty or short | No read carries `seq`, or the bases around the clip junction do not repeat. The result is a mirror-match heuristic on the read, not true microhomology | Supply `seq` for the supporting reads; use the reference sequence for a real microhomology call |
| `strand` is always `'+'` | No supporting read reports `is_reverse`, or the reads are evenly split | Add `is_reverse` to the read records |
| CNV detection is slow | CBS grows with the square of the bins per chromosome (about a minute at 256,000 bins) | Use a coarser `window_size` (see [PERFORMANCE.md](PERFORMANCE.md)) |
| `variance_explained` does not add up to 1 | Each value is a share of the total variance of all components, so only the full set of components sums to 1 | Request `n_components >= min(samples, SVs)` if you need the full set |
| No INFO log lines | Logging is not configured; only warnings are printed by default | Call `configure_logging_from_env()` and set `CORE_LOG_LEVEL` |

## Diagnostic checklist

1. Check which `SV_` variables are set; they override function arguments:
   ```bash
   env | grep '^SV_'
   ```

2. Check the keys your alignment records carry and how many evidence items they give:
   ```python
   from metainformant.structural_variants.detection.sv_calling import (
       InsertSizeStats,
       detect_discordant_pairs,
       detect_split_reads,
   )

   print(sorted(set().union(*(record.keys() for record in alignments))))
   print(len(detect_split_reads(alignments, min_clip=20)))
   print(len(detect_discordant_pairs(alignments, InsertSizeStats(mean=400.0, std=100.0))))
   ```

3. Run the pipeline on one chromosome or region first to see how time and memory behave.

## Getting help

Provide: OS, Python version, metainformant version, the call with its arguments, the full
traceback or log lines, and a small set of records (a few dozen alignment dicts or a short depth
array) that reproduces the problem.
