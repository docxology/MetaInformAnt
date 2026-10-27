# Agent Directives: structural_variants

## Role
Documentation agent for the structural variants module covering CNV/SV detection and annotation.

## Module Scope
- CNV detection from read depth and segmentation
- SV calling (split-read, discordant-pair methods), with translocations carried as `chrom`/`start` plus `chrom2`/`pos2`
- Breakpoint refinement: strand, microhomology (a mirror-match heuristic) and inserted sequence per breakpoint
- Gene and regulatory overlap annotation
- Functional impact prediction and pathogenicity scoring
- TAD disruption and dosage sensitivity analysis (HI, TS, pLI, LOEUF)
- Quality filtering and blacklist application
- Multi-caller merging (reciprocal overlap, SURVIVOR-style breakpoint distance)
- Population genotyping, allele frequency, association, PCA, LD
- Circos plots and SV visualization

## Key Source Files
- `src/metainformant/structural_variants/detection/` - CNV, SV calling, breakpoints
- `src/metainformant/structural_variants/annotation/` - Overlap and functional impact
- `src/metainformant/structural_variants/filtering/` - Quality filtering and merging
- `src/metainformant/structural_variants/population/` - Cohort genotyping and statistics
- `src/metainformant/structural_variants/visualization/` - Circos and coverage plots

## Documentation Rules
- State only what the code does: check each claim against `src/metainformant/structural_variants/` and the tests in `tests/structural_variants/`.
- Examples bind to the real signatures; `uv run python scripts/verify_documentation_code.py --docs-dir docs/structural_variants` checks imports, attributes, calls and links.
- Timings carry their conditions (input size, CPU time, machine) or are replaced by a statement of complexity.

## External Dependencies
- NumPy (CNV detection, population functions) and matplotlib (visualization) are imported by the module; SciPy is imported but not called
- pysam is optional: `survivor_merge` uses it to read VCF paths (a basic text parser otherwise); BAM reading is left to the caller
- The module calls no external program (no samtools, SURVIVOR or bedtools); the merge and interval algorithms are implemented in Python
