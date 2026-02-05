# Agent Directives: structural_variants

## Role
Documentation agent for the structural variants module covering CNV/SV detection and annotation.

## Module Scope
- CNV detection from read depth and segmentation
- SV calling (split-read, discordant-pair methods)
- Breakpoint refinement and microhomology detection
- Gene and regulatory overlap annotation
- Functional impact prediction and pathogenicity scoring
- TAD disruption and dosage sensitivity analysis
- Quality filtering and blacklist application
- Multi-caller merging (reciprocal overlap, SURVIVOR)
- Circos plots and SV visualization

## Key Source Files
- `src/metainformant/structural_variants/detection/` - CNV, SV calling, breakpoints
- `src/metainformant/structural_variants/annotation/` - Overlap and functional impact
- `src/metainformant/structural_variants/filtering/` - Quality filtering and merging
- `src/metainformant/structural_variants/visualization/` - Circos and coverage plots

## External Dependencies
- samtools for BAM processing
- SURVIVOR for multi-caller SV merging
- bedtools for interval operations
