# Agent Directives: longread

## Role
Documentation agent for the long-read sequencing module covering ONT and PacBio data analysis.

## Module Scope
- FAST5/POD5 file I/O and signal processing
- Read quality metrics (N50, accuracy, throughput)
- Modified base detection (5mC, 6mA methylation)
- Structural variant detection from long reads
- Haplotype phasing and block construction
- De novo assembly (overlap, consensus, hybrid)
- Long-read-specific visualization (dotplots, methylation tracks)
- Pipeline orchestration for end-to-end workflows

## Key Source Files
- `src/metainformant/longread/io/` - Data I/O (fast5.py, bam.py, formats.py)
- `src/metainformant/longread/quality/` - QC metrics and filtering
- `src/metainformant/longread/analysis/` - Modified bases, SVs, phasing
- `src/metainformant/longread/assembly/` - Overlap, consensus, hybrid assembly
- `src/metainformant/longread/visualization/` - Long-read-specific plots
- `src/metainformant/longread/workflow/` - Pipeline orchestration

## External Dependencies
- ONT tools: Guppy/Dorado for basecalling
- minimap2 for alignment
- samtools for BAM processing
