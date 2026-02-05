# Long-Read Sequencing Module

Analysis tools for Oxford Nanopore (ONT) and PacBio long-read sequencing data.

## Submodules

### I/O (`io/`)
- **fast5.py** - FAST5 file reading and signal extraction
- **bam.py** - Long-read BAM parsing and methylation tag extraction
- **formats.py** - Format conversion (FAST5/FASTQ, POD5/FAST5, PAF)

### Quality (`quality/`)
- **metrics.py** - N50/NX statistics, read length distributions, accuracy estimation
- **filtering.py** - Length/quality filtering, adapter trimming, chimeric read splitting

### Analysis (`analysis/`)
- **modified_bases.py** - Methylation detection (5mC, 6mA), differential methylation
- **structural.py** - SV detection from long reads (insertions, inversions)
- **phasing.py** - Haplotype phasing, block building, read tagging

### Assembly (`assembly/`)
- **overlap.py** - Overlap finding, minimizer sketching, overlap graphs
- **consensus.py** - Consensus generation, polishing, quality metrics
- **hybrid.py** - Hybrid assembly with short reads, scaffolding

### Visualization (`visualization/`)
- **plots.py** - Read length histograms, dotplots, alignment views, methylation tracks

### Workflow (`workflow/`)
- **orchestrator.py** - Pipeline orchestration (LongReadOrchestrator)
- **pipelines.py** - Pre-configured analysis pipelines
- **reporting.py** - QC report generation

## Usage

```python
from metainformant.longread.io import fast5, bam
from metainformant.longread.analysis import modified_bases, phasing
from metainformant.longread.assembly import consensus, hybrid
from metainformant.longread.workflow import orchestrator
```
