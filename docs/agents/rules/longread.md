# Long-Read Sequencing Module Rules

## Purpose
Analysis tools for Oxford Nanopore (ONT) and PacBio long-read sequencing data, including signal processing, quality control, modified base detection, structural variants, haplotype phasing, and de novo assembly.

## Source Structure
```
src/metainformant/longread/
├── io/
│   ├── fast5.py          # FAST5 file reading, signal extraction
│   ├── bam.py            # Long-read BAM parsing, methylation tags
│   └── formats.py        # Format conversion (FAST5/FASTQ, POD5/FAST5, PAF)
├── quality/
│   ├── metrics.py        # N50/NX, read length stats, accuracy, throughput
│   └── filtering.py      # Length/quality filtering, adapter trimming, chimeric splitting
├── analysis/
│   ├── modified_bases.py # 5mC/6mA methylation detection, differential methylation
│   ├── structural.py     # SV detection from long reads (insertions, inversions)
│   └── phasing.py        # Haplotype phasing, block building, read tagging
├── assembly/
│   ├── overlap.py        # Overlap finding, minimizer sketching, overlap graphs
│   ├── consensus.py      # Consensus generation, polishing, quality calculation
│   └── hybrid.py         # Hybrid assembly with short reads, scaffolding
├── visualization/
│   └── plots.py          # Read length histograms, dotplots, methylation tracks
├── utils/
│   ├── batch.py          # Batch processing, filtering, metrics
│   └── summary.py        # QC/assembly/methylation summaries, run comparisons
└── workflow/
    ├── orchestrator.py   # Pipeline orchestration (LongReadOrchestrator)
    ├── pipelines.py      # Pre-configured pipelines
    └── reporting.py      # QC report generation
```

## Dependencies
- **Required**: numpy, pysam
- **Optional**: h5py (FAST5), pod5 (POD5), minimap2 (alignment)

## Import Patterns
```python
from metainformant.longread.io import fast5, bam, formats
from metainformant.longread.quality import metrics, filtering
from metainformant.longread.analysis import modified_bases, phasing, structural
from metainformant.longread.assembly import overlap, consensus, hybrid
from metainformant.longread.workflow import orchestrator
```

## Configuration
- Environment prefix: `LR_` (e.g., `LR_THREADS`, `LR_WORK_DIR`)
- Output path: `output/longread/<analysis_type>/` (e.g., `basecalling/`, `assembly/`, `methylation/`)

## Integration
- **Longread → DNA**: Genomic coordinates, variant calling
- **Longread → Epigenome**: Methylation from modified base detection
- **Longread → Structural Variants**: SV detection complements short-read methods

## Testing
- Use `@pytest.mark.external_tool` for tests requiring minimap2, samtools
- Generate test FAST5/BAM data programmatically, never mock file I/O
- All test outputs to `tmp_path`
