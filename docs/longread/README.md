# Long-Read Sequencing

## Overview

PacBio and Oxford Nanopore long-read sequencing analysis module for METAINFORMANT. Covers signal I/O, quality assessment, assembly, methylation calling, haplotype phasing, and structural variant detection.

## Contents

- **io/** - FAST5, POD5, BAM reading/writing and format conversion
- **quality/** - Read metrics (N50, Nx, accuracy), filtering by length/quality
- **analysis/** - Modified base detection, structural variants, phasing
- **assembly/** - Minimizer overlaps, POA consensus, hybrid assembly
- **methylation/** - Signal-level methylation calling and DMR detection
- **phasing/** - Read-based haplotyping and switch error computation
- **workflow/** - Pipeline orchestration with dependency resolution
- **visualization/** - Read length, quality, alignment, and methylation plots
- **utils/** - Batch processing and run summary generation

## Architecture

```mermaid
graph TD
    subgraph "Long-Read Module"
        IO[io/] --> |fast5.py, bam.py| F[FAST5/POD5/BAM I/O]
        IO --> |formats.py| FMT[Format Conversion]

        Q[quality/] --> |metrics.py| QM[N50, Accuracy, QScore]
        Q --> |filtering.py| QF[Length/Quality Filtering]

        AN[analysis/] --> |modified_bases.py| MB[5mC/6mA Detection]
        AN --> |structural.py| SV[SV Calling]
        AN --> |phasing.py| PH[Haplotype Phasing]

        AS[assembly/] --> |overlap.py| OV[Minimizer Overlaps]
        AS --> |consensus.py| CO[POA Consensus]
        AS --> |hybrid.py| HY[Hybrid Assembly]

        ME[methylation/] --> |calling.py| MC[Signal-Level Calling]

        PHA[phasing/] --> |haplotyping.py| HT[Phase Block Construction]

        W[workflow/] --> |orchestrator.py| OR[Pipeline Orchestration]
        W --> |pipelines.py| PI[Pre-defined Pipelines]
        W --> |reporting.py| RE[QC Reports]

        V[visualization/] --> |plots.py| VP[Read/Alignment Plots]

        U[utils/] --> |batch.py, summary.py| UT[Batch Processing]
    end

    IO --> Q
    Q --> AN
    AN --> W
    AS --> W
    ME --> W
    V --> W
```

## Usage

```python
from metainformant.longread.workflow.orchestrator import LongReadOrchestrator
from metainformant.longread.io import fast5, bam
from metainformant.longread.quality import metrics, filtering
from metainformant.longread.analysis import modified_bases, structural
from metainformant.longread.assembly import overlap, consensus, hybrid
```
