# Contamination Detection

Algorithms for detecting various types of contamination in biological
sequencing data, including microbial, cross-species, adapter, vector,
duplication, RNA, mycoplasma, and rRNA contamination.

## Key Concepts

**Contamination detection** identifies unwanted sequences in sequencing data
that may arise from laboratory handling, library preparation, or biological
sources. Early detection prevents downstream analysis artefacts.

**k-mer based detection** underpins microbial and cross-species checks. The
detector builds k-mer sets (default k=20 for microbial, k=25 for cross-species)
from known contaminant sequences and computes Jaccard similarity against input
reads.

**Severity classification** aggregates results from all analyses into an overall
severity score:

| Severity | Score Range | Recommendation |
|----------|------------|----------------|
| High     | >= 50      | Re-sequence with improved protocols |
| Moderate | 20 -- 49   | Review protocols, filter reads |
| Low      | 5 -- 19    | Monitor in downstream analysis |
| None     | < 5        | Data quality appears good |

## Class Reference

### ContaminationDetector

```python
class ContaminationDetector:
    def __init__(self, reference_genomes: Optional[Dict[str, str]] = None)
```

Stateful detector that holds reference genomes and known contaminant sequences.
Ships with built-in signatures for E. coli, yeast, and human.

**Methods:**

```python
def detect_microbial_contamination(
    self, sequences: List[str], threshold: float = 0.01,
) -> Dict[str, Any]
```

k-mer Jaccard similarity against known contaminant signatures. Returns
`detected`, `contaminants` (list of matches with contamination rate and average
similarity), and `total_sequences_analyzed`.

```python
def detect_cross_species_contamination(
    self, sequences: List[str], target_species: str,
    other_species: List[str],
) -> Dict[str, Any]
```

Compares k-mer overlap with target vs potential contaminant reference genomes.
A read is flagged when its similarity to the contaminant exceeds 1.5x its
similarity to the target.

```python
def detect_adapter_contamination(
    self, sequences: List[str],
    adapters: Optional[List[str]] = None,
) -> Dict[str, Any]
```

Scans for exact and partial adapter matches. Defaults to Illumina TruSeq and
Nextera adapters. Reports per-adapter contamination rates and match positions.

```python
def detect_duplication_contamination(
    self, sequences: List[str], max_duplicates: int = 10,
) -> Dict[str, Any]
```

Counts exact sequence duplicates. Reports duplication rate, unique count, and
the most duplicated sequences with percentages.

```python
def comprehensive_contamination_analysis(
    self, sequences: List[str],
    target_species: Optional[str] = None,
) -> Dict[str, Any]
```

Runs all detection methods and produces a unified report with per-analysis
results, an overall severity score, and severity classification.

## Standalone Functions

### detect_rna_contamination

```python
def detect_rna_contamination(dna_sequences: List[str]) -> Dict[str, Any]
```

Checks for uracil (U) bases in DNA sequencing data. Threshold: 0.1%.

### detect_vector_contamination

```python
def detect_vector_contamination(
    sequences: List[str],
    vector_sequences: Optional[List[str]] = None,
) -> Dict[str, Any]
```

Scans for restriction enzyme sites (BamHI, EcoRI, XhoI, SalI, SmaI) and
custom vector sequences. Threshold: 1%.

### detect_mycoplasma_contamination

```python
def detect_mycoplasma_contamination(sequences: List[str]) -> Dict[str, Any]
```

Screens for mycoplasma-specific signatures using the microbial detection
pipeline.

### detect_rrna_contamination

```python
def detect_rrna_contamination(sequences: List[str]) -> Dict[str, Any]
```

Screens for conserved bacterial 16S, archaeal 16S, and eukaryotic 18S rRNA
sequences.

### generate_contamination_report

```python
def generate_contamination_report(
    contamination_results: Dict[str, Any],
    output_path: Optional[str | Path] = None,
) -> str
```

Generates a formatted text report with per-analysis breakdown, severity
classification, and actionable recommendations. Optionally saves to file.

## Usage Example

```python
from metainformant.quality.analysis.contamination import (
    ContaminationDetector,
    detect_rna_contamination,
    generate_contamination_report,
)

# Comprehensive analysis
detector = ContaminationDetector()
results = detector.comprehensive_contamination_analysis(sequences)
report = generate_contamination_report(results, "output/quality/contamination.txt")

# Standalone RNA check
rna_result = detect_rna_contamination(dna_sequences)
if rna_result["detected"]:
    print(f"RNA contamination: {rna_result['contamination_rate']:.2%}")
```

## Related Modules

- `metainformant.quality.analysis.metrics` -- quality scoring and outlier detection
- `metainformant.quality.io.fastq` -- FASTQ file parsing and QC
- `metainformant.quality.reporting` -- QC report generation
