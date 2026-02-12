# Integration

Cross-omics utilities that connect DNA sequence analysis with RNA expression data, including ORF finding, splice site prediction, and regulatory element detection.

## Contents

| File | Purpose |
|------|---------|
| `rna.py` | DNA-RNA integration: gene structure, regulatory elements, and expression correlation |

## Key Functions

| Function | Description |
|----------|-------------|
| `find_open_reading_frames()` | Locate ORFs in a DNA sequence |
| `predict_transcription_start_sites()` | Score potential TSS locations using promoter motifs |
| `find_transcription_factor_binding_sites()` | Scan DNA for known TF motif matches |
| `calculate_codon_usage_bias()` | Compute codon bias statistics for a coding sequence |
| `analyze_gene_structure()` | Identify exons, introns, and UTR regions |
| `correlate_dna_with_rna_expression()` | Correlate sequence features with expression levels |
| `predict_splice_sites()` | Detect canonical donor/acceptor splice signals |
| `analyze_regulatory_elements()` | Find promoters, enhancers, and CpG islands |
| `integrate_dna_rna_data()` | Merge DNA and RNA datasets into a unified result |

## Usage

```python
from metainformant.dna.integration.rna import (
    find_open_reading_frames,
    correlate_dna_with_rna_expression,
    integrate_dna_rna_data,
)

orfs = find_open_reading_frames("ATGGCCATTGTAATGGGCCGCTGAAAGGGT")
merged = integrate_dna_rna_data(dna_data, rna_data)
```
