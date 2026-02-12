# Protein Workflow

High-level workflow orchestration composing multiple protein analysis steps into complete pipelines for sequence analysis, structure analysis, comparative studies, and batch processing.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `orchestration` module |
| `orchestration.py` | End-to-end protein analysis pipelines |

## Key Functions

| Function | Description |
|----------|-------------|
| `analyze_protein_sequence()` | Complete sequence analysis (properties, SS prediction, motifs) |
| `analyze_protein_structure()` | Structure-level analysis pipeline |
| `batch_analyze_sequences()` | Analyze multiple protein sequences in batch |
| `comparative_analysis()` | Compare two or more protein sequences |
| `analyze_from_fasta()` | Load and analyze proteins from a FASTA file |
| `compare_structures()` | Compare two protein structures |
| `assess_alphafold_quality()` | Evaluate AlphaFold prediction confidence |
| `full_protein_analysis()` | Combined sequence and structure analysis |
| `batch_compare_structures()` | Pairwise structural comparison across a set |

## Usage

```python
from metainformant.protein.workflow import orchestration

result = orchestration.analyze_protein_sequence(sequence, name="my_protein")
batch = orchestration.batch_analyze_sequences(sequences)
```
