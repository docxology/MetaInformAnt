### DNA: Transcription

Functions: `transcribe_dna_to_rna`, `reverse_transcribe_rna_to_dna`

```mermaid
flowchart LR
  A[DNA] --> BtranscribeDnaToRna[transcribe_dna_to_rna] --> C[RNA]
  C --> DreverseTranscribeRnaToDna[reverse_transcribe_rna_to_dna] --> A
```

Example

```python
from metainformant.dna import transcription

# DNA to RNA transcription (T -> U)
dna = "ATCGATCG"
rna = transcription.transcribe_dna_to_rna(dna)  # "AUCGAUCG"

# RNA to DNA reverse transcription (U -> T)  
rna = "AUCGAUCG"
dna = transcription.reverse_transcribe_rna_to_dna(rna)  # "ATCGATCG"
```

Features:
- **Bidirectional conversion**: DNA ↔ RNA transcription
- **Case preservation**: Maintains original case of input sequence
- **Non-standard bases**: Passes through non-ATGCU characters unchanged
- **Complete sequences**: Processes entire sequence without position restrictions

Transcription rules:
- **DNA to RNA**: T → U (thymine to uracil)
- **RNA to DNA**: U → T (uracil to thymine)  
- **Preserved**: A, G, C remain unchanged
- **Case sensitivity**: Lowercase and uppercase handled separately

Related: Often used with [translation](./translation.md) for gene expression analysis.
