### DNA: Entrez / NCBI

Functions: `get_genome_from_ncbi` (Entrez), CLI validation via `is_valid_assembly_accession`.

```mermaid
flowchart LR
  A[Accession/ID] --> B[Entrez efetch]
  B --> C[SeqRecord (FASTA)]
```

Example

```python
from metainformant.dna import entrez

rec = entrez.get_genome_from_ncbi("NC_001422.1", email="you@example.com")
```

Note: Entrez requires a valid email and network access.
