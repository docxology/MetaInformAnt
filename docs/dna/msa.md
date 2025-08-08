### DNA: Multiple Sequence Alignment (MSA)

Functions: `align_msa` (progressive, no external deps), `align_with_cli` (MUSCLE/Clustal if available).

```mermaid
flowchart TD
  A[id->seq map] --> B[align_msa]
  A --> C[align_with_cli]
  B --> D[aligned id->seq]
  C --> D
```

Example

```python
from metainformant.dna import msa

aligned = msa.align_msa({"A": "ACGT", "B": "AG-T"})
```


