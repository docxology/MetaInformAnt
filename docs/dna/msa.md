### DNA: Multiple Sequence Alignment (MSA)

Functions: `align_msa` (progressive, no external deps), `align_with_cli` (MUSCLE/Clustal if available).

```mermaid
flowchart TD
  Aid->seqMap[id->seq map] --> BalignMsa[align_msa]
  A --> CalignWithCli[align_with_cli]
  B --> DalignedId->seq[aligned id->seq]
  C --> D
```

Example

```python
from metainformant.dna import msa

aligned = msa.align_msa({"A": "ACGT", "B": "AG-T"})

# External tools
# - If MUSCLE is in PATH, align_with_cli will use it.
# - The setup script reports whether a real MUSCLE binary is on PATH; install MUSCLE v5 (Homebrew or bioconda)
#   before running align_with_cli or the external MSA tests.
# msa.align_with_cli(id_to_seq, tool="muscle")  # or "clustalo"
```

**Parameters:**

| Parameter   | Type             | Default  | Description                                          |
|-------------|------------------|----------|------------------------------------------------------|
| `sequences` | `dict[str, str]` | required | Sequence ID to DNA sequence mapping                  |
| `method`    | `str`            | `auto`   | `auto`, `muscle`, `mafft`, `clustalo`, or `clustalw` |
| `tool`      | `str`            | `muscle` | External CLI tool name for `align_with_cli`          |
