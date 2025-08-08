### DNA: Consensus

Function: `consensus_from_alignment`

```python
from metainformant.dna import consensus

cons = consensus.consensus_from_alignment({
  "A": "ACG-T",
  "B": "ACGGT",
  "C": "ACG-T",
})
```

Assumes gapped, equal-length alignment; ignores gaps when voting.


