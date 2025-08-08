### DNA: FASTQ

Function: `average_phred_by_position`

```python
from pathlib import Path
from metainformant.dna import fastq

avg = fastq.average_phred_by_position(Path("/path/to/reads.fastq"))
```

Computes average Phred+33 score per position across reads; truncates to the shortest read length for safety.


