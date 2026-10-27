# Genome Setup Quick Guide

This checkout has no command-line wrapper for genome preparation. The former
`scripts/rna/setup_genome.py` was retired because the package function it
wrapped no longer honours its dry-run, skip and k-mer options. Call
`metainformant.rna.amalgkit.genome_prep` directly; see
[genome_preparation.md](genome_preparation.md) for the full contract.

## One species

```python
from pathlib import Path

from metainformant.core import io
from metainformant.rna.amalgkit.genome_prep import orchestrate_genome_setup, verify_genome_status

config = io.load_yaml("config/amalgkit/amalgkit_camponotus_floridanus.yaml")
work_dir = Path(config["work_dir"])
species = "camponotus_floridanus"

# Inspect existing assets without changing files
print(verify_genome_status(work_dir, species))

# Download the configured RNA FASTA and build the Kallisto index
result = orchestrate_genome_setup(config, species, work_dir)
print(result["success"], result["messages"])
```

## All configured species

There is no all-species genome wrapper with a safe external-volume contract.
Iterate explicitly so that each configuration, destination and result is
visible, and verify before building anything:

```python
from pathlib import Path

from metainformant.core import io
from metainformant.rna.amalgkit.genome_prep import verify_genome_status

for path in sorted(Path("config/amalgkit").glob("amalgkit_*.yaml")):
    if any(tag in path.stem for tag in ("template", "test", "cross_species")):
        continue
    config = io.load_yaml(path)
    species = path.stem.removeprefix("amalgkit_")
    print(path.name, verify_genome_status(Path(config["work_dir"]), species))
```

Build only after reviewing the verification output and the resource budget
for the selected data root. Verification is diagnostic; it does not establish
downstream quantification or cross-species completion.

## Kallisto index policy

Use the transcriptome FASTA named by the species configuration and retain the
Kallisto index beside the configured genome assets. Do not substitute a
different species index, silently rebuild an index with a different k-mer, or
claim quantification reproducibility without recording the index source and
parameters.

See [genome_preparation.md](genome_preparation.md) for the full current
contract and [the RNA validation guide](../VALIDATION.md) for downstream
checks.
