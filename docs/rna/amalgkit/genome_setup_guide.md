# Genome Setup Quick Guide

Use the single executable wrapper described in
[genome_preparation.md](genome_preparation.md). The current checkout does
not contain the older split commands `download_missing_genomes.py`,
`prepare_transcriptomes.py`, `build_kallisto_indexes.py`,
`verify_genomes_and_indexes.py`, or `orchestrate_genome_setup.py`.

## One species

```bash
CONFIG=config/amalgkit/amalgkit_camponotus_floridanus.yaml

# Inspect without changing files
uv run python scripts/rna/setup_genome.py --config "$CONFIG" --dry-run

# Check existing assets
uv run python scripts/rna/setup_genome.py --config "$CONFIG" --verify-only

# Build missing assets
uv run python scripts/rna/setup_genome.py --config "$CONFIG"
```

## All configured species

There is no current all-species genome wrapper with a safe external-volume
contract. Iterate explicitly so that each configuration, destination, and
failure is visible:

```bash
for config in config/amalgkit/amalgkit_*.yaml; do
  case "$config" in
    *template*|*test*|*cross_species*) continue ;;
  esac
  uv run python scripts/rna/setup_genome.py --config "$config" --verify-only \
    || echo "verification failed: $config" >&2
done
```

Execute only after reviewing the verification output and resource budget for
the selected data root. The verification lane is diagnostic; it does not
establish downstream quantification or cross-species completion.

## Kallisto index policy

Use the transcriptome FASTA named by the species configuration and retain the
Kallisto index beside the configured genome assets. Do not substitute a
different species index, silently rebuild an index with a different k-mer, or
claim quantification reproducibility without recording the index source and
parameters.

See [genome_preparation.md](genome_preparation.md) for the full current
contract and [the RNA validation guide](../VALIDATION.md) for downstream
checks.
