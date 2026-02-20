# Resuming Amalgkit on Linux

> **Snapshot**: 2026-02-19 10:36 PST — Mac processes killed, all results synced and committed.

## Prerequisites

- Linux machine with Python 3.11+ and `uv` installed
- The `blue/` external drive connected and mounted
- Git clone of MetaInformAnt

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/docxology/MetaInformAnt.git
cd MetaInformAnt

# 2. Update the blue/ symlink to point to your mount point
rm blue
ln -s /path/to/your/blue/mount/data blue
# The drive layout is: blue/ -> <mount>/data/amalgkit/<species>/...

# 3. Install dependencies
uv sync

# 4. Install amalgkit
pip install amalgkit  # or: uv pip install amalgkit

# 5. Verify the data is accessible
uv run python scripts/maintenance/count_amalgkit_progress.py
```

## Pipeline Progress at Transfer Time

| Species | GetFastQ | Quant | Status |
|---|---|---|---|
| acromyrmex_echinatior | 24 | 8 | Partial |
| amellifera | 34 | 3 | Partial |
| anoplolepis_gracilipes | 14 | 7 | Partial |
| apis_mellifera_all | 0 | 0 | Metadata only |
| atta_cephalotes | 354 | 148 | Partial |
| camponotus_floridanus | 412 | 153 | Partial |
| cardiocondyla_obscurior | 148 | 51 | Partial |
| dinoponera_quadriceps | 28 | 2 | Partial |
| formica_exsecta | 42 | 14 | Partial |
| harpegnathos_saltator | 14 | 178 | **Needs getfastq** |
| linepithema_humile | 306 | 109 | Partial |
| monomorium_pharaonis | 310 | **360** | ✅ **Complete** |
| nylanderia_fulva | 62 | 31 | Partial |
| odontomachus_brunneus | 38 | 19 | Partial |
| ooceraea_biroi | 412 | 206 | Partial |
| pbarbatus | 36 | 18 | Partial |
| solenopsis_invicta | 430 | 119 | Partial |
| temnothorax_americanus | 64 | 32 | Partial |
| temnothorax_curvispinosus | 86 | 43 | Partial |
| temnothorax_longispinosus | 548 | 214 | **Was actively downloading** |
| temnothorax_nylanderi | 60 | 30 | Partial |
| vollenhovia_emeryi | 30 | 15 | Partial |
| wasmannia_auropunctata | 52 | 26 | Partial |

## Resume Commands

### Resume getfastq for a species (download FASTQ files)

```bash
# Generic pattern:
amalgkit getfastq \
  --out_dir blue/amalgkit/<SPECIES>/fastq \
  --metadata blue/amalgkit/<SPECIES>/work/metadata/metadata.tsv \
  --threads 4 --redo no --aws yes --max_bp 50000000

# Example for Harpegnathos (highest priority — most quant done but few fastq):
nohup amalgkit getfastq \
  --out_dir blue/amalgkit/harpegnathos_saltator/fastq \
  --metadata blue/amalgkit/harpegnathos_saltator/work/metadata/metadata.tsv \
  --threads 4 --redo no --aws yes --max_bp 50000000 \
  > blue/amalgkit/harpegnathos_saltator/getfastq.log 2>&1 &
```

### Resume quant for a species (quantify with Kallisto)

```bash
# Generic pattern:
amalgkit quant \
  --out_dir blue/amalgkit/<SPECIES>/work \
  --metadata blue/amalgkit/<SPECIES>/work/metadata/metadata.tsv \
  --threads 4 --redo no --clean_fastq yes \
  --index_dir /path/to/blue/data/amalgkit/shared/genome/<Species_name>/index \
  --fasta_dir /path/to/blue/data/amalgkit/shared/genome/<Species_name>

# Example for Temnothorax longispinosus:
nohup amalgkit quant \
  --out_dir blue/amalgkit/temnothorax_longispinosus/work \
  --metadata blue/amalgkit/temnothorax_longispinosus/work/metadata/metadata.tsv \
  --threads 4 --redo no --clean_fastq yes \
  > blue/amalgkit/temnothorax_longispinosus/quant.log 2>&1 &
```

### Run all species in batch

```bash
# Process each species sequentially (getfastq then quant):
for species in harpegnathos_saltator temnothorax_longispinosus solenopsis_invicta; do
  echo "=== Processing $species ==="
  
  # getfastq
  amalgkit getfastq \
    --out_dir blue/amalgkit/$species/fastq \
    --metadata blue/amalgkit/$species/work/metadata/metadata.tsv \
    --threads 4 --redo no --aws yes --max_bp 50000000
  
  # quant (update index_dir/fasta_dir per species config)
  # See config/amalgkit/amalgkit_${species}.yaml for correct paths
done
```

## Config Files

All species configs are in `config/amalgkit/amalgkit_<species>.yaml`. Each config contains:

- `genome.dest_dir`: Path to shared genome directory
- `steps.quant.index_dir`: Path to Kallisto index
- `steps.quant.fasta_dir`: Path to genome FASTA
- `steps.getfastq.max_bp`: Max base pairs per sample (50M)

**Important**: Update paths in config YAML files if your mount point differs from `/Volumes/blue/data`.

## Syncing Results Back to Git

After processing completes on Linux, sync results back:

```bash
uv run python scripts/rna/sync_quant_results.py
git add output/amalgkit_results/
git commit -m "feat(rna): sync results from Linux processing"
git push
```

## Key Directories on blue/

```
blue/ -> /path/to/mount/data/
  amalgkit/
    shared/genome/<Species>/       # Shared genome references + Kallisto indices
    <species>/
      fastq/getfastq/<SRR_ID>/     # Downloaded FASTQ files
      work/
        metadata/metadata.tsv      # Sample metadata (from SRA)
        quant/<SRR_ID>/             # Kallisto quantification results
        logs/                       # Per-sample quant logs
        curate/                     # Curated expression tables
        sanity/                     # Sanity check outputs
```

## Known Issues

1. **Harpegnathos index**: Had duplicate `.idx` files — cleaned to single `Harpegnathos_saltator_transcripts.idx`
2. **Monomorium DRR029611**: Had extra single-end fastq alongside paired — removed `DRR029611.fastq.gz`
3. **`--redo no`**: Always use this to skip already-completed samples
4. **`--clean_fastq yes`**: On `quant` step, deletes fastq after quantification to save disk space
