# Genome Setup Scripts Reference

Commands to download reference genomes, prepare transcriptomes, and build kallisto indices for all configured species.

## Quick Reference

### Option 1: Automated (Recommended)

```bash
# Run complete genome setup for all species
python3 scripts/rna/orchestrate_genome_setup.py

# Or via bash script
bash scripts/rna/run_genome_setup.sh
```

### Option 2: Step-by-Step

```bash
# Step 1: Check current status
python3 scripts/rna/verify_genomes_and_indexes.py

# Step 2: Download missing genomes
python3 scripts/rna/download_missing_genomes.py

# Step 3: Prepare transcriptomes (extract RNA FASTA)
python3 scripts/rna/prepare_transcriptomes.py

# Step 4: Build kallisto indexes
python3 scripts/rna/build_kallisto_indexes.py

# Step 5: Verify final status
python3 scripts/rna/verify_genomes_and_indexes.py
```

## Detailed Command Reference

### Status Verification

```bash
# Check all species
python3 scripts/rna/verify_genomes_and_indexes.py

# Check specific species
python3 scripts/rna/verify_genomes_and_indexes.py --species camponotus_floridanus

# Custom output location
python3 scripts/rna/verify_genomes_and_indexes.py --output output/my_status.json
```

### Download Genomes

```bash
# All missing genomes
python3 scripts/rna/download_missing_genomes.py

# Specific species
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus

# Dry run (preview without downloading)
python3 scripts/rna/download_missing_genomes.py --dry-run
```

### Prepare Transcriptomes

```bash
# All transcriptomes
python3 scripts/rna/prepare_transcriptomes.py

# Specific species
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
```

### Build Kallisto Indexes

```bash
# All indexes
python3 scripts/rna/build_kallisto_indexes.py

# Specific species
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus

# Custom k-mer size (for short reads <100bp)
python3 scripts/rna/build_kallisto_indexes.py --kmer-size 23
```

### Master Orchestrator

```bash
# Complete pipeline, all species
python3 scripts/rna/orchestrate_genome_setup.py

# Single species
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus

# Skip specific steps
python3 scripts/rna/orchestrate_genome_setup.py --skip-download
python3 scripts/rna/orchestrate_genome_setup.py --skip-prepare
python3 scripts/rna/orchestrate_genome_setup.py --skip-build
```

## Example Workflows

### Complete Setup for All Species

```bash
# Option A: Automated pipeline
bash scripts/rna/run_genome_setup.sh

# Option B: Master orchestrator
python3 scripts/rna/orchestrate_genome_setup.py

# Option C: Manual step-by-step
python3 scripts/rna/verify_genomes_and_indexes.py
python3 scripts/rna/download_missing_genomes.py
python3 scripts/rna/prepare_transcriptomes.py
python3 scripts/rna/build_kallisto_indexes.py
python3 scripts/rna/verify_genomes_and_indexes.py
```

### Setup for Single Species

```bash
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus

# Or step-by-step:
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus
```

### Resuming After Interruption

```bash
# Check what's missing
python3 scripts/rna/verify_genomes_and_indexes.py

# Run only missing steps
python3 scripts/rna/download_missing_genomes.py  # only downloads missing
python3 scripts/rna/prepare_transcriptomes.py     # only prepares missing
python3 scripts/rna/build_kallisto_indexes.py     # only builds missing
```

## Output Files

- `output/genome_index_status.json` — verification results
- `output/genome_download_results.json` — download results
- `output/transcriptome_preparation_results.json` — preparation results
- `output/kallisto_index_build_results.json` — index build results

Genome files are stored under: `output/amalgkit/<species>/genome/`  
Kallisto indexes are stored under: `output/amalgkit/<species>/work/index/`

## Kallisto Index Naming Convention

Index files must be named `{Scientific_Name}_transcripts.idx`:

| Species | Expected Index Filename |
|---------|------------------------|
| `Apis mellifera` | `Apis_mellifera_transcripts.idx` |
| `Pogonomyrmex barbatus` | `Pogonomyrmex_barbatus_transcripts.idx` |

## Troubleshooting

```bash
# Check progress during long downloads
python3 scripts/rna/verify_genomes_and_indexes.py
cat output/amalgkit/*/genome/download.progress.txt
```

## See Also

- [genome_setup_guide.md](genome_setup_guide.md) — Detailed setup guide
- [genome_preparation.md](genome_preparation.md) — Technical API reference
- [steps/06_quant.md](steps/06_quant.md) — Using kallisto indexes for quantification
