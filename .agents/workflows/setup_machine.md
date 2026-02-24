---
description: Replicate the MetaInformAnt environment on a new machine (Linux/macOS)
---
# Move Environment to a New Computer

// turbo-all

## 1. Install `uv` (Python Package Manager)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv --version
```

## 2. Install System Dependencies

**Linux (Debian/Ubuntu/Parrot)**:

```bash
bash scripts/package/install_linux_deps.sh
```

This installs: `kallisto`, `sra-toolkit` (fasterq-dump, prefetch), `fastp`, `seqkit`, `samtools`, `pigz`, `parallel`, `r-base`, `ncbi-datasets-cli`, and R packages for the amalgkit curate step.

**macOS**:

```bash
bash scripts/package/setup.sh --with-deps
```

## 3. Run the Python Setup Script

```bash
bash scripts/package/setup.sh --skip-tests
```

*This script auto-detects FAT/exFAT filesystems and routes the venv to `/tmp/metainformant_venv` if needed.*

## 4. Activate the Virtual Environment

Standard (ext4/NTFS):

```bash
source .venv/bin/activate
```

FAT/exFAT filesystem:

```bash
source /tmp/metainformant_venv/bin/activate
```

## 5. Set Required Environment Variables

```bash
# NCBI_EMAIL is required for RNA-seq metadata downloads
export NCBI_EMAIL="YourEmail@example.com"
echo 'export NCBI_EMAIL="YourEmail@example.com"' >> ~/.bashrc

# R packages user library (set so amalgkit curate can find them)
export R_LIBS_USER="$HOME/R/library"
echo 'export R_LIBS_USER="$HOME/R/library"' >> ~/.bashrc
```

## 6. Create RNA-seq Output Directory Structure

```bash
bash scripts/rna/setup_genome_dirs.sh
```

Creates incremental storage directories for all species under `output/amalgkit/`.
**Key rule**: Never delete `output/amalgkit/*/work/quant/` — this is where processed quant data accumulates incrementally.

## 7. Verify Installation

```bash
# Python environment
python -c "import metainformant; print('✓ metainformant', metainformant.__version__)"
python -c "import amalgkit; print('✓ amalgkit')"

# External CLI tools
kallisto version && echo "✓ kallisto"
fasterq-dump --version 2>&1 | head -1 && echo "✓ fasterq-dump"
fastp --version 2>&1 | head -1 && echo "✓ fastp"
datasets version 2>/dev/null && echo "✓ ncbi-datasets-cli"
R --version | head -1 && echo "✓ R"

# Run fast tests (no network required)
uv run pytest tests/ -m "not network and not external" -q --tb=short
```

## 8. Run First RNA-seq Workflow

```bash
# Check status of a species workflow
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status

# Start workflow (incremental – resumes from where it left off)
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml \
  > output/amalgkit/pbarbatus_workflow.log 2>&1 &
echo "Workflow PID: $!"
```

Monitor progress:

```bash
tail -f output/amalgkit/pbarbatus_workflow.log
```
