# MetaInformAnt RNA-seq Workflow - Quick Setup Commands

## Current Status
- ✓ Python 3.11.2
- ✓ uv 0.8.17
- ✓ Virtual environment (.venv)
- ✓ metainformant installed
- ✓ amalgkit installed
- ✗ SRA Toolkit (fasterq-dump, prefetch) - **MISSING**
- ✗ kallisto - **MISSING**

## Commands to Run Now

### 1. Install Missing System Tools

```bash
# Install SRA Toolkit and kallisto
sudo apt-get update
sudo apt-get install -y sra-toolkit kallisto

# Verify installations
fasterq-dump --version
kallisto version
```

### 2. Set Environment Variable (Optional but Recommended)

```bash
# Add to ~/.bashrc (permanent)
echo 'export NCBI_EMAIL="your.email@example.com"' >> ~/.bashrc
source ~/.bashrc

# OR set for current session only
export NCBI_EMAIL="your.email@example.com"
```

### 3. Run the Workflow

```bash
cd /home/q/Documents/GitHub/MetaInformAnt
source .venv/bin/activate
python3 scripts/rna/run_multi_species.py
```

## Complete One-Liner Setup (if starting fresh)

```bash
cd /home/q/Documents/GitHub/MetaInformAnt && \
sudo apt-get update && \
sudo apt-get install -y sra-toolkit kallisto && \
source .venv/bin/activate && \
python3 scripts/rna/run_multi_species.py
```

## Verification Commands

```bash
# Check all tools are available
source .venv/bin/activate
command -v amalgkit && echo "✓ amalgkit"
command -v fasterq-dump && echo "✓ fasterq-dump" 
command -v kallisto && echo "✓ kallisto"
```

## What Each Tool Does

- **amalgkit**: RNA-seq metadata management and workflow orchestration
- **fasterq-dump**: Downloads FASTQ files from NCBI SRA
- **prefetch**: Pre-downloads SRA data (bundled with sra-toolkit)
- **kallisto**: Fast RNA-seq quantification

