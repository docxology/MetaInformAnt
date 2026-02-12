# R Installation Guide for Amalgkit Workflows

## Quick Install Commands

### Debian/Ubuntu (Recommended)
```bash
sudo apt-get update
sudo apt-get install -y r-base r-base-dev
```

### Verify Installation
```bash
Rscript --version
which Rscript
```

### Test R Functionality
```bash
# Test basic R execution
Rscript -e "print('R is working')"

# Test required capabilities
Rscript -e "capabilities()"
```

### Configure R Temp Directory (IMPORTANT)

**Problem**: R package installation may fail with "No space left on device" if `/tmp` (tmpfs) is full.

**Solution**: Configure R to use the repository's `.tmp/R` directory instead of system `/tmp`:

```bash
# Set environment variables for current session
export TMPDIR="$(pwd)/.tmp/R"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p .tmp/R

# Verify configuration
Rscript -e "cat('TMPDIR:', Sys.getenv('TMPDIR'), '\n')"
```

**Permanent Configuration**: Create `.Rprofile` in repository root:

```bash
cat > .Rprofile << 'EOF'
# Set R temp directory to repository .tmp/R instead of system /tmp
if (Sys.getenv("TMPDIR") == "") {
  repo_root <- normalizePath(".")
  tmp_dir <- file.path(repo_root, ".tmp", "R")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMPDIR = tmp_dir)
  Sys.setenv(TEMP = tmp_dir)
  Sys.setenv(TMP = tmp_dir)
}
EOF
```

**Note**: The `.Rprofile` file is automatically loaded by R when starting in the repository directory.

## Alternative Installation Methods

### Using apt-get (with specific version)
```bash
# Add CRAN repository for latest R version
sudo apt-get install -y software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt-get install -y r-base r-base-dev
```

### Fedora/RHEL/CentOS
```bash
sudo dnf install R
# or on older systems
sudo yum install R
```

### macOS (with Homebrew)
```bash
brew install r
```

### From Source (Advanced)
```bash
# Install build dependencies
sudo apt-get install -y build-essential gfortran libreadline-dev \
    libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev \
    xvfb libbz2-dev libzstd-dev liblzma-dev libcurl4-openssl-dev

# Download and compile R
cd /tmp
wget https://cran.r-project.org/src/base/R-4/R-4.3.2.tar.gz
tar -xzf R-4.3.2.tar.gz
cd R-4.3.2
./configure --enable-R-shlib --with-blas --with-lapack
make
sudo make install
```

## What R is Used For in Amalgkit

### Merge Step
- **Optional**: Visualization plots (PCA, correlation heatmaps)
- **Core functionality**: Still works without R (expression matrices created)
- **Impact if missing**: Plots not generated, but data analysis proceeds

### Curate Step  
- **Required**: Statistical analysis and batch effect correction
- **Functions**: Outlier detection, normalization, QC metrics
- **Impact if missing**: Curate step fails completely

## Post-Installation Verification

Run the automated check:
```bash
python3 scripts/rna/check_r_dependencies.py
```

Expected output when R is installed:
```
================================================================================
CHECKING R DEPENDENCIES
================================================================================
✅ Rscript found: /usr/bin/Rscript
   Version: R scripting front-end version 4.x.x
================================================================================
✅ R DEPENDENCIES: SATISFIED
================================================================================
```

## Testing End-to-End with R

Once R is installed, test the complete workflow:

```bash
# Run merge + curate + sanity for a species
bash scripts/rna/amalgkit/run_amalgkit.sh \
    --config config/amalgkit/amalgkit_pbarbatus.yaml \
    --steps merge,curate,sanity
```

### Expected Behavior

#### With R Installed ✅
- **Merge**: Creates expression matrices + generates plots
- **Curate**: Performs statistical QC and batch correction
- **Sanity**: Validates all outputs

#### Without R ⚠️
- **Merge**: Creates expression matrices, skips plots (warning)
- **Curate**: Fails with error about missing Rscript
- **Sanity**: Still validates core outputs

## Current Workflow Status (November 2025)

- **Batch 1**: Processing 3,820 samples across 10 ant species
- **Batch 2**: Queued - 728 samples across 10 ant species
- **R-dependent steps**: Required for `curate` and `merge` plotting
- **R installation**: See installation methods above

## Workaround: Skip R-Dependent Steps

If R cannot be installed, you can:

1. **Skip curate entirely**:
   ```bash
   # Run only merge and sanity
   bash scripts/rna/amalgkit/run_amalgkit.sh \
       --config config/amalgkit/amalgkit_pbarbatus.yaml \
       --steps merge,sanity
   ```

2. **Use Python-based QC** (alternative to curate):
   ```python
   # Use scikit-learn or pandas for basic QC
   import pandas as pd
   tpm = pd.read_csv("output/.../Solenopsis_invicta_tpm.tsv", sep="\t")
   # Apply log transformation, normalization, etc.
   ```

3. **Export data for R analysis elsewhere**:
   - Expression matrices are standard TSV format
   - Can be analyzed in R Studio or other environments

## Installation Troubleshooting

### Issue: sudo not available
```bash
# Check if you're in a restricted environment
sudo -l

# If no sudo access, contact system administrator
# Or use conda/mamba for user-level R installation (system tool):
conda install -c conda-forge r-base

# Note: For Python packages, use uv pip install (primary method)
```

### Issue: Package conflicts
```bash
# Remove existing R installations
sudo apt-get remove r-base r-base-core
sudo apt-get autoremove
sudo apt-get autoclean

# Reinstall
sudo apt-get update
sudo apt-get install -y r-base
```

### Issue: R installed but Rscript not in PATH
```bash
# Find Rscript location
find /usr -name Rscript 2>/dev/null

# Add to PATH (add to ~/.bashrc for persistence)
export PATH="/usr/bin:$PATH"
```

## Conda/Mamba Alternative (No sudo required)

```bash
# Install conda/mamba if not available
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
source $HOME/mambaforge/bin/activate

# Install R
mamba install -c conda-forge r-base

# Verify
Rscript --version
```

## Summary

**Recommended command for most users**:
```bash
sudo apt-get update && sudo apt-get install -y r-base r-base-dev
```

**Verify**:
```bash
python3 scripts/rna/check_r_dependencies.py
```

**Test end-to-end**:
```bash
# Run complete workflow with R steps
bash scripts/rna/amalgkit/run_amalgkit.sh \
    --config config/amalgkit/amalgkit_pbarbatus.yaml \
    --steps merge,curate,sanity
```

---

**Note**: The core quantification and merge functionality works without R. Only advanced QC (curate) and optional visualizations require R.

