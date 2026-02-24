#!/usr/bin/env bash
# install_linux_deps.sh – Install all system-level dependencies for METAINFORMANT on Linux/Debian/Ubuntu/Parrot
# Usage: bash scripts/package/install_linux_deps.sh [--skip-r-packages] [--no-sudo]
#
# Tools installed via apt:
#   kallisto      – RNA-seq quantification (preferred over salmon)
#   sra-toolkit   – NCBI SRA downloads (fasterq-dump, prefetch)
#   fastp         – FASTQ quality control
#   seqkit        – General FASTA/FASTQ manipulation
#   samtools      – SAM/BAM processing
#   pigz          – Parallel gzip
#   parallel      – GNU parallel (for batch operations)
#   r-base        – R language (for amalgkit curate step)
#   wget, curl    – Download tools
#
# Tools installed from binary releases:
#   ncbi-datasets-cli – NCBI Datasets CLI (not in apt)
#
# R packages installed (requires r-base):
#   BiocManager, edgeR, limma, tximport, vegan, amap, Rtsne
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info()  { echo -e "${BLUE}ℹ️${NC}  $*"; }
ok()    { echo -e "${GREEN}✅${NC} $*"; }
warn()  { echo -e "${YELLOW}⚠️${NC}  $*"; }
error() { echo -e "${RED}❌${NC} $*"; }

SKIP_R_PACKAGES=0
NO_SUDO=0

for arg in "$@"; do
  case "$arg" in
    --skip-r-packages) SKIP_R_PACKAGES=1 ;;
    --no-sudo) NO_SUDO=1 ;;
    -h|--help)
      echo "Usage: bash $0 [--skip-r-packages] [--no-sudo]"
      exit 0
      ;;
  esac
done

SUDO="sudo"
if [[ "$NO_SUDO" -eq 1 ]]; then
  SUDO=""
  warn "Running without sudo - some installations may fail if you lack write permissions"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 1. Detect Linux package manager
# ─────────────────────────────────────────────────────────────────────────────
if ! command -v apt-get >/dev/null 2>&1; then
  error "This script requires apt-get (Debian/Ubuntu/Parrot). On other distros, install manually."
  error "Required tools: kallisto, sra-toolkit, fastp, seqkit, samtools, pigz, parallel, r-base"
  exit 1
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " METAINFORMANT – Linux Dependency Installer"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# 2. Install apt packages
# ─────────────────────────────────────────────────────────────────────────────
APT_PACKAGES=(
  wget
  curl
  pigz
  parallel
  seqkit
  samtools
  kallisto
  fastp
  sra-toolkit
  r-base
  r-base-dev
)

info "Updating apt package lists..."
$SUDO apt-get update -qq

info "Installing apt packages: ${APT_PACKAGES[*]}"
for pkg in "${APT_PACKAGES[@]}"; do
  if dpkg -l "$pkg" 2>/dev/null | grep -q "^ii"; then
    ok "$pkg already installed"
  else
    echo "  Installing $pkg..."
    $SUDO apt-get install -y "$pkg" 2>/dev/null && ok "$pkg installed" || warn "Failed to install $pkg (non-fatal)"
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 3. Install ncbi-datasets-cli binary (not in apt)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
info "Checking ncbi-datasets-cli..."
if command -v datasets >/dev/null 2>&1; then
  ok "ncbi-datasets-cli already installed: $(datasets version 2>/dev/null || echo 'version unknown')"
else
  info "Installing ncbi-datasets-cli from NCBI FTP..."
  DATASETS_URL="https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"
  DATAFORMATS_URL="https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat"

  # Prefer ~/.local/bin (always user-writable); fall back to /usr/local/bin with sudo
  INSTALL_DIR="$HOME/.local/bin"
  mkdir -p "$INSTALL_DIR"

  DATASETS_TMP=$(mktemp)
  DATAFORMATS_TMP=$(mktemp)

  if curl -fsSL "$DATASETS_URL" -o "$DATASETS_TMP" && \
     curl -fsSL "$DATAFORMATS_URL" -o "$DATAFORMATS_TMP"; then
    chmod +x "$DATASETS_TMP" "$DATAFORMATS_TMP"
    mv "$DATASETS_TMP" "$INSTALL_DIR/datasets"
    mv "$DATAFORMATS_TMP" "$INSTALL_DIR/dataformat"
    ok "ncbi-datasets-cli installed to $INSTALL_DIR"
    # Ensure ~/.local/bin is on PATH
    if [[ ":$PATH:" != *":$INSTALL_DIR:"* ]]; then
      warn "Add $INSTALL_DIR to your PATH:"
      warn "  echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.bashrc && source ~/.bashrc"
    fi
  else
    rm -f "$DATASETS_TMP" "$DATAFORMATS_TMP"
    warn "ncbi-datasets-cli install failed – download manually from:"
    warn "  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 4. Install R packages for amalgkit curate step
# ─────────────────────────────────────────────────────────────────────────────
echo ""
if [[ "$SKIP_R_PACKAGES" -eq 1 ]]; then
  warn "Skipping R package installation (--skip-r-packages)"
else
  if ! command -v R >/dev/null 2>&1; then
    warn "R not found – skipping R package installation"
  else
    info "Installing R packages for amalgkit curate step..."
    # Install to user library to avoid needing sudo for R packages
    R_LIB_DIR="$HOME/R/library"
    mkdir -p "$R_LIB_DIR"

    R --quiet --no-save << 'RSCRIPT'
# Set user library
lib_dir <- Sys.getenv("HOME")
lib_dir <- file.path(lib_dir, "R", "library")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))

cat("R library path:", lib_dir, "\n")

# CRAN packages
cran_pkgs <- c("vegan", "amap", "Rtsne", "ggplot2", "reshape2", "pheatmap")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    tryCatch(
      install.packages(pkg, repos = "https://cloud.r-project.org/", lib = lib_dir, quiet = TRUE),
      error = function(e) cat("  Failed to install", pkg, ":", conditionMessage(e), "\n")
    )
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# Bioconductor packages (for amalgkit curate)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/", lib = lib_dir, quiet = TRUE)
}

bioc_pkgs <- c("edgeR", "limma", "tximport")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing Bioconductor package:", pkg, "...\n")
    tryCatch(
      BiocManager::install(pkg, lib = lib_dir, ask = FALSE, update = FALSE, quiet = TRUE),
      error = function(e) cat("  Failed to install", pkg, ":", conditionMessage(e), "\n")
    )
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

cat("R package installation complete.\n")
RSCRIPT

    ok "R packages installed (user library: $R_LIB_DIR)"
    info "Add to ~/.bashrc: export R_LIBS_USER=\"\$HOME/R/library\""
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 5. Summary
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " Dependency Check Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

check_tool() {
  local name="$1"
  local cmd="${2:-$1}"
  if command -v "$cmd" >/dev/null 2>&1; then
    ok "$name: $(command -v "$cmd")"
  else
    warn "$name: NOT FOUND"
  fi
}

check_tool "kallisto"
check_tool "fasterq-dump"
check_tool "prefetch"
check_tool "fastp"
check_tool "seqkit"
check_tool "samtools"
check_tool "R"
check_tool "ncbi-datasets-cli" "datasets"
check_tool "wget"
check_tool "curl"
check_tool "pigz"
check_tool "parallel"

echo ""
info "Run 'bash scripts/package/setup.sh' to install Python packages and amalgkit."
echo ""
