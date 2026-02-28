#!/usr/bin/env Rscript
# Comprehensive R package installer for amalgkit
# Must be run with permissions to write to R library path
#
# SYSTEM DEPENDENCIES:
# Ensure you have run: `sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev`
# BEFORE running this script, or Bioconductor packages (sva, RUVSeq) will fail to compile.
#
# Required packages:
#   CRAN:         ggplot2, Rtsne
#   Bioconductor: edgeR, RUVSeq, sva

cat("═══ Amalgkit R Package Installer ═══\n\n")

# 1. BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2. CRAN packages
cran_pkgs <- c("ggplot2", "Rtsne")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from CRAN...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# 3. Bioconductor packages
bioc_pkgs <- c("edgeR", "RUVSeq", "sva")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from Bioconductor...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# 4. Verify
cat("\n═══ Verification ═══\n")
all_ok <- TRUE
for (pkg in c(cran_pkgs, bioc_pkgs)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s: %s\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("  ✗ %s: MISSING\n", pkg))
    all_ok <- FALSE
  }
}

if (all_ok) {
  cat("\n✅ All amalgkit R packages installed.\n")
} else {
  cat("\n❌ Some packages failed to install. Check errors above.\n")
  quit(status = 1)
}
