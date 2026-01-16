# Install R dependencies for amalgkit
cran_packages <- c(
  "amap", "RColorBrewer", "colorspace", "dendextend", "NMF", 
  "MASS", "pvclust", "Rtsne", "ggplot2", "patchwork", "optparse",
  "BiocManager", "reshape2", "gridExtra"
)

# Install CRAN packages
install_packages <- cran_packages[!cran_packages %in% installed.packages()[,"Package"]]
if(length(install_packages) > 0) {
  install.packages(install_packages, repos="https://cloud.r-project.org")
}

# Install Bioconductor packages
bioc_packages <- c("pcaMethods", "edgeR", "RUVSeq", "sva")
install_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
if(length(install_bioc) > 0) {
  BiocManager::install(install_bioc, ask=FALSE, update=FALSE)
}

# Final check
missing <- c(cran_packages, bioc_packages)
missing <- missing[!missing %in% installed.packages()[,"Package"]]
if(length(missing) > 0) {
  cat("FAILED to install packages:", paste(missing, collapse=", "), "\n")
  quit(status=1)
} else {
  cat("SUCCESS: All R packages installed.\n")
}
