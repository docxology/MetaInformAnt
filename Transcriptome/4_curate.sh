#!/bin/bash

# This script performs the merging and curation of expression data using amalgkit.

# Run amalgkit merge to combine expression data from multiple samples
amalgkit merge

# Set the normalization method for curation
normalization_method="log2p1-fpkm"

# Run amalgkit curate with the specified normalization method
amalgkit curate --norm $normalization_method
