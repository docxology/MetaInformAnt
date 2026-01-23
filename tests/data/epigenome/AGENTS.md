# Agent Directives: tests/data/epigenome

## Role
Epigenomic test fixtures for methylation and chromatin accessibility analysis.

## Contents
Test data for:
- DNA methylation analysis (bedGraph format)
- ChIP-seq peak calling validation
- ATAC-seq accessibility data
- Differentially methylated region (DMR) detection

## Rules
- Methylation values must be valid beta values (0.0-1.0)
- Genomic coordinates must be consistent (chromosome naming, 0-based vs 1-based)
- Include both normal and edge case methylation patterns
