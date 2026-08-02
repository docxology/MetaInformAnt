# Comprehensive RNA analysis guide

This guide connects pipeline execution to a defensible transcriptome
meta-analysis.

## Pipeline lane

1. Declare species, tissue rules, genome/index assets, threads, and data root.
2. Retrieve metadata and record the query, timestamp, and source columns.
3. Select samples with explicit inclusion/exclusion criteria.
4. Retrieve and validate reads; retain failures and retry records.
5. Quantify each sample with a recorded index and tool version.
6. Merge abundance tables after sample-ID and feature-schema checks.
7. Apply `wsfilter`, then `finalize`, and validate each resulting table.
8. Run `sanity` and generate a compact evidence manifest.

## Metadata normalization

Tissue mappings are declared in configuration and applied through
`normalize_tissue_metadata.py`. The raw label, normalized label, mapping rule,
and unresolved values must remain traceable. Do not overwrite source metadata
without saving the original and a transformation record.

## Statistical analysis lane

Pipeline tables are inputs to statistics, not statistical conclusions. A full
analysis should specify:

- unit of replication and study-level grouping;
- count/abundance scale and filtering thresholds;
- normalization and batch model;
- within-study model and contrasts;
- effect-size scale and standard errors;
- random/fixed-effects meta-analysis model when combining studies;
- heterogeneity, sensitivity analysis, and leave-one-study-out checks;
- false-discovery control and the family of tests;
- orthology mapping and feature intersection for cross-species work.

The current project cross-species script is descriptive. It reports pairwise
expression fingerprints and visualizations, but it does not produce p-values,
confidence intervals, phylogenetic corrections, or multiplicity-adjusted
significance without an explicitly validated extension.

## Reproducible figures

Every figure must be generated from a manifest-linked table and include:

- a result-derived title and caption;
- species/sample counts after filtering;
- scale and transformation;
- missing-data handling;
- color and ordering rules;
- uncertainty or a statement that the visualization is descriptive;
- the script, tool versions, and input checksum.

## Release checklist

Run the deterministic tests, configuration/schema validators, documentation
link audit, evidence-manifest generator, report generator, and figure
rendering checks. Compare regenerated reports with committed artifacts and
explain any difference before publication.
