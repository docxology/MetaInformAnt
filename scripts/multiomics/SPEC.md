# SPEC: MultiOmics Scripts

Scripts for cross-platform data harmonization and joint multi-omic analysis.

## Workflows

- `harmonize_omics_datasets.py`: Aligns DNA, RNA, and protein datasets by sample ID and feature names.
- `run_joint_analysis.py`: Performs multi-view clustering and correlation analysis.

## Standards

- **Consistency**: Standardizes all omics data to common nomenclature (e.g., Gene Symbols) before integration.
- **Validation**: Verifies sample overlap and variance across platforms.
