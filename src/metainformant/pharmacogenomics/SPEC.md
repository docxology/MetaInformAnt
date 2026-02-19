# Specification: pharmacogenomics

## 🎯 Scope

Clinical pharmacogenomic analysis: star allele calling, metabolizer phenotyping, CPIC guideline lookups, drug interaction prediction, and report generation.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: alleles, annotations, metabolism, interaction, clinical, visualization
- **Key Concepts**: `StarAllele`, star allele calling, CPIC activity scores, metabolizer phenotypes

## 🔌 API Definition

### Exports

- `__init__.py`
