# Agent Directives: pharmacogenomics

## Role
Documentation agent for the pharmacogenomics module covering clinical variant interpretation and drug-gene interactions.

## Module Scope
- Star allele calling and diplotype determination
- Metabolizer phenotype prediction (CYP enzymes)
- CPIC guideline integration and dosing recommendations
- PharmGKB clinical annotations and evidence levels
- FDA drug label parsing and biomarker extraction
- ACMG variant classification and pathogenicity scoring
- Drug-gene interaction analysis and polypharmacy checks
- Clinical report generation with disclaimers

## Key Source Files
- `src/metainformant/pharmacogenomics/alleles/` - Star allele, diplotype, phenotype
- `src/metainformant/pharmacogenomics/annotations/` - CPIC, PharmGKB, drug labels
- `src/metainformant/pharmacogenomics/clinical/` - Pathogenicity, interactions, reporting
- `src/metainformant/pharmacogenomics/visualization/` - Clinical visualization

## External Dependencies
- ClinVar database for variant annotations
- gnomAD for population frequencies
- CPIC guidelines database
- PharmGKB API
