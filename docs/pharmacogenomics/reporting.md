# Clinical Report Generation

The reporting module generates comprehensive clinical pharmacogenomic reports that integrate genotype results, phenotype classifications, drug-specific recommendations, and clinical action items. Reports can be exported in plain text, HTML, or JSON formats.

## Key Concepts

### Report Structure

A clinical report contains the following sections:

1. **Header**: Report metadata (type, version, timestamp, generator)
2. **Patient Info**: Sanitized demographics (ID, sex, ethnicity, ordering provider)
3. **Genotype Results**: Per-gene diplotype, activity score, phenotype, clinical significance
4. **Drug Recommendations**: Per-drug clinical recommendations with severity and alternatives
5. **Interaction Summary**: Overall risk assessment with severity breakdown
6. **Clinical Actions**: Prioritized action items requiring clinical attention
7. **Disclaimer**: Clinical disclaimer text

### Export Formats

- **Text**: Plain text report with sections separated by horizontal rules
- **HTML**: Styled HTML document with color-coded severity indicators
- **JSON**: Machine-readable JSON with full report structure

### Clinical Disclaimer

Every report includes a standard clinical disclaimer stating the report is for informational purposes, is not a substitute for professional judgment, and should be interpreted in the context of the complete clinical picture.

## Function Reference

### generate_clinical_report

```python
def generate_clinical_report(
    patient_data: dict[str, Any],
    genotypes: dict[str, dict[str, Any]],
    drugs: list[str] | None = None,
) -> dict[str, Any]
```

Generate a comprehensive clinical report. `patient_data` includes optional keys like `patient_id`, `name`, `dob`, `sex`, `ethnicity`, `ordering_provider`. `genotypes` maps gene symbol to `{"diplotype": "*1/*4", ...}`. If `drugs` is provided, generates drug-specific recommendations. Returns the full report dictionary.

### format_recommendation

```python
def format_recommendation(
    drug: str,
    gene: str,
    phenotype: str,
    guideline: dict[str, Any],
) -> dict[str, Any]
```

Format a single drug-gene recommendation for the report. Adds urgency indicators: "ACTION REQUIRED" for Major severity, "ATTENTION RECOMMENDED" for Moderate, "FOR INFORMATION" otherwise.

### generate_summary_table

```python
def generate_summary_table(
    results: list[dict[str, Any]],
) -> list[dict[str, str]]
```

Create a condensed tabular view of pharmacogenomic findings. Returns rows with standardized columns: `Gene`, `Diplotype`, `Phenotype`, `Drug`, `Recommendation`, `Severity`.

### export_report

```python
def export_report(
    report: dict[str, Any],
    format: str = "text",
    output_path: str | Path | None = None,
) -> str
```

Export a report to text, HTML, or JSON format. If `output_path` is provided, writes the report to file (creating directories as needed). Returns the formatted report string.

### add_disclaimer

```python
def add_disclaimer(
    report: dict[str, Any],
    custom_disclaimer: str | None = None,
) -> dict[str, Any]
```

Add or update the clinical disclaimer in a report. If `custom_disclaimer` is None, uses the standard METAINFORMANT disclaimer.

## Usage Examples

```python
from metainformant.pharmacogenomics import (
    generate_clinical_report, export_report,
    generate_summary_table, add_disclaimer,
)

# Patient and genotype data
patient = {
    "patient_id": "PGX-001",
    "sex": "Female",
    "ethnicity": "European",
    "ordering_provider": "Dr. Smith",
}

genotypes = {
    "CYP2D6": {"diplotype": "*1/*4"},
    "CYP2C19": {"diplotype": "*1/*17"},
    "DPYD": {"diplotype": "*1/*1"},
}

# Generate the full clinical report
report = generate_clinical_report(
    patient, genotypes,
    drugs=["codeine", "clopidogrel", "fluorouracil"],
)

# Export as plain text
text_output = export_report(report, format="text")
print(text_output)

# Export as HTML to file
export_report(report, format="html", output_path="output/pgx_report.html")

# Export as JSON
json_output = export_report(report, format="json")

# Generate a summary table
summary = generate_summary_table(report["genotype_results"])
for row in summary:
    print(f"{row['Gene']}: {row['Diplotype']} -> {row['Phenotype']}")

# Add a custom disclaimer
report = add_disclaimer(report, custom_disclaimer="Custom institutional disclaimer text.")
```

## Configuration

- **Environment prefix**: `PHARMA_`
- HTML reports include responsive CSS styling with color-coded severity classes
- JSON export handles custom types (MetabolizerPhenotype, InteractionSeverity, DrugRecommendation)
- Output paths are auto-created if they do not exist

## Related Modules

- `pharmacogenomics.alleles.phenotype` -- Phenotype classification feeds into genotype results
- `pharmacogenomics.clinical.drug_interaction` -- Interaction analysis populates drug recommendations
- `pharmacogenomics.annotations.cpic` -- CPIC guidelines referenced in recommendations
