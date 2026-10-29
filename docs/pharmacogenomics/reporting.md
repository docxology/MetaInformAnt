# Clinical Reporting

`generate_clinical_report` builds a report dictionary from a patient's diplotypes (or
phenotypes) and the built-in CPIC guidelines; `export_report` renders it as plain text, HTML
or JSON for EHR integration, patient portals or research write-ups.

## Building a report

```python
from metainformant.pharmacogenomics.clinical.reporting import (
    generate_clinical_report,
    export_report,
)

patient = {'patient_id': 'PT001', 'sex': 'F'}
genotypes = {
    'CYP2D6': {'diplotype': '*4/*10'},
    'CYP2C19': {'diplotype': '*2/*2'},
}

report = generate_clinical_report(patient, genotypes)                    # every CPIC drug for these genes
report = generate_clinical_report(patient, genotypes, drugs=['codeine'])  # only codeine
report = generate_clinical_report(patient, genotypes, drugs=[])          # genotype results only
```

`generate_clinical_report(patient_data, genotypes, drugs=None)`:

- `patient_data` is a dict; the keys `patient_id`, `sex`, `ethnicity` and `ordering_provider`
  are used.
- `genotypes` maps a gene symbol to `{'diplotype': ..., 'phenotype': ..., 'activity_score': ...}`;
  give a diplotype, a phenotype, or both.
- `drugs=None` (the default) analyses every drug that has a built-in CPIC guideline against the
  genotyped genes: codeine, tramadol, tamoxifen, clopidogrel, voriconazole, warfarin,
  fluorouracil, azathioprine and simvastatin. A list names exactly the drugs to analyse, and an
  empty list requests no drug recommendations.

### Where each gene's phenotype comes from

| Input for a gene | Phenotype in the report |
|------------------|-------------------------|
| a `diplotype` that `classify_phenotype` accepts | classified from the diplotype, even when a `phenotype` key disagrees |
| no `diplotype`, a `phenotype` | the phenotype as supplied, in any form `MetabolizerPhenotype.from_string` accepts (`'Poor Metabolizer'`, `'PM'`, `'poor'`); free text such as `'Likely Poor Metabolizer'` keeps its category abbreviation |
| a malformed diplotype or a gene without a supported phenotype definition (including CYP1A2) | supplied phenotype when present, otherwise `Indeterminate` (`IND`); `clinical_significance` records the classification failure |
| a supported gene's diplotype with an uncharacterized allele | `Indeterminate`, score `None` and the allele named in `unknown_alleles` |
| neither | `Indeterminate` |

Drug recommendations are built from the phenotypes shown in `genotype_results`, so a gene needs
only a diplotype or a phenotype, not both. A drug produces an entry for a gene only when a
built-in CPIC guideline exists for that drug-gene pair and has a recommendation for the
phenotype; a gene without a guideline (CYP3A5, for example) is reported but adds none.

## Report structure

The report is a dict with these keys:

| Key | Content |
|-----|---------|
| `header` | `report_type`, `report_version` (`'1.0'`), `generated_at` (UTC ISO 8601 timestamp) and `generator` |
| `patient_info` | `patient_id` (`'UNKNOWN'` if absent), `sex`, `ethnicity`, `ordering_provider` (`'Not specified'` if absent). `name` and `dob` are not copied |
| `genotype_results` | one dict per gene: `gene`, `diplotype`, `activity_score`, `phenotype`, `phenotype_abbreviation`, `clinical_significance`; successful classification also includes `unknown_alleles` |
| `drug_recommendations` | one dict per drug-gene pair: `drug`, `gene`, `phenotype`, `recommendation`, `evidence_level`, `source`, `severity`, `urgency`, `alternatives` |
| `interaction_summary` | `overall_severity`, `major_count`, `moderate_count`, `minor_count`, `total_interactions`, `risk_level`, `summary` |
| `clinical_actions` | the Major and Moderate recommendations, Major first, each with `priority`, `drug`, `gene`, `severity`, `action`, `alternatives` |
| `disclaimer` | the clinical decision-support disclaimer (`add_disclaimer(report, custom_disclaimer=None)` replaces it) |

`urgency` is `'ACTION REQUIRED'` for a Major finding, `'ATTENTION RECOMMENDED'` for Moderate and
`'FOR INFORMATION'` otherwise. `generate_summary_table(results)` flattens genotype results or
recommendations into rows of `Gene`, `Diplotype`, `Phenotype`, `Drug`, `Recommendation` and
`Severity`.

## Export formats

```python
text = export_report(report, format='text')   # plain text
json_text = export_report(report, format='json')
html = export_report(report, format='html', output_path='output/pharmacogenomics/report.html')
```

Unknown activity scores render as `N/A` in text and HTML and as `null` in JSON.

`export_report(report, format='text', output_path=None)` returns the rendered string and, when
`output_path` is given, also writes it there (creating parent directories). Any `format` other
than `'text'`, `'html'` or `'json'` raises `ValueError`.

- **text** prints, in order: the header, patient information, genotype results, drug
  recommendations with their CPIC level (when there are any), prioritized clinical actions
  (when there are any), the interaction summary and the disclaimer.
- **html** is one self-contained document with inline CSS: patient information, genotype
  results, drug recommendations, clinical action items and the disclaimer. It has no
  interaction summary and uses no template engine or external template; the free-text fields
  (patient details, gene, diplotype, phenotype, recommendation, alternatives) are HTML-escaped.
- **json** is `json.dumps(report, indent=2)` with phenotype and severity enums written as their
  string values.

## PHI handling

`patient_info` carries only `patient_id`, `sex`, `ethnicity` and `ordering_provider`; a patient's
name and date of birth never reach the report. `patient_id` is copied as given, so pass a
pseudonymous identifier when the report must be de-identified.

## Audit fields

The header records `report_version` and the `generated_at` timestamp. It does not record a
report UUID or the versions of the package or of the CPIC data, so store those yourself with the
report (for example `metainformant.__version__`) if your audit trail needs them.

## Known limitations

- Recommendations come from the 9 built-in level-A CPIC guidelines; `generate_clinical_report`
  takes no `guidelines` argument.
- There is no PDF output; convert the HTML with an external tool such as `weasyprint`.
- Reports are in English only.
