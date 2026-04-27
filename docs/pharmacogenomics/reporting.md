# Clinical Reporting

Generate human‑readable and machine‑parseable pharmacogenomic reports for EHR
integration, patient portals, or research dissemination.

## Report sections (default order)

1. **Header** — patient demographic block, report UUID, generation timestamp,
   metainformant version, CPIC version used.
2. **Genotype summary** — per‑gene diplotype, activity score, phenotype with
   evidence breakdown (which alleles contributed what).
3. **Drug recommendations** — for each actionable gene–drug pair, CPIC level,
   recommendation text, and classification (strong/moderate/optional).
4. **Drug–Drug Interaction review** — any contraindicated/major DDI based on the
   patient's medication list.
5. **Disclaimer** — mandatory statement that this is decision‑support, not a
   prescription; clinician review required.
6. **Appendix** (optional) — PharmGKB annotations linked to observed variants.

## Export formats

```python
from metainformant.pharmacogenomics.clinical.reporting import (
    generate_clinical_report,
    export_report,
)

report = generate_clinical_report(patient)

text  = export_report(report, format='text')   # plain monospaced
json  = export_report(report, format='json')   # schema v1.0
html  = export_report(report, format='html', template='default')
```

HTML uses a Jinja2 template located in `src/metainformant/pharmacogenomics/clinical/templates/`.

### Custom template

```python
html = export_report(report, format='html', template='my_clinic.html')
```

Place your template under `~/.hermes/templates/pharmacogenomics/`. The context
passed to Jinja2 is:

```python
{
  'patient': patient_dict,
  'genotype_table': list of per‑gene dicts,
  'drug_recommendations': list of {drug,gene,phenotype,recommendation,level},
  'ddi_alerts': list of {drug_a,drug_b,severity,recommendation},
  'report_meta': {report_id, generated_at, version, cpic_version},
}
```

## Section control

```python
report = generate_clinical_report(
    patient,
    sections={
        'ddi': True,               # include interaction section
        'pharmgkb': False,         # omit annotation footnotes
        'disclaimer': True,        # mandatory for clinical use
        'appendix': False,
    }
)
```

## PHI handling

Report includes patient name / DOB if supplied. To anonymise:

```python
anon = patient.copy()
anon.pop('name', None)
anon['patient_id'] = 'REDACTED'
report = generate_clinical_report(anon, include_phi=False)
```

## Audit fields (JSON / HTML)

Each report embeds:

```json
{
  "report_id": "uuid4",
  "generated_at": "2025-06-14T10:22:15Z",
  "metainformant_version": "0.8.0",
  "cpic_version": "3.0",
  "pharmgkb_cache_date": "2025-05-15"
}
```

Store these with patient ID in your audit database for CLIA/CAP traceability.

## Known limitations

- Report timing: generation ~2 ms; negligible compared to phenotype computation.
- PDF output not directly supported; convert HTML with `weasyprint` or similar.
- Multi‑language localisation not yet available (English only).
