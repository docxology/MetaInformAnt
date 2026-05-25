# Integration Guide: Pharmacogenomics

Connecting PGx predictions to external systems — EHR, CDS hooks, databases, and
micro‑services.

## HL7 FHIR Observation resource

```python
from metainformant.pharmacogenomics.clinical.reporting import generate_clinical_report
report = generate_clinical_report(patient)
fhir_obs = report.to_fhir()   # fhir.resources.Observation
```

Extension URL: `http://hl7.org/fhir/StructureDefinition/pharmacogenetics`. Component
`valueCodeableConcept` encodes phenotype string per gene.

## CCDA / C-CDA narrative block

```html
<!-- Inside <section> ... -->
<narrative>
  <table>
    <tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Recommendation</th></tr>
    <tr><td>CYP2D6</td><td>*1/*4</td><td>Intermediate Metabolizer</td>
        <td>Avoid codeine; consider alternative analgesic (CPIC Level A)</td></tr>
  </table>
</narrative>
```

Generate HTML fragment via `export_report(report, format='html')` and embed.

## REST API pattern (FastAPI)

```python-snippet
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

app = FastAPI(title='PGx API')

class PredictInput(BaseModel):
    gene: str
    rsids: list[str]

@app.post('/pgx/predict')
async def predict(payload: PredictInput):
    try:
        alleles = call_star_alleles(set(payload.rsids), gene=payload.gene)
        dip = determine_diplotype(alleles[0].name, alleles[1].name, gene=payload.gene)
        ph  = classify_phenotype(dip, gene=payload.gene)
        return ph
    except KeyError:
        raise HTTPException(404, f'Gene {payload.gene} not supported')

@app.post('/pgx/report')
async def report(patient: dict):
    rpt = generate_clinical_report(patient)
    return export_report(rpt, format='json')
```

Deploy behind API gateway; enable caching (Redis) on `/pgx/predict` keyed by sorted
rsID list.

## Database schema (PostgreSQL example)

```sql
CREATE TABLE patient_pgx (
  patient_id TEXT NOT NULL,
  gene TEXT NOT NULL,
  diplotype TEXT,
  phenotype TEXT,
  activity_score REAL,
  generated_at TIMESTAMPTZ DEFAULT NOW(),
  report_id UUID,
  PRIMARY KEY (patient_id, gene)
);
CREATE INDEX idx_pgx_gene ON patient_pgx(gene);
CREATE INDEX idx_pgx_patient ON patient_pgx(patient_id);
```

Bulk insert from Pandas:

```python-snippet
import psycopg2, pandas as pd
df = pd.read_csv('cohort_phenotypes.tsv', sep='	')
conn = psycopg2.connect(dsn)
cur = conn.cursor()
cur.executemany(
    'INSERT INTO patient_pgx(patient_id,gene,diplotype,phenotype,activity_score) VALUES (%s,%s,%s,%s,%s)',
    df[['patient','gene','diplotype','phenotype','activity_score']].itertuples(index=False, name=None)
)
conn.commit()
```

## Snakemake / Nextflow pipeline embedding

```snakemake
rule pgx_phenotype:
    input:
        vcf='data/{sample}.vcf.gz',
        gene_list='config/genes.txt'
    output:
        pheno='results/{sample}_pgx.tsv',
        report='reports/{sample}_pgx.pdf'
    threads: 4
    run:
        import vcfpy, pandas as pd
        from metainformant.pharmacogenomics import predict_metabolizer
        # parse VCF → dict{sample:{gene: set(rsids)}}
        # loop genes, call predict_metabolizer, write TSV + PDF report
```

## Interoperability with other tools

| Tool | Bridge method |
|------|---------------|
| `PharmCAT` | Export JSON via `report.to_pharmcat()` (roadmap) |
| `OneOme` | Read phenotype table CSV → OneOme compatible columns |
| `Genelex` | Not directly compatible; use phenotype table as intermediate |

## Auditing & provenance

Every report embeds immutable fields:

```json
{
  "report_id": "a1b2c3d4-…",
  "generated_at": "2025-06-14T10:22:15Z",
  "metainformant_version": "0.8.0",
  "cpic_version": "3.0",
  "pharmgkb_cache_date": "2025-05-15"
}
```

Store alongside patient ID in your audit log for CLIA compliance.

## Security & privacy

- All computation local; no PII leaves the host unless explicitly sent.
- PharmGKB download over HTTPS; cache permissions 0600.
- Reports may contain PHI (name, DOB) — protect per HIPAA / GDPR / local law.
