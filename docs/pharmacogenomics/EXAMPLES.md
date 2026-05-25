# Examples: Pharmacogenomics

### 1. Single-patient end-to-end

```python
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation

variants = {'rs1065852', 'rs3892097', 'rs5030655'}
alleles = call_star_alleles(variants, gene='CYP2D6')
dip = determine_diplotype(alleles[0].name, alleles[1].name, gene='CYP2D6')
phenotype = classify_phenotype(dip, gene='CYP2D6')
print(phenotype['phenotype'].value)  # Intermediate Metabolizer

rec = get_dosing_recommendation('codeine', phenotype['phenotype'].value, gene='CYP2D6')
print(rec['recommendation'])
```

### 2. Batch processing of a large cohort

```python-snippet
import pandas as pd
from metainformant import parallel, pharmacogenomics

df = pd.read_csv('cohort.csv')  # column 'rsids' = ';'-separated rsIDs per patient

def predict_metabolizer(row):
    observed = set(str(row['rsids']).split(';')) if pd.notna(row['rsids']) else set()
    return pharmacogenomics.predict_metabolizer(observed, gene='CYP2D6',
                                                algorithm='balanced')

results = parallel.map(predict_metabolizer,
                        df.to_dict('records'),
                        workers=8,
                        chunk_size=10000,
                        desc='CYP2D6 metabolizer')
df['cyp2d6_phenotype'] = [r['phenotype'] for r in results]
df.to_csv('cohort_with_phenotypes.tsv', sep='\t', index=False)
```

> **Notes** — `chunk_size` dictates memory; `workers` should match CPU cores; consider `parallel.set_backend('loky')` for reliability.

### 3. Cross-gene CPIC actionable lookup

```python
from metainformant.pharmacogenomics.annotations.cpic import lookup_drug_gene, get_dosing_recommendation

actionable = lookup_drug_gene(level='A')  # Level A = 'strong' – action required
print(f"CPIC actionable drug–gene pairs: {len(actionable)}")

for entry in actionable:
    drug, gene = entry['drug'], entry['gene']
    phenotypes = entry['phenotypes']
    # show a sample recommendation
    example = list(phenotypes.values())[0]
    print(f"  {drug} ({gene}): {example['recommendation'][:60]}…")
```

### 4. Drug–drug and drug–gene interaction analysis

```python-snippet
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    analyze_drug_gene_interactions,
    polypharmacy_risk,
    cyp_inhibition_prediction,
)

meds = ['warfarin', 'amiodarone', 'simvastatin', 'fluconazole']
# PM on CYP2C9 — high bleeding risk
interactions = analyze_drug_gene_interactions(
    gene='CYP2C9', phenotype='Poor Metabolizer', drugs=meds)
for ix in interactions:
    print(f"{ix['drug']}: {ix['recommendation']} (severity {ix['severity']})")

# Composite polypharmacy score
risk = polypharmacy_risk(meds,
                         patient_phenotype='Normal Metabolizer',
                         gene='CYP2C9')
print(f"Aggregate risk = {risk['risk_score']:.2f} ({risk['high_severity_count']} major)")

# CYP inhibition net score
score = cyp_inhibition_prediction(['ketoconazole', 'rifampin'], enzyme='CYP3A4')
print(f"Net CYP3A4 inhibition category: {score['category']}")
```

### 5. ACMG classification of a pharmacogene variant

```python
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    classify_variant_acmg,
    apply_acmg_criteria,
    ACMGClassification,
    ACMGCriteria,
)

variant = {
    'rsid': 'rs4244285',
    'gene': 'CYP2C19',
    'gnomad_af': 0.0002,
    'clinvar_sig': 'Pathogenic',
    'sift': 0.0,
    'polyphen': 1.0,
    # optional extra evidence
    'cadd_score': 25.0,
    'provean_score': -5.0,
}

result = classify_variant_acmg(
    variant,
    gene='CYP2C19',
    criteria_applied=[ACMGCriteria.PVS1, ACMGCriteria.PM2, ACMGCriteria.PP3],
)
print(f"Classification: {result['classification'].value}")
print('Evidence:', result['evidence_summary'])
```

### 6. Full clinical report generation (CLIA-ready structure)

```python-snippet
from metainformant.pharmacogenomics.clinical.reporting import (
    generate_clinical_report,
    export_report,
    ClinicalReportBuilder,
)

patient = {
    'patient_id': 'PT001',
    'name': 'Sample Patient',
    'dob': '1970-01-01',
    'sex': 'F',
    'medications': ['codeine', 'clopidogrel', 'atorvastatin'],
    'genotype': {
        'CYP2D6': {'diplotype': '*1/*4', 'phenotype': 'Intermediate Metabolizer'},
        'CYP2C19': {'diplotype': '*1/*1', 'phenotype': 'Normal Metabolizer'},
        'CYP2C9':  {'diplotype': '*1/*3', 'phenotype': 'Intermediate Metabolizer'},
        'SLCO1B1': {'diplotype': '*1/*1', 'phenotype': 'Normal Function'},
    }
}

report = generate_clinical_report(patient)
print(export_report(report, format='text'))
```

Report sections: header → genotype summary per gene → drug recommendations by CPIC → DDI warnings → disclaimer & UUID.

### 7. Custom allele table loading for a non-coverage gene

```python
import json
from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions

cyp2e1 = {
    'gene': 'CYP2E1',
    'star_alleles': [
        {
            'name': '*1',
            'variants': ['rs2031920'],
            'function': 'normal',
            'activity': 1.0,
        },
        {
            'name': '*5B',
            'variants': ['rs3813867'],
            'function': 'decreased',
            'activity': 0.5,
        },
    ]
}
with open('/tmp/cyp2e1.json', 'w') as fh:
    json.dump(cyp2e1, fh, indent=2)

load_allele_definitions(gene='CYP2E1', filepath='/tmp/cyp2e1.json')
variants = {'rs2031920'}
alleles = call_star_alleles(variants, gene='CYP2E1')
print([a.name for a in alleles])
```

### 8. Population phenotype frequency tabulation

```python
from metainformant.pharmacogenomics.alleles.phenotype import population_phenotype_frequencies
import pandas as pd

df = pd.read_csv('large_pharmacogenomic_cohort.csv')

def parse_variants(row):
    raw = row.get('CYP2D6_rsids', '')
    return set(str(raw).split(';')) if pd.notna(raw) else set()

freqs = population_phenotype_frequencies(
    genotype_matrix=df,
    parse_variants=parse_variants,
    gene='CYP2D6',
    sample_weight_col=None,  # optional column for survey weights
)
print(freqs)  # {'PM':0.07, 'IM':0.34, 'NM':0.50, 'RM':0.07, 'UM':0.02}
```

### 9. Matplotlib bar-chart visualization of metabolizer distribution

```python
from metainformant.pharmacogenomics.visualization.plots import plot_metabolizer_status
counts = {'PM':70, 'IM':320, 'NM':500, 'RM':30, 'UM':10}
plot_metabolizer_status(
    counts,
    gene='CYP2D6',
    title='CYP2D6 phenotype distribution — N=1030',
    save='cyp2d6_distribution.png',
    show=False,
    style='darkgrid',
)
```

### 10. PharmGKB annotation enrichment

```python-snippet
from metainformant.pharmacogenomics.annotations.pharmgkb import (
    query_pharmgkb_annotations,
    get_evidence_level,
    get_drug_phenotypes,
)

ann = query_pharmgkb_annotations(gene='CYP2D6', drug='codeine')
for a in ann:
    lvl = get_evidence_level(a)
    print(f"{a['variant']} – {a['drug']} (Level {lvl['level']}): {a['annotation_text'][:80]}")

pheno_guides = get_drug_phenotypes('codeine', gene='CYP2D6')
print(pheno_guides)
```

### 11. Error-handling wrapper for production inference service

```python
def safe_predict(variants, gene, fallback='NM'):
    try:
        return pharmacogenomics.predict_metabolizer(variants, gene=gene)
    except ValueError as e:
        # Known failure modes: unknown gene, empty variant set, allele mismatch
        return {'phenotype': fallback, 'error': str(e), 'gene': gene}
    except KeyError as e:
        # allele activity missing
        return {'phenotype': fallback, 'error': f'activity table: {e}', 'gene': gene}

result = safe_predict({'rsxxxx'}, 'CYP2E1')
if result.get('error'):
    log.error('PGx failed for gene %s: %s', result['gene'], result['error'])
```

### 12. Multi-gene and multi-drug recommendation matrix

```python
import pandas as pd

genes  = ['CYP2D6','CYP2C19','CYP2C9','SLCO1B1','UGT1A1']
drugs  = ['codeine','tamoxifen','warfarin','atorvastatin','irinotecan']

matrix = {}
for drug in drugs:
    matrix[drug] = {}
    for gene in genes:
        ph = predict_metabolizer(patient_variants[gene], gene=gene)['phenotype']
        rec = get_dosing_recommendation(drug, ph, gene=gene)
        matrix[drug][gene] = rec['recommendation']

df = pd.DataFrame(matrix).T  # drugs × genes
df.to_csv('drug_gene_matrix.csv')
```

### 13. Parallel MapReduce via Dask for >1M records

```python-snippet
import dask.dataframe as dd
from dask.distributed import Client
client = Client(n_workers=8, threads_per_worker=1)

ddf = dd.read_csv('massive_cohort_*.csv')
def predict_partition(df):
    out = []
    for _, row in df.iterrows():
        out.append(pharmacogenomics.predict_metabolizer(
            set(str(row['rsids']).split(';')), gene='CYP2D6'))
    return pd.Series([o['phenotype'] for o in out], index=df.index)

ddf['phenotype'] = ddf.map_partitions(predict_partition, meta=('phenotype', 'object'))
ddf.to_csv('output/phenotypes-*.csv', index=False)
```

### 14. Build a dosing decision-tool CLI

```python
import argparse

parser = argparse.ArgumentParser(description='Drug–gene dosing lookup')
parser.add_argument('--drug', required=True)
parser.add_argument('--gene', required=True)
parser.add_argument('--phenotype', required=True)
parser.add_argument('--format', default='text', choices=['text','json'])
args = parser.parse_args()

rec = get_dosing_recommendation(args.drug, args.phenotype, args.gene)
if args.format == 'json':
    import json; print(json.dumps(rec, indent=2))
else:
    print(f"{args.drug} ({args.gene}, {args.phenotype}): {rec['recommendation']}")
```

### 15. PharmGKB annotation caching across runs

```python-snippet
from metainformant.pharmacogenomics.annotations.pharmgkb import get_cached_pharmgkb_path
import os

cache_dir = get_cached_pharmgkb_path()  # ~/.cache/metainformant/pharmgkb/
print(f"PharmGKB cache lives at {cache_dir}")
# For offline machines, pre-populate it before deployment:
#   download latest annotations.tsv manually into that folder
```

### 16. Debug logging for allele mismatch investigation

```python
import os
os.environ['PG_LOGGING'] = 'DEBUG'   # must be set before importing module
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles

variants = {'rs1065852', 'rs5030655', 'rs3892097'}
alleles = call_star_alleles(variants, gene='CYP2D6')
# console shows: attempted matches, score, why winner chosen
```

### 17. Activity score table export for regulatory submission

```python-snippet
from metainformant.pharmacogenomics.alleles.activity import _ACTIVITY_SCORE_TABLES
import json
gene = 'CYP2D6'
table = _ACTIVITY_SCORE_TABLES[gene]
with open(f'{gene}_activity_scores.json', 'w') as fh:
    json.dump(table, fh, indent=2)
print(f"Exported {len(table)} allele entries")
```

### 18. Multi-threaded phenotype classification pipeline

```python
from concurrent.futures import ThreadPoolExecutor

def pheno(row):
    alleles = call_star_alleles(set(row['rsids'].split(';')), gene='CYP2D6', algorithm='balanced')
    dip = determine_diplotype(alleles[0].name, alleles[1].name, gene='CYP2D6')
    return classify_phenotype(dip, gene='CYP2D6')['phenotype'].value

with ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(pheno, df.to_dict('records')))
df['phenotype'] = results
```

### 19. HTML clinical report generation (patient-friendly)

```python
report = generate_clinical_report(patient)
html = export_report(report, format='html')
with open('report.html', 'w') as fh:
    fh.write(html)
# Open in browser or embed in patient portal
```

### 20. Integration into FastAPI micro-service

```python-snippet
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

app = FastAPI(title='Pharmacogenomics API')

class GenotypeInput(BaseModel):
    gene: str
    rsids: list[str]

@app.post('/predict/{gene}')
async def predict_endpoint(gene: str, payload: GenotypeInput):
    try:
        alleles = call_star_alleles(set(payload.rsids), gene=gene)
        dip = determine_diplotype(alleles[0].name, alleles[1].name, gene=gene)
        pheno = classify_phenotype(dip, gene=gene)
        return pheno
    except KeyError:
        raise HTTPException(404, f'Gene {gene} not supported')

if __name__ == '__main__':
    import uvicorn; uvicorn.run(app, host='0.0.0.0', port=8080)
```

### 21. Prescription-safety check list

```python
def check_patient_meds(patient_genotype: dict, meds: list[str]) -> list[dict]:
    findings = []
    for med in meds:
        for gene, geno in patient_genotype.items():
            ph = geno['phenotype']
            rec = get_dosing_recommendation(med, ph, gene=gene)
            if rec.get('action') != 'no_action':
                findings.append({
                    'drug': med,
                    'gene': gene,
                    'phenotype': ph,
                    'recommendation': rec['recommendation'],
                    'level': rec.get('cpic_level'),
                })
    return findings

findings = check_patient_meds(patient['genotype'], patient['medications'])
for f in findings:
    print(f"[ALERT] {f['drug']} — {f['gene']} ({f['phenotype']}): {f['recommendation']}")
```

### 22. Export phenotype to VCF INFO field (for downstream genomics pipelines)

```python-snippet
import vcfpy  # pip install vcfpy

reader = vcf.Reader(open('cohort.vcf.gz'))
writer = vcf.Writer(open('phenotypes.vcf.gz','wb'), reader.header)

for record in reader:
    gene = record.INFO.get('gene', ['NA'])[0]
    if gene in ['CYP2D6','CYP2C19']:
        rsids = {sid.split('rs')[1] for sid in record.ID if sid.startswith('rs')}
        ph = pharmacogenomics.predict_metabolizer(rsids, gene=gene)['phenotype'].value
        record.INFO['PGX_PHENO'] = ph
    writer.write_record(record)
```
