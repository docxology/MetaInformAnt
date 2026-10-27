# Examples: Pharmacogenomics

### 1. Single-patient end-to-end

```python
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation

variants = {'rs3892097', 'rs1065852'}  # CYP2D6 *4 and *10
alleles = call_star_alleles(variants, gene='CYP2D6')  # [*10, *4]: every match, sorted by name
dip = determine_diplotype(alleles[0], alleles[1], gene='CYP2D6')
print(dip.diplotype_string)  # *10/*4 (alleles are ordered as strings)
phenotype = classify_phenotype(dip, gene='CYP2D6')
print(phenotype['phenotype'].value)  # Intermediate Metabolizer

rec = get_dosing_recommendation('codeine', phenotype['phenotype'].value)
print(rec['recommendation'])  # Use codeine per standard of care. Monitor for reduced efficacy.
```

`call_star_alleles` returns every allele it matched, so a variant set can give one allele (or just `*1`), or more than two; `get_dosing_recommendation(drug, phenotype)` takes no gene and returns `None` when no guideline matches. Examples 2, 11–13, 18, 20 and 22 reuse this one-call helper, which pads with `*1` and uses the first two alleles:

```python
def phenotype_from_variants(variants, gene):
    names = [a.name for a in call_star_alleles(variants, gene)] + ['*1', '*1']
    return classify_phenotype(determine_diplotype(names[0], names[1], gene), gene)
```

### 2. Batch processing of a large cohort

```python-snippet
import pandas as pd
from metainformant.core.execution import parallel

df = pd.read_csv('cohort.csv')  # column 'rsids' = ';'-separated rsIDs per patient

def cyp2d6_phenotype(rsids):
    if pd.isna(rsids):
        return None  # missing genotype data is not "no variants" (that would read as *1/*1)
    return phenotype_from_variants(set(str(rsids).split(';')), 'CYP2D6')['phenotype'].value

df['cyp2d6_phenotype'] = parallel.thread_map(cyp2d6_phenotype, df['rsids'], max_workers=8)
df.to_csv('cohort_with_phenotypes.tsv', sep='\t', index=False)
```

> **Notes** — `thread_map(func, items, max_workers=8)` runs a thread pool and has no `desc` argument; the matching is plain Python, so for CPU-bound cohorts use `parallel.process_map(func, items, max_workers=...)` with a top-level function.

### 3. Cross-gene CPIC actionable lookup

```python
from metainformant.pharmacogenomics.annotations.cpic import list_actionable_genes, lookup_drug_gene

actionable = list_actionable_genes(min_level='A')  # [{'gene', 'drug', 'cpic_level'}, ...]
print(f"CPIC actionable drug–gene pairs: {len(actionable)}")

for entry in actionable:
    guideline = lookup_drug_gene(entry['drug'], entry['gene'])
    # show a sample recommendation (the first listed phenotype)
    example = next(iter(guideline['recommendations'].values()))
    print(f"  {entry['drug']} ({entry['gene']}): {example['recommendation'][:60]}…")
```

### 4. Drug–drug and drug–gene interaction analysis

```python-snippet
from metainformant.pharmacogenomics.clinical.drug_interaction import analyze_drug_gene_interactions
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    polypharmacy_risk,
    cyp_inhibition_prediction,
)

meds = ['warfarin', 'amiodarone', 'simvastatin', 'fluconazole']
# PM on CYP2C9 — high bleeding risk; only drug–gene pairs with a built-in CPIC
# guideline (here warfarin/CYP2C9) produce an entry
interactions = analyze_drug_gene_interactions(
    meds, {'CYP2C9': {'phenotype': 'Poor Metabolizer'}})
for ix in interactions:
    print(f"{ix.drug}: {ix.recommendation} (severity {ix.severity.value})")

# Composite polypharmacy score
risk = polypharmacy_risk(meds)
print(f"Aggregate risk = {risk['risk_score']:.2f} ({risk['severity_counts']['major']} major)")

# CYP inhibition/induction profile of a single drug
profile = cyp_inhibition_prediction('ketoconazole')
print(profile['summary'])
```

### 5. ACMG classification of a pharmacogene variant

```python
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    classify_variant_acmg,
    apply_acmg_criteria,
    ACMGClassification,
)

# Illustrative evidence, not a real gnomAD/ClinVar record
variant = {
    'rsid': 'rs0000000',
    'gene': 'CYP2C19',
    'consequence': 'frameshift',
    'allele_frequency': 0.00002,  # a float, or {population: frequency}
    'computational_predictions': {'sift': 'deleterious', 'polyphen': 'probably_damaging', 'cadd': 'damaging'},
}

criteria = apply_acmg_criteria(variant)  # {criterion: bool} for every criterion; here PVS1, PM2, PP3 are True
result = classify_variant_acmg(variant)  # evaluates the criteria from the variant data
print(f"Classification: {result['classification'].value}")  # Pathogenic
print('Criteria met:', result['criteria_met'])  # ['PM2', 'PP3', 'PVS1']

# Or supply pre-evaluated criteria ({criterion: bool}) instead of auto-evaluating
result = classify_variant_acmg(variant, evidence={'PVS1': True, 'PM2': True})
print(result['classification'] is ACMGClassification.LIKELY_PATHOGENIC)  # True
```

### 6. Full clinical report generation

```python-snippet
from metainformant.pharmacogenomics.clinical.reporting import (
    generate_clinical_report,
    export_report,
)

patient = {
    'patient_id': 'PT001',
    'name': 'Sample Patient',
    'dob': '1970-01-01',
    'sex': 'F',
}
genotypes = {
    'CYP2D6': {'diplotype': '*4/*10', 'phenotype': 'Intermediate Metabolizer'},
    'CYP2C19': {'diplotype': '*2/*2', 'phenotype': 'Poor Metabolizer'},
    'CYP2C9':  {'diplotype': '*1/*3', 'phenotype': 'Intermediate Metabolizer'},
    'SLCO1B1': {'diplotype': '*1/*1', 'phenotype': 'Normal Metabolizer'},
}
medications = ['codeine', 'clopidogrel', 'simvastatin']

report = generate_clinical_report(patient, genotypes, drugs=medications)
print(export_report(report, format='text'))
```

The report is a dict with the keys `header`, `patient_info`, `genotype_results`, `drug_recommendations`, `interaction_summary`, `clinical_actions` and `disclaimer`; the text export prints them as header → patient information → genotype results per gene → drug recommendations (CPIC level) → clinical actions → interaction summary → disclaimer. `patient_info` keeps only `patient_id`, `sex`, `ethnicity` and `ordering_provider`, so `name` and `dob` do not appear. Only drug–gene pairs with a built-in CPIC guideline produce a recommendation (here codeine/CYP2D6, clopidogrel/CYP2C19 and simvastatin/SLCO1B1). Drug recommendations use the `phenotype` you supply in `genotypes`, while the genotype table recomputes the phenotype from `diplotype`, so keep the two consistent (for example by filling `phenotype` from `classify_phenotype`).

### 7. Custom allele table loading for a non-coverage gene

```python
import json
from pathlib import Path
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles, load_allele_definitions

cyp2e1 = {
    '*1': {
        'defining_variants': [],
        'function': 'Normal function',
        'activity_value': 1.0,
    },
    '*5B': {
        'defining_variants': ['rs3813867'],
        'function': 'Decreased function',
        'activity_value': 0.5,
    },
}
path = Path('output/pharmacogenomics/cyp2e1.json')
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(cyp2e1, indent=2))

definitions = load_allele_definitions('CYP2E1', source='file', filepath=path)
variants = {'rs3813867'}
alleles = call_star_alleles(variants, 'CYP2E1', allele_definitions=definitions)
print([a.name for a in alleles])  # ['*5B']
```

### 8. Population phenotype frequency tabulation

```python
from metainformant.pharmacogenomics.alleles.phenotype import population_phenotype_frequencies

# allele name -> population frequency (illustrative values; they should sum to about 1)
allele_freqs = {'*1': 0.35, '*2': 0.30, '*4': 0.20, '*10': 0.08, '*41': 0.07}

freqs = population_phenotype_frequencies(allele_freqs, gene='CYP2D6')
print(freqs)  # phenotype abbreviation -> expected frequency under Hardy–Weinberg (empty categories omitted)
```

For observed frequencies in a cohort, classify each sample (example 2) and count: `df['cyp2d6_phenotype'].value_counts(normalize=True)`.

### 9. Matplotlib bar-chart visualization of metabolizer distribution

```python
from metainformant.pharmacogenomics.visualization.plots import plot_metabolizer_status

counts = {'PM': 70, 'IM': 320, 'NM': 500, 'RM': 30, 'UM': 10}
path = plot_metabolizer_status(
    counts,
    'output/pharmacogenomics/cyp2d6_distribution.png',
    title=f'CYP2D6 phenotype distribution — N={sum(counts.values())}',
    figsize=(10, 6),
)
print(path)  # the saved file (a pathlib.Path); the figure is written to disk, never shown
```

### 10. PharmGKB annotation enrichment

```python-snippet
from metainformant.pharmacogenomics.annotations.cpic import lookup_drug_gene
from metainformant.pharmacogenomics.annotations.pharmgkb import (
    query_pharmgkb_annotations,
    get_evidence_level,
)

ann = query_pharmgkb_annotations(gene='CYP2D6', drug='codeine')
for a in ann:
    lvl = get_evidence_level(a)
    print(f"{a['gene']} – {a['drug']} (Level {lvl['level']}): {a['annotation_text'][:80]}")

pheno_guides = lookup_drug_gene('codeine', 'CYP2D6')['recommendations']
print(pheno_guides)
```

### 11. Error-handling wrapper for production inference service

```python
import logging

log = logging.getLogger(__name__)

def safe_phenotype(variants, gene, fallback='Indeterminate'):
    try:
        return phenotype_from_variants(variants, gene)
    except ValueError as e:
        # unknown gene: there is no built-in allele table for it
        return {'phenotype': fallback, 'error': str(e), 'gene': gene}  # never fall back to a normal phenotype

result = safe_phenotype({'rs3892097'}, 'CYP2E1')
if result.get('error'):
    log.error('PGx failed for gene %s: %s', result['gene'], result['error'])
```

### 12. Multi-gene and multi-drug recommendation matrix

```python
import pandas as pd
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation, lookup_drug_gene

genes = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'SLCO1B1', 'TPMT']
drugs = ['codeine', 'clopidogrel', 'warfarin', 'simvastatin', 'azathioprine']
patient_variants = {  # observed variants per gene
    'CYP2D6': {'rs3892097', 'rs1065852'}, 'CYP2C19': {'rs4244285', 'rs4986893'},
    'CYP2C9': {'rs1057910'}, 'SLCO1B1': set(), 'TPMT': set(),
}

matrix = {}
for drug in drugs:
    matrix[drug] = {}
    for gene in genes:
        guideline = lookup_drug_gene(drug, gene)
        if guideline is None:
            continue  # no built-in CPIC guideline for this drug–gene pair
        phenotype = phenotype_from_variants(patient_variants[gene], gene)['phenotype'].value
        rec = get_dosing_recommendation(drug, phenotype, guidelines=[guideline])
        matrix[drug][gene] = rec['recommendation'] if rec else f'no guidance for {phenotype}'

df = pd.DataFrame(matrix).T  # drugs × genes; NaN where there is no guideline
df.to_csv('drug_gene_matrix.csv')
```

### 13. Parallel MapReduce via Dask for >1M records

```python-snippet
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client

client = Client(n_workers=8, threads_per_worker=1)

ddf = dd.read_csv('massive_cohort_*.csv')

def predict_partition(df):
    return pd.Series(
        [None if pd.isna(rsids)  # missing genotype data, not "no variants"
         else phenotype_from_variants(set(rsids.split(';')), 'CYP2D6')['phenotype'].value
         for rsids in df['rsids']],
        index=df.index,
    )

ddf['phenotype'] = ddf.map_partitions(predict_partition, meta=('phenotype', 'object'))
ddf.to_csv('output/phenotypes-*.csv', index=False)
```

### 14. Build a dosing decision-tool CLI

```python
import argparse
import json

from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation, lookup_drug_gene

parser = argparse.ArgumentParser(description='Drug–gene dosing lookup')
parser.add_argument('--drug', required=True)
parser.add_argument('--gene', required=True)
parser.add_argument('--phenotype', required=True)  # e.g. "Poor Metabolizer" or PM
parser.add_argument('--format', default='text', choices=['text', 'json'])
args = parser.parse_args()

guideline = lookup_drug_gene(args.drug, args.gene)  # get_dosing_recommendation has no gene filter
rec = get_dosing_recommendation(args.drug, args.phenotype, guidelines=[guideline]) if guideline else None
if rec is None:
    raise SystemExit(f'No built-in CPIC recommendation for {args.drug}/{args.gene} ({args.phenotype})')
if args.format == 'json':
    print(json.dumps(rec, indent=2))
else:
    print(f"{args.drug} ({args.gene}, {args.phenotype}): {rec['recommendation']}")
```

### 15. PharmGKB annotations from a local export

```python-snippet
from metainformant.pharmacogenomics.annotations.pharmgkb import parse_clinical_annotations

# The package does not download or cache PharmGKB data; parse an export you downloaded yourself
annotations = parse_clinical_annotations('clinical_annotations.tsv')
print(f"Parsed {len(annotations)} clinical annotations")
```

### 16. Debug logging for allele mismatch investigation

```python
import logging
logging.getLogger('metainformant.pharmacogenomics').setLevel(logging.DEBUG)  # any time before the call
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles

variants = {'rs3892097', 'rs1065852'}
alleles = call_star_alleles(variants, gene='CYP2D6')
# console shows one DEBUG line per matched allele ("Matched allele *4 for CYP2D6 ...")

call_star_alleles({'rs0000000'}, gene='CYP2D6')
# console shows "No variant alleles matched for CYP2D6, returning reference *1"
```

### 17. Activity score table export for regulatory submission

```python-snippet
from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions
import json
gene = 'CYP2D6'
table = {a.name: a.activity_value for a in load_allele_definitions(gene)}
with open(f'{gene}_activity_scores.json', 'w') as fh:
    json.dump(table, fh, indent=2)
print(f"Exported {len(table)} allele entries")
```

### 18. Multi-threaded phenotype classification pipeline

```python
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

def pheno(row):
    if pd.isna(row['rsids']):
        return None  # missing genotype data, not "no variants"
    return phenotype_from_variants(set(row['rsids'].split(';')), 'CYP2D6')['phenotype'].value

with ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(pheno, df.to_dict('records')))  # df as in example 2
df['phenotype'] = results
```

### 19. HTML clinical report generation

```python
from metainformant.pharmacogenomics.clinical.reporting import generate_clinical_report, export_report

# patient, genotypes and medications as in example 6
report = generate_clinical_report(patient, genotypes, drugs=medications)
html = export_report(report, format='html', output_path='output/pharmacogenomics/report.html')
# the HTML string is returned and also written to output_path; open it in a browser or embed it in a portal
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
        return phenotype_from_variants(set(payload.rsids), gene)
    except ValueError:  # unknown gene: no built-in allele table
        raise HTTPException(404, f'Gene {gene} not supported')

if __name__ == '__main__':
    import uvicorn; uvicorn.run(app, host='0.0.0.0', port=8080)
```

### 21. Prescription-safety check list

```python
from metainformant.pharmacogenomics.clinical.drug_interaction import analyze_drug_gene_interactions

def check_patient_meds(genotypes: dict, meds: list[str]) -> list:
    """Drug–gene findings of Major or Moderate severity, i.e. those that call for action."""
    return [ix for ix in analyze_drug_gene_interactions(meds, genotypes) if ix.severity.requires_action]

findings = check_patient_meds(genotypes, medications)  # as in example 6
for f in findings:
    print(f"[ALERT] {f.drug} — {f.gene} ({f.phenotype}): {f.recommendation}")
```

### 22. Export phenotype to VCF INFO field (for downstream genomics pipelines)

```python-snippet
import vcfpy  # uv pip install vcfpy (not a metainformant dependency)

reader = vcfpy.Reader.from_path('cohort.vcf.gz')
reader.header.add_info_line({'ID': 'PGX_PHENO', 'Number': '1', 'Type': 'String',
                             'Description': 'Metabolizer phenotype from metainformant'})
writer = vcfpy.Writer.from_path('phenotypes.vcf.gz', reader.header)

for record in reader:
    gene = record.INFO.get('gene', ['NA'])[0]
    if gene in ['CYP2D6', 'CYP2C19']:
        rsids = {sid for sid in record.ID if sid.startswith('rs')}
        record.INFO['PGX_PHENO'] = phenotype_from_variants(rsids, gene)['phenotype'].value
    writer.write_record(record)
writer.close()
```
