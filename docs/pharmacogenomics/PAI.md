# Pre‑Implementation Guide: Pharmacogenomics

## Problem statement

Translate raw genotypes (variants) into drug‑dosing guidance, following CPIC
standards. Requires: variant→star‑allele mapping, diplotype formation, activity
scoring, phenotype classification, drug‑specific recommendation lookup.

## Prior art & gaps

| Solution | Coverage | Limitation |
|----------|----------|------------|
| `pgx-js` (Node) | allele lookup only | no CPIC, slower |
| R `pharmvar` | definitions only | no phenotype→recommend |
| Commercial CDS (OneOme) | end‑to‑end | closed, expensive, update lag |
| **metainformant** | **full pipeline** | open, fast (<10 ms), reproducible |

## Data scale

| Scale | Inputs | Strategy |
|-------|--------|----------|
| 1 patient | 10–40 rsIDs | direct call (<10 ms) |
| 1–10 k cohort | 10–40 rsIDs each | `parallel.map` workers 4–8 |
| 100 k–1 M | same | Dask distributed + Parquet output |

Memory: static tables ~4 MB shared; per‑sample transient ~2 KB.

## User journeys

- **Clinician**: upload VCF → get PDF/JSON with dosing alerts → import to EHR.
- **Bioinformatician**: batch 100 k genomes → phenotype table → downstream GWAS of drug response.

## Non‑goals

- Phasing from BAM (assume already phased/imputed).
- Somatic pharmacogenes (future oncology module may cover).
- Sub‑ms realtime API (typical 7–12 ms; fine for batch/clinic-day-ahead).

## Performance highlights

| Function | Latency (median) |
|----------|------------------|
| `call_star_alleles()` | 6.2 ms |
| `determine_diplotype()` | 0.3 ms |
| `classify_phenotype()` | 0.08 ms |
| `predict_metabolizer()` | 7.1 ms |
| CPIC lookup | 0.2 ms |



