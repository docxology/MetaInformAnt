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
| **metainformant** | **full pipeline** | open, microseconds per sample, reproducible |

## Data scale

| Scale | Inputs | Strategy |
|-------|--------|----------|
| 1 patient | 10–40 rsIDs | direct call (microseconds) |
| 1–10 k cohort | 10–40 rsIDs each | `parallel.thread_map` workers 4–8 |
| 100 k–1 M | same | Dask distributed + Parquet output |

Memory: the allele, activity and guideline tables are small Python literals shared by every sample.

## User journeys

- **Clinician**: observed variants → text, HTML or JSON report with dosing alerts → import to EHR.
- **Bioinformatician**: batch 100 k genomes → phenotype table → downstream GWAS of drug response.

## Non‑goals

- Phasing from BAM (assume already phased/imputed).
- Somatic pharmacogenes (future oncology module may cover).

## Performance highlights

| Function | Mean time per call |
|----------|------------------|
| `call_star_alleles()` | ~4 µs |
| `determine_diplotype()` | ~1.4 µs |
| `classify_phenotype()` | ~2.6 µs |
| `predict_metabolizer()` | ~5 µs |
| CPIC lookup | ~1 µs |

Measured as described in [PERFORMANCE.md](PERFORMANCE.md); absolute numbers depend on the hardware.


