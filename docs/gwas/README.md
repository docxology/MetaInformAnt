# GWAS Module

Genome-Wide Association Studies analysis, fine-mapping, and variant-expression integration.

## 📊 Architecture

```mermaid
graph TD
    subgraph "GWAS Module"
        A[analysis/] --> |association.py| AS[Association Testing]
        A --> |mixed_model.py| LM[EMMA Mixed Models]
        A --> |quality.py| QC[VCF Parsing & QC]
        A --> |structure.py| ST[PCA & Kinship]

        F[finemapping/] --> |credible_sets.py| CS[Credible Sets & SuSiE]
        F --> |colocalization.py| CO[Colocalization]
        F --> |eqtl.py| EQ[eQTL Analysis]

        D[data/] --> |download.py| DL[Reference Genome & SRA Downloads]
        D --> |genome.py| GN[Chromosome Mapping & GFF3]
        D --> |metadata.py| MD[Sample Metadata]

        H[heritability/] --> |estimation.py| HE[LDSC, GREML, HE Regression]

        V[visualization/] --> VIS[Manhattan, QQ, Interactive Plots]

        W[workflow/] --> WF[Pipeline Orchestration]
    end
```

## 🔑 Key Capabilities

### Fine-Mapping & Colocalization

```python
from metainformant.gwas.finemapping.colocalization import (
    eqtl_coloc,        # GWAS-eQTL colocalization
    multi_trait_coloc, # Multi-trait analysis
    compute_clpp,      # CLPP method
    regional_coloc,    # Regional analysis
)

# Test if GWAS signal colocalizes with expression QTL
result = eqtl_coloc(
    gwas_z=[1.2, 2.5, 3.1],  # GWAS Z-scores per variant
    eqtl_z=[1.1, 2.3, 2.9],  # eQTL Z-scores
    gene_id="LOC12345"
)
# Returns: PP_H4 (shared causal), interpretation, credible sets
```

### Association Testing

| Function | Purpose |
|----------|---------|
| `run_gwas()` | Full GWAS pipeline |
| `logistic_regression()` | Case-control associations |
| `linear_regression()` | Quantitative trait analysis |
| `lmm_association()` | Mixed models for population structure |

### Data I/O

| Function | Purpose |
|----------|---------|
| `read_vcf()` | Parse VCF variant files |
| `read_plink()` | Load PLINK bed/bim/fam |
| `write_sumstats()` | Export summary statistics |

## 📦 Submodules

| Module | Purpose |
|--------|---------|
| [`analysis/`](analysis/) | Association testing |
| [`data/`](data/) | VCF, PLINK I/O |
| [`finemapping/`](finemapping/) | Fine-mapping, colocalization |
| [`heritability/`](heritability/) | Heritability estimation |
| [`visualization/`](visualization/) | Manhattan, QQ, LocusZoom |
| [`workflow/`](workflow/) | Pipeline orchestration |

## 🧬 Integration with Expression Data

The GWAS module integrates with RNA-seq via the multiomics module:

```python
from metainformant.multiomics.analysis import integration
from metainformant.gwas.finemapping.colocalization import eqtl_coloc

# Convert GWAS and expression data to common format
dna_data = integration.from_dna_variants(gwas_sumstats)
rna_data = integration.from_rna_expression(expression_matrix)

# Run integrated analysis
integrated = integration.integrate_omics_data(
    dna_data=dna_data,
    rna_data=rna_data,
)
```

## 🔗 Related

- [metainformant.multiomics](../multiomics/) - Multi-omic integration
- [metainformant.rna](../rna/) - RNA-seq analysis
- [config/gwas/](../../../config/gwas/) - Configuration files

## Apis mellifera Real Data

- [Real genome and read organization](../../projects/apis_gwas/doc/current/amellifera_real_genome_reads.md)
- [Raw reads to genomic estimators runbook](../../projects/apis_gwas/doc/current/beewas_raw_reads_to_genomic_estimators.md)
- [BeeWAS synthetic-phenotype validation](../../projects/apis_gwas/doc/current/beewas_synthetic_validation.md)
- Real-cohort QC reporter: `scripts/gwas/pipelines/analyze_beewas_2026_real.py`
- Genomic-estimator validator: `scripts/gwas/pipelines/validate_beewas_genomic_estimators.py`
- Shared BeeWAS reporting helpers: `scripts/gwas/pipelines/beewas_reporting.py`
- Real guarded config: `projects/apis_gwas/config/beewas_2026_full_guarded.yaml`
- Generated/demo config: `config/gwas/gwas_amellifera.yaml`
