# Variation

VCF parsing and variant analysis including variant calling from pileups, mutation classification, and variant effect prediction.

## Contents

| File | Purpose |
|------|---------|
| `calling.py` | Variant calling from pileup data, genotyping, and filtering |
| `mutations.py` | Mutation rate calculation, substitution matrices, and hotspot detection |
| `variants.py` | VCF parsing/writing, allele frequency, Ti/Tv ratio, and effect prediction |

## Key Functions

| Function | Description |
|----------|-------------|
| `parse_vcf()` | Parse a VCF file into a structured dict with header and variants |
| `write_vcf()` | Write variant data back to VCF format |
| `filter_variants_by_quality()` | Remove variants below a QUAL threshold |
| `filter_variants_by_maf()` | Keep variants above a minor allele frequency cutoff |
| `calculate_variant_statistics()` | Counts by type (SNP, indel), Ti/Tv ratio, per-chromosome summary |
| `detect_variant_type()` | Classify a variant as SNP, insertion, deletion, or MNP |
| `predict_variant_effect()` | Predict coding impact (synonymous, missense, nonsense) |
| `calculate_ti_tv_ratio()` | Transition-to-transversion ratio across a VCF dataset |
| `call_variants_pileup()` | Call variants from read pileup data |
| `genotype_variants()` | Assign genotypes to called variant sites |
| `calculate_mutation_rate()` | Substitution rate between ancestral and derived sequences |
| `classify_mutations()` | Count transitions, transversions, and mutation types |
| `simulate_sequence_evolution()` | Generate evolved sequences over multiple generations |
| `analyze_mutation_spectrum()` | Summarize mutation types between a sequence and reference |

## Usage

```python
from metainformant.dna.variation.variants import parse_vcf, filter_variants_by_quality
from metainformant.dna.variation.mutations import classify_mutations, calculate_mutation_rate
from metainformant.dna.variation.calling import call_variants_pileup

vcf_data = parse_vcf("variants.vcf")
filtered = filter_variants_by_quality(vcf_data, min_qual=30.0)
mutations = classify_mutations("ATCGATCG", "ATCAATCG")
```
