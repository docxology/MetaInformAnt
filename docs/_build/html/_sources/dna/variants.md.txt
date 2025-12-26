### DNA: Variants

Function: `parse_vcf`

```mermaid
flowchart LR
  A[VCF File] --> B[parse_vcf] --> C[Sample & Variant Counts]
```

Example

```python
from metainformant.dna import variants

# Parse VCF file for basic statistics
vcf_info = variants.parse_vcf("variants.vcf")
print(f"Samples: {vcf_info['samples']}")
print(f"Variants: {vcf_info['num_variants']}")
```

Features:
- **Basic parsing**: Extracts sample names and variant counts from VCF files
- **Header processing**: Reads VCF header to identify samples
- **Variant counting**: Counts total number of variant records
- **Simple interface**: Returns dictionary with key statistics

VCF Format Support:
- **Standard VCF**: Variant Call Format specification compliance
- **Sample extraction**: Identifies samples from header #CHROM line
- **Variant records**: Counts lines that aren't headers or comments

Note: This is a lightweight parser for basic statistics. For comprehensive VCF analysis, consider specialized libraries like pyvcf or cyvcf2.

Related: Part of variant analysis workflow, complements [population genetics](./population.md) analyses.
