### DNA: Accessions

Tools for validating and working with NCBI accession identifiers. NCBI uses specific
accession formats for genome assemblies, chromosomes, contigs, and scaffolds. The
`metainformant.dna.external.genomes` module provides validation functions that check
accession strings against these patterns.

---

## Accession Formats

NCBI uses standardized accession formats for different types of genomic data:

| Format Pattern             | Example                | Description                          |
|----------------------------|------------------------|--------------------------------------|
| `GCF_NNNNNNNNN.V`         | `GCF_000001405.39`     | RefSeq genome assembly               |
| `GCA_NNNNNNNNN.V`         | `GCA_000001405.28`     | GenBank genome assembly              |
| `NC_NNNNNN`               | `NC_000001`            | RefSeq chromosome                    |
| `NT_NNNNNN`               | `NT_000001`            | RefSeq contig                        |
| `NW_NNNNNNNN`             | `NW_001838827`         | RefSeq scaffold (WGS)               |
| `ASMNNNN`                 | `ASM2732`              | Older assembly format                |
| `chrN`                     | `chr1`                 | Chromosome shorthand                 |

The `GCF_` prefix denotes RefSeq (curated) assemblies, while `GCA_` denotes GenBank
(submitted) assemblies. The version suffix (`.39`) indicates the assembly version.

---

## Validation Functions

### `is_valid_assembly_accession`

Validates strictly against the NCBI assembly accession pattern (`GCF_` or `GCA_`
prefix followed by at least 9 digits, with an optional version suffix).

```python
from metainformant.dna.genomes import is_valid_assembly_accession

# Valid assembly accessions
is_valid_assembly_accession("GCF_000001405.39")  # True  (human GRCh38, RefSeq)
is_valid_assembly_accession("GCA_000001405")      # True  (GenBank, no version)
is_valid_assembly_accession("GCA_000001405.28")   # True  (GenBank, versioned)

# Invalid formats
is_valid_assembly_accession("NC_000001")           # False (chromosome, not assembly)
is_valid_assembly_accession("invalid")             # False
is_valid_assembly_accession("")                    # False
```

**Regex pattern used:** `^(GCF|GCA)_[0-9]{9,}(\.[0-9]+)?$`

### `validate_accession`

A broader validation function that recognizes multiple NCBI accession types including
assemblies, chromosomes, contigs, and scaffolds.

```python
from metainformant.dna.external.genomes import validate_accession

# Assembly accessions
validate_accession("GCF_000001405.39")  # True
validate_accession("GCA_000001405.39")  # True

# Chromosome and contig accessions
validate_accession("NC_000001")          # True
validate_accession("NT_000001")          # True
validate_accession("NW_001838827")       # True

# Chromosome shorthand
validate_accession("chr1")               # True

# Older assembly format
validate_accession("ASM2732")            # True

# Invalid
validate_accession("random_string")      # False
validate_accession("")                   # False
```

**Supported patterns:**

| Pattern                    | Matches                                          |
|----------------------------|--------------------------------------------------|
| `^GC[FA]_\d{9}\.\d+$`     | Versioned RefSeq/GenBank assembly accessions     |
| `^GCA_\d{9}\.\d+$`        | Versioned GenBank assembly accessions            |
| `^ASM\d+$`                | Older-format assembly identifiers                |
| `^chr\d+$`                | Chromosome shorthand (chr1, chr2, ...)           |
| `^NC_\d{6}$`              | RefSeq chromosome accessions                    |
| `^NT_\d{6}$`              | RefSeq contig accessions                        |
| `^NW_\d{8}$`              | RefSeq WGS scaffold accessions                  |

---

## Usage in Genome Retrieval

These validation functions are called automatically by the genome download functions.
If an invalid accession is passed, a `ValueError` is raised before any network
request is made.

```python
from metainformant.dna.external.genomes import download_genome_package, get_genome_metadata

# Validation happens automatically
try:
    metadata = get_genome_metadata("GCF_000001405.39")
except ValueError as e:
    print(f"Invalid accession: {e}")

# Download genome assembly (validates before downloading)
# package_dir = download_genome_package("GCF_000001405.39", "output/genomes/")
```

---

## See Also

- **[NCBI Integration](ncbi.md)** -- Full NCBI API integration
- **[Entrez Utilities](ncbi.md#entrez-client)** -- Entrez search and fetch
- `metainformant.dna.external.genomes` -- Source module
