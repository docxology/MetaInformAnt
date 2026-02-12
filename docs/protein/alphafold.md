# AlphaFold Integration

Tools for accessing and validating AlphaFold predicted protein structures from
the AlphaFold Protein Structure Database (AFDB).

## Key Concepts

**AlphaFold DB** provides over 200 million predicted protein structures. Each
structure is identified by its UniProt accession and an AlphaFold version
number. Structures are available in PDB and mmCIF formats.

**pLDDT confidence scores** are stored in the B-factor column of AlphaFold PDB
files. Values range from 0 to 100:

| Range  | Interpretation                       |
|--------|--------------------------------------|
| 90-100 | Very high confidence (well-modelled)  |
| 70-89  | Confident (backbone well-modelled)    |
| 50-69  | Low confidence (caution advised)      |
| 0-49   | Very low confidence (disordered)      |

## Function Reference

### build_alphafold_url

```python
def build_alphafold_url(
    uniprot_acc: str,
    *,
    version: int = 4,
    fmt: str = "pdb",
) -> str
```

Constructs the download URL for an AlphaFold model. The `fmt` parameter
accepts `"pdb"` or `"cif"`.

### fetch_alphafold_model

```python
def fetch_alphafold_model(
    uniprot_acc: str,
    out_dir: Path,
    *,
    version: int = 4,
    fmt: str = "pdb",
) -> Path
```

Downloads an AlphaFold structure file to `out_dir`. Returns the path to the
downloaded file. Raises `requests.RequestException` on network failure.

### batch_download_alphafold_models

```python
def batch_download_alphafold_models(
    uniprot_accessions: List[str],
    out_dir: Path,
    max_workers: int = 4,
) -> Dict[str, Path]
```

Downloads multiple structures. Returns a mapping from accession to file path
(or `None` for failed downloads).

### get_alphafold_metadata

```python
def get_alphafold_metadata(uniprot_acc: str) -> Dict[str, Any]
```

Returns metadata including accession, source, method, and download URL for a
given UniProt accession.

### parse_alphafold_confidence

```python
def parse_alphafold_confidence(pdb_path: Path) -> List[float]
```

Parses per-atom pLDDT confidence scores from the B-factor column of an
AlphaFold PDB file.

### validate_alphafold_structure

```python
def validate_alphafold_structure(pdb_path: Path) -> Dict[str, Any]
```

Validates a downloaded PDB file. Returns a dictionary with:
- `is_valid`: Whether the file contains ATOM records.
- `n_atoms`: Total atom count.
- `has_confidence_scores`: Whether B-factors are present.
- `avg_confidence`: Mean pLDDT score.
- `issues`: List of detected problems.

### get_alphafold_structure_quality

```python
def get_alphafold_structure_quality(pdb_path: Path) -> Dict[str, float]
```

Detailed quality assessment: mean pLDDT score, counts of high/medium/low
confidence residues, and a four-bin confidence distribution (very_low, low,
medium, high).

### find_alphafold_models_by_sequence

```python
def find_alphafold_models_by_sequence(
    sequence: str,
    identity_threshold: float = 0.9,
) -> List[Dict[str, Any]]
```

Searches for AlphaFold models matching a protein sequence above the identity
threshold. (API integration placeholder.)

### search_alphafold_by_keyword

```python
def search_alphafold_by_keyword(
    keyword: str,
    max_results: int = 100,
) -> List[Dict[str, Any]]
```

Keyword search against the AlphaFold database. (API integration placeholder.)

### get_alphafold_coverage

```python
def get_alphafold_coverage() -> Dict[str, Any]
```

Returns approximate database coverage statistics (total structures, unique
proteins, coverage percentage, last update date).

## Usage Example

```python
from pathlib import Path
from metainformant.protein.structure.alphafold import (
    fetch_alphafold_model,
    validate_alphafold_structure,
    get_alphafold_structure_quality,
)

# Download and validate
pdb_path = fetch_alphafold_model("P04637", Path("output/structures"))
validation = validate_alphafold_structure(pdb_path)
print(f"Valid: {validation['is_valid']}, Atoms: {validation['n_atoms']}")

quality = get_alphafold_structure_quality(pdb_path)
print(f"Mean pLDDT: {quality['plddt_score']:.1f}")
print(f"High confidence residues: {quality['high_confidence_residues']}")
```

## Configuration

| Environment Variable | Description               | Default |
|---------------------|---------------------------|---------|
| `PROT_TIMEOUT`      | HTTP request timeout (s)  | 30      |

## Related Modules

- `metainformant.protein.database.uniprot` -- UniProt record retrieval
- `metainformant.protein.structure.contacts` -- contact map analysis
- `metainformant.protein.structure.general` -- general structure utilities
- `metainformant.protein.structure.pdb` -- PDB file parsing
