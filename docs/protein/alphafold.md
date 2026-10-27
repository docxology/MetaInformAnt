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

```python-snippet
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

```python-snippet
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

```python-snippet
def batch_download_alphafold_models(
    uniprot_accessions: List[str],
    out_dir: Path,
    max_workers: int = 4,
) -> Dict[str, Path]
```

Downloads multiple structures. Returns a mapping from accession to file path
(or `None` for failed downloads).

### get_alphafold_metadata

```python-snippet
def get_alphafold_metadata(
    uniprot_acc: str,
    *,
    version: int = 4,
    pdb_path: Optional[Path] = None,
    pae_path: Optional[Path] = None,
) -> Dict[str, Any]
```

Nothing is looked up online. The accession and version give the download `url`;
everything else is read from the local model files that are passed, and is `None`
when the file it comes from is not given. Facts that only the AlphaFold service
holds (organism, gene, sequence length, release date) are not returned.

Returns `uniprot_accession`, `version`, `source`, `method`, `url` and:
- `confidence_score`: mean pLDDT over the residues of `pdb_path`.
- `n_residues`: number of residues of `pdb_path`, or the size of the PAE matrix
  when only `pae_path` is given.
- `plddt_bands`: residues per AlphaFold confidence band: `very_high` (>= 90),
  `confident` (>= 70), `low` (>= 50) and `very_low` (below 50).
- `mean_pae` and `max_pae`: the off-diagonal mean and the largest predicted
  aligned error of `pae_path`.

Raises `FileNotFoundError` for a file that does not exist and `ValueError` for a
file that cannot be parsed, or a PAE matrix and a model with different residue
counts.

### parse_alphafold_pae

```python-snippet
def parse_alphafold_pae(pae_path: Path) -> Dict[str, Any]
```

Parses the predicted aligned error (PAE) JSON file of an AlphaFold model, in
either layout of the AlphaFold database (a full `predicted_aligned_error`
matrix, or the older sparse `residue1`/`residue2`/`distance` lists), inside a
one-element list or as a bare object. Returns `pae_matrix` (N x N, row = aligned
residue, column = scored residue), `n_residues`, `mean_pae` (over the
off-diagonal pairs; `None` for a single residue), `max_pae` (the largest value of
the matrix), `reported_max_pae` (the file's `max_predicted_aligned_error`, the
cap AlphaFold reports for the model, or `None`) and `format` (`"matrix"` or
`"sparse"`). Raises `FileNotFoundError` for a missing file and `ValueError` when
the file is not valid JSON or does not hold a complete square matrix of finite,
non-negative errors.

### parse_alphafold_confidence

```python-snippet
def parse_alphafold_confidence(pdb_path: Path) -> List[float]
```

Returns one pLDDT per residue. AlphaFold writes a residue's pLDDT into the
B-factor column of every one of its atoms, so the score of a residue is the mean
B-factor of its ATOM records. Residues are told apart by chain, residue number
and insertion code and listed in file order; only the first model is read and
HETATM records are ignored. A missing file, or an ATOM record without a numeric
B-factor, is logged and gives an empty list.

### validate_alphafold_structure

```python-snippet
def validate_alphafold_structure(pdb_path: Path) -> Dict[str, Any]
```

Validates a downloaded PDB file. Returns a dictionary with:
- `is_valid`: Whether the file contains ATOM records.
- `n_atoms`: Number of ATOM records.
- `n_residues`: Number of distinct residues (chain, residue number and insertion code).
- `has_confidence_scores`: Whether B-factors are present.
- `avg_confidence`: Mean pLDDT over the residues, not over the atoms.
- `issues`: List of detected problems.

### get_alphafold_structure_quality

```python-snippet
def get_alphafold_structure_quality(pdb_path: Path) -> Dict[str, float]
```

Quality assessment; every figure is per residue. Returns `plddt_score` (mean
pLDDT over the residues), the counts `high_confidence_residues` (>= 90),
`medium_confidence_residues` (>= 70) and `low_confidence_residues` (below 70),
and `confidence_distribution`, a four-bin histogram of the residues: `very_low`
(0-50), `low` (50-70), `medium` (70-90) and `high` (90-100). It does not return a
residue count: `n_residues` comes from `validate_alphafold_structure` and
`get_alphafold_metadata`.

### find_alphafold_models_by_sequence

```python-snippet
def find_alphafold_models_by_sequence(
    sequence: str,
    identity_threshold: float = 0.9,
) -> List[Dict[str, Any]]
```

Not implemented: it raises `NotImplementedError`. A sequence-similarity search
needs a search service that this module cannot call, and an empty list would
claim that no similar model exists. Map the sequence to a UniProt accession with
a sequence-search tool, then call `fetch_alphafold_model`.

### search_alphafold_by_keyword

```python-snippet
def search_alphafold_by_keyword(
    keyword: str,
    max_results: int = 100,
) -> List[Dict[str, Any]]
```

Not implemented: it raises `NotImplementedError`. Find accessions with a UniProt
or AlphaFold DB search, then pass them to `fetch_alphafold_model`.

### get_alphafold_coverage

```python-snippet
def get_alphafold_coverage() -> Dict[str, Any]
```

Not implemented: it raises `NotImplementedError`. Database coverage statistics
would have to be queried from the AlphaFold Protein Structure Database; the
function holds no figures of its own.

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
print(f"Valid: {validation['is_valid']}, Residues: {validation['n_residues']}")

quality = get_alphafold_structure_quality(pdb_path)
print(f"Mean pLDDT: {quality['plddt_score']:.1f}")
print(f"High confidence residues: {quality['high_confidence_residues']}")
```

Confidence and predicted aligned error of a model you already have on disk are read without any request:

```python
from pathlib import Path
from metainformant.protein.structure.alphafold import get_alphafold_metadata

metadata = get_alphafold_metadata(
    "P04637",
    pdb_path=Path("output/structures/AF-P04637-F1-model_v4.pdb"),
    pae_path=Path("output/structures/AF-P04637-F1-predicted_aligned_error_v4.json"),
)
print(metadata["confidence_score"], metadata["plddt_bands"], metadata["mean_pae"])
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
