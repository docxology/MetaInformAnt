# Protein Contact Analysis

Functions for computing and analysing protein residue-residue contacts,
hydrogen bonds, salt bridges, hydrophobic interactions, disulfide bonds,
and contact network properties.

## Key Concepts

**Contact maps** are binary symmetric matrices where entry (i, j) = 1 if the
minimum distance between any atom of residue i and residue j is below a
distance threshold (default 8.0 Angstroms).

**Interaction types** classified by this module:

| Type               | Criterion                          | Default Threshold |
|-------------------|------------------------------------|-------------------|
| Residue contact    | Any atom pair within threshold     | 8.0 A             |
| Hydrogen bond      | N-O donor-acceptor distance        | 3.5 A             |
| Salt bridge        | Charged group centre distance      | 4.0 A             |
| Hydrophobic        | Non-polar sidechain centre distance| 5.0 A             |
| Disulfide bond     | Cysteine SG-SG distance            | 2.5 A             |

**C-alpha contacts** provide a simplified contact representation using only
backbone C-alpha atom positions, useful for coarse-grained analysis and
phylogenetic comparison.

## Function Reference

### calculate_residue_contacts

```python
def calculate_residue_contacts(
    coords: np.ndarray,
    residue_ranges: List[Tuple[int, int]],
    threshold: float = 8.0,
) -> np.ndarray
```

Computes a binary contact map from atomic coordinates. `residue_ranges` maps
each residue to its (start_atom, end_atom) index range. Returns an
(n_residues x n_residues) numpy array.

### identify_hydrogen_bonds

```python
def identify_hydrogen_bonds(
    atoms: List[Dict[str, Any]],
    coords: np.ndarray,
    distance_threshold: float = 3.5,
    angle_threshold: float = 120.0,
) -> List[Dict[str, Any]]
```

Identifies potential hydrogen bonds between nitrogen (donor) and oxygen
(acceptor) atoms. Each result contains `donor_atom`, `acceptor_atom`,
`distance`, `angle`, and `strength` (inverse distance normalised to [0, 1]).

### identify_salt_bridges

```python
def identify_salt_bridges(
    atoms: List[Dict[str, Any]],
    coords: np.ndarray,
    distance_threshold: float = 4.0,
) -> List[Dict[str, Any]]
```

Detects electrostatic interactions between positively charged residues (ARG,
LYS, HIS) and negatively charged residues (ASP, GLU). Returns pairs with
residue identifiers and centre-of-charge distances.

### identify_hydrophobic_contacts

```python
def identify_hydrophobic_contacts(
    atoms: List[Dict[str, Any]],
    coords: np.ndarray,
    distance_threshold: float = 5.0,
) -> List[Dict[str, Any]]
```

Finds contacts between non-polar residues (ALA, VAL, LEU, ILE, MET, PHE, TRP,
PRO) based on sidechain carbon atom centres.

### identify_disulfide_bonds

```python
def identify_disulfide_bonds(
    atoms: List[Dict[str, Any]],
    coords: np.ndarray,
    distance_threshold: float = 2.5,
) -> List[Dict[str, Any]]
```

Identifies covalent S-S bonds between cysteine SG atoms.

### classify_contact_types

```python
def classify_contact_types(
    atoms: List[Dict[str, Any]],
    coords: np.ndarray,
) -> Dict[str, List[Dict[str, Any]]]
```

Runs all contact analysis functions and returns a unified dictionary with
keys: `residue_contacts`, `hydrogen_bonds`, `salt_bridges`,
`hydrophobic_contacts`, `disulfide_bonds`, and a `summary` with total counts.

### analyze_contact_network

```python
def analyze_contact_network(
    contact_map: np.ndarray,
) -> Dict[str, Any]
```

Graph analysis of the contact map. Returns:
- `n_residues`, `n_contacts` (undirected)
- `average_degree`, `max_degree`
- `clustering_coefficient` (triangle-based)
- `n_components`, `component_sizes`, `largest_component_size`

### calculate_contact_persistence

```python
def calculate_contact_persistence(
    contact_maps: List[np.ndarray],
) -> np.ndarray
```

Computes the fraction of input structures in which each contact is present.
Useful for comparing homologous structures or MD trajectory snapshots.

### compute_ca_contact_pairs

```python
def compute_ca_contact_pairs(
    coords: list,
    threshold: float = 8.0,
) -> list
```

Returns all (i, j) pairs (i < j) of C-alpha atoms within the distance
threshold. Uses scipy `pdist` when available, falls back to numpy
broadcasting.

```python
>>> coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
>>> compute_ca_contact_pairs(coords, threshold=4.0)
[(0, 1)]
```

## Usage Example

```python
import numpy as np
from metainformant.protein.structure.contacts import (
    calculate_residue_contacts,
    analyze_contact_network,
    compute_ca_contact_pairs,
)

# Compute contact map from coordinates
coords = np.random.rand(100, 3) * 50.0
residue_ranges = [(i*5, (i+1)*5) for i in range(20)]
contact_map = calculate_residue_contacts(coords, residue_ranges, threshold=8.0)

# Analyse the contact network
network = analyze_contact_network(contact_map)
print(f"Contacts: {network['n_contacts']}")
print(f"Avg degree: {network['average_degree']:.1f}")
print(f"Clustering: {network['clustering_coefficient']:.3f}")

# C-alpha contact pairs
ca_coords = [(float(i), 0.0, 0.0) for i in range(10)]
pairs = compute_ca_contact_pairs(ca_coords, threshold=2.0)
```

## Related Modules

- `metainformant.protein.structure.alphafold` -- predicted structures
- `metainformant.protein.structure.analysis` -- structure analysis utilities
- `metainformant.protein.structure.general` -- general structure functions
- `metainformant.protein.structure.secondary` -- secondary structure
