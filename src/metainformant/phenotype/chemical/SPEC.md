# Chemical Module Technical Specification

## Architecture
The `chemical` module provides a structured way to handle metabolomic and chemically-defined phenotypic data.

## Components

### `Compound`
- **Attributes**:
    - `name`: Common name (e.g., "n-C25").
    - `formula`: Chemical formula.
    - `retention_time`: Standardized retention index or time.
    - `identifiers`: Dict of DB IDs (e.g., PubChem, CAS).

### `ChemicalProfile`
- **Attributes**:
    - `sample_id`: Unique identifier.
    - `peaks`: List/Dict of `(Compound, abundance)` tuples.
    - `metadata`: Experimental conditions.
- **Methods**:
    - `normalize()`: Normalizes abundances (e.g., to total ion current).
    - `distance(other)`: Calculates chemical distance (e.g., Bray-Curtis).

## Data Models
- Supports integration with common GC-MS output formats (CSV, simple tabular).
