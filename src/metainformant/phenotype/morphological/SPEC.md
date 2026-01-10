# Morphological Module Technical Specification

## Architecture
The `morphological` module provides strict typing and unit awareness for morphometric data.

## Components

### `Measurement`
- **Attributes**:
    - `name`: Standardized name (e.g., "Head Width").
    - `value`: Float value.
    - `unit`: Unit string (default "mm").
    - `uncertainty`: Optional error margin.

### `MorphometricProfile`
- **Attributes**:
    - `specimen_id`: Unique identifier.
    - `measurements`: Dict of name->Measurement.
- **Methods**:
    - `get(name)`: Retrieve value.
    - `calculate_index(name, numerator, denominator)`: Returns ratio * 100.

## Integration
- Should map readily to the structure used in `antwiki.py`.
