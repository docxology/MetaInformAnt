# SPEC: Test Infrastructure

The MetaInformAnt test suite is designed to ensure mathematical and biological accuracy through rigorous, mock-free verification.

## Testing Strategy

### 1. No-Mocking Policy
We do not use `unittest.mock` or similar libraries. All tests must use real data and live methods. If a method requires an external resource (e.g., a database or API), a local instance or a file-based cache must be provided.

### 2. Composition Over Integration
Tests are categorized by their scope:
- **Unit**: Single function verification with minimal dependencies.
- **Comprehensive**: Domain-level testing involving multiple functions within a module.
- **Integration**: Cross-domain workflows (e.g., DNA + RNA + Multi-omics).

## Test Data Management

Test data is located in `tests/data/`.
- **Static Assets**: Small FASTA/FASTQ files, configuration templates.
- **Dynamic Caching**: Downloads from NCBI/SRA are cached in `temp` or `tests/data/cache` during test runs to avoid redundant network calls.

## Verification Requirements

- **Accuracy**: Mathematical results must match known benchmarks or analytical solutions.
- **Reproducibility**: Tests must pass consistently across different environments (macOS/Linux).
- **Compliance**: All new features must include tests that verify the "Triple Play" documentation exists.

## Usage

Tests are executed using `pytest`:
```bash
uv run pytest tests/
```
Specific domain tests:
```bash
uv run pytest tests/rna/
```
