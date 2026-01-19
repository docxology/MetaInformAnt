# Test Examples Runner

A comprehensive test runner for discovering, executing, and reporting on MetaInformAnt example scripts.

## Purpose

This module provides infrastructure to:
- Discover example scripts across all domains
- Run them sequentially or in parallel
- Generate reports in JSON, HTML, or JUnit XML formats

## Key Components

| File | Description |
|------|-------------|
| [main.py](main.py) | Entry point with `ExampleTester` class |
| [discover.py](discover.py) | `discover_examples()` for finding scripts |
| [runner.py](runner.py) | Sequential/parallel execution logic |
| [validator.py](validator.py) | Verifies script outputs |
| [reporting.py](reporting.py) | JSON, HTML, JUnit report generation |

## Usage

```bash
python -m scripts.test_examples.main --verbose --parallel
```

### CLI Options

| Flag | Description |
|------|-------------|
| `--verbose`, `-v` | Enable verbose output |
| `--domain`, `-d` | Filter by domain (e.g., `dna`, `rna`) |
| `--parallel`, `-p` | Run examples in parallel |
| `--html` | Generate HTML report |
| `--junit-xml` | Generate JUnit XML for CI/CD |

## Related Documentation

- **Parent**: [scripts/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md) - Technical specification.
- **AGENTS**: [AGENTS.md](AGENTS.md) - AI contributions.
- **Test Suite**: [tests/README.md](../../tests/README.md)
