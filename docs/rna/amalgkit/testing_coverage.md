# Amalgkit integration testing

The test suite follows a real-implementation policy: command builders and
workflow validators use the actual wrapper code, while network and large
external-tool runs are marked separately.

## Current coverage lanes

- command availability and version reporting;
- command-specific option filtering for the installed CLI;
- step registry and default plan order;
- configuration parsing and output-root isolation;
- prerequisite validation for metadata, reads, quantification, merge,
  filtering, finalization, and ortholog inputs;
- completion detection based on output content;
- manifest and log creation;
- mounted-data report generation and cross-species input validation.

## Run the deterministic suite

```bash
env -u VIRTUAL_ENV uv run pytest -q \
  -m "not network and not external_tool" tests/rna
```

## Run environment-backed checks

```bash
env -u VIRTUAL_ENV uv run python scripts/rna/validate_all_species_workflow.py
env -u VIRTUAL_ENV uv run python scripts/rna/check_environment.py
```

These checks validate planning and environment contracts; they do not claim
that all configured species have materialized data or that a complete
biological analysis has been run.
