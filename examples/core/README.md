# CORE

## Overview
Runnable examples for Core workflows. Keep these examples small, deterministic, and aligned with the current public APIs.

## Contents
- [example_config.py](example_config.py)
- [example_io.py](example_io.py)
- [example_logging.py](example_logging.py)
- [example_paths.py](example_paths.py)
- [example_workflow.py](example_workflow.py)

## Usage
```bash
uv run python examples/core/example_config.py
```

The [workflow example](example_workflow.py) reads the checked-out project metadata
and example dependency manifest, joins both results through the workflow engine,
and writes its configuration and results under `output/examples/core/`. It runs
without network access and raises an error if validation or execution fails.
