# Agent Directives: information

**Context**: Information theory analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `integration/`
- `metrics/` — exports: `advanced`, `analysis`, `core`
- `network_info/` — exports: `information_flow`
- `workflow/`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
