# Agent Directives: menu

**Context**: Interactive menu and CLI interface for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `core/` — exports: `discovery`, `executor`
- `ui/` — exports: `display`, `navigation`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
