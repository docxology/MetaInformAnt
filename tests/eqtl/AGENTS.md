# Agent Directives: eqtl

## 🤖 Role
Specialized agent context for the `eqtl` component.

## 🛠️ Tools & Capabilities
- **Context**: Zero-mocks test suite for the `metainformant.eqtl` pipeline — parameter resolution, synthetic cohorts, variant calling, and variant statistics.
- **Pattern**: Test Suite Pattern

## ⚠️ Rules & Constraints
- **Imports**: Prefer absolute imports from `metainformant`.
- **I/O**: Use `metainformant.core.io` for all file operations.
- **Logging**: Use `metainformant.core.utils.logging`.
