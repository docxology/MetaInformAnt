# Ecology Module Technical Specification

## Architecture
The `ecology` module is a component of the `metainformant` package, designed to encapsulate ecological metadata and community analysis.

## Components
*   **`__init__.py`**: Exposes public API.
*   **Submodules**: Specialized components for domain logic.

## Dependencies
*   **Internal**: `metainformant.core` and related domain modules.
*   **External**: Standard scientific stack (numpy, pandas, scipy) and domain-specific tools.

## Standards
*   **Code Style**: PEP 8
*   **Docstrings**: Google Style
*   **Type Hints**: Full coverage required
*   **Testing**: Pytest with 100% coverage target
