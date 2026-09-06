# Specification: ontology

## Scope
Gene ontology and functional annotation module for METAINFORMANT.

## Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures
- **Modules**: 12 Python modules across 6 subpackages (annotation, core, pathway_enrichment, query, visualization, workflow)
- **Key Concepts**: Refer to Pydantic models in source.

## API Definition
### Exports
- `__init__.py`
- `annotation/annotate.py`
- `core/go.py`
- `core/go_api.py`
- `core/hpo.py`
- `core/obo.py`
- `core/types.py`
- `pathway_enrichment/enrichment.py`
- `query/query.py`
- `query/serialize.py`
- `visualization/plots.py`
- `visualization/visualization.py`
- `workflow/run_ontology.py`
