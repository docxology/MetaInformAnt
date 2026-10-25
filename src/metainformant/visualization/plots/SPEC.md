# Specification: plots

## 🎯 Scope
Visualization plots subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: Listed below; current source is authoritative.
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `animations.py`
- `basic.py`
- `cross_species.py`
- `general.py`
- `multidim.py`
- `rna_campaign.py`: `plot_acquisition_observation(payload, output_prefix, title=...)`
  validates the scope partition and exports PNG/PDF with the shared palette and
  timestamp-independent PDF metadata. It plots operational counts, not replication.
- `specialized.py`
