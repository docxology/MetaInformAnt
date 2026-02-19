# Specification: structural_variants

## 🎯 Scope

Structural variant analysis: CNV detection via circular binary segmentation, SV calling from split/discordant reads, annotation, filtering, population genotyping, and visualization.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: detection, annotation, filtering, population, visualization
- **Key Concepts**: `CNVSegment`, `StructuralVariant`, `SVType`, circular binary segmentation

## 🔌 API Definition

### Exports

- `__init__.py`
