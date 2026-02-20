# Specification: phenotype

## 🎯 Scope

Phenotype module for MetaInformAnt. Multi-modal phenotyping across morphological,
behavioral, chemical, acoustic, and electronic tracking domains.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: analysis, behavior, chemical, data, electronic, gwas_integration, integration, morphological, sonic, visualization, workflow
- **Key Concepts**: Morphology, behavior, chemical profiles, acoustic, electronic tracking

## 🔌 API Definition

### Exports — `analysis/life_course.py`

- `extract_phenotypes_from_events` — Extract phenotype features from event sequences (life_events or local)
- `aggregate_temporal_phenotypes` — Time-windowed phenotype aggregation
- `map_events_to_traits` — Map events to trait definitions with counts/timestamps

### Exports — `morphological/`

- `Measurement` — Dataclass for morphometric measurements with unit conversion (mm↔cm↔m)
- `MorphometricProfile` — Collection of measurements for a specimen
