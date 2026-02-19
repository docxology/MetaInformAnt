# Specification: longread

## 🎯 Scope

PacBio and Oxford Nanopore long-read analysis: signal I/O, quality assessment, assembly, methylation calling, haplotype phasing, and structural variant detection.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: io, quality, analysis, assembly, methylation, phasing, workflow, visualization, utils
- **Key Concepts**: `LongReadOrchestrator`, `PipelineStep`, `PipelineResult`

## 🔌 API Definition

### Exports

- `__init__.py`
