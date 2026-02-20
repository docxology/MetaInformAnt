# Specification: protein

## 🎯 Scope

Protein sequence and structure analysis module for METAINFORMANT. Covers proteome
retrieval, domain analysis, structure prediction, and workflow orchestration.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: analysis, data, database, domains, function, sequence, structure, visualization, workflow
- **Key Concepts**: Protein structure, domains, function prediction, sequence analysis

## 🔌 API Definition

### Exports — `analysis/orchestration.py`

- `ProteinPipeline` — Full protein analysis pipeline orchestrator

### Exports — `data/proteomes.py`

- `read_taxon_ids` — Read taxon IDs from file (returns `List[str]`)
- `download_proteome` — Download proteome from UniProt
- `parse_fasta` — Parse FASTA format proteome files
