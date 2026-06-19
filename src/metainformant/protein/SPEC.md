# Specification: protein

## Scope

Protein sequence and structure analysis module for METAINFORMANT. Covers proteome
retrieval, domain analysis, structure prediction, and workflow orchestration.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: database, domains, function, sequence, structure, visualization, workflow
- **Key Concepts**: Protein structure, domains, function prediction, sequence analysis

## API Definition

### Exports — `workflow/orchestration.py`

- `analyze_protein_sequence` — Run sequence metrics, composition, hydropathy, motifs, and optional secondary-structure prediction
- `analyze_protein_structure` — Parse a structure file and compute summary structure statistics
- `batch_analyze_sequences` — Analyze a mapping of protein names to sequences
- `comparative_analysis` — Run pairwise and multiple-sequence comparison for related proteins

### Exports — `sequence/proteomes.py`

- `read_taxon_ids` — Read taxon IDs from file
- `download_proteome_fasta` — Download proteome FASTA from UniProt; honors `PROT_TIMEOUT`
- `proteome_statistics` — Calculate summary statistics for a proteome FASTA file
- `compare_proteomes` — Compare two proteome FASTA files

### Exports — `database/uniprot.py`

- `fetch_uniprot_record`, `fetch_uniprot_fasta`, `search_uniprot_proteins` — UniProt REST retrieval helpers; honor `PROT_TIMEOUT`
- `get_uniprot_annotations` — Return normalized GO, keyword, domain, PTM, and location annotations
- `validate_uniprot_accession` — Validate current 6-character and 10-character UniProtKB accession formats
- `map_ids_uniprot` — Map protein identifiers with explicit source database selection

### Exports — `database/interpro.py`

- `fetch_interpro_domains`, `fetch_interpro_by_accession`, `search_interpro_entries` — InterPro REST helpers; honor `PROT_TIMEOUT`
- `parse_interpro_xml` — Parse InterProScan-style XML domain annotations
- `get_interpro_hierarchy`, `get_interpro_statistics`, `find_similar_interpro_entries` — Explicitly unsupported and raise `NotImplementedError`

### Exports — `sequence/sequences.py`

- `read_fasta` — Read protein FASTA files into an ID-to-sequence mapping
- `write_fasta` — Write protein sequences to FASTA
- `validate_protein_sequence` — Validate amino-acid alphabets
- `amino_acid_composition` / `calculate_aa_composition` — Calculate amino-acid composition
- `molecular_weight`, `isoelectric_point`, `gravy`, `hydropathy_score` — Common physicochemical metrics
