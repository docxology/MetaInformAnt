# PROTEIN

## Overview
Protein sequence and structure analysis module for METAINFORMANT.

## 📦 Contents
- **[database/](database/)**
- **[sequence/](sequence/)**
- **[structure/](structure/)**
- **[visualization/](visualization/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Protein Module"
        DB[database/] --> |uniprot.py| UP[UniProt Record Fetching]
        DB --> |interpro.py| IP[InterPro Domain Lookup]

        ST[structure/] --> |analysis.py| AN[Contact Maps, Domain ID]
        ST --> |alphafold.py| AF[AlphaFold Model Retrieval]
        ST --> |pdb.py| PD[PDB Parsing]
        ST --> |contacts.py| CO[Residue Contacts]
        ST --> |secondary.py| SS[Secondary Structure]

        SQ[sequence/] --> |sequences.py| FA[FASTA Read/Write]
        SQ --> |alignment.py| AL[Sequence Alignment]
        SQ --> |proteomes.py| PR[Proteome Downloads]

        FN[function/] --> |prediction.py| FP[Function Prediction]
        DO[domains/] --> |detection.py| DD[Domain Detection]
    end
```

## Usage
Import module:
```python
from metainformant.protein import ...
```
