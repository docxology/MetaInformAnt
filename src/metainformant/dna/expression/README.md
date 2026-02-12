# Expression

DNA transcription, translation, and codon usage analysis including codon optimization for heterologous expression studies.

## Contents

| File | Purpose |
|------|---------|
| `codon.py` | Codon usage bias analysis, CAI calculation, and codon optimization |
| `transcription.py` | DNA-to-RNA transcription, reverse complement, and TSS prediction |
| `translation.py` | RNA/DNA-to-protein translation, ORF finding, and six-frame translation |

## Key Functions

| Function | Description |
|----------|-------------|
| `transcribe()` | Convert DNA sequence to RNA (T to U replacement) |
| `transcribe_reverse_complement()` | Transcribe the antisense strand to RNA |
| `transcribe_with_introns()` | Transcribe DNA after splicing out intron regions |
| `find_transcription_start_sites()` | Locate promoter motifs (e.g., TATA box) in DNA |
| `translate()` | Translate RNA sequence to amino acid sequence |
| `translate_dna()` | Direct DNA-to-protein translation |
| `find_orfs()` | Find open reading frames above a minimum length |
| `six_frame_translation()` | Translate all six reading frames of a DNA sequence |
| `codon_usage()` | Calculate relative synonymous codon usage frequencies |
| `cai()` | Compute Codon Adaptation Index for expression prediction |
| `optimize_codons()` | Recode sequence to match a target codon usage table |
| `calculate_enc()` | Effective number of codons (ENC) as a bias metric |
| `back_translate()` | Convert protein sequence back to DNA using codon preferences |

## Usage

```python
from metainformant.dna.expression.transcription import transcribe
from metainformant.dna.expression.translation import translate, find_orfs
from metainformant.dna.expression.codon import cai, codon_usage

rna = transcribe("ATGGCCATTGTAATG")
protein = translate(rna)
orfs = find_orfs(rna, min_length=30)
usage = codon_usage("ATGGCCATTGTAATGGCCATT")
```
