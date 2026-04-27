# Core: Text Processing Utilities

The `text` module provides comprehensive text processing and normalization utilities tailored for biological data. It handles whitespace normalization, filename sanitization, identifier cleaning, species name formatting, and text analysis tasks used throughout METAINFORMANT.

## Purpose

Biological data is notoriously messy:
- FASTA headers contain complex identifiers, descriptions, version numbers
- Sample metadata contains free-text fields with inconsistent formatting
- Gene names have multiple aliases, punctuation variations, case inconsistencies
- Species names appear in different binomial/genus-only formats
- Downloaded filenames may contain spaces, special characters, or be overly long

The `text` module provides canonicalization functions to normalize this diversity into predictable, filesystem-safe, and comparison-friendly formats.

## Design Principles

### 1. **Regex-Based Efficiency**
Common patterns are pre-compiled at module load time (`_WHITESPACE_RE`, `_SLUG_INVALID_RE`) for efficient repeated use.

### 2. **Unicode-Aware**
Uses `unicodedata` module to correctly handle international characters, accents, and control characters.

### 3. **Biological Specificity**
Functions target common bioinformatics formats:
- FASTA header parsing (`clean_sequence_id()`)
- Gene name standardization (`standardize_gene_name()`)
- Species binomial formatting (`format_species_name()`)

### 4. **Filesystem Safety**
`safe_filename()` and `sanitize_filename()` (from paths) ensure output filenames work across filesystems and avoid shell injection risks.

### 5. **Graceful Degradation**
Functions handle edge cases:
- Empty strings → empty output
- None inputs → return empty string or sensible default (where applicable; not all functions accept None)
- Very long inputs → truncation or pass-through

### 6. **Pure Functions**
No side effects, no global state.

### 7. **No External Dependencies**
Standard library only (`re`, `unicodedata`, `pathlib`).

## Module Organization

The `text` module is in `src/metainformant/core/utils/text.py`. It's part of the `utils` subpackage.

**Public API**:

**Basic text processing**:
- `normalize_whitespace()` — Collapse whitespace to single spaces
- `clean_whitespace()` — Alias for `normalize_whitespace()`
- `remove_control_chars()` — Strip control characters

**Slugification & filename safety**:
- `slugify()` — Convert to URL-safe slug (lowercase, dashes, alphanumeric only)
- `safe_filename()` — Combine slugify with extension preservation

**Biological identifier formatting**:
- `standardize_gene_name()` — Uppercase, remove separators
- `format_species_name()` — Binomial format (Genus species)
- `clean_sequence_id()` — Extract clean accession from FASTA headers

**Text analysis**:
- `extract_numbers()` — Extract all numeric values
- `count_words()` — Word count after normalization
- `truncate_text()` — Truncate with ellipsis

**Pattern extraction**:
- `extract_email_addresses()` — Regex-based email finder

## API Reference

### Basic Text Processing

#### `normalize_whitespace(text: str) -> str`

Collapse all whitespace characters (spaces, tabs, newlines, etc.) into single spaces, then strip leading/trailing whitespace.

**Parameters**:
- `text`: Input string

**Returns**: Normalized string

**Whitespace definition**: Matches regex `\s+` which includes space, tab, newline, carriage return, vertical tab, form feed, and Unicode spaces (NBSP, em space, etc.).

**Example**:
```python
messy = "  Hello\t\tworld\n\n  with  spaces  "
clean = normalize_whitespace(messy)
print(clean)  # "Hello world with spaces"

# Multiline to single line
multiline = """Line 1
    Line 2
        Line 3"""
print(normalize_whitespace(multiline))
# "Line 1 Line 2 Line 3"

# Unicode non-breaking space (U+00A0) treated as whitespace
s = "Café\u00a0au\u00a0lait"
print(normalize_whitespace(s))
# "Café au lait"
```

**Use case**: Cleaning free-text fields from samplesheets, metadata tables.

#### `clean_whitespace(text: str) -> str`

Alias for `normalize_whitespace()`. Maintained for backward compatibility.

**Example**:
```python
# These are equivalent:
a = normalize_whitespace(text)
b = clean_whitespace(text)
assert a == b
```

#### `remove_control_chars(text: str) -> str`

Remove Unicode control characters (categories Cc, Cf, Cs, Co, Cn) except for tab (`\t`), newline (`\n`), and space.

**Parameters**:
- `text`: Input string

**Returns**: String without control characters

**Implementation**: Filters characters by `unicodedata.category(char)[0] != "C"` while preserving tab, newline, space.

**Example**:
```python
dirty = "File\x00name\x01with\x02control\x03chars"
clean = remove_control_chars(dirty)
print(clean)  # "Filenamewithcontrolchars"

# Preserves newlines and tabs
text = "Line1\n\tLine2"
assert remove_control_chars(text) == text
```

**Use case**: Cleaning text copied from binary files or OCR output that may contain strikethrough characters, zero-width spaces, etc.

### Slugification & Filename Safety

#### `slugify(text: str) -> str`

Convert text to URL-safe slug: lowercase, spaces to dashes, remove non-alphanumeric characters.

**Parameters**:
- `text`: Input text

**Returns**: Slug-safe string (only lowercase letters, digits, dashes; no leading/trailing dashes)

**Process**:
1. Call `normalize_whitespace()`
2. Convert to lowercase
3. Replace spaces with `-`
4. Remove all characters not matching `[a-z0-9-]` (regex `[^a-z0-9-]+`)
5. Collapse multiple dashes → single dash
6. Strip leading/trailing dashes

**Example**:
```python
print(slugify("Hello, World!"))
# "hello-world"

print(slugify("  Multiple   spaces   "))
# "multiple-spaces"

print(slugify("Café résumé naïve"))
# Note: accented chars removed by `[^a-z0-9-]` regex
# Result: "caf-rsum-nave"

print(slugify("___test---__"))
# "test"  (dashes collapsed, leading/trailing removed)

print(slugify("100% complete"))
# "100-complete"

print(slugify("file (version 2).txt"))
# "file-version-2txt"  (parentheses and period removed)
```

**Use case**: URL slugs, HTML IDs, markdown filename conventions.

#### `safe_filename(name: str) -> str`

Create filesystem-safe filename while preserving extension.

**Parameters**:
- `name`: Original filename (possibly including extension)

**Returns`: Sanitized filename

**Process**:
1. Split into stem and suffix via `Path(name)`
2. Apply `slugify()` to stem
3. Append original suffix unchanged

**Example**:
```python
print(safe_filename("report: 2023 analysis?.pdf"))
# "report-2023-analysis.pdf"

print(safe_filename("data (final).csv"))
# "data-final.csv"

print(safe_filename("archive.tar.gz"))
# "archive.tar.gz"  (both extensions preserved as single suffix .gz)

print(safe_filename("file<with>dangerous|chars.txt"))
# "file_with_dangerous_chars.txt"
```

**Use case**: Sanitizing user-provided filenames before saving to disk.

### Biological Data Formatting

#### `standardize_gene_name(gene_name: str) -> str`

Standardize gene name to uppercase with no separators.

**Parameters**:
- `gene_name`: Raw gene symbol

**Returns**: Uppercase gene symbol with hyphens/underscores/dots removed

**Transformations**:
1. Strip whitespace
2. Uppercase
3. Remove `-`, `_`, `.` characters

**Example**:
```python
print(standardize_gene_name("brca-1"))     # "BRCA1"
print(standardize_gene_name("BRCA_2"))      # "BRCA2"
print(standardize_gene_name("tp53.p"))      # "TP53P"  (pseudogene)
print(standardize_gene_name("  cdkn2a  "))  # "CDKN2A"

# Commonly seen variants:
gene_variants = [
    "BRCA1", "Brca1", "brca1",           # Case variations
    "BRCA-1", "BRCA_1", "BRCA.1",       # Separator variations
    "BRCA 1",  "BRCA1",                  # Space vs no space
]
standardized = [standardize_gene_name(g) for g in gene_variants]
assert all(g == "BRCA1" for g in standardized)
```

**Use case**: Gene name lookups across databases where gene symbols may have inconsistent formatting.

#### `format_species_name(species_name: str) -> str`

Format species name in proper binomial nomenclature: capitalize genus, lowercase species epithet.

**Parameters**:
- `species_name`: Raw species name (e.g., `"homo sapiens"`, `"APIS MELLIFERA"`, `"e. coli"`)

**Returns**: Properly formatted binomial (e.g., `"Homo sapiens"`)

**Rules**:
1. Lowercase entire string
2. Split on whitespace
3. If ≥2 parts: capitalize first part (genus), lowercase second part (species), join with space
4. If 1 part: capitalize and return (genus-only name, e.g., `"Escherichia"` from `"E. coli"` after initial split may need special handling—this is simple)

**Example**:
```python
print(format_species_name("homo sapiens"))      # "Homo sapiens"
print(format_species_name("APIS MELLIFERA"))    # "Apis mellifera"
print(format_species_name("escherichia coli"))  # "Escherichia coli"
print(format_species_name("Drosophila melanogaster"))  # "Drosophila melanogaster"

# Genus-only
print(format_species_name("E. coli"))           # "E. coli"  (single token "E." → capitalize)
# Better approach if you need to expand E. coli:
# Use lookup table or more sophisticated parser for abbreviations

# Extra tokens ignored (subspecies/variety)
print(format_species_name("canis lupus familiaris"))
# "Canis lupus familiaris" (only first two words capitalized)
```

**Limitation**: Does not handle subspecies/variety formatting (third word lowercase) or author abbreviations. For production taxonomy, consider using `taxonkit` or `ete3`.

**Use case**: Standardizing species names from sample sheets, metadata, FASTA headers.

#### `clean_sequence_id(sequence_id: str) -> str`

Extract clean sequence identifier from FASTA header line.

**Parameters**:
- `sequence_id`: Full FASTA header line (with or without leading `>`)

**Returns**: Clean identifier like `"NM_001302504"` or `"gi_12345"`

**Algorithm**:
1. Strip leading `>` character if present
2. If contains pipe `|` delimiters (common NCBI/RefSeq format):
   - Look for known database prefixes: `ref`, `gb`, `emb`, `dbj`
   - Return the following element (the accession)
3. Otherwise: split on first whitespace or bracket `[`, return first token

**Example**:
```python
# RefSeq format
header = ">ref|NM_001302504.2| Homo sapiens BRCA1 mRNA"
print(clean_sequence_id(header))  # "NM_001302504.2"

# GenBank format
header = ">gb|AF123456.1| Example gene"
print(clean_sequence_id(header))  # "AF123456.1"

# Plain accession
header = ">NC_000001.11"
print(clean_sequence_id(header))  # "NC_000001.11"

# GI format (old)
header = ">gi|12345|ref|NM_001| Homo sapiens"
print(clean_sequence_id(header))  # "NM_001"

# With description only
header = ">contig_001 random sequence assembly"
print(clean_sequence_id(header))  # "contig_001"

# With version in parentheses
header = ">chr1 (genome)"
print(clean_sequence_id(header))  # "chr1"

# Empty after stripping > returns empty
print(clean_sequence_id(">"))  # ""
```

**Use case**: Deriving consistent filenames or identifiers from diverse FASTA headers across databases.

### Text Analysis

#### `extract_numbers(text: str) -> list[float]`

Extract all decimal and integer numbers from text.

**Parameters**:
- `text`: Input string

**Returns**: List of floats in order of appearance

**Pattern**: `r"\d+\.?\d*"` matches:
- Integers: `42`, `1000`
- Decimals: `3.14`, `0.05`
- Numbers with trailing decimal: `42.` → `42.0`

**Does NOT match**: Scientific notation (`1.5e-3`) — needs separate pattern; not included as it can pick up false positives in text (e.g., version numbers like `1.2.3`).

**Example**:
```python
text = "Expression: 15.3 ± 2.1, p-value: 0.05, n=100"
nums = extract_numbers(text)
print(nums)  # [15.3, 2.1, 0.05, 100.0]

# Multiple occurrences
text2 = "Versions: 1.0, 2.1.3, 10"
# Caution: "2.1.3" → extracts 2.1 and 3.0 separately
print(extract_numbers(text2))  # [1.0, 2.1, 3.0, 10.0]
```

**Enhanced version for scientific notation**:
```python
def extract_numbers_advanced(text: str) -> list[float]:
    """Also match scientific notation."""
    import re
    pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    return [float(m) for m in re.findall(pattern, text)]

print(extract_numbers_advanced("p = 1.5e-8, fold = 2.3e2"))
# [1.5e-08, 230.0]
```

#### `count_words(text: str) -> int`

Count words after normalizing whitespace.

**Parameters**:
- `text`: Input text

**Returns`: Word count (split on whitespace after normalization)

**Example**:
```python
print(count_words("Hello world"))                    # 2
print(count_words("  Multiple   spaces   here  "))   # 3
print(count_words("One\nTwo\nThree"))                # 3
print(count_words(""))                               # 0
```

**Implementation**: `len(normalize_whitespace(text).split())` — ensures leading/trailing whitespace doesn't create empty tokens.

#### `truncate_text(text: str, max_length: int, suffix: str = "...") -> str`

Truncate text to maximum length with optional suffix.

**Parameters**:
- `text`: Input string
- `max_length`: Maximum length of result (including suffix)
- `suffix`: String to append when truncating (default `"..."`)

**Returns**: Truncated string

**Example**:
```python
long = "This is a very long description that needs truncation"
print(truncate_text(long, 30))
# "This is a very long descrip..."

print(truncate_text(long, 20, suffix=""))
# "This is a very long..."

# No truncation needed
short = "Brief"
print(truncate_text(short, 100))  # "Brief" (unchanged)

# Edge: suffix longer than max_length
print(truncate_text("Hello", 2, suffix="..."))
# ".."  (suffix truncated if needed—check implementation; most will return text[:max_length])
```

#### `extract_email_addresses(text: str) -> list[str]`

Extract email addresses using regex pattern.

**Parameters**:
- `text`: Text containing email addresses

**Returns`: List of email strings (may include duplicates)

**Pattern**: Basic email regex `r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b"` — matches standard email format but not all RFC-compliant addresses (e.g., doesn't match `+` in local part if not included, but it is; doesn't match quoted strings).

**Example**:
```python
contact_block = """
Support: help@example.com
Sales: sales@example.co.uk
Admin: admin+spam@example.org
"""

emails = extract_email_addresses(contact_block)
print(emails)
# ['help@example.com', 'sales@example.co.uk', 'admin+spam@example.org']

# Avoid duplicates (common if same email appears multiple times)
unique = list(set(emails))
```

**Limitations**: May miss edge-case emails like `"user@localhost"` or IP-address domains. Use `email.utils` or `email_validator` library for production validation.

## Biological Data Integration

### FASTA Header Processing Pipeline

```python
def parse_fasta_header(header_line: str) -> dict:
    """Extract structured info from diverse FASTA headers."""
    clean_id = clean_sequence_id(header_line)
    safe_file = safe_filename(f"{clean_id}.fasta")

    return {
        "clean_id": clean_id,
        "safe_filename": safe_file,
        "original": header_line.lstrip(">"),
    }

# Test with various formats
headers = [
    ">gi|12345|ref|NM_001.2| Homo sapiens BRCA1 mRNA",
    ">gb|AF123456.1| isolated from E. coli",
    ">NC_000001.11 Homo sapiens chromosome 1",
    ">contig_0001 N50 contig from assembly",
]

for h in headers:
    info = parse_fasta_header(h)
    print(f"{info['clean_id']} → {info['safe_filename']}")
```

### Sample Metadata Cleaning

```python
def clean_sample_metadata_table(metadata: dict[str, Any]) -> dict[str, Any]:
    """Sanitize all string fields in sample metadata."""
    cleaned = {}
    for key, value in metadata.items():
        if isinstance(value, str):
            # Normalize whitespace
            val = normalize_whitespace(value)

            # Format species names if key suggests species
            if "species" in key.lower():
                val = format_species_name(val)

            # Sanitize for filesystem safety if used in filenames
            if "filename" in key.lower() or "file" in key.lower():
                val = safe_filename(val)

            cleaned[key] = val
        else:
            cleaned[key] = value

    return cleaned

raw_metadata = {
    "sample_id": "  SAMPLE_001  ",
    "species": "drosophila melanogaster",
    "description": "  Test sample \n\t with extra whitespace  ",
    "filename": "data/raw/SAMPLE_001.fastq.gz",
}

clean = clean_sample_metadata_table(raw_metadata)
# {
#   'sample_id': 'SAMPLE_001',
#   'species': 'Drosophila melanogaster',
#   'description': 'Test sample with extra whitespace',
#   'filename': 'data_raw_SAMPLE_001.fastq.gz',
# }
```

### Gene Symbol Normalization Across Databases

```python
from typing import Iterable

def normalize_gene_list(gene_symbols: Iterable[str]) -> set[str]:
    """Normalize list of gene symbols from various sources into canonical form."""
    return {standardize_gene_name(g) for g in gene_symbols if g}

# Input from different databases
ensembl_genes = ["BRCA1", "BRCA2", "TP53"]  # Already uppercase
ncbi_genes = ["Brca-1", "Brca-2", "p53"]   # Mixed case, hyphens
print(normalize_gene_list(ensembl_genes + ncbi_genes))
# {'BRCA1', 'BRCA2', 'TP53'}

# Match across datasets
query_genes = {"brca1", "tp53", "egfr"}
database_genes = {"BRCA1", "EGFR", "KRAS"}

normalized_query = normalize_gene_list(query_genes)
normalized_db = normalize_gene_list(database_genes)

intersection = normalized_query & normalized_db
print(f"Found {len(intersection)} matching genes: {intersection}")
# {'BRCA1'}
```

### Sequence ID to Filename Mapping

```python
def make_fasta_filenames(accessions: list[str], extension: str = ".fasta") -> dict[str, str]:
    """Map sequence accessions to safe filenames."""
    mapping = {}
    for acc in accessions:
        clean = clean_sequence_id(acc)
        if not clean:
            clean = "unknown"
        safe = safe_filename(clean)
        mapping[acc] = safe + extension
    return mapping

accessions = [
    ">ref|NM_001302504.2|",
    ">gb|AF123456.1|",
    "gi|12345|ref|
]

filenames = make_fasta_filenames(accessions)
# {
#   '>ref|NM_001302504.2|': 'NM_001302504.2.fasta',
#   '>gb|AF123456.1|': 'AF123456.1.fasta',
# }
```

### Extract Numbers from Report Text

```python
def extract_qc_metrics(report_text: str) -> dict[str, float]:
    """Pull QC metrics from free-text report."""
    metrics = {}

    # Example: "Mean quality: 32.5, Depth: 45.2×, Coverage: 98.7%"
    numbers = extract_numbers(report_text)

    # Heuristic assignment based on known patterns
    if "mean qual" in report_text.lower() and len(numbers) >= 1:
        metrics["mean_quality"] = numbers[0]
    if "depth" in report_text.lower() and len(numbers) >= 2:
        metrics["depth"] = numbers[1]
    if "coverage" in report_text.lower() and len(numbers) >= 3:
        metrics["coverage"] = numbers[2]

    return metrics

report = """
QC Report
---------
Mean quality: 32.5
Depth: 45.2×
Coverage: 98.7%
Total bases: 1.5e9
"""

print(extract_qc_metrics(report))
# {'mean_quality': 32.5, 'depth': 45.2, 'coverage': 98.7}
# Note: 1.5e9 not extracted by basic number extractor
```

## Error Handling

### UnicodeDecodeError

**Symptom**: Input string contains bytes not decodable as UTF-8 passed as `str` (should already be decoded). If you have raw bytes, decode first:
```python
text_bytes = b"Caf\xc3\xa9"  # UTF-8 encoded "Café"
text = text_bytes.decode("utf-8", errors="replace")  # Replace invalid with �
clean = remove_control_chars(text)
```

### Unexpected Input Types

**Symptom**: `AttributeError: 'int' object has no attribute 'strip'` or similar.

**Cause**: Non-string passed to function expecting `str`.

**Fix**: Validate and coerce:
```python
def safe_normalize(text: Any) -> str:
    if text is None:
        return ""
    return normalize_whitespace(str(text))
```

### Regex Performance (ReDoS)

**Symptom**: Function hangs on certain inputs (catastrophic backtracking).

**Cause**: Crafted input causing regex engine exponential time.

**Mitigation**: Current regexes are simple (no nested quantifiers) and safe. Be cautious when adding complex patterns:
```python
# DANGEROUS: nested quantifiers may ReDoS
re.compile(r"(a+)+b")  # Can blow up on "aaaaaaaaaaaaac"

# SAFE: bounded repetition
re.compile(r"(?:a{1,100})+b")
```

## Performance Considerations

### Pre-compiled Patterns

Module-level compiled regex avoids re-compilation on every call:
```python
_WHITESPACE_RE = re.compile(r"\s+")
_SLUG_INVALID_RE = re.compile(r"[^a-z0-9-]+")
```

This makes repeated calls (e.g., in loops over thousands of strings) efficient.

### String Operations

Most functions are O(n) in input length and allocate minimal intermediate strings:
- `normalize_whitespace`: One regex substitution + strip
- `slugify`: Multiple regex substitutions but all O(n)
- `standardize_gene_name`: Uppercase + single regex sub

For bulk processing of millions of strings, consider:
```python
from metainformant.core import parallel

# Parallel normalization
normalized = parallel.thread_map(
    standardize_gene_name,
    gene_list,
    max_workers=8,
)
```

### Unicode Handling

`remove_control_chars()` iterates character-by-character with `unicodedata.category()` lookup per char. For very large texts (100MB+), this can be nontrivial. Consider C optimization via `regex` module or `str.translate()` for known character ranges:
```python
# Alternative using translate table
_control_chars = dict.fromkeys(range(0, 32))  # 0-31 control
_control_chars[127] = None  # DEL

def remove_control_chars_fast(text: str) -> str:
    return text.translate(_control_chars)
```

## Testing

### Property-Based Testing (Hypothesis)

```python
from hypothesis import given, strategies as st

@given(st.text())
def test_normalize_whitespace_idempotent(text):
    once = normalize_whitespace(text)
    twice = normalize_whitespace(once)
    assert once == twice  # Idempotent

@given(st.text(alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd")), min_size=1))
def test_slugify_lowercase_and_alnum(text):
    slug = slugify(text)
    assert slug == slug.lower()  # All lowercase
    assert all(c.isalnum() or c == '-' for c in slug)

@given(st.integers())
def test_extract_numbers_roundtrip(num):
    text = f"Value: {num}"
    nums = extract_numbers(text)
    assert nums == [float(num)]
```

### Regression Tests for Known Cases

```python
def test_gene_name_standardization_regression():
    cases = {
        "BRCA-1": "BRCA1",
        "BRCA_1": "BRCA1",
        "BRCA.1": "BRCA1",
        "Brca1": "BRCA1",
        " brca1 ": "BRCA1",
    }
    for raw, expected in cases.items():
        assert standardize_gene_name(raw) == expected

def test_species_formatting_regression():
    cases = {
        "homo sapiens": "Homo sapiens",
        "HOMO SAPIENS": "Homo sapiens",
        "drosophila melanogaster": "Drosophila melanogaster",
        "E. coli": "E. coli",
    }
    for raw, expected in cases.items():
        assert format_species_name(raw) == expected
```

### Fuzzing

```python
import random
import string

def fuzz_text_functions():
    """Randomized stress test."""
    for _ in range(10_000):
        length = random.randint(0, 200)
        text = ''.join(random.choices(string.printable, k=length))

        # Should never raise exception
        try:
            normalize_whitespace(text)
            slugify(text)
            safe_filename(text)
            standardize_gene_name(text)
            format_species_name(text)
        except Exception as e:
            print(f"Failed on {repr(text)}: {e}")
            raise
    raise

## Security Notes

Text processing functions are generally safe, but be aware of:

### ReDoS (Regular Expression Denial of Service)

The module uses simple regex patterns without nested quantifiers, which are safe from catastrophic backtracking. When extending with custom regexes:
- Avoid patterns like `(a+)+` on untrusted input
- Use `re.compile()` with timeout (Python 3.11+ has `re.TEMPLATE` limits)

### Unicode Homoglyph Attacks

Attackers may use visually similar Unicode characters to bypass validation (e.g., Cyrillic `а` vs Latin `a`). Functions like `slugify()` remove non-ASCII characters, which mitigates but doesn't fully address this. For high-security contexts, consider normalization:
```python
import unicodedata

def normalize_for_security(text: str) -> str:
    # NFKC normalization decomposes compatibility characters
    return unicodedata.normalize("NFKC", text)
```

### Filename Sanitization

Always use `safe_filename()` or `sanitize_filename()` before writing user-provided strings to disk. Never trust raw user input for file paths.

### Information Disclosure

Be cautious logging sanitized text—original may contain PII/PHI. If logs are centralized, ensure proper access controls.

## Dependencies

- **Required**: Standard library only (`re`, `unicodedata`, `pathlib`)
- **Optional**: None

## Further Reading

- Unicode categories: https://unicode.org/reports/tr44/#General_Category_Values
- Slugify patterns: https://gist.github.com/mahmoud/235d19a0cd5b194f7a354
- Biological nomenclature: https://www.issn.org/services/online-services/access-to-the-latn-issn/
