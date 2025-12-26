# Core: Text Processing Utilities

The `text` module provides comprehensive text processing and normalization utilities for biological data, including filename safety, identifier cleaning, and format standardization.

## Functions

### Basic Text Processing
- **`normalize_whitespace(text)`** → `str`
  - Collapse all whitespace to single spaces and strip ends
  - Handles tabs, newlines, and multiple spaces

- **`clean_whitespace(text)`** → `str`
  - Alias for normalize_whitespace (legacy compatibility)

- **`remove_control_chars(text)`** → `str`
  - Remove control characters while preserving tabs and newlines
  - Uses Unicode categories for proper character classification

### Identifier Formatting
- **`slugify(text)`** → `str`
  - Convert text to URL-safe slugs (lowercase, dashes, alphanumeric)
  - Remove special characters and normalize separators

- **`safe_filename(name)`** → `str`
  - Create filesystem-safe filenames while preserving extensions
  - Combines slugify with proper file extension handling

### Biological Data Formatting
- **`standardize_gene_name(gene_name)`** → `str`
  - Standardize gene names (uppercase, remove separators)
  - Normalize common gene identifier formats

- **`format_species_name(species_name)`** → `str`
  - Format species names in proper binomial nomenclature
  - Capitalize genus, lowercase species

- **`clean_sequence_id(sequence_id)`** → `str`
  - Extract clean sequence identifiers from FASTA headers
  - Handle NCBI accession formats and complex headers

### Text Analysis
- **`extract_numbers(text)`** → `List[float]`
  - Extract all numeric values from text
  - Handle decimal numbers and scientific notation

- **`count_words(text)`** → `int`
  - Count words after whitespace normalization

- **`truncate_text(text, max_length, suffix="...")`** → `str`
  - Truncate text to specified length with optional suffix

### Pattern Extraction
- **`extract_email_addresses(text)`** → `List[str]`
  - Extract valid email addresses using regex pattern matching

## Usage Examples

### Basic Text Cleaning
```python
from metainformant.core import text

# Normalize messy text
messy = "  Hello   world \n\t with\t\ttabs  "
clean = text.normalize_whitespace(messy)
print(clean)  # "Hello world with tabs"

# Create URL-safe slugs
title = "My Research: COVID-19 Analysis!"
slug = text.slugify(title)
print(slug)  # "my-research-covid-19-analysis"

# Safe filenames
unsafe = "Report: 2023/Genomic Analysis?.pdf"
safe = text.safe_filename(unsafe)
print(safe)  # "report-2023-genomic-analysis.pdf"
```

### Biological Identifier Processing
```python
from metainformant.core import text

# Gene names
raw_gene = "BRCA-1"
standardized = text.standardize_gene_name(raw_gene)
print(standardized)  # "BRCA1"

# Species names
species = "apis mellifera"
formatted = text.format_species_name(species)
print(formatted)  # "Apis mellifera"

# Sequence IDs from FASTA headers
header = ">gi|12345|ref|NM_001| Homo sapiens BRCA1"
clean_id = text.clean_sequence_id(header)
print(clean_id)  # "NM_001"
```

### Text Analysis and Extraction
```python
from metainformant.core import text

# Extract numbers from text
data_text = "Expression: 15.3 ± 2.1, p-value: 0.05"
numbers = text.extract_numbers(data_text)
print(numbers)  # [15.3, 2.1, 0.05]

# Count words
article = "This is a sample article with several words."
word_count = text.count_words(article)
print(word_count)  # 8

# Truncate long text
long_text = "This is a very long description that needs to be shortened for display"
short = text.truncate_text(long_text, 30)
print(short)  # "This is a very long descript..."
```

### Email and Pattern Extraction
```python
from metainformant.core import text

# Extract emails from text
contact_info = "Contact: john.doe@example.com or jane.smith@university.edu"
emails = text.extract_email_addresses(contact_info)
print(emails)  # ['john.doe@example.com', 'jane.smith@university.edu']
```

## Biological Data Integration

### FASTA Processing
```python
from metainformant.core import text

def process_fasta_header(header_line):
    """Process FASTA header for consistent IDs."""
    sequence_id = text.clean_sequence_id(header_line)
    safe_filename = text.safe_filename(f"{sequence_id}.fasta")
    return sequence_id, safe_filename

# Example usage
header = ">ref|NM_001302504.2| Homo sapiens BRCA1"
seq_id, filename = process_fasta_header(header)
print(f"ID: {seq_id}, File: {filename}")
```

### Metadata Cleaning
```python
from metainformant.core import text

def clean_sample_metadata(metadata):
    """Clean and standardize sample metadata."""
    cleaned = {}
    for key, value in metadata.items():
        if isinstance(value, str):
            # Clean text fields
            cleaned[key] = text.normalize_whitespace(value)
            # Format species names
            if 'species' in key.lower():
                cleaned[key] = text.format_species_name(cleaned[key])
        else:
            cleaned[key] = value
    return cleaned

# Usage
raw_metadata = {
    "species": "drosophila melanogaster",
    "description": "  Test sample  \n\t with extra whitespace  ",
    "expression": 15.7
}
clean_metadata = clean_sample_metadata(raw_metadata)
```

## Error Handling

Text processing functions are designed to be robust:
- Handle None inputs gracefully (where appropriate)
- Preserve empty strings and edge cases
- Use safe regex patterns to prevent ReDoS attacks
- Provide sensible defaults for malformed input

## Performance Considerations

- **Regex Compilation**: Common patterns are pre-compiled for efficiency
- **Unicode Handling**: Proper Unicode support for international text
- **Memory Efficiency**: Process text in-place where possible
- **Streaming Compatible**: Functions work with individual strings

## Dependencies

- **Required**: Standard library `re`, `unicodedata`, `pathlib`
- **Optional**: None (all functionality uses stdlib)
