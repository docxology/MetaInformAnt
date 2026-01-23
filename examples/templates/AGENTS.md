# Agent Directives: examples/templates

## Role
Template files for creating new examples.

## Contents
- `base_template.py` - Basic example template
- `dna_template.py` - DNA module example template
- `ml_template.py` - Machine learning example template

## Usage
Copy the appropriate template when creating a new example:

```bash
cp examples/templates/base_template.py examples/{module}/example_{feature}.py
```

Then customize with module-specific imports and functionality.

## Template Structure
Templates include:
- Module docstring
- Standard imports
- Main function scaffold
- Output handling
- Error handling patterns
