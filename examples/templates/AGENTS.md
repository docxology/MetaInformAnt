# Agent Directives: examples/templates

## Role
Template files for creating new examples.

## Contents
- `base_template.py.j2` - Basic example template
- `dna_template.py.j2` - DNA module example template
- `ml_template.py.j2` - Machine learning example template

## Usage
Copy the appropriate template when creating a new example:

```bash
python scripts/test_examples/generate_example.py {module} {feature}
```

Then customize with module-specific imports and functionality.

## Template Structure
Templates include:
- Module docstring
- Standard imports
- Main function scaffold
- Output handling
- Error handling patterns
