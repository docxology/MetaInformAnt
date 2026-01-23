# Agent Directives: tests/data/ontology

## Role
Gene Ontology and biological ontology test fixtures.

## Directory Structure
- `GO_v3/` - Gene Ontology version 3 format test files

## Contents
Test data for:
- OBO file format parsing
- GO term hierarchy traversal
- Semantic similarity calculations
- Annotation propagation

## Rules
- OBO files must follow official OBO format specification
- Include representative subset of GO terms (not full ontology)
- Test data should cover all three GO aspects (BP, MF, CC)
- Include obsolete terms for handling edge cases
