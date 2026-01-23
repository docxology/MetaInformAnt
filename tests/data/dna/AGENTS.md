# Agent Directives: tests/data/dna

## Role
DNA sequence test fixtures for validating sequence analysis functions.

## Contents
Test data files for:
- FASTA sequence parsing and validation
- FASTQ quality score handling
- Multiple sequence alignment inputs
- Variant calling test cases

## Rules
- All sequences must be valid DNA (ACGT only, or with standard ambiguity codes)
- Include edge cases: empty sequences, single nucleotide, very long sequences
- Quality scores in FASTQ must be valid Phred scores
