# Agent Directives: config/config_base

## Role
Base configuration files used by amalgkit for sample filtering and grouping.

## Contents
- `control_term.config` - Terms identifying control samples
- `exclude_keyword.config` - Keywords to exclude samples
- `group_attribute.config` - Attributes for sample grouping

## Usage
These files are referenced by amalgkit during the metadata and curation steps to:
- Filter out unwanted samples (controls, problematic conditions)
- Group samples by biological attributes
- Standardize sample classification

## Format
Plain text files with one entry per line. Comments start with `#`.
