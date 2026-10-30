# Specification: core

## Scope
Core utilities for METAINFORMANT bioinformatics toolkit.

## Architecture
- **Dependency Level**: Core
- **Component Type**: Source Code

## Data Structures
- **Sub-packages**: data, engine, execution, io, ui, utils
- **Key Concepts**: I/O, logging, config, execution engine, UI

## API Definition
### Exports
- `__init__.py`

## Verification observations

`io.verification.BoundReads` hashes each physical absolute path once per
operation and rechecks its device/inode/mode/ownership/size/mtime/ctime at
reuse and successful context exit. Every new instance reads bytes again.
Changes or redirected/nonregular paths refuse. This observation does not
replace locks, backup verification, or deployment approval.
