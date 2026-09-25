# Specification: detection

## 🎯 Scope
Structural variant detection subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `breakpoints.py`
- `cnv.py`
- `sv_calling.py`

### Function contracts (behavior pinned by tests)

- `sv_calling.classify_sv_type(evidence)` makes a call only when strand and
  direction evidence support it: TRA for inter-chromosomal pairs; INS for
  split reads with nearly coincident breakpoints; INV for same-strand
  pairs with actual strand characters; DUP for everted (upstream '-',
  downstream '+') pairs after normalising to genomic order; DEL for
  head-to-head ('+', '-') pairs with split-read support or a span of at
  least `_DEL_MIN_BREAKPOINT_DISTANCE`. Weak evidence (short-span FR pairs
  without split support, unknown strand characters) returns
  `SVType.UNKNOWN`; the classifier never manufactures a deletion.
