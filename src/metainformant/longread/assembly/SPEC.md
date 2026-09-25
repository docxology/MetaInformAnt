# Specification: assembly

## 🎯 Scope
Long-read assembly module for overlap computation, consensus, and hybrid assembly.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `consensus.py`
- `hybrid.py`
- `overlap.py`

### Function contracts (behavior pinned by tests)

- `overlap.find_overlaps(reads, ...)` chains shared canonical minimizers in
  two frames per candidate pair: a constant `qpos - tpos` diagonal for
  same-strand overlaps (`strand = "+"`) and, after mapping target k-mer
  starts into the target's reverse-complement frame
  (`t' = target_length - k - tpos`), a constant diagonal for
  antiparallel overlaps (`strand = "-"`, with target coordinates mapped
  back to the target's forward orientation). The frame with the most
  chained matches wins; ties prefer the forward strand.
