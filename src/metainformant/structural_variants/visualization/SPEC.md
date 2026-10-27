# Specification: visualization

## 🎯 Scope
Structural variant visualization subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `plots.py`

### Function contracts (behavior pinned by tests)

- `plots.plot_circos(variants, chromosomes, output_path)` draws one chord per variant from
  `chrom:start` to its second breakpoint: `chrom2:pos2` when `pos2` is given (translocations from
  `call_structural_variants`), else `chrom2:end` for a hand-built dictionary whose `chrom2` is
  another chromosome (the VCF `END`-on-`CHR2` layout), else `chrom:end`. `chrom2` of `""` or `None`
  means `chrom`.
