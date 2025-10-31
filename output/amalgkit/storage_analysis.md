# AmalgKit Storage Analysis

**Date:** 2024-10-31

## Overall Storage Usage

| Species | Total Size | FASTQ Directory |
|---------|-----------|----------------|
| pbarbatus | 55G | 54G |
| cfloridanus | 35G | 34G |
| mpharaonis | 825M | 0B |
| sinvicta | 908M | 0B |

**Total:** ~91 GB

## P. barbatus (55G)

### Quantified Samples (Can Delete SRA)
**18 samples** with completed quantification still retaining SRA files:

| Sample | Size | Status |
|--------|------|--------|
| SRR14740487 | 7.1G | ✓ Quantified |
| SRR14740535 | 7.0G | ✓ Quantified |
| SRR14740551 | 6.5G | ✓ Quantified |
| SRR14740519 | 4.4G | ✓ Quantified |
| SRR1817186 | 2.0G | ✓ Quantified |
| SRR14740520 | 1.9G | ✓ Quantified |
| SRR1817187 | 1.9G | ✓ Quantified |
| SRR1817188 | 1.9G | ✓ Quantified |
| SRR1817185 | 1.7G | ✓ Quantified |
| SRR1817189 | 1.4G | ✓ Quantified |
| SRR14740543 | 1.3G | ✓ Quantified |
| SRR14740495 | 934M | ✓ Quantified |
| SRR1817183 | 666M | ✓ Quantified |
| SRR14740511 | 582M | ✓ Quantified |
| SRR14740503 | 390M | ✓ Quantified |
| SRR14740488 | 313M | ✓ Quantified |
| SRR14740536 | 279M | ✓ Quantified |
| SRR14740527 | 8.4M | ✓ Quantified |

**Reclaimable space: ~19.83 GB**

### Unquantified Samples (Need Processing)
**2 samples** with SRA files not yet quantified:

| Sample | Size | Status |
|--------|------|--------|
| DRR029869 | 363M | ✗ Not quantified |
| DRR048431 | 226M | ✗ Not quantified |

**Should be quantified before deletion: ~0.57 GB**

## C. floridanus (35G)

### All Samples (Need Quantification)
**18 samples** with SRA files, **NONE** quantified yet:

| Sample | Size | Status |
|--------|------|--------|
| SRR22031377 | 7.4G | ✗ Not quantified |
| SRR32143976 | 4.8G | ✗ Not quantified |
| SRR25496822 | 4.0G | ✗ Not quantified |
| SRR32143984 | 3.8G | ✗ Not quantified |
| SRR25496830 | 3.4G | ✗ Not quantified |
| SRR5931476 | 1.8G | ✗ Not quantified |
| SRR10960077 | 1.6G | ✗ Not quantified |
| SRR5931475 | 1.1G | ✗ Not quantified |
| SRR2060722 | 1.0G | ✗ Not quantified |
| SRR5931468 | 1.0G | ✗ Not quantified |
| SRR32143992 | 1.0G | ✗ Not quantified |
| SRR32143986 | 795M | ✗ Not quantified |
| SRR25496798 | 789M | ✗ Not quantified |
| SRR25496806 | 553M | ✗ Not quantified |
| SRR32143985 | 513M | ✗ Not quantified |
| SRR25496814 | 326M | ✗ Not quantified |
| SRR25496790 | 257M | ✗ Not quantified |
| SRR22031369 | 12M | ✗ Not quantified |

**Total needs quantification: ~32.12 GB**

Note: Only SRR32143984 has extracted FASTQ files visible (paired-end), but no quantification output.

## M. pharaonis (825M) & S. invicta (908M)

Both species have minimal storage footprint with no FASTQ/SRA files in the output directory.

## Recommendations

### Immediate Actions (P. barbatus)

1. **Delete 18 quantified SRA files** to reclaim ~19.83 GB:
   ```bash
   cd /Users/4d/Documents/GitHub/metainformant/output/amalgkit/pbarbatus/fastq/getfastq
   rm -f SRR14740487/SRR14740487.sra \
         SRR14740488/SRR14740488.sra \
         SRR14740495/SRR14740495.sra \
         SRR14740503/SRR14740503.sra \
         SRR14740511/SRR14740511.sra \
         SRR14740519/SRR14740519.sra \
         SRR14740520/SRR14740520.sra \
         SRR14740527/SRR14740527.sra \
         SRR14740535/SRR14740535.sra \
         SRR14740536/SRR14740536.sra \
         SRR14740543/SRR14740543.sra \
         SRR14740551/SRR14740551.sra \
         SRR1817183/SRR1817183.sra \
         SRR1817185/SRR1817185.sra \
         SRR1817186/SRR1817186.sra \
         SRR1817187/SRR1817187.sra \
         SRR1817188/SRR1817188.sra \
         SRR1817189/SRR1817189.sra
   ```

2. **Quantify remaining 2 P. barbatus samples** (DRR029869, DRR048431) to enable deletion of additional 0.57 GB

### Next Steps (C. floridanus)

1. **Quantify all 18 C. floridanus samples** (~32 GB of SRA data)
2. **Delete SRA files after successful quantification** to reclaim space

### Workflow Optimization

Consider implementing automatic SRA cleanup after successful quantification to prevent accumulation of redundant data:
- Keep FASTQ files temporarily during quantification
- Delete SRA files immediately after FASTQ extraction
- Delete FASTQ files after successful quantification (retaining only quant outputs)

## Total Reclaimable Space

- **Immediate (P. barbatus quantified):** 19.83 GB
- **After P. barbatus completion:** 20.40 GB
- **After C. floridanus quantification:** 52.52 GB

**Grand Total Potential Savings: ~52.5 GB**

