# Ant Species RNA-seq Discovery System

**Last Updated**: November 3, 2025

## Overview

Automated system for discovering all ant species with publicly available RNA-seq data and generating production-ready `amalgkit` YAML configurations.

## Discovery Results

### Summary Statistics

- **Total Species Discovered**: 55 ant species
- **Total RNA-seq Samples**: ~5,800 samples
- **Validated Genomes**: 20 species (manually curated)
- **Configs Generated**: 20 production-ready YAML files
- **Currently Processing**: Batch 1 (10 species, 3,820 samples)

### Species Tiers

**Batch 1 (Top 10 - Currently Running)**:
1. Harpegnathos saltator - 695 samples
2. Temnothorax longispinosus - 557 samples  
3. Solenopsis invicta - 451 samples
4. Monomorium pharaonis - 370 samples
5. Camponotus floridanus - 359 samples
6. Temnothorax rugatulus - 316 samples
7. Ooceraea biroi - 278 samples
8. Atta cephalotes - 239 samples
9. Cardiocondyla obscurior - 191 samples
10. Lasius niger - 191 samples

**Batch 2 (Queued - 10 species, 728 samples)**:
11-20. Additional species ready for processing

**Future Candidates**: 35 species with genomes requiring validation

## Scripts

### Discovery
- `scripts/rna/discover_ant_rnaseq_by_genus.py` - Find species with RNA-seq
- `scripts/rna/generate_ant_configs_with_genomes.py` - Generate configs

### Execution
- `scripts/rna/run_top10_ant_species.sh` - Launch Batch 1
- `scripts/rna/run_batch2_ant_species.sh` - Launch Batch 2

### Monitoring  
- `scripts/rna/get_current_status.py` - Check workflow status

## Discovery Data

All discovery data stored in `output/ant_discovery/`:
- `ant_species_rnaseq_data.json` - Complete discovery results
- `ant_species_with_genomes.json` - 20 validated species
- `ant_species_without_genomes.json` - 35 pending validation

## Generated Configurations

26 YAML configs in `config/amalgkit/`:
- 20 newly discovered species (validated genomes)
- 6 pre-existing configs (sinvicta, mpharaonis, cfloridanus, pbarbatus, test, template)

## Current Processing Status

**Batch 1**: Running (started Nov 3, 2025 16:18 PST)
- 10 workflows active
- 100+ processes running (orchestrators + workers)
- Auto-cleanup enabled (FASTQs deleted after quantification)
- ETA: 24-48 hours

**Batch 2**: Ready to launch after Batch 1 completes

## Usage

See parent directory documentation:
- `../ANT_SPECIES_DISCOVERY.md` - Comprehensive guide
- `../ANT_DISCOVERY_QUICK_REF.md` - Quick reference
- `../README.md` - RNA analysis overview

## Scientific Impact

### Phylogenetic Diversity
- 13 genera represented
- 3 subfamilies: Formicinae, Myrmicinae, Ponerinae
- Diverse ecological niches

### Research Applications
- Caste differentiation genomics
- Social behavior evolution
- Invasive species biology
- Fungus-growing symbiosis
- Clonal reproduction
- Aging and longevity studies

---

For detailed workflow documentation, see `../amalgkit/README.md`

