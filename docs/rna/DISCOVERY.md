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

## Quick Reference

### One-Line Commands

```bash
# Basic discovery (all ants with RNA-seq data)
python3 scripts/rna/discover_species.py

# Discovery results stored in output/ant_discovery/
# See below for filtering and configuration options

# Deploy all generated configs
cp output/ant_discovery/configs/*.yaml config/amalgkit/

# Deploy specific species
cp output/ant_discovery/configs/amalgkit_camponotus_floridanus.yaml config/amalgkit/

# Run workflow for discovered species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_SPECIES.yaml
```

### Prerequisites

```bash
# Using uv (recommended)
uv pip install biopython ncbi-datasets-pylib
export NCBI_EMAIL="your.email@example.com"
```

### Typical Discovery Results

Expected to find 30-100 ant species with RNA-seq data, including:

**High Priority** (>100 samples, chromosome genomes):
- Camponotus floridanus (Florida carpenter ant) - ~300 samples
- Solenopsis invicta (Red fire ant) - ~350 samples

**Medium Priority** (10-100 samples, good genomes):
- Pogonomyrmex barbatus (Red harvester ant) - ~80 samples
- Monomorium pharaonis (Pharaoh ant) - ~100 samples
- Acromyrmex echinatior (Leafcutter ant)
- Atta cephalotes (Leafcutter ant)

**Emerging** (1-10 samples):
- Many additional species with initial studies

## Discovery Workflow

### Step 1: RNA-seq Data Discovery

Searches NCBI SRA for:
- **Taxonomy**: All species in Formicidae (TaxID: 7389)
- **Strategy**: RNA-Seq experiments
- **Platform**: Illumina sequencing
- **Access**: Public data only

For each species, collects:
- Scientific name and taxonomy ID
- Number of RNA-seq samples
- SRA run IDs (SRR*)
- Study IDs (SRP*, ERP*, DRP*)

### Step 2: Genome Assembly Retrieval

For each discovered species:
- Queries NCBI Datasets API for genome assemblies
- Retrieves assembly metadata:
  - Accession (GCF/GCA)
  - Assembly name and level (chromosome/scaffold/contig)
  - Sequencing technology
  - Contig/scaffold N50 statistics
  - Annotation release version
- Selects best assembly using prioritization:
  1. RefSeq (GCF) > GenBank (GCA)
  2. Chromosome-level > Scaffold > Contig
  3. Latest release date
  4. Highest contig N50

### Step 3: YAML Configuration Generation

For each species, generates a complete `amalgkit_SPECIES.yaml` with:
- **Species metadata**: Scientific name, taxonomy ID, sample counts
- **Genome configuration**: 
  - Assembly accession and name
  - FTP URLs for genome files
  - Specific file paths (genomic, transcriptome, CDS, protein, GFF)
- **Workflow steps**: Pre-configured for all 11 amalgkit steps
- **Optimized settings**:
  - 12 threads for parallel processing
  - Cloud acceleration (AWS, GCP, NCBI)
  - Auto-delete FASTQs after quantification
  - Auto-build Kallisto index

### Step 4: Comprehensive Reporting

Generates:
- **DISCOVERY_REPORT.md**: Summary with tables of all species
- **ant_species_rnaseq_data.json**: Raw metadata for programmatic access
- **configs/**: Directory with all generated YAML files

## Output Structure

```
output/ant_discovery/
├── DISCOVERY_REPORT.md          # Human-readable summary
├── ant_species_rnaseq_data.json # Machine-readable data
└── configs/
    ├── amalgkit_camponotus_floridanus.yaml
    ├── amalgkit_solenopsis_invicta.yaml
    ├── amalgkit_pogonomyrmex_barbatus.yaml
    ├── amalgkit_monomorium_pharaonis.yaml
    └── ... (all discovered species)
```

## Example Generated YAML

Here's what a generated configuration looks like:

```yaml
# METAINFORMANT Amalgkit Configuration
# Species: Camponotus floridanus
# NCBI Taxonomy ID: 104421
# Assembly: GCF_003227725.1_Cflo_v7.5 (Cflo_v7.5)
# Assembly Level: Chromosome
# Sequencing: PacBio
# Annotation Release: 101
# RNA-seq Samples Available: 307
# Generated: 2025-11-03

work_dir: /home/q/Documents/GitHub/MetaInformAnt/output/amalgkit/camponotus_floridanus/work
log_dir: /home/q/Documents/GitHub/MetaInformAnt/output/amalgkit/camponotus_floridanus/logs
threads: 12

auto_install_amalgkit: true

filters:
  require_tissue: false

species_list:
  - Camponotus_floridanus

taxon_id: 104421

genome:
  accession: GCF_003227725.1
  assembly_name: Cflo_v7.5
  annotation_release: 101
  dest_dir: output/amalgkit/camponotus_floridanus/genome
  
  include:
    - genome
    - gff3
    - rna
    - cds
    - protein
    - seq-report
    - feature-table
    - gene-ontology
  
  ftp_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/GCF_003227725.1_Cflo_v7.5/
  
  files:
    genomic_fasta: GCF_003227725.1_Cflo_v7.5_genomic.fna.gz
    transcriptome_fasta: GCF_003227725.1_Cflo_v7.5_rna_from_genomic.fna.gz
    cds_fasta: GCF_003227725.1_Cflo_v7.5_cds_from_genomic.fna.gz
    protein_fasta: GCF_003227725.1_Cflo_v7.5_protein.faa.gz
    annotation_gff: GCF_003227725.1_Cflo_v7.5_genomic.gff.gz
    annotation_gtf: GCF_003227725.1_Cflo_v7.5_genomic.gtf.gz

steps:
  metadata:
    out_dir: output/amalgkit/camponotus_floridanus/work
    search_string: '"Camponotus floridanus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
    redo: yes
  # ... (complete step configuration)
```

## Using Generated Configurations

### 1. Review and Validate

```bash
# Check the discovery report (program-generated output)
cat output/ant_discovery/DISCOVERY_REPORT.md

# Review a specific configuration
cat output/ant_discovery/configs/amalgkit_solenopsis_invicta.yaml
```

### 2. Deploy Configurations

```bash
# Copy desired configurations to main config directory
cp output/ant_discovery/configs/amalgkit_*.yaml config/amalgkit/

# Or copy selectively
cp output/ant_discovery/configs/amalgkit_camponotus_floridanus.yaml config/amalgkit/
```

### 3. Run Workflows

```bash
# Single species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_camponotus_floridanus.yaml

# Multiple species in parallel
for species in camponotus_floridanus solenopsis_invicta pogonomyrmex_barbatus; do
    nohup python3 scripts/rna/run_workflow.py \
        --config config/amalgkit/amalgkit_${species}.yaml \
        > output/${species}_workflow.log 2>&1 &
done
```

## Batch Processing

### Immediate Per-Sample Processing Strategy

Processing 20 ant species with immediate per-sample processing (download → immediately quantify → immediately delete FASTQs) using 24 total threads distributed across all species.

**Rationale:**
- ✅ **Maximum Disk Efficiency**: Only one sample's FASTQs exist at any time  
- ✅ **Immediate Processing**: Download → immediately quantify → immediately delete  
- ✅ **Total Thread Allocation**: 24 threads distributed across all species (not per species)  
- ✅ **Dynamic Redistribution**: Threads redistribute as species complete  
- ✅ **Scalability**: Can process all species simultaneously without disk space issues  

### Processing Strategy

- **All Species**: 20 species processed simultaneously
- **Total Samples**: 1,578 samples across all species
- **Thread Allocation**: 24 threads TOTAL distributed evenly (minimum 1 per species)
- **Processing Mode**: Each sample: download → immediately quantify → immediately delete FASTQs
- **Script**: `scripts/rna/run_workflow.py` (configure `num_download_workers` in each species config file)

### Disk Management

**Per Sample** (immediate processing):
1. Download FASTQ (~2-10 GB)
2. **Immediately** quantify with Kallisto
3. ✅ **Immediately delete FASTQ** after quantification
4. Keep results only (~10-50 MB)
5. Move to next sample

**Peak Usage**: ~2-10 GB (only one sample's FASTQs exist at a time)  
**Final Usage**: ~40-55 GB (all results across all species)

### Timeline

| Phase | Duration | Batch 1 | Batch 2 |
|-------|----------|---------|---------|
| Metadata | 1-2 hours | ✅ | ⏳ |
| Downloads | 12-36 hours | 🔄 | ⏳ |
| Quantification | 6-18 hours | ⏳ | ⏳ |
| Merge + QC | 2-4 hours | ⏳ | ⏳ |
| **Total** | **24-48 hours** | | **12-24 hours** |

**Complete Project**: 36-72 hours for all 20 species

### Launch Commands

```bash
# Immediate processing with 24 threads TOTAL distributed across all species (recommended)
# Run separately for each species (recommended)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species2.yaml

# Or run in parallel (background)
nohup python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species1.yaml > logs/species1.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species2.yaml > logs/species2.log 2>&1 &
```

### Monitoring

```bash
# Check status
python3 scripts/rna/orchestrate_workflows.py --status

# View logs
tail -f output/top10_*.log

# Check disk
df -h output/

# Active processes
ps aux | grep run_workflow | wc -l
```

## Advanced Usage

### Filter by Sample Count

Only generate configurations for well-studied species:

```bash
# Species with ≥50 RNA-seq samples
python3 scripts/rna/discover_species.py
```

### Prioritize Species

The discovery report automatically ranks species by:
1. Genome availability (with genome > without)
2. Sample count (more samples = better powered analyses)
3. Assembly quality (chromosome > scaffold > contig)

Example report excerpt:

```markdown
## Species with Genomes (Recommended for Analysis)

| Species | Samples | Genome Accession | Assembly Level | Annotation |
|---------|---------|------------------|----------------|------------|
| Camponotus floridanus | 307 | GCF_003227725.1 | Chromosome | 101 |
| Solenopsis invicta | 354 | GCF_016802725.1 | Chromosome | 101 |
| Pogonomyrmex barbatus | 83 | GCF_000187915.1 | Scaffold | 101 |
```

### Programmatic Access

Use the JSON output for custom analysis:

```python
import json

with open('output/ant_discovery/ant_species_rnaseq_data.json') as f:
    species_data = json.load(f)

# Find species with >100 samples and chromosome-level genomes
high_quality = [
    name for name, data in species_data.items()
    if data['sample_count'] > 100
    and data.get('genome', {}).get('level') == 'Chromosome'
]

print(f"High-quality species: {high_quality}")
```

## Genome Assembly Selection Logic

The script uses intelligent ranking to select the best assembly:

### Scoring System

| Criterion | Points | Rationale |
|-----------|--------|-----------|
| RefSeq (GCF_*) | +1000 | More curated, higher quality |
| GenBank (GCA_*) | +0 | Original submissions |
| Chromosome-level | +500 | Most complete |
| Scaffold-level | +100 | Good contiguity |
| Contig-level | +10 | Basic assembly |
| High contig N50 | +0 to +100 | Better contiguity (scaled) |

### Example Selection

Given assemblies:
```
1. GCA_XXXXXX.1 (GenBank, Scaffold, N50=1Mb) → Score: 100
2. GCF_YYYYYY.1 (RefSeq, Contig, N50=100kb) → Score: 1010
3. GCF_ZZZZZZ.1 (RefSeq, Chromosome, N50=10Mb) → Score: 1600 ✅ SELECTED
```

## Scripts

### Discovery
- `scripts/rna/discover_species.py` - Find species with RNA-seq data
- Discovery results stored in `output/ant_discovery/`

### Execution
- `scripts/rna/run_workflow.py` - Launch workflow for each species (recommended)
  - Run separately for each species config
  - Parallel downloads configured via `num_download_workers` in each config file

### Monitoring  
- `scripts/rna/orchestrate_workflows.py --status` - Check workflow status

## Current Processing Status

**Batch 1**: Running (started Nov 3, 2025 16:18 PST)
- 10 workflows active
- 100+ processes running (orchestrators + workers)
- Auto-cleanup enabled (FASTQs deleted after quantification)
- ETA: 24-48 hours

**Batch 2**: Ready to launch after Batch 1 completes

## Troubleshooting

### No Species Found

**Problem**: Script reports 0 species found

**Solutions**:
- Check NCBI_EMAIL is set: `echo $NCBI_EMAIL`
- Verify internet connection
- Try again later (NCBI API rate limits)

### Genome Lookup Fails

**Problem**: "ncbi-datasets-pylib not installed"

**Solution**:
```bash
# Using uv (recommended)
uv pip install ncbi-datasets-pylib

# Or with conda (alternative for system-level installation)
conda install -c conda-forge ncbi-datasets-pylib
```

### FTP URLs Don't Work

**Problem**: Generated FTP URLs return 404

**Solutions**:
- Manually verify assembly exists on NCBI FTP
- Assembly may have been deprecated
- Try alternative assembly for that species
- Update `ftp_url` field manually in YAML

## Expected Ant Species

The script should discover at least:

### Well-Studied Species (≥100 samples)
- **Camponotus floridanus** - Florida carpenter ant
- **Solenopsis invicta** - Red fire ant  
- **Acromyrmex echinatior** - Leafcutter ant
- **Atta cephalotes** - Leafcutter ant

### Moderately-Studied (10-100 samples)
- **Pogonomyrmex barbatus** - Red harvester ant
- **Monomorium pharaonis** - Pharaoh ant
- **Harpegnathos saltator** - Indian jumping ant
- **Vollenhovia emeryi** - Asian ant

### Emerging Models (1-10 samples)
- Many species with initial transcriptome studies

## Performance

Typical runtime:
- **SRA search**: 2-5 minutes (for ~10,000 SRA records)
- **Genome lookup**: ~5-10 seconds per species
- **YAML generation**: <1 second per species
- **Total**: 10-15 minutes for ~50 species

## Integration with Existing Workflows

Generated configurations are fully compatible with:

1. **Manual workflows**: `python3 scripts/rna/run_workflow.py`
2. **Orchestrator**: `python3 scripts/rna/orchestrate_workflows.py`
3. **Status monitoring**: `python3 scripts/rna/orchestrate_workflows.py --status`
4. **Batch processing**: Works with existing batch scripts

## Updating Configurations

Rerun discovery periodically to:
- Capture newly published RNA-seq datasets
- Get updated genome assemblies
- Incorporate improved annotations

```bash
# Monthly update
python3 scripts/rna/discover_species.py \
    --output-dir output/ant_discovery_$(date +%Y%m) \
    --min-samples 1
```

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

## Related Documentation

- **[workflow.md](workflow.md)**: Workflow planning and execution
- **[amalgkit/README.md](amalgkit/README.md)**: Amalgkit guide
- **[ORCHESTRATION.md](ORCHESTRATION.md)**: Orchestrator overview
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

## See Also

- **Configuration Template**: `config/amalgkit/amalgkit_template.yaml`
- **Discovery Data**: `output/ant_discovery/ant_species_rnaseq_data.json`
- **Generated Configs**: `output/ant_discovery/configs/`

---

**Last Updated**: November 3, 2025  
**Status**: ✅ Production-ready  
**Tested with**: Biopython 1.81, ncbi-datasets-pylib 14.22.0

