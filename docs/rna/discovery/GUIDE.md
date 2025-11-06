# Ant Species RNA-seq Discovery System

## Overview

This system automatically discovers all ant species with publicly available RNA-seq data and generates complete, accurate amalgkit YAML configurations with the latest genome assemblies from NCBI.

## Key Features

✅ **Comprehensive Discovery**: Searches all Formicidae (ant family) species in NCBI SRA  
✅ **Genome Integration**: Retrieves latest genome assemblies with annotations  
✅ **Accurate Metadata**: Real NCBI taxonomy IDs, assembly accessions, and FTP URLs  
✅ **Auto-Configuration**: Generates production-ready YAML files  
✅ **Quality Ranking**: Prioritizes RefSeq, chromosome-level assemblies  
✅ **Complete Reports**: Detailed markdown reports with statistics  

## Quick Start

### Prerequisites

```bash
# Install required dependencies (using uv)
uv pip install biopython ncbi-datasets-pylib

# Set NCBI email (required for Entrez API)
export NCBI_EMAIL="your.email@example.com"
```

### Basic Usage

```bash
# Discover all ant species with RNA-seq data (≥1 sample)
# Current method: Genus-based search (recommended)
python3 scripts/rna/discover_ant_rnaseq_by_genus.py

# Generate YAML configs for discovered species
python3 scripts/rna/generate_ant_configs_with_genomes.py

# Alternative: Direct species search (legacy)
python3 scripts/rna/discover_ant_species_with_rnaseq.py --min-samples 10
```

## What It Does

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
# Check the discovery report
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
bash scripts/rna/amalgkit/run_amalgkit.sh \
    --config config/amalgkit/amalgkit_camponotus_floridanus.yaml \
    --steps metadata,select,getfastq,quant,merge,curate,sanity

# Multiple species in parallel
for species in camponotus_floridanus solenopsis_invicta pogonomyrmex_barbatus; do
    nohup bash scripts/rna/amalgkit/run_amalgkit.sh \
        --config config/amalgkit/amalgkit_${species}.yaml \
        > output/${species}_workflow.log 2>&1 &
done
```

## Advanced Usage

### Filter by Sample Count

Only generate configurations for well-studied species:

```bash
# Species with ≥50 RNA-seq samples
python3 scripts/rna/discover_ant_species_with_rnaseq.py --min-samples 50
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

1. **Manual workflows**: `bash scripts/rna/amalgkit/run_amalgkit.sh`
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
python3 scripts/rna/discover_ant_species_with_rnaseq.py \
    --output-dir output/ant_discovery_$(date +%Y%m) \
    --min-samples 1
```

## Documentation

Generated configurations include:
- **Inline comments**: Explain each parameter
- **Assembly metadata**: Version, level, sequencing tech
- **Sample counts**: Number of available RNA-seq runs
- **FTP URLs**: Direct links to genome files

## Future Enhancements

Planned improvements:
- Validate genome FTP URLs automatically
- Check for genome annotation quality
- Estimate disk space requirements
- Generate species phylogeny for multi-species analyses
- Add functional annotation summaries

## Related Documentation

- **[WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution
- **[amalgkit/README.md](../amalgkit/README.md)**: Amalgkit guide
- **Configuration Template**: `config/amalgkit/amalgkit_template.yaml`
- **[MULTI_SPECIES_QUICK_START.md](../MULTI_SPECIES_QUICK_START.md)**: Multi-species quick start
- **[GETTING_STARTED.md](../GETTING_STARTED.md)**: Setup and installation
- **[ORCHESTRATION/README.md](../ORCHESTRATION/README.md)**: Orchestrator overview

## See Also

- **[QUICK_REF.md](QUICK_REF.md)**: Quick reference for discovery
- **[README.md](README.md)**: Discovery results and current status

---

**Last Updated**: November 3, 2025  
**Status**: ✅ Production-ready  
**Generated by**: AI-assisted development (grok-code-fast-1)  
**Tested with**: Biopython 1.81, ncbi-datasets-pylib 14.22.0

