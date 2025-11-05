#!/usr/bin/env python3
"""
Discover ant species with RNA-seq data and generate amalgkit YAML configurations.

This script:
1. Searches NCBI SRA for all ant species with RNA-seq data
2. Retrieves genome assembly information from NCBI
3. Validates genome completeness and annotation
4. Generates species-specific amalgkit YAML configurations
5. Creates comprehensive reports of available ant RNA-seq resources

Usage:
    python3 scripts/rna/discover_ant_species_with_rnaseq.py [--output-dir OUTPUT_DIR] [--min-samples N]
"""

import argparse
import json
import logging
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError
from urllib.request import urlopen

try:
    from Bio import Entrez
except ImportError:
    print("ERROR: Biopython not installed. Install with: pip install biopython")
    sys.exit(1)

try:
    from ncbi.datasets import GenomeApi
    from ncbi.datasets.openapi import ApiClient
    NCBI_DATASETS_AVAILABLE = True
except ImportError:
    NCBI_DATASETS_AVAILABLE = False
    print("WARNING: ncbi-datasets-pylib not installed. Genome data will be limited.")
    print("Install with: pip install ncbi-datasets-pylib")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)
logger = logging.getLogger(__name__)

# Set your email for NCBI Entrez
Entrez.email = "your.email@example.com"  # Override from environment


class AntSpeciesDiscoverer:
    """Discover ant species with RNA-seq data and generate configurations."""
    
    def __init__(self, output_dir: Path, min_samples: int = 1):
        """Initialize the discoverer.
        
        Args:
            output_dir: Directory for output files
            min_samples: Minimum RNA-seq samples required
        """
        self.output_dir = Path(output_dir)
        self.min_samples = min_samples
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get email from environment
        import os
        email = os.environ.get('NCBI_EMAIL') or os.environ.get('ENTREZ_EMAIL')
        if email:
            Entrez.email = email
    
    def search_ants_with_rnaseq(self) -> Dict[str, Dict]:
        """Search NCBI SRA for all ant species with RNA-seq data.
        
        Returns:
            Dictionary mapping species names to metadata
        """
        logger.info("Searching NCBI SRA for ant species with RNA-seq data...")
        
        # Formicidae is the family that includes all ants
        # NCBI Taxonomy ID: 7389
        search_query = (
            'txid7389[Organism:exp] AND '
            '"RNA-Seq"[Strategy] AND '
            '"Illumina"[Platform] AND '
            '"public"[Access]'
        )
        
        logger.info(f"Search query: {search_query}")
        
        try:
            # Search SRA
            handle = Entrez.esearch(
                db="sra",
                term=search_query,
                retmax=10000,  # Get up to 10,000 results
                usehistory="y"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            total_records = int(search_results["Count"])
            logger.info(f"Found {total_records} RNA-seq records for Formicidae")
            
            if total_records == 0:
                logger.warning("No RNA-seq data found for ants")
                return {}
            
            # Fetch details in batches
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]
            
            species_data = defaultdict(lambda: {
                'sample_count': 0,
                'run_ids': [],
                'study_ids': set(),
                'taxonomy_id': None,
                'scientific_name': None
            })
            
            batch_size = 500
            for start in range(0, total_records, batch_size):
                end = min(start + batch_size, total_records)
                logger.info(f"Fetching records {start+1} to {end} of {total_records}...")
                
                try:
                    fetch_handle = Entrez.efetch(
                        db="sra",
                        rettype="xml",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key
                    )
                    records = Entrez.read(fetch_handle)
                    fetch_handle.close()
                    
                    # Parse records
                    for record in records:
                        try:
                            # Extract organism information
                            organism_info = self._extract_organism_info(record)
                            if organism_info:
                                species_name = organism_info['scientific_name']
                                species_data[species_name]['sample_count'] += 1
                                species_data[species_name]['taxonomy_id'] = organism_info['taxonomy_id']
                                species_data[species_name]['scientific_name'] = species_name
                                
                                # Extract run and study IDs
                                run_ids = self._extract_run_ids(record)
                                study_ids = self._extract_study_ids(record)
                                
                                species_data[species_name]['run_ids'].extend(run_ids)
                                species_data[species_name]['study_ids'].update(study_ids)
                        except Exception as e:
                            logger.debug(f"Error parsing record: {e}")
                            continue
                
                except Exception as e:
                    logger.error(f"Error fetching batch {start}-{end}: {e}")
                    continue
            
            # Filter by minimum samples
            filtered_species = {
                name: data for name, data in species_data.items()
                if data['sample_count'] >= self.min_samples
            }
            
            # Convert sets to lists for JSON serialization
            for species_name in filtered_species:
                filtered_species[species_name]['study_ids'] = list(
                    filtered_species[species_name]['study_ids']
                )
            
            logger.info(f"Found {len(filtered_species)} ant species with ≥{self.min_samples} RNA-seq samples")
            
            return dict(filtered_species)
            
        except Exception as e:
            logger.error(f"Error searching SRA: {e}")
            return {}
    
    def _extract_organism_info(self, record: Dict) -> Optional[Dict]:
        """Extract organism information from SRA record.
        
        Args:
            record: SRA record from Entrez
            
        Returns:
            Dictionary with organism info or None
        """
        try:
            # Navigate SRA XML structure
            pool = record.get('Pool', {})
            if 'Member' in pool:
                members = pool['Member']
                if not isinstance(members, list):
                    members = [members]
                
                for member in members:
                    organism = member.get('organism', {})
                    if organism:
                        return {
                            'scientific_name': organism.get('ScientificName', ''),
                            'taxonomy_id': organism.get('taxid', '')
                        }
            
            # Try alternative structure
            if 'EXPERIMENT_PACKAGE_SET' in record:
                exp_pkg = record['EXPERIMENT_PACKAGE_SET']
                if 'EXPERIMENT_PACKAGE' in exp_pkg:
                    exp_pkg = exp_pkg['EXPERIMENT_PACKAGE']
                    if not isinstance(exp_pkg, list):
                        exp_pkg = [exp_pkg]
                    
                    for pkg in exp_pkg:
                        sample = pkg.get('SAMPLE', {})
                        sample_name = sample.get('SAMPLE_NAME', {})
                        if 'SCIENTIFIC_NAME' in sample_name:
                            taxon_id = sample_name.get('TAXON_ID', '')
                            return {
                                'scientific_name': sample_name['SCIENTIFIC_NAME'],
                                'taxonomy_id': str(taxon_id)
                            }
        except Exception as e:
            logger.debug(f"Error extracting organism info: {e}")
        
        return None
    
    def _extract_run_ids(self, record: Dict) -> List[str]:
        """Extract SRA run IDs from record."""
        run_ids = []
        try:
            pool = record.get('Pool', {})
            if 'Member' in pool:
                members = pool['Member']
                if not isinstance(members, list):
                    members = [members]
                
                for member in members:
                    run_id = member.get('@accession', '')
                    if run_id and run_id.startswith('SRR'):
                        run_ids.append(run_id)
        except Exception:
            pass
        return run_ids
    
    def _extract_study_ids(self, record: Dict) -> set:
        """Extract study IDs from record."""
        study_ids = set()
        try:
            if 'EXPERIMENT_PACKAGE_SET' in record:
                exp_pkg = record['EXPERIMENT_PACKAGE_SET']
                if 'EXPERIMENT_PACKAGE' in exp_pkg:
                    exp_pkg = exp_pkg['EXPERIMENT_PACKAGE']
                    if not isinstance(exp_pkg, list):
                        exp_pkg = [exp_pkg]
                    
                    for pkg in exp_pkg:
                        study = pkg.get('STUDY', {})
                        study_id = study.get('@accession', '')
                        if study_id:
                            study_ids.add(study_id)
        except Exception:
            pass
        return study_ids
    
    def get_genome_info(self, taxonomy_id: str, species_name: str) -> Optional[Dict]:
        """Get genome assembly information for a species.
        
        Args:
            taxonomy_id: NCBI taxonomy ID
            species_name: Scientific name
            
        Returns:
            Genome information dictionary or None
        """
        if not NCBI_DATASETS_AVAILABLE:
            logger.warning(f"Skipping genome lookup for {species_name} (ncbi-datasets-pylib not installed)")
            return None
        
        logger.info(f"Fetching genome info for {species_name} (TaxID: {taxonomy_id})...")
        
        try:
            with ApiClient() as api_client:
                genome_api = GenomeApi(api_client)
                
                # Get assemblies for this taxon
                genome_metadata = genome_api.assembly_descriptors_by_taxon(
                    taxon=taxonomy_id,
                    filters_reference_only=False,  # Include non-reference assemblies
                    filters_assembly_source="RefSeq",  # Prefer RefSeq assemblies
                    page_size=10
                )
                
                assemblies = genome_metadata.get("assemblies", [])
                
                if not assemblies:
                    # Try GenBank if RefSeq not available
                    genome_metadata = genome_api.assembly_descriptors_by_taxon(
                        taxon=taxonomy_id,
                        filters_reference_only=False,
                        filters_assembly_source="GenBank",
                        page_size=10
                    )
                    assemblies = genome_metadata.get("assemblies", [])
                
                if not assemblies:
                    logger.warning(f"No genome assemblies found for {species_name}")
                    return None
                
                # Select best assembly (prefer RefSeq, then latest, then most complete)
                best_assembly = self._select_best_assembly(assemblies)
                
                if not best_assembly:
                    return None
                
                assembly_info = best_assembly.get('assembly', {})
                assembly_stats = best_assembly.get('assembly_stats', {})
                annotation_info = best_assembly.get('annotation_info', {})
                
                genome_info = {
                    'accession': assembly_info.get('assembly_accession', ''),
                    'assembly_name': assembly_info.get('assembly_name', ''),
                    'level': assembly_info.get('assembly_level', ''),
                    'release_date': assembly_info.get('submission_date', ''),
                    'contig_n50': assembly_stats.get('contig_n50', 0),
                    'scaffold_n50': assembly_stats.get('scaffold_n50', 0),
                    'total_sequence_length': assembly_stats.get('total_sequence_length', 0),
                    'number_of_contigs': assembly_stats.get('number_of_contigs', 0),
                    'annotation_release': annotation_info.get('release_version', ''),
                    'annotation_date': annotation_info.get('release_date', ''),
                    'sequencing_tech': assembly_info.get('sequencing_tech', ''),
                }
                
                return genome_info
                
        except Exception as e:
            logger.error(f"Error fetching genome info for {species_name}: {e}")
            return None
    
    def _select_best_assembly(self, assemblies: List[Dict]) -> Optional[Dict]:
        """Select the best genome assembly from available options.
        
        Preference order:
        1. RefSeq assemblies
        2. Chromosome-level assemblies
        3. Latest release date
        4. Highest contig N50
        
        Args:
            assemblies: List of assembly metadata dictionaries
            
        Returns:
            Best assembly dictionary or None
        """
        if not assemblies:
            return None
        
        # Score each assembly
        scored_assemblies = []
        for assembly in assemblies:
            assembly_info = assembly.get('assembly', {})
            assembly_stats = assembly.get('assembly_stats', {})
            
            score = 0
            
            # Prefer RefSeq (GCF) over GenBank (GCA)
            accession = assembly_info.get('assembly_accession', '')
            if accession.startswith('GCF_'):
                score += 1000
            
            # Prefer chromosome-level assemblies
            level = assembly_info.get('assembly_level', '').lower()
            if level == 'chromosome':
                score += 500
            elif level == 'scaffold':
                score += 100
            elif level == 'contig':
                score += 10
            
            # Prefer higher N50
            contig_n50 = assembly_stats.get('contig_n50', 0)
            if contig_n50:
                score += min(contig_n50 / 10000, 100)  # Cap at 100 points
            
            scored_assemblies.append((score, assembly))
        
        # Return highest scoring assembly
        scored_assemblies.sort(key=lambda x: x[0], reverse=True)
        return scored_assemblies[0][1]
    
    def generate_yaml_config(self, species_name: str, species_data: Dict, 
                            genome_info: Optional[Dict]) -> str:
        """Generate amalgkit YAML configuration for a species.
        
        Args:
            species_name: Scientific name
            species_data: RNA-seq metadata
            genome_info: Genome assembly metadata
            
        Returns:
            YAML configuration string
        """
        # Convert species name to safe filename format
        species_slug = species_name.replace(' ', '_').lower()
        genus, species = species_name.split(' ', 1) if ' ' in species_name else (species_name, 'sp')
        
        # Get taxonomy ID
        taxonomy_id = species_data.get('taxonomy_id', 'UNKNOWN')
        sample_count = species_data.get('sample_count', 0)
        
        # Build YAML content
        yaml_lines = [
            "# METAINFORMANT Amalgkit Configuration",
            f"# Species: {species_name}",
            f"# NCBI Taxonomy ID: {taxonomy_id}",
        ]
        
        if genome_info:
            yaml_lines.extend([
                f"# Assembly: {genome_info['accession']} ({genome_info['assembly_name']})",
                f"# Assembly Level: {genome_info['level']}",
                f"# Sequencing: {genome_info.get('sequencing_tech', 'Not specified')}",
            ])
            if genome_info.get('annotation_release'):
                yaml_lines.append(f"# Annotation Release: {genome_info['annotation_release']}")
        
        yaml_lines.extend([
            f"# RNA-seq Samples Available: {sample_count}",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d')}",
            "",
            "",
            f"work_dir: /home/q/Documents/GitHub/MetaInformAnt/output/amalgkit/{species_slug}/work",
            f"log_dir: /home/q/Documents/GitHub/MetaInformAnt/output/amalgkit/{species_slug}/logs",
            "threads: 12",
            "",
            "",
            "auto_install_amalgkit: true",
            "",
            "",
            "# Pre-selection filtering rules",
            "filters:",
            "  require_tissue: false  # Set to true to filter by tissue metadata",
            "",
            "",
            "# Species configuration",
            "species_list:",
            f"  - {species_name.replace(' ', '_')}",
            "",
            "",
            "# Alternative taxonomy specification",
            f"taxon_id: {taxonomy_id}",
            "",
            "",
        ])
        
        # Add genome configuration if available
        if genome_info:
            accession = genome_info['accession']
            assembly_name = genome_info['assembly_name']
            assembly_base = accession.replace('.', '_')
            
            # Construct FTP URL
            accession_parts = accession.split('_')
            if len(accession_parts) >= 2:
                gcf_gca = accession_parts[0]
                numbers = accession_parts[1].split('.')[0]
                # GCF_XXXXXXXXX.Y -> GCF/XXX/XXX/XXX/GCF_XXXXXXXXX.Y_AssemblyName
                num_groups = [numbers[i:i+3] for i in range(0, min(9, len(numbers)), 3)]
                ftp_path = '/'.join([gcf_gca] + num_groups + [f"{accession}_{assembly_name}"])
                ftp_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{ftp_path}/"
            else:
                ftp_url = "# FTP URL construction failed - please verify manually"
            
            yaml_lines.extend([
                "# Reference genome configuration",
                f"# Based on NCBI assembly: {accession}",
                "genome:",
                f"  accession: {accession}",
                f"  assembly_name: {assembly_name}",
            ])
            
            if genome_info.get('annotation_release'):
                yaml_lines.append(f"  annotation_release: {genome_info['annotation_release']}")
            
            yaml_lines.extend([
                f"  dest_dir: output/amalgkit/{species_slug}/genome",
                "",
                "",
                "  # Files to download (based on available NCBI FTP resources)",
                "  include:",
                "    - genome           # Genomic sequences (FASTA)",
                "    - gff3            # Gene annotations (GFF3)",
                "    - rna             # RNA sequences (FASTA)",
                "    - cds             # CDS sequences (FASTA)",
                "    - protein         # Protein sequences (FASTA)",
                "    - seq-report      # Sequence report",
                "    - feature-table   # Feature table",
                "    - gene-ontology   # Gene Ontology annotations (if available)",
                "",
                "",
                "  # Direct FTP URL for all genome files",
                f"  ftp_url: {ftp_url}",
                "",
                "",
                "  # Specific important files for transcriptome analysis",
                "  files:",
                f"    genomic_fasta: {accession}_{assembly_name}_genomic.fna.gz",
                f"    transcriptome_fasta: {accession}_{assembly_name}_rna_from_genomic.fna.gz",
                f"    cds_fasta: {accession}_{assembly_name}_cds_from_genomic.fna.gz",
                f"    protein_fasta: {accession}_{assembly_name}_protein.faa.gz",
                f"    annotation_gff: {accession}_{assembly_name}_genomic.gff.gz",
                f"    annotation_gtf: {accession}_{assembly_name}_genomic.gtf.gz",
                "",
                "",
            ])
        else:
            yaml_lines.extend([
                "# Reference genome configuration",
                "# NOTE: No genome assembly found in NCBI for this species",
                "# You may need to provide a custom genome or use a related species",
                "# genome:",
                "#   accession: GCF_XXXXXXXXX.X",
                "#   assembly_name: AssemblyName",
                "#   dest_dir: output/amalgkit/{species_slug}/genome",
                "",
                "",
            ])
        
        # Add step configuration
        yaml_lines.extend([
            "# Per-step parameters (merged on top of common params)",
            "steps:",
            "  metadata:",
            "    # Required by amalgkit metadata",
            f"    out_dir: output/amalgkit/{species_slug}/work",
            "    # Uses NCBI_EMAIL from environment; override here if needed",
            "    # entrez_email: you@example.com",
            f"    search_string: '\"{species_name}\"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'",
            "    redo: yes",
            "  integrate:",
            f"    fastq_dir: output/amalgkit/{species_slug}/fastq",
            "  config: {}",
            "  select: {}",
            "  getfastq:",
            f"    out_dir: output/amalgkit/{species_slug}/fastq",
            "    threads: 12",
            "    # Cloud acceleration via individual mirrors",
            "    aws: yes",
            "    gcp: yes",
            "    ncbi: yes",
            "    pfd: no           # DISABLED: parallel-fastq-dump usually not needed",
            "    fastp: yes        # ENABLED: fastp QC tool",
            "  quant:",
            "    # If a Salmon/Kallisto index is available, set it; otherwise genome_dir is injected from genome step",
            f"    out_dir: output/amalgkit/{species_slug}/quant",
            "    threads: 12",
            "    redo: no",
            "    # Delete SRA/FASTQ files after quantification to save space",
            "    keep_fastq: no",
            "    # Build kallisto index if not present",
            "    build_index: yes",
            "  merge:",
            f"    out: output/amalgkit/{species_slug}/merged/merged_abundance.tsv",
            f"    out_dir: output/amalgkit/{species_slug}/merged",
            "  cstmm:",
            f"    out_dir: output/amalgkit/{species_slug}/cstmm",
            "  curate:",
            f"    out_dir: output/amalgkit/{species_slug}/curate",
            "  csca:",
            f"    out_dir: output/amalgkit/{species_slug}/csca",
            "  sanity: {}",
            "",
        ])
        
        return '\n'.join(yaml_lines)
    
    def run(self):
        """Run the complete discovery and configuration generation process."""
        logger.info("=" * 80)
        logger.info("ANT SPECIES RNA-SEQ DISCOVERY AND CONFIGURATION GENERATION")
        logger.info("=" * 80)
        logger.info("")
        
        # Step 1: Search for ant species with RNA-seq data
        species_data = self.search_ants_with_rnaseq()
        
        if not species_data:
            logger.error("No ant species found with RNA-seq data")
            return
        
        # Save raw species data
        species_json_path = self.output_dir / "ant_species_rnaseq_data.json"
        with open(species_json_path, 'w') as f:
            json.dump(species_data, f, indent=2)
        logger.info(f"Saved species data to {species_json_path}")
        
        # Step 2: Get genome information for each species
        logger.info("")
        logger.info("=" * 80)
        logger.info("FETCHING GENOME INFORMATION")
        logger.info("=" * 80)
        logger.info("")
        
        species_with_genomes = {}
        species_without_genomes = {}
        
        for species_name, data in sorted(species_data.items(), key=lambda x: x[1]['sample_count'], reverse=True):
            taxonomy_id = data.get('taxonomy_id')
            if not taxonomy_id:
                logger.warning(f"No taxonomy ID for {species_name}, skipping genome lookup")
                species_without_genomes[species_name] = data
                continue
            
            genome_info = self.get_genome_info(taxonomy_id, species_name)
            
            if genome_info:
                data['genome'] = genome_info
                species_with_genomes[species_name] = data
                logger.info(f"✅ {species_name}: {genome_info['accession']} ({genome_info['level']} level)")
            else:
                species_without_genomes[species_name] = data
                logger.warning(f"⚠️  {species_name}: No genome assembly found")
        
        # Step 3: Generate YAML configurations
        logger.info("")
        logger.info("=" * 80)
        logger.info("GENERATING YAML CONFIGURATIONS")
        logger.info("=" * 80)
        logger.info("")
        
        config_dir = self.output_dir / "configs"
        config_dir.mkdir(exist_ok=True)
        
        generated_configs = []
        
        for species_name, data in sorted(species_data.items(), key=lambda x: x[1]['sample_count'], reverse=True):
            species_slug = species_name.replace(' ', '_').lower()
            genome_info = data.get('genome')
            
            yaml_config = self.generate_yaml_config(species_name, data, genome_info)
            
            config_path = config_dir / f"amalgkit_{species_slug}.yaml"
            with open(config_path, 'w') as f:
                f.write(yaml_config)
            
            generated_configs.append({
                'species': species_name,
                'config_file': str(config_path),
                'samples': data['sample_count'],
                'has_genome': genome_info is not None
            })
            
            logger.info(f"✅ Generated config for {species_name} ({data['sample_count']} samples)")
        
        # Step 4: Generate summary report
        logger.info("")
        logger.info("=" * 80)
        logger.info("GENERATING SUMMARY REPORT")
        logger.info("=" * 80)
        logger.info("")
        
        report_lines = [
            "# Ant Species RNA-seq Discovery Report",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Summary Statistics",
            f"- **Total ant species with RNA-seq data**: {len(species_data)}",
            f"- **Species with genome assemblies**: {len(species_with_genomes)}",
            f"- **Species without genome assemblies**: {len(species_without_genomes)}",
            f"- **Total RNA-seq samples**: {sum(d['sample_count'] for d in species_data.values())}",
            f"- **Minimum samples per species**: {self.min_samples}",
            "",
            "## Species with Genomes (Recommended for Analysis)",
            ""
        ]
        
        report_lines.append("| Species | Samples | Genome Accession | Assembly Level | Annotation |")
        report_lines.append("|---------|---------|------------------|----------------|------------|")
        
        for species_name in sorted(species_with_genomes.keys(), key=lambda x: species_with_genomes[x]['sample_count'], reverse=True):
            data = species_with_genomes[species_name]
            genome = data['genome']
            annotation = genome.get('annotation_release', 'None')
            report_lines.append(
                f"| {species_name} | {data['sample_count']} | {genome['accession']} | "
                f"{genome['level']} | {annotation} |"
            )
        
        if species_without_genomes:
            report_lines.extend([
                "",
                "## Species without Genome Assemblies",
                "",
                "These species have RNA-seq data but no public genome assembly in NCBI.",
                "You may need to use a closely related species' genome or wait for assembly publication.",
                "",
                "| Species | Samples | Taxonomy ID |",
                "|---------|---------|-------------|"
            ])
            
            for species_name in sorted(species_without_genomes.keys(), key=lambda x: species_without_genomes[x]['sample_count'], reverse=True):
                data = species_without_genomes[species_name]
                report_lines.append(f"| {species_name} | {data['sample_count']} | {data['taxonomy_id']} |")
        
        report_lines.extend([
            "",
            "## Generated Configuration Files",
            "",
            f"All YAML configuration files have been generated in: `{config_dir}/`",
            "",
            "To use a configuration:",
            "```bash",
            "# Copy to main config directory",
            "cp output/ant_discovery/configs/amalgkit_SPECIES.yaml config/amalgkit/",
            "",
            "# Run amalgkit workflow",
            "bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_SPECIES.yaml",
            "```",
            "",
            "## Next Steps",
            "",
            "1. **Review configurations**: Check generated YAML files for accuracy",
            "2. **Verify genome URLs**: Test FTP URLs for genome downloads",
            "3. **Prioritize species**: Focus on species with high sample counts and chromosome-level assemblies",
            "4. **Start workflows**: Begin with species that have complete genomes and many samples",
            "",
            "## Notes",
            "",
            "- Configurations use direct ENA download for maximum reliability",
            "- Quantification uses Kallisto with automatic index building",
            "- FASTQs are automatically deleted after quantification to save space",
            "- All paths are configured for the METAINFORMANT repository structure",
            ""
        ])
        
        report_path = self.output_dir / "DISCOVERY_REPORT.md"
        with open(report_path, 'w') as f:
            f.write('\n'.join(report_lines))
        
        logger.info(f"✅ Summary report saved to {report_path}")
        
        # Final summary
        logger.info("")
        logger.info("=" * 80)
        logger.info("DISCOVERY COMPLETE")
        logger.info("=" * 80)
        logger.info(f"✅ Discovered {len(species_data)} ant species with RNA-seq data")
        logger.info(f"✅ Generated {len(generated_configs)} YAML configurations")
        logger.info(f"✅ {len(species_with_genomes)} species have genome assemblies")
        logger.info(f"📁 Output directory: {self.output_dir}")
        logger.info("=" * 80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Discover ant species with RNA-seq data and generate amalgkit configurations",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output/ant_discovery',
        help='Output directory for configurations and reports (default: output/ant_discovery)'
    )
    
    parser.add_argument(
        '--min-samples',
        type=int,
        default=1,
        help='Minimum RNA-seq samples required per species (default: 1)'
    )
    
    args = parser.parse_args()
    
    # Run discovery
    discoverer = AntSpeciesDiscoverer(
        output_dir=Path(args.output_dir),
        min_samples=args.min_samples
    )
    
    try:
        discoverer.run()
    except KeyboardInterrupt:
        logger.info("\nDiscovery interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Discovery failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()




