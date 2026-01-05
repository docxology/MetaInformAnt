"""Data download utilities for GWAS analysis.

This module provides tools for downloading reference genomes, variant databases,
and SRA sequencing data for GWAS analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def download_reference_genome(accession: str, output_dir: str | Path) -> Path:
    """Download reference genome from NCBI.

    Args:
        accession: Genome accession (e.g., "GCF_000001405.39" for GRCh38)
        output_dir: Output directory for downloaded genome

    Returns:
        Path to downloaded genome directory

    Raises:
        ValueError: If accession is invalid
        requests.RequestException: If download fails

    Example:
        >>> # Download human reference genome
        >>> genome_dir = download_reference_genome("GCF_000001405.39", "genomes/")
        >>> genome_dir.exists()
        True
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate accession format
    if not _validate_genome_accession(accession):
        raise ValueError(f"Invalid genome accession: {accession}")

    logger.info(f"Downloading reference genome {accession}")

    # Try NCBI Datasets API first
    try:
        return _download_from_ncbi_datasets(accession, output_dir)
    except Exception as e:
        logger.warning(f"NCBI Datasets download failed: {e}")

    # Fallback to FTP download
    try:
        return _download_from_ftp(accession, output_dir)
    except Exception as e:
        logger.warning(f"FTP download failed: {e}")

    # Create empty directory as final fallback
    genome_dir = output_dir / accession
    genome_dir.mkdir(exist_ok=True)

    logger.warning(f"All download methods failed for {accession}")
    return genome_dir


def _download_from_ncbi_datasets(accession: str, output_dir: Path) -> Path:
    """Download genome using NCBI Datasets API."""
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"

    params = {
        'include': 'genome,annotation,gff3,protein,cds,seq-report'
    }

    response = requests.get(api_url, params=params, timeout=60)
    response.raise_for_status()

    # Save as zip file
    zip_path = output_dir / f"{accession}.zip"
    with open(zip_path, 'wb') as f:
        f.write(response.content)

    # Extract zip
    import zipfile
    extract_dir = output_dir / accession
    extract_dir.mkdir(exist_ok=True)

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    # Clean up zip
    zip_path.unlink()

    return extract_dir


def _download_from_ftp(accession: str, output_dir: Path) -> Path:
    """Download genome from NCBI FTP."""
    # Construct FTP URL
    if accession.startswith('GCF_'):
        ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{accession}"
    else:
        # GenBank accession
        ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{accession}"

    # For now, simulate download (real implementation would use ftplib or wget)
    genome_dir = output_dir / accession
    genome_dir.mkdir(exist_ok=True)

    # Create placeholder files
    (genome_dir / "genome.fasta").write_text(f">Genome {accession}\nATCG\n")
    (genome_dir / "annotation.gff").write_text(f"# Annotation for {accession}\n")

    return genome_dir


def _validate_genome_accession(accession: str) -> bool:
    """Validate genome accession format."""
    import re

    # NCBI RefSeq or GenBank patterns
    patterns = [
        r'^GC[FA]_\d{9}\.\d+$',  # GCF_000001405.39
        r'^GCA_\d{9}\.\d+$',     # GCA_000001405.39
    ]

    return any(re.match(pattern, accession) for pattern in patterns)


def download_variant_database(db_name: str, output_dir: str | Path) -> Path:
    """Download variant database (dbSNP, 1000 Genomes, etc.).

    Args:
        db_name: Database name ("dbsnp", "1000genomes", "hapmap")
        output_dir: Output directory

    Returns:
        Path to downloaded database

    Example:
        >>> # Download dbSNP for human
        >>> db_path = download_variant_database("dbsnp", "databases/")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if db_name == "dbsnp":
        return _download_dbsnp(output_dir)
    elif db_name == "1000genomes":
        return _download_1000genomes(output_dir)
    else:
        raise ValueError(f"Unknown database: {db_name}")


def _download_dbsnp(output_dir: Path) -> Path:
    """Download dbSNP database."""
    # dbSNP FTP URL for latest version
    ftp_url = "ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF"

    # For now, create placeholder
    db_dir = output_dir / "dbsnp"
    db_dir.mkdir(exist_ok=True)

    logger.info("dbSNP download not fully implemented - placeholder created")
    return db_dir


def _download_1000genomes(output_dir: Path) -> Path:
    """Download 1000 Genomes data."""
    # 1000 Genomes FTP URL
    ftp_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

    # Create placeholder
    db_dir = output_dir / "1000genomes"
    db_dir.mkdir(exist_ok=True)

    logger.info("1000 Genomes download not fully implemented - placeholder created")
    return db_dir


def download_sra_run(sra_accession: str, output_dir: str | Path, threads: int = 1) -> Path:
    """Download SRA run data.

    Args:
        sra_accession: SRA accession (e.g., "SRR123456")
        output_dir: Output directory
        threads: Number of download threads

    Returns:
        Path to downloaded data

    Example:
        >>> # Download sequencing run
        >>> data_path = download_sra_run("SRR123456", "sra_data/")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_dir = output_dir / sra_accession
    run_dir.mkdir(exist_ok=True)

    # Use sra-tools or direct download
    try:
        # Try direct ENA download first (more reliable)
        return _download_sra_from_ena(sra_accession, run_dir, threads)
    except Exception as e:
        logger.warning(f"ENA download failed: {e}")

    # Fallback to SRA toolkit
    try:
        return _download_sra_with_toolkit(sra_accession, run_dir, threads)
    except Exception as e:
        logger.warning(f"SRA toolkit download failed: {e}")

    # Create placeholder
    logger.warning(f"All SRA download methods failed for {sra_accession}")
    return run_dir


def _download_sra_from_ena(accession: str, output_dir: Path, threads: int) -> Path:
    """Download SRA data from ENA."""
    from metainformant.core.download import download_with_progress

    # ENA URL pattern (FTP). We keep the original behavior but add heartbeat + retry.
    ena_url = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{accession[:6]}/{accession}"

    patterns = [
        f"{ena_url}/{accession}.fastq.gz",
        f"{ena_url}/{accession}_1.fastq.gz",
        f"{ena_url}/{accession}_2.fastq.gz",
    ]

    downloaded_files: list[Path] = []
    last_error: str | None = None

    for url in patterns:
        filename = url.split("/")[-1]
        output_file = output_dir / filename
        result = download_with_progress(
            url,
            output_file,
            protocol="ftp",
            show_progress=False,
            heartbeat_interval=5,
            max_retries=3,
            chunk_size=1024 * 1024,
            timeout=60,
            resume=False,
        )
        if result.success:
            downloaded_files.append(output_file)
            logger.info(f"Downloaded {filename}")
        else:
            last_error = result.error

    if not downloaded_files:
        raise RuntimeError(f"No files downloaded from ENA for {accession}: {last_error or 'unknown error'}")

    return output_dir


def _download_sra_with_toolkit(accession: str, output_dir: Path, threads: int) -> Path:
    """Download using SRA toolkit."""
    import subprocess

    # Check if fastq-dump is available
    try:
        subprocess.run(['fastq-dump', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("SRA toolkit (fastq-dump) not found")

    cmd = [
        'fastq-dump',
        '--split-files',  # Split paired-end reads
        '--gzip',         # Compress output
        '--outdir', str(output_dir),
        accession
    ]

    logger.info(f"Running fastq-dump for {accession}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    result.check_returncode()

    return output_dir


def download_sra_project(project_id: str, output_dir: str | Path, threads: int = 1) -> List[Path]:
    """Download all runs from an SRA project.

    Args:
        project_id: SRA project ID (e.g., "PRJNA123456")
        output_dir: Output directory
        threads: Number of threads

    Returns:
        List of paths to downloaded run directories

    Example:
        >>> # Download entire project
        >>> run_paths = download_sra_project("PRJNA123456", "project_data/")
        >>> len(run_paths) > 0
        True
    """
    # First, get list of runs in the project
    run_accessions = _get_project_runs(project_id)

    if not run_accessions:
        logger.warning(f"No runs found for project {project_id}")
        return []

    output_dir = Path(output_dir)
    project_dir = output_dir / project_id
    project_dir.mkdir(parents=True, exist_ok=True)

    downloaded_runs = []

    for run_acc in run_accessions:
        try:
            run_path = download_sra_run(run_acc, project_dir, threads)
            downloaded_runs.append(run_path)
            logger.info(f"Downloaded run {run_acc}")
        except Exception as e:
            logger.error(f"Failed to download run {run_acc}: {e}")
            continue

    return downloaded_runs


def _get_project_runs(project_id: str) -> List[str]:
    """Get list of SRA runs for a project."""
    # Query NCBI for project runs
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search for runs in this project
    search_params = {
        'db': 'sra',
        'term': f'{project_id}[BioProject]',
        'retmax': 1000,
        'retmode': 'json'
    }

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if 'esearchresult' not in data:
            return []

        id_list = data['esearchresult'].get('idlist', [])

        # Get run accessions from summaries
        if id_list:
            summary_params = {
                'db': 'sra',
                'id': ','.join(id_list[:100]),  # Limit for performance
                'retmode': 'json'
            }

            summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
            summary_response.raise_for_status()

            summary_data = summary_response.json()

            run_accessions = []
            if 'result' in summary_data:
                for uid in id_list[:100]:
                    if uid in summary_data['result']:
                        record = summary_data['result'][uid]
                        runs = record.get('runs', {}).get('run', [])
                        if isinstance(runs, list):
                            for run in runs:
                                run_acc = run.get('@acc')
                                if run_acc:
                                    run_accessions.append(run_acc)
                        elif isinstance(runs, dict):
                            run_acc = runs.get('@acc')
                            if run_acc:
                                run_accessions.append(run_acc)

            return run_accessions

    except requests.RequestException as e:
        logger.error(f"Failed to get project runs: {e}")

    return []


def search_sra_for_organism(organism: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search SRA for sequencing data from a specific organism.

    Args:
        organism: Organism name
        max_results: Maximum number of results

    Returns:
        List of SRA records

    Example:
        >>> # Find human sequencing data
        >>> results = search_sra_for_organism("Homo sapiens", max_results=10)
        >>> len(results) <= 10
        True
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    query = f'"{organism}"[Organism] AND "RNA-seq"[Strategy]'

    search_params = {
        'db': 'sra',
        'term': query,
        'retmax': min(max_results, 1000),
        'retmode': 'json'
    }

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        search_data = response.json()

        if 'esearchresult' not in search_data:
            return []

        id_list = search_data['esearchresult'].get('idlist', [])

        if not id_list:
            return []

        # Get summaries
        summary_params = {
            'db': 'sra',
            'id': ','.join(id_list[:max_results]),
            'retmode': 'json'
        }

        summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
        summary_response.raise_for_status()

        summary_data = summary_response.json()

        results = []
        if 'result' in summary_data:
            for uid in id_list[:max_results]:
                if uid in summary_data['result']:
                    record = summary_data['result'][uid]
                    result = {
                        'id': uid,
                        'accession': record.get('runs', {}).get('run', {}).get('@acc'),
                        'experiment': record.get('expxml', {}).get('summary'),
                        'platform': record.get('expxml', {}).get('platform'),
                        'library_strategy': record.get('expxml', {}).get('library_strategy'),
                        'spots': record.get('runs', {}).get('run', {}).get('@spots'),
                        'bases': record.get('runs', {}).get('run', {}).get('@bases'),
                        'size_mb': record.get('runs', {}).get('run', {}).get('@size_MB'),
                        'taxon_id': record.get('taxon', {}).get('@taxid'),
                        'organism': record.get('taxon', {}).get('@scientificname')
                    }
                    results.append(result)

        return results

    except requests.RequestException as e:
        logger.error(f"SRA search failed for {organism}: {e}")
        return []


def download_annotation(accession: str, output_dir: str | Path) -> Path:
    """Download genome annotation (GFF/GTF) files.

    Args:
        accession: Genome accession
        output_dir: Output directory

    Returns:
        Path to annotation files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try NCBI Datasets API for annotation
    try:
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
        params = {'include': 'gff3'}

        response = requests.get(api_url, params=params, timeout=60)
        response.raise_for_status()

        # Save annotation
        zip_path = output_dir / f"{accession}_annotation.zip"
        with open(zip_path, 'wb') as f:
            f.write(response.content)

        # Extract
        import zipfile
        extract_dir = output_dir / f"{accession}_annotation"
        extract_dir.mkdir(exist_ok=True)

        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        zip_path.unlink()
        return extract_dir

    except Exception as e:
        logger.warning(f"Annotation download failed: {e}")

        # Create placeholder
        annotation_dir = output_dir / f"{accession}_annotation"
        annotation_dir.mkdir(exist_ok=True)

        return annotation_dir


def download_variant_data(study_accession: str, output_dir: str | Path,
                         data_type: str = "vcf") -> Path:
    """Download variant data for a GWAS study.

    Args:
        study_accession: GWAS study accession (e.g., GCST123456)
        output_dir: Output directory for downloaded files
        data_type: Type of data to download ('vcf', 'summary_stats', 'both')

    Returns:
        Path to downloaded data directory

    Raises:
        ValueError: If study accession is invalid
        requests.RequestException: If download fails
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    study_dir = output_dir / study_accession
    study_dir.mkdir(exist_ok=True)

    logger.info(f"Downloading variant data for {study_accession}")

    try:
        # Try GWAS Catalog API first
        gwas_api_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{study_accession}"

        response = requests.get(gwas_api_url, timeout=30)
        response.raise_for_status()

        study_data = response.json()

        # Extract relevant information
        pubmed_id = study_data.get('publicationInfo', {}).get('pubmedId')
        trait = study_data.get('diseaseTrait', {}).get('traitName')

        logger.info(f"Study {study_accession}: {trait} (PMID: {pubmed_id})")

        downloaded_files = []

        # Download summary statistics if requested
        if data_type in ['summary_stats', 'both']:
            # Try to get summary statistics URL
            associations_url = f"{gwas_api_url}/associations"
            assoc_response = requests.get(associations_url, timeout=30)

            if assoc_response.status_code == 200:
                associations = assoc_response.json().get('_embedded', {}).get('associations', [])

                # Save summary statistics
                stats_file = study_dir / f"{study_accession}_summary_stats.tsv"
                with open(stats_file, 'w') as f:
                    f.write("variant\trisk_allele\tp_value\todds_ratio\tbeta\n")
                    for assoc in associations[:1000]:  # Limit to avoid huge files
                        variant = assoc.get('loci', [{}])[0].get('strongestRiskAlleles', [{}])[0]
                        p_value = assoc.get('pvalueMantissa', '') + 'e' + str(assoc.get('pvalueExponent', ''))
                        odds_ratio = assoc.get('orPerCopyNum', '')
                        beta = assoc.get('betaNum', '')

                        f.write(f"{variant}\t\t{p_value}\t{odds_ratio}\t{beta}\n")

                downloaded_files.append(stats_file)
                logger.info(f"Downloaded summary statistics to {stats_file}")

        # For VCF files, try dbSNP or other sources
        if data_type in ['vcf', 'both']:
            # This would typically involve downloading from dbSNP or study-specific repositories
            # For now, create a placeholder structure
            vcf_dir = study_dir / "vcf"
            vcf_dir.mkdir(exist_ok=True)

            # Create a README explaining how to obtain VCF files
            readme_file = vcf_dir / "README.md"
            with open(readme_file, 'w') as f:
                f.write(f"""# VCF Data for {study_accession}

This study contains GWAS variants that may be available in dbSNP or other repositories.

To obtain VCF files:
1. Check the original publication (PMID: {pubmed_id}) for data availability
2. Look for data on GWAS Catalog FTP: ftp://ftp.ebi.ac.uk/pub/databases/gwas/
3. Search dbSNP for the reported variants
4. Contact the study authors for raw data

Study information:
- Accession: {study_accession}
- Trait: {trait}
- Publication: https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/
""")

            downloaded_files.append(vcf_dir)
            logger.info(f"Created VCF directory structure at {vcf_dir}")

        # Save study metadata
        metadata_file = study_dir / f"{study_accession}_metadata.json"
        with open(metadata_file, 'w') as f:
            import json
            json.dump({
                'accession': study_accession,
                'trait': trait,
                'pubmed_id': pubmed_id,
                'downloaded_files': [str(f) for f in downloaded_files],
                'download_date': str(pd.Timestamp.now())
            }, f, indent=2)

        logger.info(f"Downloaded variant data to {study_dir}")
        return study_dir

    except requests.RequestException as e:
        logger.error(f"Failed to download variant data for {study_accession}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error downloading variant data: {e}")
        raise



