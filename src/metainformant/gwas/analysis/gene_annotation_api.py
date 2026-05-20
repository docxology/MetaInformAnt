"""Gene region annotation using live NCBI E-utilities for Apis mellifera.

Replaces the Ensembl-based lookup which does not serve Apis mellifera
(a non-vertebrate; hosted in Ensembl Metazoa but no public REST API for region
queries at the time of writing). Uses NCBI E-utilities to:

1. Convert chromosome accession + position to gene IDs via ``elink/esearch``.
2. Fetch gene details (symbol, description, genomic location) via ``esummary``.

Works for any organism in the NCBI Gene database.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_TIMEOUT = 15  # seconds per request


# ---------------------------------------------------------------------------
# Low-level NCBI helpers
# ---------------------------------------------------------------------------


def _ncbi_get(url: str) -> dict | list | None:
    """HTTP GET with JSON response; returns None on failure."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "metainformant/1.0"})
        with urllib.request.urlopen(req, timeout=NCBI_TIMEOUT) as resp:
            return json.load(resp)
    except urllib.error.HTTPError as exc:
        logger.debug("NCBI HTTP %d for %s", exc.code, url)
        return None
    except Exception as exc:
        logger.debug("NCBI request failed (%s): %s", type(exc).__name__, exc)
        return None


def lookup_genes_by_region_ncbi(
    chrom_accession: str,
    start: int,
    end: int,
    taxon_id: int = 7460,
) -> list[dict[str, Any]]:
    """Query NCBI Gene for genes overlapping a genomic region.

    Uses NCBI E-Utilities ``esearch`` with chromosome accession-based query
    and gene position filter.  Results are returned as minimal annotation dicts.

    Args:
        chrom_accession: Chromosome accession (e.g. ``CM009934.2``).
        start: Region start position (1-based, inclusive).
        end: Region end position (1-based, inclusive).
        taxon_id: NCBI taxon ID (default ``7460`` — *Apis mellifera*).

    Returns:
        List of gene annotation dicts::

            {
                gene_id, gene_name, description, chromosome,
                chr_start, chr_end, strand, biotype
            }

        Returns ``[]`` when no genes are found or the API is unavailable.

    Examples:
        >>> genes = lookup_genes_by_region_ncbi('CM009934.2', 4600000, 4900000)
        >>> isinstance(genes, list)
        True
    """
    # Apis mellifera (Amel_HAv3.1) mapping from GenBank CM* to RefSeq NC*
    if chrom_accession.startswith("CM0099") and chrom_accession.endswith(".2"):
        try:
            num = int(chrom_accession[6:8])
            if 31 <= num <= 46:
                # CM009931.2 is LG1 -> NC_037638.1
                chrom_accession = f"NC_0376{num + 7}.1"
        except ValueError:
            pass

    # Build query: taxon + chromosome accession + position range
    # NCBI gene does not support WGS accession + position natively through
    # a simple term, so we query by accession and organism, then filter.
    term = f"{taxon_id}[taxid] AND {chrom_accession}[accession]"
    search_url = f"{NCBI_EUTILS_BASE}/esearch.fcgi" f"?db=gene&term={urllib.parse.quote(term)}&retmax=100&retmode=json"
    result = _ncbi_get(search_url)
    if not result:
        logger.debug("NCBI esearch returned nothing for %s", term)
        return []

    gene_ids = result.get("esearchresult", {}).get("idlist", [])
    if not gene_ids:
        return []

    # Fetch summaries in one batch
    ids_joined = ",".join(gene_ids[:50])
    summary_url = f"{NCBI_EUTILS_BASE}/esummary.fcgi?db=gene&id={ids_joined}&retmode=json"
    summary = _ncbi_get(summary_url)
    if not summary:
        return []

    result_map = summary.get("result", {})
    genes: list[dict[str, Any]] = []
    for gid in gene_ids:
        g = result_map.get(gid)
        if not g or g.get("status") == "discontinued":
            continue
        # Filter by genomic range
        loc_info = g.get("genomicinfo", [{}])
        in_region = False
        for loc in loc_info:
            gs = int(loc.get("chrstart", 0))
            ge = int(loc.get("chrstop", 0))
            # Overlap check (NCBI uses 0-based half-open intervals)
            if gs <= end and ge >= start:
                in_region = True
                break
        if not in_region and loc_info and loc_info[0].get("chrstart") is not None:
            continue  # has location data but outside range

        chrinfo = loc_info[0] if loc_info else {}
        gene_mid = (int(chrinfo.get("chrstart", 0)) + int(chrinfo.get("chrstop", 0))) / 2

        genes.append(
            {
                "gene_id": gid,
                "gene_name": g.get("symbol", g.get("name", gid)),
                "description": g.get("description", ""),
                "chromosome": g.get("chromosome", ""),
                "chr_start": int(chrinfo.get("chrstart", 0)),
                "chr_end": int(chrinfo.get("chrstop", 0)),
                "strand": chrinfo.get("exongapcount", 0),
                "biotype": g.get("type", ""),
                "distance_bp": 0,  # will be set by caller
                "_gene_mid": gene_mid,
            }
        )

    logger.debug(
        "NCBI gene region %s:%d-%d → %d genes (taxon %d)",
        chrom_accession,
        start,
        end,
        len(genes),
        taxon_id,
    )
    return genes


# ---------------------------------------------------------------------------
# Top-level annotation function
# ---------------------------------------------------------------------------


def annotate_top_hits_ncbi(
    association_results: list[dict[str, Any]],
    *,
    radius_bp: int = 100_000,
    taxon_id: int = 7460,
    top_n: int = 20,
    rate_limit_s: float = 0.34,  # NCBI allows ~3 req/s without API key
) -> list[dict[str, Any]]:
    """Annotate top GWAS hits with nearby genes via NCBI Gene database.

    For each of the top ``top_n`` hits (by p-value), queries NCBI Gene for
    genes within ±``radius_bp`` of the variant position.

    Args:
        association_results: GWAS results with ``chrom``, ``pos``,
            ``p_value``, ``snp`` keys.
        radius_bp: Search radius around each hit (default 100 kbp).
        taxon_id: NCBI taxon ID.
        top_n: Number of top hits to annotate.
        rate_limit_s: Sleep between API calls.

    Returns:
        List of annotation dicts, one per hit::

            {
                snp, chrom, pos, p_value, beta, se,
                nearby_genes: [{gene_id, gene_name, description, ...}],
                nearest_gene: str | None,
                nearest_gene_dist_bp: int | None,
            }

    Examples:
        >>> annotations = annotate_top_hits_ncbi(results, taxon_id=7460)
        >>> isinstance(annotations, list)
        True
    """
    sorted_results = sorted(association_results, key=lambda r: r.get("p_value", 1.0))[:top_n]

    annotated: list[dict[str, Any]] = []
    n_genes_total = 0

    for r in sorted_results:
        chrom = r.get("chrom", "")
        pos = int(r.get("pos", 0))
        snp = r.get("snp", "")
        region_start = max(1, pos - radius_bp)
        region_end = pos + radius_bp

        genes = lookup_genes_by_region_ncbi(
            chrom_accession=chrom,
            start=region_start,
            end=region_end,
            taxon_id=taxon_id,
        )

        # Compute distance from SNP position to gene midpoint
        for g in genes:
            mid = g.pop("_gene_mid", (g["chr_start"] + g["chr_end"]) / 2)
            g["distance_bp"] = int(abs(pos - mid))
            g["overlaps"] = g["chr_start"] <= pos <= g["chr_end"]

        genes.sort(key=lambda g: g["distance_bp"])
        n_genes_total += len(genes)

        annotated.append(
            {
                "snp": snp,
                "chrom": chrom,
                "pos": pos,
                "p_value": r.get("p_value"),
                "beta": r.get("beta"),
                "se": r.get("se"),
                "nearby_genes": genes,
                "nearest_gene": genes[0]["gene_name"] if genes else None,
                "nearest_gene_dist_bp": genes[0]["distance_bp"] if genes else None,
            }
        )

        time.sleep(rate_limit_s)

    logger.info(
        "NCBI annotation: %d hits queried, %d genes found total (taxon %d)",
        len(annotated),
        n_genes_total,
        taxon_id,
    )
    return annotated


# ---------------------------------------------------------------------------
# Backwards-compatible wrapper kept for existing callers of Ensembl function
# ---------------------------------------------------------------------------


def annotate_top_hits_ensembl(
    association_results: list[dict[str, Any]],
    radius_bp: int = 100_000,
    species: str = "apis_mellifera",
    top_n: int = 10,
    rate_limit_s: float = 0.15,
) -> list[dict[str, Any]]:
    """Backwards-compatible wrapper — delegates to :func:`annotate_top_hits_ncbi`.

    The original Ensembl REST API does not serve non-vertebrate genomes
    (*Apis mellifera* is in Ensembl Metazoa, which has no public region-query
    REST endpoint).  This wrapper falls through to the NCBI implementation.

    Args:
        association_results: GWAS result dicts.
        radius_bp: Search radius in bp.
        species: Ensembl species name (used only to infer taxon_id).
        top_n: Number of top hits to annotate.
        rate_limit_s: Passed to NCBI annotation; bumped to 0.34 minimum.

    Returns:
        Same format as :func:`annotate_top_hits_ncbi`.
    """
    # Derive taxon_id from species name
    _TAXON_MAP = {
        "apis_mellifera": 7460,
        "homo_sapiens": 9606,
        "mus_musculus": 10090,
        "drosophila_melanogaster": 7227,
        "danio_rerio": 7955,
        "caenorhabditis_elegans": 6239,
    }
    taxon_id = _TAXON_MAP.get(species, 7460)
    logger.info(
        "annotate_top_hits_ensembl: delegating to NCBI (taxon=%d, species=%s)",
        taxon_id,
        species,
    )
    return annotate_top_hits_ncbi(
        association_results,
        radius_bp=radius_bp,
        taxon_id=taxon_id,
        top_n=top_n,
        rate_limit_s=max(rate_limit_s, 0.34),
    )
