#!/usr/bin/env python3
"""Stage 10: GO/HPO Ontology Analysis of GWAS Results.

Orchestrates a full functional annotation and enrichment pipeline for
GWAS hits, using the metainformant.ontology library modules:

    10a  Load GWAS summary statistics and gene annotations
    10b  Extract hit gene list and ranked gene list
    10c  Fetch GO annotations via QuickGO REST API (taxon-specific)
    10d  ORA  – over-representation analysis (hypergeometric, BH-FDR)
    10e  GSEA – gene set enrichment analysis (permutation NES)
    10f  Pathway similarity network (Jaccard)
    10g  HPO  – map GWAS phenotype label to HPO terms
    10h  Visualisation: enrichment dotplot + pathway network plot
    10i  Save all outputs (TSV / JSON / PNG)

Outputs:
    results/<phenotype>/<model>/ontology/
        ora_results.tsv
        gsea_results.tsv
        enrichment_dotplot.png
        pathway_network.png
        hpo_mapping.json
        ontology_summary.json

Usage::

    uv run scripts/10_ontology_gwas.py [--config CONFIG] [--phenotype PHENOTYPE]
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

# Set non-interactive matplotlib backend BEFORE any other matplotlib imports.
# Must happen before any module that might trigger matplotlib import.
import matplotlib

matplotlib.use("Agg")

# ── Path bootstrap ────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent.parent
for _candidate in [
    REPO_ROOT / "src",
    PROJECT_ROOT / "../../src",
]:
    _cand = _candidate.resolve()
    if _cand.exists() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

# ── Logging ───────────────────────────────────────────────────────────────────
try:
    from metainformant.core.utils.logging import configure_logging_from_env, get_logger

    configure_logging_from_env()
except ImportError:
    import logging

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    get_logger = logging.getLogger  # type: ignore

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_yaml(path: Path) -> dict:
    try:
        import yaml

        with open(path) as f:
            return yaml.safe_load(f) or {}
    except Exception as exc:
        logger.warning("Could not load YAML %s: %s", path, exc)
        return {}


def _load_tsv(path: Path) -> list[dict]:
    if not path.exists():
        logger.warning("TSV not found: %s", path)
        return []
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def _coerce_float(val: str | float, default: float = 1.0) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _save_tsv(rows: list[dict], path: Path, fields: list[str] | None = None) -> None:
    if not rows:
        logger.info("No rows to save: %s", path)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = fields or list(rows[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Saved %d rows -> %s", len(rows), path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_ontology_analysis(config: dict, phenotype: str, model: str, project_root: Path) -> None:
    """Run the full ontology analysis workflow."""
    config.get("gwas", {})

    ont_cfg = config.get("ontology", {})
    paths_cfg = config.get("paths", {})

    taxon_id: int = int(ont_cfg.get("taxon_id", 7460))
    top_n: int = int(ont_cfg.get("top_n_hits", 20))
    go_aspects: list[str] = ont_cfg.get(
        "go_aspects", ["biological_process", "molecular_function", "cellular_component"]
    )
    n_permutations: int = int(ont_cfg.get("gsea_permutations", 500))
    min_set_size: int = int(ont_cfg.get("min_go_set_size", 3))
    max_set_size: int = int(ont_cfg.get("max_go_set_size", 500))
    hpo_phenotype: str = ont_cfg.get("hpo_phenotype", phenotype)
    fdr_threshold: float = float(ont_cfg.get("fdr_threshold", 0.05))
    jaccard_threshold: float = float(ont_cfg.get("jaccard_threshold", 0.3))
    rate_limit_s: float = float(ont_cfg.get("rate_limit_s", 0.2))

    # Paths — read from config, same as all other scripts
    paths_cfg = config.get("paths", {})
    results_dir = project_root / paths_cfg.get("results_dir", "output")
    results_base = results_dir / phenotype / model
    output_dir = results_base / "ontology"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_stats_path = results_base / "summary_statistics.tsv"
    post_gwas_json = results_base / "post_gwas" / "post_gwas_results.json"

    # ── Step 10a: Load data ────────────────────────────────────────────────
    logger.info("Step 10a: Loading GWAS summary statistics and gene annotations")

    assoc_results_raw = _load_tsv(summary_stats_path)
    if not assoc_results_raw:
        logger.error("No summary statistics found at %s. Run stage 6 first.", summary_stats_path)
        sys.exit(1)

    # Parse numeric fields — handle both lowercase and uppercase column headers
    association_results: list[dict] = []
    for row in assoc_results_raw:
        try:
            # Accept both "chrom"/"CHR" etc.
            chrom = row.get("chrom") or row.get("CHR") or row.get("CHROM") or ""
            pos_raw = row.get("pos") or row.get("POS") or "0"
            snp = row.get("snp") or row.get("SNP") or row.get("variant_id") or row.get("ID") or ""
            p_raw = row.get("p_value") or row.get("P") or row.get("PVALUE") or "1"
            beta_raw = row.get("beta") or row.get("BETA") or "0"
            se_raw = row.get("se") or row.get("SE") or "1"
            maf_raw = row.get("maf") or row.get("MAF") or "0"
            association_results.append(
                {
                    "snp": snp,
                    "chrom": str(chrom),
                    "pos": int(float(pos_raw or 0)),
                    "p_value": _coerce_float(p_raw, 1.0),
                    "beta": _coerce_float(beta_raw, 0.0),
                    "se": _coerce_float(se_raw, 1.0),
                    "maf": _coerce_float(maf_raw, 0.0),
                }
            )
        except Exception:
            continue

    logger.info("  Loaded %d variants from summary statistics", len(association_results))

    # Load gene annotations from post_gwas_results.json
    gene_annotations: list[dict] = []
    if post_gwas_json.exists():
        try:
            with open(post_gwas_json) as f:
                pg_data = json.load(f)
            gene_annotations = pg_data.get("gene_annotations", [])
            logger.info("  Loaded %d gene annotations from step 8g", len(gene_annotations))
        except Exception as exc:
            logger.warning("  Could not read post_gwas_results.json: %s", exc)

    # ── Step 10b: Extract genes and ranked list ────────────────────────────
    logger.info("Step 10b: Extracting hit genes and ranked gene list")
    from metainformant.ontology.annotation.annotate import (
        build_background_from_vcf_genes,
        gwas_hits_to_genes,
        rank_genes_by_pvalue,
    )

    hit_genes = gwas_hits_to_genes(
        association_results,
        top_n=top_n,
        gene_annotations=gene_annotations,
    )
    ranked_genes = rank_genes_by_pvalue(
        association_results,
        gene_annotations=gene_annotations,
    )
    background_genes = build_background_from_vcf_genes(gene_annotations)

    logger.info(
        "  Hit genes: %d | Ranked list: %d | Background: %d", len(hit_genes), len(ranked_genes), len(background_genes)
    )

    if not hit_genes:
        logger.info(
            "  No gene annotations found — gene annotation step (8g) must be run with "
            "real NCBI/Ensembl IDs. Skipping ontology analysis for this dataset."
        )
        # Write empty summary so the pipeline does not hard-fail
        (output_dir / "ontology_summary.json").write_text(
            json.dumps({"status": "skipped", "reason": "no_gene_annotations"}, indent=2)
        )
        sys.exit(0)

    # ── Step 10c: Fetch GO gene sets (taxon-wide QuickGO reverse lookup) ──
    # Strategy: instead of per-gene-product queries (requires UniProtKB IDs),
    # we pull all GO annotations for the taxon (paginated) to build the full
    # GO-term → gene-product mapping.  Hit gene NCBI symbols are then mapped
    # to UniProt IDs and intersected with this mapping for ORA.
    logger.info(
        "Step 10c: Building taxon-wide GO gene-set reference via QuickGO " "(taxon_id=%d, aspects=%s)",
        taxon_id,
        go_aspects,
    )
    from metainformant.ontology.core.go_api import (
        build_taxon_go_gene_sets,
        map_symbols_to_uniprot,
    )

    gene_sets: dict[str, set[str]] = {}
    uniprot_hit_genes: list[str] = []
    try:
        # 1. Map NCBI gene symbols → UniProtKB accessions for hit genes
        symbol_to_uniprot = map_symbols_to_uniprot(hit_genes, taxon_id=taxon_id, rate_limit_s=rate_limit_s)
        uniprot_hit_genes = []
        for sym in hit_genes:
            mapped = symbol_to_uniprot.get(sym, [])
            if mapped:
                uniprot_hit_genes.extend(mapped)
            else:
                uniprot_hit_genes.append(sym)
        logger.info(
            "  UniProt mapping: %d/%d hit genes mapped (%d UniProt accessions)",
            len(symbol_to_uniprot),
            len(hit_genes),
            len(uniprot_hit_genes),
        )

        # 2. Fetch all-taxon GO gene sets (paged, up to max_pages)
        max_go_pages = int(ont_cfg.get("max_go_pages", 25))
        gene_sets = build_taxon_go_gene_sets(
            taxon_id=taxon_id,
            go_aspects=go_aspects,
            max_pages=max_go_pages,
            rate_limit_s=rate_limit_s,
        )
        logger.info(
            "  Taxon GO gene sets: %d GO terms, %d pages fetched",
            len(gene_sets),
            max_go_pages,
        )

        # 3. Explicitly fetch and merge annotations for our hit genes to ensure they aren't missed by pagination
        from metainformant.ontology.core.go_api import fetch_gene_go_annotations

        for uni_id in uniprot_hit_genes:
            try:
                annots = fetch_gene_go_annotations(
                    uni_id, taxon_id=taxon_id, go_aspects=go_aspects, rate_limit_s=rate_limit_s
                )
                for a in annots:
                    goid = a.get("go_id")
                    if goid:
                        if goid not in gene_sets:
                            gene_sets[goid] = set()
                        gene_sets[goid].add(uni_id)
            except Exception as e:
                logger.debug("Failed to fetch explicit annotations for %s: %s", uni_id, e)
        logger.info("  Merged explicit gene-product annotations. Total GO terms now: %d", len(gene_sets))

    except Exception as exc:
        logger.error("  GO gene-set fetch failed with unexpected error: %s", exc, exc_info=True)

    # ── Step 10d: ORA ─────────────────────────────────────────────────────
    # Use UniProt-mapped hit genes against the taxon-wide GO gene sets.
    # If UniProt mapping succeeded, query_list = uniprot_hit_genes;
    # otherwise fall back to raw NCBI symbols (which may not match QuickGO
    # gene product IDs — ORA will run but yield 0 overlaps, which is logged).
    logger.info("Step 10d: Over-Representation Analysis (ORA)")
    ora_results: list[dict] = []
    active_gene_sets: dict = {}
    try:
        from metainformant.ontology.pathway_enrichment.enrichment import over_representation_analysis

        # Use UniProt-mapped IDs for hit genes if available; otherwise raw symbols
        query_genes = uniprot_hit_genes if uniprot_hit_genes else hit_genes

        # Filter GO gene sets by size
        filtered_gs = {k: v for k, v in gene_sets.items() if min_set_size <= len(v) <= max_set_size}

        if not gene_sets:
            logger.info(
                "  ORA skipped: no GO gene sets available (QuickGO fetch returned 0 terms "
                "for taxon %d — check connectivity or run with real gene annotations).",
                taxon_id,
            )
        elif not filtered_gs:
            logger.info(
                "  ORA skipped: 0 GO terms in size range [%d, %d] (got %d total terms "
                "— increase max_go_pages in config to fetch more).",
                min_set_size,
                max_set_size,
                len(gene_sets),
            )
        else:
            active_gene_sets = {k: list(v) for k, v in filtered_gs.items()}
            ora_results = over_representation_analysis(
                gene_list=query_genes,
                gene_sets=filtered_gs,
                background=list(gene_sets.get("GO:0008150", set()) or set(query_genes)) or None,
                correction="fdr_bh",
            )
            n_sig = sum(1 for r in ora_results if r.get("adjusted_p", 1.0) < fdr_threshold)
            logger.info(
                "  ORA: %d terms tested, %d significant (FDR<%.2f), "
                "query genes=%d (UniProt), background from taxon gene sets",
                len(ora_results),
                n_sig,
                fdr_threshold,
                len(query_genes),
            )

    except Exception as exc:
        logger.error("  ORA failed: %s", exc, exc_info=True)

    # ── Step 10e: GSEA ────────────────────────────────────────────────────
    logger.info("Step 10e: Gene Set Enrichment Analysis (GSEA, n_perm=%d)", n_permutations)
    gsea_results: list[dict] = []
    try:
        from metainformant.ontology.pathway_enrichment.enrichment import gsea

        if ranked_genes and gene_sets:
            gsea_input = [(g, float(w)) for g, w in ranked_genes]
            gsea_results = gsea(
                ranked_genes=gsea_input,
                gene_sets={k: set(v) for k, v in gene_sets.items()},
                n_permutations=min(n_permutations, 200),  # cap for speed
                min_size=min_set_size,
                max_size=max_set_size,
            )
            n_sig_gsea = sum(1 for r in gsea_results if r.get("fdr", 1.0) < 0.25)
            logger.info("  GSEA: %d gene sets tested, %d significant (FDR<0.25)", len(gsea_results), n_sig_gsea)
        else:
            logger.info("  GSEA: skipped (no ranked genes or gene sets)")
    except Exception as exc:
        logger.warning("  GSEA failed: %s", exc)

    # ── Step 10f: Pathway network ─────────────────────────────────────────
    logger.info("Step 10f: Pathway similarity network (Jaccard threshold=%.2f)", jaccard_threshold)
    network: dict = {"nodes": [], "edges": [], "clusters": []}
    try:
        from metainformant.ontology.pathway_enrichment.enrichment import pathway_network

        if ora_results and active_gene_sets:
            network = pathway_network(
                enrichment_results=ora_results,
                gene_sets=active_gene_sets,
                similarity_threshold=jaccard_threshold,
            )
            logger.info(
                "  Pathway network: %d nodes, %d edges, %d clusters",
                len(network["nodes"]),
                len(network["edges"]),
                len(network["clusters"]),
            )
        else:
            logger.info("  Pathway network: no ORA results to build network from")
    except Exception as exc:
        logger.warning("  Pathway network failed: %s", exc)

    # ── Step 10g: HPO mapping ─────────────────────────────────────────────
    logger.info("Step 10g: HPO trait mapping for '%s'", hpo_phenotype)
    hpo_mapping: dict = {"phenotype": hpo_phenotype, "hp_ids": [], "terms": []}
    try:
        from metainformant.ontology.core.hpo import describe_hpo_terms, map_phenotype_to_hpo

        hp_ids = map_phenotype_to_hpo(hpo_phenotype, rate_limit_s=rate_limit_s)
        hpo_mapping["hp_ids"] = hp_ids

        if hp_ids:
            terms = describe_hpo_terms(hp_ids, rate_limit_s=rate_limit_s)
            hpo_mapping["terms"] = terms
            for t in terms[:3]:
                logger.info("  HPO: %s — %s", t.get("id", "?"), t.get("name", "?"))
        else:
            logger.info("  HPO: no terms found for '%s' (agricultural trait, limited coverage)", hpo_phenotype)
    except Exception as exc:
        logger.warning("  HPO mapping failed: %s", exc)

    # ── Step 10h: Visualisations ──────────────────────────────────────────
    logger.info("Step 10h: Generating enrichment dotplot and pathway network plot")
    try:
        import matplotlib.pyplot as plt

        from metainformant.ontology.visualization.plots import enrichment_dotplot, pathway_network_plot

        # Enrichment dotplot
        dotplot_path = output_dir / "enrichment_dotplot.png"
        if ora_results:
            # Ensure we have term_name populated
            for r in ora_results:
                if "term_name" not in r or r["term_name"] == r.get("term_id", ""):
                    r["term_name"] = r.get("term_id", "?")
            fig_dot = enrichment_dotplot(
                ora_results,
                output_path=dotplot_path,
                top_n=20,
                title=f"GO Enrichment - {phenotype.replace('_', ' ').title()}",
            )
            if fig_dot is not None:
                plt.close(fig_dot)
            logger.info("  Enrichment dotplot saved to %s", dotplot_path)
        else:
            logger.info("  Dotplot skipped: no ORA results")

        # Pathway network plot
        net_path = output_dir / "pathway_network.png"
        if network["nodes"]:
            fig_net = pathway_network_plot(
                network,
                output_path=net_path,
                top_n=30,
                title=f"Pathway Network - {phenotype.replace('_', ' ').title()}",
            )
            if fig_net is not None:
                plt.close(fig_net)
            logger.info("  Pathway network plot saved to %s", net_path)
        else:
            logger.info("  Network plot skipped: no nodes")

    except Exception as exc:
        logger.warning("  Visualisation failed: %s", exc)

    # ── Step 10i: Save all outputs ────────────────────────────────────────
    logger.info("Step 10i: Saving all outputs")

    # ORA TSV
    if ora_results:
        ora_fields = [
            "term_id",
            "term_name",
            "p_value",
            "adjusted_p",
            "odds_ratio",
            "n_genes",
            "n_overlap",
            "overlap_genes",
        ]
        ora_rows = []
        for r in ora_results:
            row = {k: r.get(k, "") for k in ora_fields}
            if isinstance(row.get("overlap_genes"), list):
                row["overlap_genes"] = ",".join(row["overlap_genes"])
            ora_rows.append(row)
        _save_tsv(ora_rows, output_dir / "ora_results.tsv", fields=ora_fields)

    # GSEA TSV
    if gsea_results:
        gsea_fields = ["term_id", "term_name", "es", "nes", "p_value", "fdr", "leading_edge", "leading_edge_genes"]
        gsea_rows = []
        for r in gsea_results:
            row = {k: r.get(k, "") for k in gsea_fields}
            if isinstance(row.get("leading_edge_genes"), list):
                row["leading_edge_genes"] = ",".join(str(g) for g in row["leading_edge_genes"])
            gsea_rows.append(row)
        _save_tsv(gsea_rows, output_dir / "gsea_results.tsv", fields=gsea_fields)

    # HPO JSON
    with open(output_dir / "hpo_mapping.json", "w") as f:
        json.dump(hpo_mapping, f, indent=2)
    logger.info("  HPO mapping saved")

    # Ontology summary JSON
    summary: dict = {
        "stage": "10_ontology_gwas",
        "phenotype": phenotype,
        "model": model,
        "taxon_id": taxon_id,
        "n_hit_genes": len(hit_genes),
        "n_ranked_genes": len(ranked_genes),
        "n_background_genes": len(background_genes),
        "n_go_terms_fetched": len(gene_sets),
        "ora": {
            "n_tested": len(ora_results),
            "n_significant_fdr05": sum(1 for r in ora_results if r.get("adjusted_p", 1.0) < fdr_threshold),
        },
        "gsea": {
            "n_tested": len(gsea_results),
            "n_significant_fdr025": sum(1 for r in gsea_results if r.get("fdr", 1.0) < 0.25),
        },
        "pathway_network": {
            "n_nodes": len(network.get("nodes", [])),
            "n_edges": len(network.get("edges", [])),
            "n_clusters": len(network.get("clusters", [])),
        },
        "hpo": {
            "phenotype": hpo_phenotype,
            "n_hp_ids": len(hpo_mapping.get("hp_ids", [])),
            "hp_ids": hpo_mapping.get("hp_ids", []),
        },
        "outputs": [str(p.relative_to(project_root)) for p in output_dir.iterdir() if p.is_file()],
    }

    with open(output_dir / "ontology_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print("=" * 60)
    print("  Stage 10: ONTOLOGY ANALYSIS COMPLETE")
    print(f"  Phenotype : {phenotype} / {model}")
    print(f"  Hit genes : {summary['n_hit_genes']}")
    print(f"  GO terms  : {summary['n_go_terms_fetched']} fetched")
    print(f"  ORA       : {summary['ora']['n_tested']} tested, {summary['ora']['n_significant_fdr05']} sig.")
    print(f"  GSEA      : {summary['gsea']['n_tested']} tested, {summary['gsea']['n_significant_fdr025']} sig.")
    print(f"  HPO terms : {summary['hpo']['n_hp_ids']}")
    print(f"  Outputs   -> {output_dir}")
    print("=" * 60)
    print()
    print(f"SUCCESS: Ontology results in {output_dir.relative_to(project_root)}/")
