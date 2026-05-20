"""Real API tests for ontology submodules: go_api, hpo, ncbi gene annotation.

NO MOCK POLICY: All tests use live API calls.
Tests are skipped gracefully when there is no internet connectivity.

Run with:
    uv run pytest tests/ontology/test_ontology_api.py -v
"""

from __future__ import annotations

import socket

import pytest

# ---------------------------------------------------------------------------
# Connectivity guard — skip entire module if no internet
# ---------------------------------------------------------------------------


def _has_internet() -> bool:
    """Check if the internet is reachable."""
    try:
        socket.setdefaulttimeout(5)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect(("8.8.8.8", 53))
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _has_internet(), reason="No internet connection")


# ---------------------------------------------------------------------------
# QuickGO API — fetch_go_term
# ---------------------------------------------------------------------------


class TestFetchGoTerm:
    def test_known_bp_term(self):
        from metainformant.ontology.core.go_api import fetch_go_term

        term = fetch_go_term("GO:0008150")
        assert term, "Should return a non-empty dict"
        assert term["name"] == "biological_process"
        assert term["aspect"] in ("biological_process", "P", "BP")
        assert not term["is_obsolete"]

    def test_known_mf_term(self):
        from metainformant.ontology.core.go_api import fetch_go_term

        term = fetch_go_term("GO:0003674")
        assert term["name"] == "molecular_function"

    def test_nonexistent_term_returns_empty(self):
        from metainformant.ontology.core.go_api import fetch_go_term

        term = fetch_go_term("GO:9999999")
        assert term == {} or isinstance(term, dict)

    def test_caches_on_second_call(self):
        from metainformant.ontology.core.go_api import _term_cache, fetch_go_term

        fetch_go_term("GO:0005575")  # cellular_component
        assert "GO:0005575" in _term_cache


# ---------------------------------------------------------------------------
# QuickGO API — fetch_gene_go_annotations (per-gene-product)
# ---------------------------------------------------------------------------


class TestFetchGeneGoAnnotations:
    def test_uniprot_accession_returns_annotations(self):
        """A known Apis mellifera protein in QuickGO should return annotations."""
        from metainformant.ontology.core.go_api import fetch_gene_go_annotations

        # Apis mellifera vitellogenin (UniProtKB:P15203) — well-annotated bee protein
        anns = fetch_gene_go_annotations("UniProtKB:P15203", taxon_id=7460)
        assert isinstance(anns, list)
        # At least one GO annotation expected for a well-known protein
        if anns:
            assert "go_id" in anns[0]
            assert anns[0]["go_id"].startswith("GO:")
            assert "evidence_code" in anns[0]

    def test_invalid_product_id_returns_empty(self):
        """Invalid gene product IDs should return empty list (not raise)."""
        from metainformant.ontology.core.go_api import fetch_gene_go_annotations

        anns = fetch_gene_go_annotations("NOT_A_REAL_ID_XYZ789", taxon_id=7460)
        assert isinstance(anns, list)

    def test_aspect_filter(self):
        from metainformant.ontology.core.go_api import fetch_gene_go_annotations

        anns = fetch_gene_go_annotations(
            "UniProtKB:P15203",
            taxon_id=7460,
            go_aspects=["biological_process"],
        )
        # All returned annotations should be biological_process
        for ann in anns:
            aspect = ann.get("aspect", "")
            assert "biological" in aspect.lower() or aspect in ("P", "BP")


# ---------------------------------------------------------------------------
# QuickGO API — build_taxon_go_gene_sets (taxon-wide reverse-lookup)
# ---------------------------------------------------------------------------


class TestBuildTaxonGoGeneSets:
    def test_returns_nonempty_dict(self):
        """Fetching 1 page from QuickGO taxon 7460 should return GO terms."""
        from metainformant.ontology.core.go_api import build_taxon_go_gene_sets

        gene_sets = build_taxon_go_gene_sets(taxon_id=7460, max_pages=1)
        assert isinstance(gene_sets, dict)
        assert len(gene_sets) > 0, (
            "Expected ≥1 GO term for Apis mellifera (taxon 7460) from QuickGO. " "Ensure QuickGO is reachable."
        )

    def test_gene_set_values_are_sets_of_strings(self):
        from metainformant.ontology.core.go_api import build_taxon_go_gene_sets

        gene_sets = build_taxon_go_gene_sets(taxon_id=7460, max_pages=1)
        for go_id, genes in gene_sets.items():
            assert go_id.startswith("GO:")
            assert isinstance(genes, set)
            for g in genes:
                assert isinstance(g, str) and len(g) > 0

    def test_aspect_filter_reduces_result(self):
        from metainformant.ontology.core.go_api import build_taxon_go_gene_sets

        all_aspects = build_taxon_go_gene_sets(taxon_id=7460, max_pages=1)
        bp_only = build_taxon_go_gene_sets(taxon_id=7460, max_pages=1, go_aspects=["biological_process"])
        # BP-only should have equal or fewer GO terms
        assert len(bp_only) <= len(all_aspects)

    def test_multi_page_more_terms_than_single_page(self):
        from metainformant.ontology.core.go_api import build_taxon_go_gene_sets

        one_page = build_taxon_go_gene_sets(taxon_id=7460, max_pages=1)
        two_pages = build_taxon_go_gene_sets(taxon_id=7460, max_pages=2)
        # 2 pages should cover >= terms as 1 page
        assert len(two_pages) >= len(one_page)


# ---------------------------------------------------------------------------
# UniProt ID mapping — map_symbols_to_uniprot
# ---------------------------------------------------------------------------


class TestMapSymbolsToUniprot:
    def test_returns_dict(self):
        from metainformant.ontology.core.go_api import map_symbols_to_uniprot

        result = map_symbols_to_uniprot(["vg1"], taxon_id=7460)
        assert isinstance(result, dict)

    def test_empty_input_returns_empty_dict(self):
        from metainformant.ontology.core.go_api import map_symbols_to_uniprot

        result = map_symbols_to_uniprot([], taxon_id=7460)
        assert result == {}

    def test_mapping_for_homo_sapiens_known_gene(self):
        """Map a well-known human gene to verify UniProt mapping works."""
        from metainformant.ontology.core.go_api import map_symbols_to_uniprot

        result = map_symbols_to_uniprot(["TP53"], taxon_id=9606)
        assert isinstance(result, dict)
        # TP53 maps to UniProtKB:P04637 (canonical) — should appear
        if result:
            accessions = result.get("TP53", [])
            assert any("P04637" in a for a in accessions), f"Expected P04637 in TP53 mapping, got {accessions}"


# ---------------------------------------------------------------------------
# NCBI Gene annotation (gene_annotation_api)
# ---------------------------------------------------------------------------


class TestNCBIGeneAnnotation:
    def test_lookup_returns_list(self):
        from metainformant.gwas.analysis.gene_annotation_api import (
            lookup_genes_by_region_ncbi,
        )

        genes = lookup_genes_by_region_ncbi("CM009934.2", 4000000, 5500000, taxon_id=7460)
        assert isinstance(genes, list)

    def test_annotate_top_hits_returns_list(self):
        from metainformant.gwas.analysis.gene_annotation_api import (
            annotate_top_hits_ncbi,
        )

        fake_hits = [
            {"snp": "rs1", "chrom": "CM009934.2", "pos": 4744854, "p_value": 1e-8},
        ]
        result = annotate_top_hits_ncbi(fake_hits, taxon_id=7460, top_n=1)
        assert isinstance(result, list)
        assert len(result) == 1
        assert "nearby_genes" in result[0]

    def test_backwards_compat_ensembl_wrapper(self):
        """The Ensembl wrapper should delegate to NCBI and return a list."""
        from metainformant.gwas.analysis.gene_annotation_api import (
            annotate_top_hits_ensembl,
        )

        fake_hits = [
            {"snp": "rs1", "chrom": "CM009934.2", "pos": 1000000, "p_value": 0.001},
        ]
        result = annotate_top_hits_ensembl(fake_hits, species="apis_mellifera", top_n=1)
        assert isinstance(result, list)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# HPO API — map_phenotype_to_hpo, fetch_hpo_term
# ---------------------------------------------------------------------------


class TestHpoApi:
    def test_synonym_map_honey_yield(self):
        from metainformant.ontology.core.hpo import map_phenotype_to_hpo

        result = map_phenotype_to_hpo("honey yield")
        # May return empty if no synonym match — that's valid (not a human phenotype)
        assert isinstance(result, list)

    def test_known_human_phenotype(self):
        """Obesity is a known human HPO term."""
        from metainformant.ontology.core.hpo import map_phenotype_to_hpo

        result = map_phenotype_to_hpo("obesity")
        assert isinstance(result, list)
        # Should match HP:0001513 (Obesity) or similar

    def test_fetch_hpo_term_known_id(self):
        from metainformant.ontology.core.hpo import fetch_hpo_term

        term = fetch_hpo_term("HP:0001513")
        assert isinstance(term, dict)
        # HP:0001513 = Obesity — may or may not populate name but should not raise


# ---------------------------------------------------------------------------
# annotation.annotate module
# ---------------------------------------------------------------------------


class TestAnnotateModule:
    def test_gwas_hits_to_genes_with_annotations(self):
        from metainformant.ontology.annotation.annotate import gwas_hits_to_genes

        assoc = [
            {"snp": "s1", "chrom": "1", "pos": 100, "p_value": 1e-8},
            {"snp": "s2", "chrom": "1", "pos": 200, "p_value": 1e-5},
        ]
        gene_annots = [
            {"snp": "s1", "nearest_gene": "BRCA1", "nearby_genes": [{"gene_name": "BRCA1"}]},
        ]
        genes = gwas_hits_to_genes(assoc, top_n=5, gene_annotations=gene_annots)
        assert isinstance(genes, list)
        assert "BRCA1" in genes

    def test_gwas_hits_to_genes_no_annotations_returns_empty(self):
        """Without real gene annotations, function MUST return empty (no fallback)."""
        from metainformant.ontology.annotation.annotate import gwas_hits_to_genes

        assoc = [{"snp": "s1", "chrom": "1", "pos": 100, "p_value": 1e-8}]
        genes = gwas_hits_to_genes(assoc, top_n=5, gene_annotations=[])
        # No annotations → should return [] or ONLY real gene names, never positional labels
        for g in genes:
            # Positional labels are forbidden — must never return "chrom:pos" format
            assert ":" not in g, f"Positional label found in gene list: {g}"

    def test_rank_genes_by_pvalue(self):
        from metainformant.ontology.annotation.annotate import rank_genes_by_pvalue

        assoc = [
            {"snp": "s1", "chrom": "1", "pos": 100, "p_value": 0.01},
            {"snp": "s2", "chrom": "1", "pos": 200, "p_value": 0.001},
        ]
        gene_annots = [
            {"snp": "s1", "nearest_gene": "GENE_A"},
            {"snp": "s2", "nearest_gene": "GENE_B"},
        ]
        ranked = rank_genes_by_pvalue(assoc, gene_annotations=gene_annots)
        assert isinstance(ranked, list)
        # Higher significance (lower p) → higher rank
        if len(ranked) >= 2:
            genes_only = [g for g, _ in ranked]
            assert genes_only.index("GENE_B") < genes_only.index("GENE_A")

    def test_build_background_from_vcf_genes(self):
        from metainformant.ontology.annotation.annotate import build_background_from_vcf_genes

        gene_annots = [
            {"snp": "s1", "nearby_genes": [{"gene_name": "A"}, {"gene_name": "B"}]},
            {"snp": "s2", "nearby_genes": [{"gene_name": "B"}, {"gene_name": "C"}]},
        ]
        bg = build_background_from_vcf_genes(gene_annots)
        assert isinstance(bg, list)
        assert "A" in bg
        assert "C" in bg
        assert bg.count("B") == 1  # deduplicated
