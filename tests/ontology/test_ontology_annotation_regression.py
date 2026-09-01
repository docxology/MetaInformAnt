"""Tests for the ontology annotation bridge (GWAS hits -> genes -> GO sets).

Round-4 test-depth lane T4: pins regression behaviour of
``metainformant.ontology.annotation.annotate`` using real files (GAF parsing)
and real computation. No test doubles; API paths are not exercised.
"""

from __future__ import annotations

import gzip
import math

import pytest

from metainformant.ontology.annotation.annotate import (
    build_background_from_vcf_genes,
    genes_to_go_annotations,
    gwas_hits_to_genes,
    rank_genes_by_pvalue,
)

GAF_SAMPLE = (
    "![GAF header line]\n"
    "UniProtKB\tP12345\tGENE1\t\tGO:0003824\tPMID:1\tIDA\t\tF\t-\n"
    "UniProtKB\tP12345\tGENE1\t\tGO:0008152\tPMID:2\tIEA\t\tP\t-\n"
    "UniProtKB\tP67890\tGENE2\t\tGO:0005634\tPMID:3\tISS\t\tC\t-\n"
    "UniProtKB\tP11111\tOUTSIDER\t\tGO:0006915\tPMID:4\tEXP\t\tP\t-\n"
)


# ---------------------------------------------------------------------------
# gwas_hits_to_genes
# ---------------------------------------------------------------------------


class TestGwasHitsToGenes:
    def _hits(self):
        return [
            {"snp": "rs1", "p_value": 1e-8, "chrom": "1", "pos": 100},
            {"snp": "rs2", "p_value": 1e-6, "chrom": "2", "pos": 200},
            {"snp": "rs3", "p_value": 1e-4, "chrom": "3", "pos": 300},
            {"snp": "rs_unannotated", "p_value": 1e-3, "chrom": "4", "pos": 400},
        ]

    def _annotations(self):
        return [
            {"snp": "rs1", "nearest_gene": "GENEA", "nearby_genes": [{"gene_name": "GENEB"}]},
            {"snp": "rs2", "nearest_gene": "GENEC", "nearby_genes": []},
        ]

    def test_top_n_limits(self):
        genes = gwas_hits_to_genes(self._hits(), top_n=2, gene_annotations=self._annotations())
        assert set(genes) == {"GENEA", "GENEB", "GENEC"}
        assert len(genes) <= 2 + 1  # nearest + nearby may expand one hit

    def test_unannotated_hits_skipped(self):
        # rs_unannotated has no annotation: no positional labels, no SNP placeholders
        genes = gwas_hits_to_genes(self._hits(), top_n=10, gene_annotations=self._annotations())
        assert "rs_unannotated" not in genes
        assert "4:400" not in genes

    def test_deduplication(self):
        annotations = [
            {"snp": "rs1", "nearest_gene": "GENEA", "nearby_genes": [{"gene_name": "GENEA"}]},
        ]
        genes = gwas_hits_to_genes(self._hits(), top_n=5, gene_annotations=annotations)
        assert genes == ["GENEA"]

    def test_sorted_by_pvalue_not_input_order(self):
        hits = [
            {"snp": "rsB", "p_value": 1e-2},
            {"snp": "rsA", "p_value": 1e-9},
        ]
        annotations = [
            {"snp": "rsA", "nearest_gene": "STRONG", "nearby_genes": []},
            {"snp": "rsB", "nearest_gene": "WEAK", "nearby_genes": []},
        ]
        genes = gwas_hits_to_genes(hits, top_n=1, gene_annotations=annotations)
        assert genes == ["STRONG"]

    def test_no_annotations_returns_empty(self):
        assert gwas_hits_to_genes(self._hits(), top_n=10) == []


# ---------------------------------------------------------------------------
# rank_genes_by_pvalue
# ---------------------------------------------------------------------------


class TestRankGenesByPvalue:
    def test_ranking_order_and_metric(self):
        hits = [
            {"snp": "rs1", "p_value": 1e-8},
            {"snp": "rs2", "p_value": 1e-4},
        ]
        annotations = [
            {"snp": "rs1", "nearest_gene": "G1"},
            {"snp": "rs2", "nearest_gene": "G2"},
        ]
        ranked = rank_genes_by_pvalue(hits, gene_annotations=annotations)
        assert ranked[0] == ("G1", pytest.approx(8.0))
        assert ranked[1] == ("G2", pytest.approx(4.0))
        assert ranked[0][1] >= ranked[-1][1]

    def test_unannotated_snps_excluded(self):
        hits = [{"snp": "rsX", "p_value": 1e-10}]
        assert rank_genes_by_pvalue(hits) == []

    def test_same_gene_takes_most_significant(self):
        hits = [
            {"snp": "rs1", "p_value": 1e-3},
            {"snp": "rs2", "p_value": 1e-9},
        ]
        annotations = [
            {"snp": "rs1", "nearest_gene": "SHARED"},
            {"snp": "rs2", "nearest_gene": "SHARED"},
        ]
        ranked = rank_genes_by_pvalue(hits, gene_annotations=annotations)
        assert len(ranked) == 1
        assert ranked[0][1] == pytest.approx(9.0)

    def test_min_neg_log_p_filter(self):
        hits = [
            {"snp": "rs1", "p_value": 1e-8},
            {"snp": "rs2", "p_value": 0.5},
        ]
        annotations = [
            {"snp": "rs1", "nearest_gene": "G1"},
            {"snp": "rs2", "nearest_gene": "G2"},
        ]
        ranked = rank_genes_by_pvalue(hits, gene_annotations=annotations, min_neg_log_p=5.0)
        assert [g for g, _ in ranked] == ["G1"]

    def test_extreme_pvalue_clamped(self):
        hits = [{"snp": "rs1", "p_value": 0.0}]
        annotations = [{"snp": "rs1", "nearest_gene": "G1"}]
        ranked = rank_genes_by_pvalue(hits, gene_annotations=annotations)
        # Clamped at 1e-300 -> -log10 = 300, not -inf
        assert ranked[0][1] == pytest.approx(300.0)
        assert math.isfinite(ranked[0][1])


# ---------------------------------------------------------------------------
# genes_to_go_annotations (GAF path — real file, no network)
# ---------------------------------------------------------------------------


class TestGenesToGoAnnotations:
    def test_gaf_parsing_plain(self, tmp_path):
        gaf = tmp_path / "annot.gaf"
        gaf.write_text(GAF_SAMPLE)
        result = genes_to_go_annotations(["GENE1", "GENE2"], source="gaf", annotation_file=str(gaf))
        assert result == {"GO:0003824": {"GENE1"}, "GO:0008152": {"GENE1"}, "GO:0005634": {"GENE2"}}

    def test_gaf_parsing_gzipped(self, tmp_path):
        gaf = tmp_path / "annot.gaf.gz"
        with gzip.open(gaf, "wt") as fh:
            fh.write(GAF_SAMPLE)
        result = genes_to_go_annotations(["GENE1", "GENE2"], source="gaf", annotation_file=str(gaf))
        assert len(result) == 3

    def test_gaf_aspect_filter(self, tmp_path):
        gaf = tmp_path / "annot.gaf"
        gaf.write_text(GAF_SAMPLE)
        result = genes_to_go_annotations(["GENE1", "GENE2"], source="gaf", annotation_file=str(gaf), go_aspects=["P"])
        assert result == {"GO:0008152": {"GENE1"}}

    def test_gaf_short_aspect_codes_accepted(self, tmp_path):
        gaf = tmp_path / "annot.gaf"
        gaf.write_text(GAF_SAMPLE)
        result = genes_to_go_annotations(
            ["GENE1", "GENE2"], source="gaf", annotation_file=str(gaf), go_aspects=["F", "C"]
        )
        assert set(result) == {"GO:0003824", "GO:0005634"}

    def test_gaf_outsider_gene_excluded(self, tmp_path):
        gaf = tmp_path / "annot.gaf"
        gaf.write_text(GAF_SAMPLE)
        result = genes_to_go_annotations(["GENE1", "GENE2"], source="gaf", annotation_file=str(gaf))
        assert all("OUTSIDER" not in genes for genes in result.values())

    def test_gaf_missing_file_returns_empty(self, tmp_path):
        result = genes_to_go_annotations(["GENE1"], source="gaf", annotation_file=str(tmp_path / "nope.gaf"))
        assert result == {}

    def test_gaf_no_file_returns_empty(self):
        assert genes_to_go_annotations(["GENE1"], source="gaf") == {}

    def test_unknown_source_raises(self):
        with pytest.raises(ValueError, match="Unknown annotation source"):
            genes_to_go_annotations(["GENE1"], source="magic")


# ---------------------------------------------------------------------------
# build_background_from_vcf_genes
# ---------------------------------------------------------------------------


class TestBuildBackgroundFromVcfGenes:
    def test_collects_nearest_and_nearby(self):
        annotations = [
            {"nearest_gene": "GENEB", "nearby_genes": [{"gene_name": "GENEA"}, {"gene_name": "GENEC"}]},
            {"nearest_gene": "GENED", "nearby_genes": []},
        ]
        result = build_background_from_vcf_genes(annotations)
        assert result == ["GENEA", "GENEB", "GENEC", "GENED"]  # sorted, deduped

    def test_deduplication(self):
        annotations = [
            {"nearest_gene": "G1", "nearby_genes": [{"gene_name": "G1"}]},
            {"nearest_gene": "G1", "nearby_genes": []},
        ]
        assert build_background_from_vcf_genes(annotations) == ["G1"]

    def test_empty_annotations(self):
        assert build_background_from_vcf_genes([]) == []

    def test_empty_gene_fields_skipped(self):
        annotations = [{"nearest_gene": "", "nearby_genes": [{"gene_name": ""}]}]
        assert build_background_from_vcf_genes(annotations) == []
