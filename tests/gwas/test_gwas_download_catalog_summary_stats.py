"""Regression tests: the GWAS Catalog association table of ``download_variant_data`` is parsed, paged and checked.

Defects of the summary-statistics branch, all reproduced from the code itself:

* the p-value was built as ``assoc.get("pvalueMantissa", "") + "e" + str(...)``, a ``TypeError`` for the JSON
  numbers the service sends (``2`` and ``-8``), and the text ``"e"`` for an association without a p-value;
* the ``variant`` column received the whole ``strongestRiskAlleles[0]`` dict (its repr), and ``risk_allele`` was
  always empty;
* the table request was only read when its status was 200: any other status was skipped, and the function
  returned a study directory without an association table;
* ``associations[:1000]`` cut the table off without a word, and only the first page of the service was ever read.

No real response could be fetched when these were fixed, so the fixtures below are hand-written JSON that
follows the shape the code itself expects (documented at ``parse_catalog_associations``); they are data, not a
stand-in for the service. The parsers are pure functions of response text, and ``AssociationsPager`` is a
sans-IO state machine, so everything here runs without a network.
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any, Dict, List

import pytest
import requests

from metainformant.gwas.data import download
from metainformant.gwas.data.download import download_variant_data

STUDY_URL = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000001"
PROXY_VARIABLES = ("HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "http_proxy", "https_proxy", "all_proxy")


@pytest.fixture(autouse=True)
def no_external_network() -> Iterator[None]:
    """Route any HTTP request of the code under test to a closed local port (nothing leaves the machine)."""
    saved = {name: os.environ.get(name) for name in (*PROXY_VARIABLES, "NO_PROXY", "no_proxy")}
    for name in PROXY_VARIABLES:
        os.environ[name] = "http://127.0.0.1:9"
    os.environ["NO_PROXY"] = os.environ["no_proxy"] = ""
    try:
        yield
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


STUDY_JSON = """
{
  "accessionId": "GCST000001",
  "publicationInfo": {"pubmedId": "17554300", "title": "A made-up study"},
  "diseaseTrait": {"traitName": "Type 2 diabetes"}
}
"""

#: First of two pages: three associations, JSON-number p-values, a two-allele haplotype and absent OR / beta.
PAGE_0_JSON = """
{
  "_embedded": {"associations": [
    {
      "pvalueMantissa": 2, "pvalueExponent": -8, "orPerCopyNum": 1.07, "betaNum": null,
      "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs7329174-A"}]}]
    },
    {
      "pvalueMantissa": 1.5, "pvalueExponent": -12, "orPerCopyNum": null, "betaNum": -0.25,
      "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs1111-G"}, {"riskAlleleName": "rs2222-T"}]}]
    },
    {
      "pvalueMantissa": "3", "pvalueExponent": "-9", "orPerCopyNum": "1.2",
      "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "chr10:114758349-?"}]}]
    }
  ]},
  "_links": {"next": {"href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000001/associations?page=1&size=3"}},
  "page": {"size": 3, "totalElements": 5, "totalPages": 2, "number": 0}
}
"""

#: Second and last page: two associations, no ``next`` link.
PAGE_1_JSON = """
{
  "_embedded": {"associations": [
    {
      "pvalueMantissa": 9.9, "pvalueExponent": -7, "betaNum": 0.1,
      "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs4444-C"}]}]
    },
    {
      "pvalueMantissa": 4, "pvalueExponent": -10, "orPerCopyNum": 0.9,
      "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs5555-T"}]}]
    }
  ]},
  "_links": {"self": {"href": "..."}},
  "page": {"size": 3, "totalElements": 5, "totalPages": 2, "number": 1}
}
"""


def association_json(index: int) -> Dict[str, Any]:
    """A minimal association record whose variant is ``rs<index>``."""
    return {
        "pvalueMantissa": 1 + index % 9,
        "pvalueExponent": -8,
        "loci": [{"strongestRiskAlleles": [{"riskAlleleName": f"rs{index}-A"}]}],
    }


def page_json(indexes: List[int], number: int, size: int, total: int) -> str:
    """One page of ``total`` associations numbered ``0 .. total - 1`` holding ``indexes``."""
    pages = max(1, math.ceil(total / size))
    document: Dict[str, Any] = {
        "_embedded": {"associations": [association_json(i) for i in indexes]},
        "page": {"size": size, "totalElements": total, "totalPages": pages, "number": number},
    }
    return json.dumps(document)


# ---------------------------------------------------------------------------
# p-value formatting: numeric fields
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mantissa, exponent, expected",
    [
        (2, -8, "2e-8"),  # JSON numbers: the old string concatenation raised TypeError
        (1.5, -12, "1.5e-12"),
        (2.0, -8.0, "2e-8"),  # whole floats are not written as 2.0e-8.0
        ("3", "-9", "3e-9"),  # numeric text still works
        (" 4 ", " -10 ", "4e-10"),
    ],
)
def test_p_value_is_written_from_numeric_or_text_fields(mantissa: Any, exponent: Any, expected: str) -> None:
    assert download.format_catalog_p_value(mantissa, exponent) == expected


def test_a_p_value_text_reads_back_as_the_number_it_came_from() -> None:
    assert float(download.format_catalog_p_value(1.5, -12)) == 1.5e-12
    assert float(download.format_catalog_p_value(9.9, -7)) == 9.9e-7


@pytest.mark.parametrize("mantissa, exponent", [(None, None), ("", ""), (None, ""), ("  ", None)])
def test_no_p_value_is_an_empty_cell_not_the_text_e(mantissa: Any, exponent: Any) -> None:
    """The old code turned an association without a p-value into the string ``"e"`` (or crashed)."""
    assert download.format_catalog_p_value(mantissa, exponent) == ""


@pytest.mark.parametrize(
    "mantissa, exponent, message",
    [
        (2, None, "pair"),
        (None, -8, "pair"),
        ("abc", -8, "pvalueMantissa must be a number"),
        (2, "x", "pvalueExponent must be a number"),
        (2, -8.5, "whole number"),
        (float("nan"), -8, "finite"),
        (True, -8, "must be a number"),
    ],
)
def test_a_malformed_p_value_is_an_error_naming_the_field(mantissa: Any, exponent: Any, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        download.format_catalog_p_value(mantissa, exponent)


# ---------------------------------------------------------------------------
# the variant column
# ---------------------------------------------------------------------------


def test_the_variant_column_holds_the_variant_not_a_dict_repr() -> None:
    page = download.parse_catalog_associations(PAGE_0_JSON)

    assert [row["variant"] for row in page.rows] == ["rs7329174", "rs1111;rs2222", "chr10:114758349"]
    assert all("{" not in row["variant"] and "riskAlleleName" not in row["variant"] for row in page.rows)


def test_the_risk_allele_column_is_filled() -> None:
    page = download.parse_catalog_associations(PAGE_0_JSON)

    assert [row["risk_allele"] for row in page.rows] == ["A", "G;T", "?"]


def test_the_rows_carry_every_column_and_nothing_else() -> None:
    page = download.parse_catalog_associations(PAGE_0_JSON)

    assert all(tuple(row) == download.SUMMARY_STATS_COLUMNS for row in page.rows)
    assert page.rows[0] == {
        "variant": "rs7329174",
        "risk_allele": "A",
        "p_value": "2e-8",
        "odds_ratio": "1.07",
        "beta": "",
    }
    assert page.rows[1]["p_value"] == "1.5e-12" and page.rows[1]["odds_ratio"] == "" and page.rows[1]["beta"] == "-0.25"
    assert page.rows[2]["p_value"] == "3e-9" and page.rows[2]["odds_ratio"] == "1.2"


def test_a_name_without_an_allele_is_a_variant_with_no_stated_allele() -> None:
    payload = {"_embedded": {"associations": [{"loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs99"}]}]}]}}

    row = download.parse_catalog_associations(payload).rows[0]

    assert (row["variant"], row["risk_allele"]) == ("rs99", "")


def test_an_association_without_loci_has_empty_variant_cells_and_no_invented_value() -> None:
    row = download.parse_catalog_associations(
        {"_embedded": {"associations": [{"pvalueMantissa": 1, "pvalueExponent": -9}]}}
    ).rows[0]

    assert row == {"variant": "", "risk_allele": "", "p_value": "1e-9", "odds_ratio": "", "beta": ""}


def test_a_study_with_no_associations_gives_no_rows() -> None:
    page = download.parse_catalog_associations(
        '{"_links": {}, "page": {"size": 20, "totalElements": 0, "totalPages": 0, "number": 0}}'
    )

    assert page.rows == []
    assert (page.total_elements, page.total_pages, page.has_next) == (0, 0, False)


@pytest.mark.parametrize(
    "payload, message",
    [
        ("[]", "must be a JSON object"),
        ("not json", "not valid JSON"),
        ('{"_embedded": []}', "_embedded must be a JSON object"),
        ('{"_embedded": {"associations": {}}}', "associations must be a list"),
        ('{"_embedded": {"associations": [3]}}', "association 0 must be a JSON object"),
        ('{"_embedded": {"associations": [{"loci": {}}]}}', "loci must be a list"),
        (
            '{"_embedded": {"associations": [{"loci": [{"strongestRiskAlleles": [{"riskAlleleName": 7}]}]}]}}',
            "not text",
        ),
        (
            '{"_embedded": {"associations": [{"pvalueMantissa": "x", "pvalueExponent": -8}]}}',
            "association 0: pvalueMantissa",
        ),
        ('{"page": {"number": -1}}', "page.number"),
        ('{"page": 3}', "page must be a JSON object"),
    ],
)
def test_a_response_of_the_wrong_shape_is_an_error(payload: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        download.parse_catalog_associations(payload)


def test_the_error_names_the_association_that_is_malformed() -> None:
    bad = {"_embedded": {"associations": [association_json(0), association_json(1), {"orPerCopyNum": "n/a"}]}}

    with pytest.raises(ValueError, match="association 2: orPerCopyNum must be a number"):
        download.parse_catalog_associations(bad)


def test_a_decoded_mapping_is_accepted_as_well_as_text() -> None:
    assert (
        download.parse_catalog_associations(json.loads(PAGE_0_JSON)).rows
        == download.parse_catalog_associations(PAGE_0_JSON).rows
    )


# ---------------------------------------------------------------------------
# the study block
# ---------------------------------------------------------------------------


def test_the_study_trait_and_pubmed_id_are_read() -> None:
    assert download.parse_catalog_study(STUDY_JSON) == {"trait": "Type 2 diabetes", "pubmed_id": "17554300"}


@pytest.mark.parametrize(
    "payload",
    [
        '{"publicationInfo": null, "diseaseTrait": null}',
        "{}",
        '{"publicationInfo": {}, "diseaseTrait": {"traitName": ""}}',
    ],
)
def test_a_study_without_trait_or_publication_gives_none_instead_of_crashing(payload: str) -> None:
    """``publicationInfo: null`` used to raise ``AttributeError`` on ``None.get``."""
    assert download.parse_catalog_study(payload) == {"trait": None, "pubmed_id": None}


def test_a_numeric_pubmed_id_is_kept_as_text() -> None:
    assert download.parse_catalog_study('{"publicationInfo": {"pubmedId": 17554300}}')["pubmed_id"] == "17554300"


# ---------------------------------------------------------------------------
# the status check: a non-200 response is an error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("status", [201, 204, 301, 400, 403, 404, 429, 500, 502, 503])
def test_any_status_but_200_raises_and_names_the_status_and_url(status: int) -> None:
    with pytest.raises(requests.HTTPError, match=rf"{STUDY_URL}/associations returned HTTP {status}"):
        download.decode_catalog_response(status, '{"_embedded": {"associations": []}}', f"{STUDY_URL}/associations")


def test_the_body_of_a_failed_response_is_quoted_in_the_error() -> None:
    with pytest.raises(requests.HTTPError, match="Service Unavailable"):
        download.decode_catalog_response(503, "<html>\n  Service   Unavailable\n</html>", STUDY_URL)


def test_a_200_response_is_decoded() -> None:
    assert download.decode_catalog_response(200, STUDY_JSON, STUDY_URL)["accessionId"] == "GCST000001"


@pytest.mark.parametrize(
    "body, message", [("", "not valid JSON"), ("<html></html>", "not valid JSON"), ("[1]", "JSON object")]
)
def test_a_200_response_that_is_not_a_json_object_is_an_error(body: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        download.decode_catalog_response(200, body, STUDY_URL)


# ---------------------------------------------------------------------------
# paging: the whole table, or an explicit report
# ---------------------------------------------------------------------------


def drive(pager: download.AssociationsPager, pages: Dict[int, str]) -> List[Dict[str, int]]:
    """Run the pager against response texts keyed by page number; return the queries it made."""
    queries: List[Dict[str, int]] = []
    while (query := pager.request()) is not None:
        queries.append(query)
        pager.accept(download.parse_catalog_associations(pages[query["page"]]))
    pager.check_complete()
    return queries


def test_every_page_is_read_and_the_rows_concatenate() -> None:
    pager = download.AssociationsPager(page_size=3)

    queries = drive(pager, {0: PAGE_0_JSON, 1: PAGE_1_JSON})

    assert queries == [{"page": 0, "size": 3}, {"page": 1, "size": 3}]
    assert [row["variant"] for row in pager.rows] == [
        "rs7329174",
        "rs1111;rs2222",
        "chr10:114758349",
        "rs4444",
        "rs5555",
    ]
    assert (pager.available, pager.truncated, pager.pages_read) == (5, False, 2)


def test_a_table_longer_than_the_old_1000_row_limit_is_downloaded_whole() -> None:
    total, size = 2345, 500
    pages = {
        n: page_json(list(range(n * size, min(total, (n + 1) * size))), n, size, total)
        for n in range(math.ceil(total / size))
    }
    pager = download.AssociationsPager(page_size=size)

    queries = drive(pager, pages)

    assert len(queries) == 5
    assert len(pager.rows) == total == pager.available
    assert pager.truncated is False
    assert pager.rows[-1]["variant"] == f"rs{total - 1}"
    assert len({row["variant"] for row in pager.rows}) == total  # no association twice, none skipped


def test_the_default_page_size_is_requested_when_none_is_given() -> None:
    assert download.AssociationsPager().request() == {"page": 0, "size": download.DEFAULT_ASSOCIATIONS_PAGE_SIZE}


def test_without_a_page_block_the_next_link_decides() -> None:
    first = json.dumps({"_embedded": {"associations": [association_json(0)]}, "_links": {"next": {"href": "x"}}})
    second = json.dumps({"_embedded": {"associations": [association_json(1)]}})
    pager = download.AssociationsPager(page_size=1)

    queries = drive(pager, {0: first, 1: second})

    assert [q["page"] for q in queries] == [0, 1]
    assert [row["variant"] for row in pager.rows] == ["rs0", "rs1"]
    assert pager.available is None  # nothing to check the count against


def test_a_single_page_is_one_request() -> None:
    only = page_json([0, 1], number=0, size=20, total=2)
    pager = download.AssociationsPager()

    assert drive(pager, {0: only}) == [{"page": 0, "size": download.DEFAULT_ASSOCIATIONS_PAGE_SIZE}]
    assert len(pager.rows) == 2


def test_a_service_that_ignores_the_page_parameter_is_an_error() -> None:
    pager = download.AssociationsPager(page_size=3)
    pager.request()
    pager.accept(download.parse_catalog_associations(PAGE_0_JSON))
    pager.request()

    with pytest.raises(ValueError, match="requested association page 1 but the service returned page 0"):
        pager.accept(download.parse_catalog_associations(PAGE_0_JSON))  # page 0 served again


def test_fewer_rows_than_the_service_reports_is_an_error() -> None:
    short = page_json([0, 1], number=0, size=20, total=3)  # claims three associations, carries two
    pager = download.AssociationsPager()
    pager.request()
    pager.accept(download.parse_catalog_associations(short))

    with pytest.raises(RuntimeError, match="reports 3 associations but 2 were read"):
        pager.check_complete()


def test_a_service_that_never_stops_announcing_pages_is_stopped() -> None:
    endless = json.dumps({"_embedded": {"associations": []}, "_links": {"next": {"href": "x"}}})
    pager = download.AssociationsPager()

    with pytest.raises(RuntimeError, match=f"more than {download.MAX_ASSOCIATION_PAGES} pages"):
        for _ in range(download.MAX_ASSOCIATION_PAGES + 1):
            assert pager.request() is not None
            pager.accept(download.parse_catalog_associations(endless))


def test_accepting_a_page_when_none_is_outstanding_is_an_error() -> None:
    pager = download.AssociationsPager()
    pager.request()
    pager.accept(download.parse_catalog_associations(page_json([0], 0, 20, 1)))

    assert pager.request() is None
    with pytest.raises(RuntimeError, match="every page has been read"):
        pager.accept(download.parse_catalog_associations(PAGE_1_JSON))


@pytest.mark.parametrize("size", [0, -1, 1.5, True, "20", None])
def test_a_bad_page_size_is_rejected(size: Any) -> None:
    with pytest.raises(ValueError, match="page_size"):
        download.AssociationsPager(page_size=size)


@pytest.mark.parametrize("cap", [0, -3, 2.5, True, "5"])
def test_a_bad_cap_is_rejected(cap: Any) -> None:
    with pytest.raises(ValueError, match="max_associations"):
        download.AssociationsPager(max_associations=cap)


# ---------------------------------------------------------------------------
# an explicit cap is reported, never silent
# ---------------------------------------------------------------------------


def test_a_cap_below_the_table_cuts_it_and_says_so() -> None:
    pager = download.AssociationsPager(page_size=3, max_associations=4)

    queries = drive(pager, {0: PAGE_0_JSON, 1: PAGE_1_JSON})

    assert len(queries) == 2  # the second page is needed to reach the fourth row
    assert [row["variant"] for row in pager.rows] == ["rs7329174", "rs1111;rs2222", "chr10:114758349", "rs4444"]
    assert pager.truncated is True
    assert pager.available == 5  # the study has more than was kept


def test_a_cap_that_ends_on_a_page_boundary_stops_requesting_and_still_says_truncated() -> None:
    pager = download.AssociationsPager(page_size=3, max_associations=3)

    queries = drive(pager, {0: PAGE_0_JSON})

    assert len(queries) == 1
    assert len(pager.rows) == 3
    assert pager.truncated is True  # the service announced a second page


def test_a_cap_equal_to_the_table_size_is_not_a_truncation() -> None:
    pager = download.AssociationsPager(page_size=3, max_associations=5)

    drive(pager, {0: PAGE_0_JSON, 1: PAGE_1_JSON})

    assert len(pager.rows) == 5
    assert pager.truncated is False


def test_a_cap_above_the_table_size_changes_nothing() -> None:
    pager = download.AssociationsPager(page_size=3, max_associations=100)

    drive(pager, {0: PAGE_0_JSON, 1: PAGE_1_JSON})

    assert len(pager.rows) == 5
    assert pager.truncated is False


# ---------------------------------------------------------------------------
# download_variant_data: argument checks still come before any request or write
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("accession", ["GCST", "GCST12a", "gcst1", "GCST٣٤٥", "GCST1/../x", "../GCST1", "GCST1\n"])
def test_an_accession_must_be_ascii_digits_after_gcst(tmp_path: Path, accession: str) -> None:
    """The pattern used ``\\d``, which also matches the digits of other scripts."""
    with pytest.raises(ValueError, match="Invalid GWAS Catalog study accession"):
        download_variant_data(accession, tmp_path / "out", data_type="summary_stats")

    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "kwargs, name", [({"page_size": 0}, "page_size"), ({"max_associations": 0}, "max_associations")]
)
def test_paging_arguments_are_checked_before_anything_is_requested_or_written(
    tmp_path: Path, kwargs: Dict[str, Any], name: str
) -> None:
    with pytest.raises(ValueError, match=name):
        download_variant_data("GCST000001", tmp_path / "out", data_type="summary_stats", **kwargs)

    assert not (tmp_path / "out").exists()


def test_the_summary_columns_and_service_url_are_the_documented_ones() -> None:
    assert download.SUMMARY_STATS_COLUMNS == ("variant", "risk_allele", "p_value", "odds_ratio", "beta")
    assert download.GWAS_CATALOG_API_URL == "https://www.ebi.ac.uk/gwas/rest/api"


def test_the_table_is_written_with_the_documented_header_and_one_line_per_association(tmp_path: Path) -> None:
    rows = download.parse_catalog_associations(PAGE_0_JSON).rows

    path = download.write_summary_stats_table(tmp_path / "t.tsv", rows)

    assert path == tmp_path / "t.tsv"
    assert path.read_text(encoding="utf-8") == (
        "variant\trisk_allele\tp_value\todds_ratio\tbeta\n"
        "rs7329174\tA\t2e-8\t1.07\t\n"
        "rs1111;rs2222\tG;T\t1.5e-12\t\t-0.25\n"
        "chr10:114758349\t?\t3e-9\t1.2\t\n"
    )


def test_a_table_with_no_rows_is_just_the_header(tmp_path: Path) -> None:
    path = download.write_summary_stats_table(tmp_path / "empty.tsv", [])

    assert path.read_text(encoding="utf-8") == "variant\trisk_allele\tp_value\todds_ratio\tbeta\n"


def test_a_tab_quote_or_newline_inside_a_cell_cannot_add_columns_or_rows(tmp_path: Path) -> None:
    rows = [{"variant": "rs1\twith tab", "risk_allele": 'A"B', "p_value": "1e-8", "odds_ratio": "", "beta": "x\ny"}]

    path = download.write_summary_stats_table(tmp_path / "t.tsv", rows)

    with open(path, newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == rows


def test_the_1000_genomes_download_has_no_population_argument_it_would_ignore() -> None:
    """The release VCFs hold every sample, so there is no per-population file; the argument was accepted and unused."""
    import inspect

    assert list(inspect.signature(download._download_1000genomes).parameters) == ["output_dir", "chromosome"]
