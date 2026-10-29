# NCBI Datasets v2 fixtures

Responses of the NCBI Datasets v2 REST API (`https://api.ncbi.nlm.nih.gov/datasets/v2`), saved so the
request builders and response parsers in `src/metainformant/dna/external/ncbi_datasets.py` can be tested
offline on what the service really sends. Used by `tests/dna/test_dna_ncbi_datasets_v2.py`.

## Provenance

All files were captured on **2026-10-28** with a plain `GET` through the WebFetch tool of the agent session
`r4-ncbi`. The tool hands back the response body but not the status line or headers, so the body was written
to disk as returned (NCBI's `{"key":value, "key":value}` serialization is kept; one trailing newline was
added). The tool shows no body for HTTP 4xx responses, so there are no fixtures of error *bodies*; the status
codes below were read from the tool's error message.

The tool passes a page through a language model before returning it, so every file was checked after saving:
all of them validate against the schemas of the OpenAPI specification (`v2reportsAssemblyDataReportPage`,
`v2SequenceReportPage`, `v2DownloadSummary`: every key is a schema field, with the schema's type and enum
values; the specification itself is `https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml`,
554 KB, not stored here).

Cross-check: the Apis mellifera report was fetched through three different endpoints (accession lookup,
`reference_only` taxon lookup, `assembly_source=refseq` taxon lookup); the three `reports[0]` objects are
identical, and `test_honey_bee_report_is_identical_across_endpoints` asserts it.

| File | Request (`GET https://api.ncbi.nlm.nih.gov/datasets/v2` + ...) | Shows |
|------|------------------------------------------------------------------|-------|
| `dataset_report_accession_GCF_003254395.2.json` | `/genome/accession/GCF_003254395.2/dataset_report` | Complete report, *Apis mellifera* Amel_HAv3.1; `uint64` fields as strings; `total_count`, no `next_page_token` |
| `dataset_report_accession_GCF_000819615.1.json` | `/genome/accession/GCF_000819615.1/dataset_report` | Sparse report (phage phiX174): no `submitter`, `assembly_method`, `refseq_category`, `common_name` |
| `dataset_report_accession_GCF_999999999.9_not_found.json` | `/genome/accession/GCF_999999999.9/dataset_report` | Unknown accession: HTTP 200 with `{}` |
| `dataset_report_taxon_unknown_name.json` | `/genome/taxon/Notaspecies%20nonexistens/dataset_report` | Unknown taxon: HTTP 200 with `{}` |
| `dataset_report_taxon_2886930_reference_only_empty.json` | `/genome/taxon/2886930/dataset_report?filters.reference_only=true&returned_content=ASSM_ACC` | Valid taxon without a reference genome: `{"content_type":"ASSM_ACC"}`, no `reports`, no `total_count` |
| `dataset_report_taxon_Apis_mellifera_reference_only.json` | `/genome/taxon/Apis%20mellifera/dataset_report?filters.reference_only=true` | One reference genome, `refseq_category` "reference genome" |
| `dataset_report_taxon_Apis_mellifera_refseq_page1.json` | `/genome/taxon/Apis%20mellifera/dataset_report?page_size=1&filters.assembly_source=refseq` | A page that holds the only result and still carries `next_page_token` |
| `dataset_report_taxon_Apis_mellifera_refseq_page2_empty.json` | same query plus `&page_token=` the token of the file above | The follow-up page: `{"total_count":1}` |
| `dataset_report_taxon_Apis_mellifera_assm_acc_page1.json` | `/genome/taxon/Apis%20mellifera/dataset_report?page_size=1&returned_content=ASSM_ACC` | Abbreviated report, page 1 of 13 with token |
| `dataset_report_taxon_Apis_mellifera_assm_acc_page2.json` | same query plus `&page_token=` the token of page 1 | Abbreviated report, page 2 of 13 with token |
| `dataset_report_taxon_two_names_joined_by_comma.json` | `/genome/taxon/Apis%20mellifera,Pogonomyrmex%20barbatus/dataset_report?returned_content=ASSM_ACC&filters.reference_only=true` (identical answer with the comma sent as `%2C`) | One path segment, two taxa: the service splits on commas even when they are percent-encoded |
| `dataset_report_taxon_2886930_assm_acc_page1.json` | `/genome/taxon/2886930/dataset_report?returned_content=ASSM_ACC` | Default page of 20 abbreviated reports out of 150; mixed GenBank/RefSeq |
| `dataset_report_bioproject_PRJNA000000_error_message.json` | `/genome/bioproject/PRJNA000000/dataset_report` | HTTP 200 whose only content is `messages[].error` (same page schema) |
| `sequence_reports_GCF_003254395.2_page1.json` | `/genome/accession/GCF_003254395.2/sequence_reports?page_size=3` | Sequence rows (`chr_name`, `refseq_accession`, `length`, `role`), 3 of 177, with token |
| `download_summary_GCF_003254395.2_default.json` | `/genome/accession/GCF_003254395.2/download_summary` | NCBI's own description of the default package and its `download` URL |
| `download_summary_GCF_003254395.2_legacy_include_ignored.json` | `/genome/accession/GCF_003254395.2/download_summary?include=genome,gff3` | Byte-for-byte the default summary: the `include` parameter is ignored |
| `download_summary_GCF_003254395.2_genome_gff3.json` | `/genome/accession/GCF_003254395.2/download_summary?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF` | NCBI echoes the repeated-parameter URL and `--include gff3,genome` |
| `download_summary_GCF_003254395.2_rna_protein_cds_gtf_gbff_seqreport.json` | `/genome/accession/GCF_003254395.2/download_summary?include_annotation_type=RNA_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=GENOME_GTF&include_annotation_type=GENOME_GBFF&include_annotation_type=SEQUENCE_REPORT` | NCBI's CLI name for each API file type: `--include gbff,gtf,rna,cds,protein,seq-report` |

Page tokens in the `.json` files are NCBI cursors; they are opaque and may stop working, which does not
matter offline.

## Status codes seen (no fixture, the tool shows no body)

| Request | Result |
|---------|--------|
| `/genome` (what `list_genome_assemblies` called before the port) | HTTP 404 |
| `/genome/accession/GCF_003254395.2` (what `get_genome_metadata` called) | HTTP 404 |
| `/genome/taxon/Apis%20mellifera/dataset_report?filters.assembly_source=bogus&page_size=1` | HTTP 400 |
| `/genome/accession/GCF_003254395.2/download_summary?include_annotation_type=NOT_A_TYPE` | HTTP 400 |
| `/genome/accession/GCF_003254395/check` (no version) | HTTP 200, `{"invalid_assemblies":["GCF_003254395"], ...}`: versionless accessions are invalid |

Live expectations used by the network-marked tests (2026-10-28): `filters.reference_only=true` returns
GCF_003254395.2 for *Apis mellifera*, GCF_055798465.1 for *Pogonomyrmex barbatus*, GCF_040556925.1 for
*Mycoplasma genitalium* (a 580 kb genome, small enough to download in a test), and nothing for phage phiX174
or SARS-CoV-2.

## Refreshing

Re-fetch a file by sending the request in the table and replacing the body. If the service changes shape,
the tests that read these files fail first, which is the point.
