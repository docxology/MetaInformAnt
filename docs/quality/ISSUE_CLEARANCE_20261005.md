# Historical issue and PR clearance — 2026-10-05

The repository owner requested clearance of all open issues and pull requests.
The inventory contained **83 issues and one PR**. The nested Hymenoptera
repository had no open issues or PRs. Closure preserves the original GitHub
history and does not certify implementation of archived research requirements.

## Disposition

- 33 sync-test artifacts: retired as not planned.
- 13 imported priority/date metadata records: retired as not planned.
- 12 unrelated imported records: retired as not planned; personal details are not repeated here.
- 20 historical roadmap entries: three existing test/documentation requests are closed as completed; the remaining entries are archived as not planned.
- Five historical research requests, including a duplicated genome-ingestion request: archived as not planned, with the unfinished scope retained below.

## Unfinished historical research scope

- [#4](https://github.com/docxology/MetaInformAnt/issues/4): genome/RNA/protein functional annotation, locus associations, persisted Postgres outputs and reproducible reannotation. No completed database-scale EnTAP/alternative annotation campaign is asserted.
- [#5](https://github.com/docxology/MetaInformAnt/issues/5): incremental OMA hierarchical orthology groups, preserved pairwise intermediates, proteome/genome handling and database execution records. No completed OMA campaign is asserted.
- [#6](https://github.com/docxology/MetaInformAnt/issues/6): Amalgkit across the entire historical Postgres genome database. The active 27-species frozen campaign is a bounded scope; it is unfinished and does not establish completion of that larger request.
- [#7](https://github.com/docxology/MetaInformAnt/issues/7) and [#36](https://github.com/docxology/MetaInformAnt/issues/36): taxonomy-driven recent-genome acquisition and durable database ingestion. The duplicated old database acceptance criteria remain unverified.

Historical v1.0 release, PyPI/Conda publication, generalized optimization,
notebook execution and complete multi-omics integration acceptance remain
unverified. Existing software capabilities and tests do not establish those
broad release or scientific acceptance criteria. Archived issue links retain
requirements for any future explicitly scoped campaign.

## CI pull request

[#193](https://github.com/docxology/MetaInformAnt/pull/193) is closed without
merging. Its private-submodule authentication patch requires an absent
`SUBMODULES_TOKEN` secret; the hosted checks produced zero tests. It is not a
standalone working repair. A future replacement must supply independently
reviewed least-privilege read access, ensure the parent checkout can authenticate
as well as both private submodules, and demonstrate actual hosted test execution.
The patch and its discussion remain retrievable at the closed PR.

The current main build and security checks passed at `af91ef7fa`; hosted tests
remain blocked at private-submodule checkout. Closing the PR does not resolve
that credential limitation.

## Issue inventory

Numbers link to preserved original records. Imported personal/test content is
intentionally represented only by its disposition category.

| Issue | Category | Closure reason |
|---|---|---|
| [#4](https://github.com/docxology/MetaInformAnt/issues/4) | historical research request | not planned |
| [#5](https://github.com/docxology/MetaInformAnt/issues/5) | historical research request | not planned |
| [#6](https://github.com/docxology/MetaInformAnt/issues/6) | historical research request | not planned |
| [#7](https://github.com/docxology/MetaInformAnt/issues/7) | historical research request | not planned |
| [#8](https://github.com/docxology/MetaInformAnt/issues/8) | sync-test artifact | not planned |
| [#15](https://github.com/docxology/MetaInformAnt/issues/15) | sync-test artifact | not planned |
| [#16](https://github.com/docxology/MetaInformAnt/issues/16) | sync-test artifact | not planned |
| [#17](https://github.com/docxology/MetaInformAnt/issues/17) | sync-test artifact | not planned |
| [#36](https://github.com/docxology/MetaInformAnt/issues/36) | historical research request | not planned |
| [#37](https://github.com/docxology/MetaInformAnt/issues/37) | out-of-scope imported record | not planned |
| [#38](https://github.com/docxology/MetaInformAnt/issues/38) | out-of-scope imported record | not planned |
| [#39](https://github.com/docxology/MetaInformAnt/issues/39) | out-of-scope imported record | not planned |
| [#40](https://github.com/docxology/MetaInformAnt/issues/40) | out-of-scope imported record | not planned |
| [#41](https://github.com/docxology/MetaInformAnt/issues/41) | out-of-scope imported record | not planned |
| [#42](https://github.com/docxology/MetaInformAnt/issues/42) | out-of-scope imported record | not planned |
| [#43](https://github.com/docxology/MetaInformAnt/issues/43) | out-of-scope imported record | not planned |
| [#44](https://github.com/docxology/MetaInformAnt/issues/44) | out-of-scope imported record | not planned |
| [#45](https://github.com/docxology/MetaInformAnt/issues/45) | out-of-scope imported record | not planned |
| [#46](https://github.com/docxology/MetaInformAnt/issues/46) | out-of-scope imported record | not planned |
| [#47](https://github.com/docxology/MetaInformAnt/issues/47) | out-of-scope imported record | not planned |
| [#48](https://github.com/docxology/MetaInformAnt/issues/48) | out-of-scope imported record | not planned |
| [#60](https://github.com/docxology/MetaInformAnt/issues/60) | sync-test artifact | not planned |
| [#61](https://github.com/docxology/MetaInformAnt/issues/61) | sync-test artifact | not planned |
| [#62](https://github.com/docxology/MetaInformAnt/issues/62) | sync-test artifact | not planned |
| [#63](https://github.com/docxology/MetaInformAnt/issues/63) | sync-test artifact | not planned |
| [#78](https://github.com/docxology/MetaInformAnt/issues/78) | sync-test artifact | not planned |
| [#79](https://github.com/docxology/MetaInformAnt/issues/79) | sync-test artifact | not planned |
| [#80](https://github.com/docxology/MetaInformAnt/issues/80) | sync-test artifact | not planned |
| [#89](https://github.com/docxology/MetaInformAnt/issues/89) | historical roadmap entry | not planned |
| [#90](https://github.com/docxology/MetaInformAnt/issues/90) | imported metadata | not planned |
| [#91](https://github.com/docxology/MetaInformAnt/issues/91) | imported metadata | not planned |
| [#92](https://github.com/docxology/MetaInformAnt/issues/92) | historical roadmap entry | not planned |
| [#93](https://github.com/docxology/MetaInformAnt/issues/93) | historical roadmap entry | not planned |
| [#94](https://github.com/docxology/MetaInformAnt/issues/94) | imported metadata | not planned |
| [#95](https://github.com/docxology/MetaInformAnt/issues/95) | historical roadmap entry | not planned |
| [#96](https://github.com/docxology/MetaInformAnt/issues/96) | imported metadata | not planned |
| [#97](https://github.com/docxology/MetaInformAnt/issues/97) | historical roadmap entry | not planned |
| [#98](https://github.com/docxology/MetaInformAnt/issues/98) | imported metadata | not planned |
| [#99](https://github.com/docxology/MetaInformAnt/issues/99) | imported metadata | not planned |
| [#100](https://github.com/docxology/MetaInformAnt/issues/100) | historical roadmap entry | not planned |
| [#101](https://github.com/docxology/MetaInformAnt/issues/101) | historical roadmap entry | completed |
| [#102](https://github.com/docxology/MetaInformAnt/issues/102) | historical roadmap entry | not planned |
| [#103](https://github.com/docxology/MetaInformAnt/issues/103) | imported metadata | not planned |
| [#104](https://github.com/docxology/MetaInformAnt/issues/104) | historical roadmap entry | not planned |
| [#105](https://github.com/docxology/MetaInformAnt/issues/105) | historical roadmap entry | not planned |
| [#106](https://github.com/docxology/MetaInformAnt/issues/106) | historical roadmap entry | not planned |
| [#107](https://github.com/docxology/MetaInformAnt/issues/107) | imported metadata | not planned |
| [#108](https://github.com/docxology/MetaInformAnt/issues/108) | imported metadata | not planned |
| [#109](https://github.com/docxology/MetaInformAnt/issues/109) | historical roadmap entry | not planned |
| [#110](https://github.com/docxology/MetaInformAnt/issues/110) | historical roadmap entry | not planned |
| [#111](https://github.com/docxology/MetaInformAnt/issues/111) | imported metadata | not planned |
| [#112](https://github.com/docxology/MetaInformAnt/issues/112) | historical roadmap entry | completed |
| [#114](https://github.com/docxology/MetaInformAnt/issues/114) | sync-test artifact | not planned |
| [#115](https://github.com/docxology/MetaInformAnt/issues/115) | sync-test artifact | not planned |
| [#116](https://github.com/docxology/MetaInformAnt/issues/116) | sync-test artifact | not planned |
| [#125](https://github.com/docxology/MetaInformAnt/issues/125) | historical roadmap entry | not planned |
| [#126](https://github.com/docxology/MetaInformAnt/issues/126) | imported metadata | not planned |
| [#127](https://github.com/docxology/MetaInformAnt/issues/127) | imported metadata | not planned |
| [#128](https://github.com/docxology/MetaInformAnt/issues/128) | historical roadmap entry | not planned |
| [#129](https://github.com/docxology/MetaInformAnt/issues/129) | historical roadmap entry | completed |
| [#130](https://github.com/docxology/MetaInformAnt/issues/130) | historical roadmap entry | not planned |
| [#131](https://github.com/docxology/MetaInformAnt/issues/131) | imported metadata | not planned |
| [#132](https://github.com/docxology/MetaInformAnt/issues/132) | historical roadmap entry | not planned |
| [#133](https://github.com/docxology/MetaInformAnt/issues/133) | historical roadmap entry | not planned |
| [#135](https://github.com/docxology/MetaInformAnt/issues/135) | sync-test artifact | not planned |
| [#136](https://github.com/docxology/MetaInformAnt/issues/136) | sync-test artifact | not planned |
| [#137](https://github.com/docxology/MetaInformAnt/issues/137) | sync-test artifact | not planned |
| [#147](https://github.com/docxology/MetaInformAnt/issues/147) | sync-test artifact | not planned |
| [#148](https://github.com/docxology/MetaInformAnt/issues/148) | sync-test artifact | not planned |
| [#149](https://github.com/docxology/MetaInformAnt/issues/149) | sync-test artifact | not planned |
| [#151](https://github.com/docxology/MetaInformAnt/issues/151) | sync-test artifact | not planned |
| [#152](https://github.com/docxology/MetaInformAnt/issues/152) | sync-test artifact | not planned |
| [#153](https://github.com/docxology/MetaInformAnt/issues/153) | sync-test artifact | not planned |
| [#163](https://github.com/docxology/MetaInformAnt/issues/163) | sync-test artifact | not planned |
| [#168](https://github.com/docxology/MetaInformAnt/issues/168) | sync-test artifact | not planned |
| [#169](https://github.com/docxology/MetaInformAnt/issues/169) | sync-test artifact | not planned |
| [#170](https://github.com/docxology/MetaInformAnt/issues/170) | sync-test artifact | not planned |
| [#180](https://github.com/docxology/MetaInformAnt/issues/180) | sync-test artifact | not planned |
| [#181](https://github.com/docxology/MetaInformAnt/issues/181) | sync-test artifact | not planned |
| [#182](https://github.com/docxology/MetaInformAnt/issues/182) | sync-test artifact | not planned |
| [#185](https://github.com/docxology/MetaInformAnt/issues/185) | sync-test artifact | not planned |
| [#187](https://github.com/docxology/MetaInformAnt/issues/187) | sync-test artifact | not planned |
| [#188](https://github.com/docxology/MetaInformAnt/issues/188) | sync-test artifact | not planned |
