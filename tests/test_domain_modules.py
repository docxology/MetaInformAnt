from pathlib import Path

from metainformant.ontology.go import count_go_scripts
from metainformant.phenotype.data.antwiki import load_antwiki_json
from metainformant.protein.sequence.proteomes import read_taxon_ids
from metainformant.rna.engine.pipeline import summarize_curate_tables


def test_rna_curate_summary_apis_mellifera():
    repo_root = Path(__file__).resolve().parents[1]
    curate_dir = repo_root / "tests/data/rna/curate/Apis_mellifera"
    counts = summarize_curate_tables(curate_dir)
    assert any(name.endswith("metadata.tsv") for name in counts)
    assert any(name.endswith("tc.tsv") for name in counts)


def test_protein_taxon_ids_readable():
    repo_root = Path(__file__).resolve().parents[1]
    path = repo_root / "tests/data/protein/taxon_id_list.txt"
    ids = read_taxon_ids(path)
    assert len(ids) >= 5
    assert all(isinstance(x, int) for x in ids)


def test_phenotype_antwiki_json_loads():
    repo_root = Path(__file__).resolve().parents[1]
    path = repo_root / "tests/data/phenotype/antwiki_dataset_sorted_final_01.json"
    entries = load_antwiki_json(path)
    assert isinstance(entries, list)
    assert len(entries) >= 1


def test_ontology_go_dir_counts():
    repo_root = Path(__file__).resolve().parents[1]
    go_dir = repo_root / "tests/data/ontology/GO_v3"
    n = count_go_scripts(go_dir)
    assert n >= 4
