import os


def test_expected_package_layout_exists():
    repo_root = os.path.dirname(os.path.dirname(__file__))
    src_dir = os.path.join(repo_root, "src")
    assert os.path.isdir(src_dir), "Missing src/ directory"

    pkg_dir = os.path.join(src_dir, "metainformant")
    assert os.path.isdir(pkg_dir), "Missing src/metainformant package"

    for subpkg in [
        "core",
        "dna",
        "rna",
        "protein",
        "epigenome",
        "ontology",
        "phenotype",
        "ecology",
        "visualization",
    ]:
        assert os.path.isdir(os.path.join(pkg_dir, subpkg)), f"Missing subpackage: {subpkg}"


def test_pyproject_exists():
    repo_root = os.path.dirname(os.path.dirname(__file__))
    assert os.path.isfile(os.path.join(repo_root, "pyproject.toml"))


