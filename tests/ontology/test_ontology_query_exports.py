"""The ``metainformant.ontology.query`` package re-exports what its docstring says it re-exports.

``query`` used to be a single module; when it became a package, ``__init__`` documented that it re-exports every
public symbol of the ``query`` and ``serialize`` modules "for backward compatibility" but only bound the two submodules,
so ``from metainformant.ontology.query import ancestors`` raised ImportError.
"""

from __future__ import annotations

import inspect
from types import ModuleType

import pytest

import metainformant.ontology.query as query_package
from metainformant.ontology.query import query as query_module, serialize as serialize_module


def _public_functions(module: ModuleType) -> dict[str, object]:
    return {
        name: obj
        for name, obj in vars(module).items()
        if inspect.isfunction(obj) and obj.__module__ == module.__name__ and not name.startswith("_")
    }


@pytest.mark.parametrize("module", [query_module, serialize_module], ids=lambda m: m.__name__)
def test_every_public_function_is_reexported_and_listed(module: ModuleType) -> None:
    functions = _public_functions(module)
    assert functions, "the module under test should define public functions"
    for name, func in functions.items():
        assert getattr(query_package, name, None) is func, f"{name} is not re-exported"
        assert name in query_package.__all__, f"{name} is missing from __all__"


def test_all_names_resolve() -> None:
    assert len(query_package.__all__) == len(set(query_package.__all__)), "__all__ has duplicates"
    for name in query_package.__all__:
        assert hasattr(query_package, name), f"__all__ lists {name!r}, which the package does not define"


def test_named_imports_work() -> None:
    from metainformant.ontology.query import ancestors, descendants, load_ontology, save_ontology  # noqa: F401

    assert ancestors is query_module.ancestors
    assert save_ontology is serialize_module.save_ontology


def test_star_import_exposes_the_documented_names() -> None:
    namespace: dict[str, object] = {}
    exec("from metainformant.ontology.query import *", namespace)  # noqa: S102 - exercising __all__ itself
    for name in ("ancestors", "descendants", "common_ancestors", "save_ontology", "load_ontology"):
        assert name in namespace


def test_submodules_stay_importable_from_the_package() -> None:
    assert query_package.query is query_module
    assert query_package.serialize is serialize_module
    assert {"query", "serialize"} <= set(query_package.__all__)
