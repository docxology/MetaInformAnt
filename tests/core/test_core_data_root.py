"""Tests for metainformant.core.io.data_root (canonical data-root resolution)."""

from __future__ import annotations

from metainformant.core.io.data_root import (
    AMALGKIT_DATA_ROOT_ENV,
    remap_amalgkit_prefix,
    resolve_data_root,
)


class TestResolveDataRoot:
    def test_explicit_arg_wins_over_environ(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, "/env/root")
        assert resolve_data_root(tmp_path) == tmp_path.resolve()

    def test_environ_used_when_no_arg(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, str(tmp_path))
        assert resolve_data_root() == tmp_path.resolve()

    def test_default_fallback_without_environ(self, monkeypatch) -> None:
        monkeypatch.delenv(AMALGKIT_DATA_ROOT_ENV, raising=False)
        resolved = resolve_data_root()
        assert resolved.name == "amalgkit"
        assert resolved.is_absolute()

    def test_environ_injection_parameter(self, tmp_path) -> None:
        # No monkeypatch: pass a mapping directly (non-process context).
        assert resolve_data_root(environ={AMALGKIT_DATA_ROOT_ENV: str(tmp_path)}) == tmp_path.resolve()

    def test_tilde_expansion(self, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, "~/amalgkit-data")
        resolved = resolve_data_root()
        assert not str(resolved).startswith("~")
        assert resolved.is_absolute()

    def test_relative_env_value_resolves_to_absolute(self, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, "relative/root")
        resolved = resolve_data_root()
        assert resolved.is_absolute()
        assert resolved.parts[-2:] == ("relative", "root")


class TestRemapAmalgkitPrefix:
    def test_exact_prefix_maps_to_root(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, str(tmp_path))
        assert remap_amalgkit_prefix("output/amalgkit") == str(tmp_path.resolve())

    def test_prefixed_path_maps_under_root(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, str(tmp_path))
        result = remap_amalgkit_prefix("output/amalgkit/Apis_mellifera/quant")
        assert result == str(tmp_path.resolve() / "Apis_mellifera" / "quant")

    def test_dot_slash_prefix_also_maps(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, str(tmp_path))
        result = remap_amalgkit_prefix("./output/amalgkit/species/work")
        assert result == str(tmp_path.resolve() / "species" / "work")

    def test_non_prefixed_path_untouched(self, monkeypatch) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, "/data/root")
        assert remap_amalgkit_prefix("config/config_base/settings.yaml") == "config/config_base/settings.yaml"

    def test_explicit_data_root_overrides_environ(self, monkeypatch, tmp_path) -> None:
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, "/env/root")
        result = remap_amalgkit_prefix("output/amalgkit/x", tmp_path)
        assert result == str(tmp_path.resolve() / "x")

    def test_similar_prefix_not_partial_matched(self, monkeypatch, tmp_path) -> None:
        # "output/amalgkit2/..." must NOT be remapped.
        monkeypatch.setenv(AMALGKIT_DATA_ROOT_ENV, str(tmp_path))
        assert remap_amalgkit_prefix("output/amalgkit2/x") == "output/amalgkit2/x"
