"""MetaInformAnt's Amalgkit option tables match the installed (pinned) parser exactly.

The tables in ``metainformant.rna.amalgkit.commands`` are hand-maintained. These
tests enumerate the real argparse parser of the installed Amalgkit build, so a
pin bump that renames, removes or adds an option fails here instead of turning
into an argparse error (or a silently ignored setting) during a campaign.
"""

from __future__ import annotations

import argparse
import logging
from uuid import uuid4

import pytest

from metainformant.rna.amalgkit.amalgkit import REQUIRED_AMALGKIT_VERSION, build_cli_args
from metainformant.rna.amalgkit.commands import (
    CURRENT_AMALGKIT_STEPS,
    METAINFORMANT_ONLY_STEP_KEYS,
    SUPPORTED_CLI_OPTIONS,
    YES_NO_CLI_OPTIONS,
)

amalgkit_main = pytest.importorskip("amalgkit.main")


@pytest.fixture(scope="module")
def subparsers() -> dict[str, argparse.ArgumentParser]:
    parser = amalgkit_main.build_main_parser()
    (action,) = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    return {name: subparser for name, subparser in action.choices.items() if name != "help"}


def _options(subparser: argparse.ArgumentParser) -> list[argparse.Action]:
    """Public options of one command: no help action, no hidden ``obsolete_*`` aliases."""

    return [
        action
        for action in subparser._actions
        if not isinstance(action, argparse._HelpAction) and not action.dest.startswith("obsolete_")
    ]


def test_command_registry_matches_installed_commands(subparsers: dict[str, argparse.ArgumentParser]) -> None:
    assert set(subparsers) == set(CURRENT_AMALGKIT_STEPS) == set(SUPPORTED_CLI_OPTIONS) == set(YES_NO_CLI_OPTIONS)


def test_allow_list_equals_installed_parser_options(subparsers: dict[str, argparse.ArgumentParser]) -> None:
    missing_upstream: dict[str, list[str]] = {}
    not_allow_listed: dict[str, list[str]] = {}
    for command, subparser in subparsers.items():
        parsed = {action.dest for action in _options(subparser)}
        allowed = SUPPORTED_CLI_OPTIONS[command]
        if allowed - parsed:
            missing_upstream[command] = sorted(allowed - parsed)
        if parsed - allowed:
            not_allow_listed[command] = sorted(parsed - allowed)
    # An allow-listed option the pinned parser lacks would be an argparse error at run time.
    assert missing_upstream == {}
    # An upstream option missing from the allow-list would be dropped from configurations.
    assert not_allow_listed == {}


def test_every_allow_listed_option_is_spelled_as_its_dest(subparsers: dict[str, argparse.ArgumentParser]) -> None:
    for command, subparser in subparsers.items():
        for action in _options(subparser):
            assert f"--{action.dest}" in action.option_strings, (command, action.dest, action.option_strings)


def test_yes_no_table_matches_parser_value_types(subparsers: dict[str, argparse.ArgumentParser]) -> None:
    for command, subparser in subparsers.items():
        yes_no = {action.dest for action in _options(subparser) if getattr(action.type, "__name__", "") == "strtobool"}
        assert YES_NO_CLI_OPTIONS[command] == yes_no, command


def test_boolean_settings_build_arguments_the_parser_accepts(subparsers: dict[str, argparse.ArgumentParser]) -> None:
    for command, subparser in subparsers.items():
        # Required options (busco --lineage) get a placeholder from their choices.
        required = {
            action.dest: str(next(iter(action.choices))) if action.choices else "placeholder"
            for action in _options(subparser)
            if action.required
        }
        for action in _options(subparser):
            is_switch = isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction))
            if action.dest not in YES_NO_CLI_OPTIONS[command] and not is_switch:
                continue
            for value in (True, False):
                args = build_cli_args({"out_dir": "work", **required, action.dest: value}, subcommand=command)
                namespace = subparser.parse_args(args)
                if action.dest in YES_NO_CLI_OPTIONS[command]:
                    assert args[args.index(f"--{action.dest}") + 1] == ("yes" if value else "no")
                    assert bool(getattr(namespace, action.dest)) is value, (command, action.dest)
                else:
                    assert (f"--{action.dest}" in args) is value, (command, action.dest)


def test_unknown_keys_are_reported_and_known_orchestration_keys_are_not(caplog: pytest.LogCaptureFixture) -> None:
    unknown = f"typo_{uuid4().hex[:8]}"
    params = {
        "out_dir": "work",
        "metadata": "work/metadata/metadata.tsv",
        unknown: "1",
        "tissue": ["brain"],
        "show_progress": True,
        "norm": "log2p1-fpkm",  # a finalize/wsfilter option in a shared map
    }
    assert {"tissue", "show_progress"} <= METAINFORMANT_ONLY_STEP_KEYS

    with caplog.at_level(logging.DEBUG, logger="metainformant"):
        args = build_cli_args(params, subcommand="quant")
        build_cli_args(params, subcommand="quant")

    assert "--metadata" in args
    assert not any(token.lstrip("-") in {unknown, "tissue", "show_progress", "norm"} for token in args)
    warnings = [record for record in caplog.records if record.levelno == logging.WARNING]
    # Reported once per (command, key), however often the map is rebuilt.
    assert [record.getMessage() for record in warnings if unknown in record.getMessage()] == [
        f"amalgkit quant: dropping configuration key '{unknown}': "
        f"no command of the pinned Amalgkit {REQUIRED_AMALGKIT_VERSION} accepts it"
    ]
    for expected in ("tissue", "show_progress", "norm"):
        assert not [record for record in warnings if repr(expected) in record.getMessage()], expected
