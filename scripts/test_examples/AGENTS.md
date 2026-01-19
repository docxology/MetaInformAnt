# AGENTS: Test Examples Runner

AI-assisted development of the example test runner for MetaInformAnt.

## AI Contributions

The **Code Assistant Agent** developed:
- `ExampleTester` class for orchestrating discovery and execution.
- Parallel execution using `concurrent.futures`.
- Multi-format reporting (JSON, HTML, JUnit XML).
- Graceful error handling and summary statistics.

## Function Index

| Function/Class | Signature |
|----------------|-----------|
| `ExampleTester` | `(verbose, domain_filter, continue_on_error, parallel, max_workers)` |
| `ExampleTester.run_all_examples` | `() -> dict[str, Any]` |
| `discover_examples` | `(examples_dir, domain_filter) -> list[Path]` |
| `run_examples_sequential` | `(examples, repo_root, verbose) -> list[dict]` |
| `run_examples_parallel` | `(examples, repo_root, verbose, max_workers) -> list[dict]` |
| `generate_html_report` | `(results, output_dir) -> None` |
| `generate_junit_xml` | `(results, output_dir) -> None` |
| `save_json_report` | `(results, output_dir) -> None` |
| `print_summary` | `(results) -> None` |

## Related Documentation

- **README**: [README.md](README.md) - Usage guide.
- **SPEC**: [SPEC.md](SPEC.md) - Technical specification.
- **Parent**: [scripts/AGENTS.md](../AGENTS.md)
