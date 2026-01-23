# Agent Directives: tests/rna

## Role
RNA module-specific test organization for complex amalgkit workflow testing.

## Purpose
Separates RNA-specific test infrastructure from main test directory when tests require:
- Complex fixtures specific to RNA workflows
- Multi-step workflow validation
- Integration with amalgkit CLI tool

## Rules
- Tests here extend the main `test_rna_*.py` tests in parent directory
- Follow same NO MOCKING policy as all other tests
- Use real amalgkit commands when `@pytest.mark.external_tool` is appropriate
