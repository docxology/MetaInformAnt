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
- Follow same REAL IMPLEMENTATION policy as all other tests
- Use real amalgkit commands when `@pytest.mark.external_tool` is appropriate


## Operator-local qualification prerequisites

The EHD runtime/admission, HD checkpoint, local CLI/slot/storage-profile and
campaign cold-start cases qualify ignored campaign operators rather than
portable package code. They carry the external_tool marker. Their default
operator source is output/hymenoptera_local_20261019; set
RNA_OPERATOR_SOURCE_DIR to an explicitly supplied qualification bundle
containing the scripts each case requests. This binding does not copy or
version production operators.

EHD cases additionally require RNA_EHD_DATA_NAMESPACE,
RNA_EHD_DATA_UUID and UV_PROJECT_ENVIRONMENT to be declared by the
qualification owner; their values, mounted-root identity, sandbox and native
libraries remain subject to the real checks. Cross-volume checkpoint cases
require AMALGKIT_STAGING_TEST_BULK_ROOT. A missing prerequisite is an explicit
skip. A present malformed asset, invalid binding value, import error or failed
subprocess remains a failure, never a reason to skip after execution.

The package acquisition-worker cold-start, lazy public-module identity and
unknown-attribute cases remain portable and unmarked. A portable test result
does not qualify the unavailable operator/native cases; report their skips
and retain earlier failed qualification receipts separately.
