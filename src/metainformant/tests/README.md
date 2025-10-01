# Tests Source Directory

This directory contains test implementations for METAINFORMANT's internal functionality and integration testing.

## Overview

The tests module provides comprehensive testing for core functionality, module integration, and system validation.

## Files

### Test Infrastructure
- **`__init__.py`**: Test module initialization
- **`runner.py`**: Test execution and reporting utilities

## Usage in Development

These files support the testing infrastructure:
- Test discovery and execution
- Test result reporting and analysis
- Integration testing coordination
- Performance benchmarking

## Integration

Tests module integrates with:
- **pytest** for test framework
- **Core utilities** for I/O and configuration
- **All domain modules** for functionality testing
- **CI/CD systems** for automated testing

## Testing Philosophy

Follows METAINFORMANT's NO_MOCKING_POLICY:
- Real implementations for all functionality
- Real external API calls where applicable
- Comprehensive error condition testing
- Performance and scalability validation

## Contributing

When adding new test functionality:
1. Follow established testing patterns
2. Include comprehensive edge case testing
3. Ensure compatibility with CI/CD workflows
4. Update test documentation and coverage

## Related Documentation

- See `tests/README_tests.md` for comprehensive test suite documentation
- See `docs/testing.md` for testing philosophy and guidelines
- See `scripts/run_tests.sh` for test execution tooling

This module provides essential testing infrastructure for METAINFORMANT development.
