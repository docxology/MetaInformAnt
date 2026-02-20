---
description: Running the zero-mock test suite
---
# Run the Test Suite

Our tests mandate a strict Zero-Mock policy, testing real I/O and outputs.

1. **Run Full Test Suite**

   ```bash
   uv run pytest tests/ -v
   ```

2. **Run Module-Specific Tests**

   ```bash
   uv run pytest tests/test_dna_sequence.py -v
   ```

3. **Run with Coverage**

   ```bash
   uv run pytest tests/ --cov=src/metainformant
   ```
