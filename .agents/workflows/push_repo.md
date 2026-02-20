---
description: How to push the MetaInformAnt repository securely
---
# Prepare and Push Repository

This workflow outlines the precise steps for pushing the `metainformant` repository, ensuring all tests pass, the environment is clean, and documentation is updated.

1. **Verify Test Suite**
   Ensure that the zero-mock test suite passes locally.

   ```bash
   uv run pytest tests/ -v
   ```

2. **Check Git Status**
   Ensure working tree is clean and there are no lingering ephemeral files in `output/` tracked by git, to avoid large unneeded blobs.

   ```bash
   git status
   ```

3. **Stage Changes**

   ```bash
   git add .
   ```

4. **Commit with Conventions**
   Use substantive commit messages, referring to the module or feature.

   ```bash
   git commit -m "chore: Finalize documentation for robust cross-machine movement"
   ```

5. **Push to Remote**

   ```bash
   git push origin main
   ```
