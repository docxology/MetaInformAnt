# Validation Scripts

Validation scripts run repository-level audits and signposting checks.

## Commands

```bash
uv run python scripts/validate/project_completeness_audit.py --root .
```

Add new validation scripts here when they inspect repository structure rather
than package runtime behavior.

## Output

Write reports under `output/` or print summaries to stdout. Do not place audit
cache files under `src/`.
