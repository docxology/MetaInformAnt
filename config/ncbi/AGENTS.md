# Agent Directives: config/ncbi

## Role

NCBI API and data retrieval configuration (email, API keys, rate limiting).

## Contents

| File | Description |
| :--- | :--- |
| `ncbi.yaml` | NCBI API credentials and download settings |

## Configuration Structure

```yaml
# NCBI configuration
email: user@example.com
api_key: null  # Optional, increases rate limit
max_retries: 3
timeout: 30
```

## Rules

- Validate with schema before committing new configs
- Follow NO MOCKING policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
