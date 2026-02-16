# RNA Test Suite

Tests for the RNA-seq pipeline engine, retrieval, and amalgkit integration.

## 📦 Test Categories

| Area | What's Tested |
|------|---------------|
| Streaming Orchestrator | Multi-species pipeline, parallel downloads, timeout handling |
| ENA Downloader | FTP/HTTP download, MD5 verification, size-based skip |
| Amalgkit Wrapper | CLI command construction, config loading, step execution |
| Workflow Engine | Config parsing, step sequencing, progress tracking |

## 🚀 Running Tests

```bash
# Run all RNA tests
uv run pytest tests/rna/ -v

# Run specific test file
uv run pytest tests/rna/test_ena_downloader.py -v

# Run with coverage
uv run pytest tests/rna/ --cov=metainformant.rna -v
```

## 🔗 Related

- [src/metainformant/rna/](../../src/metainformant/rna/) - Source module
- [scripts/rna/](../../scripts/rna/) - Pipeline scripts
- [config/amalgkit/](../../config/amalgkit/) - Workflow configurations
