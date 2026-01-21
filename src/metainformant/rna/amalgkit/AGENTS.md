# Amalgkit Agents API

## Purpose

This document describes the agentic interfaces and capabilities exposed by the `amalgkit` module. This module serves as the primary interface for the Amalgkit RNA-seq toolkit within the MetaInformAnt ecosystem.

## Integration

* **Tooling**: Functions in this module are designed to be used as tools by high-level agents (like the RNA Workflow Engine).
* **Observability**: Agents can monitor the execution and state of components in this module via return codes and process outputs.

## Capabilities

### 1. RNA-seq Workflow Components

Agents can inspect and execute the following atomic functions, which correspond to Amalgkit CLI subcommands:

* **`metadata()`**: Fetches and curates metadata for selected species.
* **`config()`**: helper to generate configuration files.
* **`select()`**: Filters samples based on metadata criteria.
* **`getfastq()`**: **Robust** SRA downloader. Handles retries, network glitches, and parallel execution. Supports AWS, NCBI, and PFD backends.
* **`integrate()`**: Verifies and prepares FASTQ files for quantification.
* **`quant()`**: Quantifies gene expression using Kallisto (or Salmon). Handles index building automatically.
* **`merge()`**: Aggregates quantification results into gene count/TPM matrices.
* **`cstmm()`**: Normalizes expression data across species/samples.
* **`curate()`**: Finalizes dataset structure and metadata.
* **`csca()`**: Performs Cross-Species Correlational Analysis (QC).
* **`sanity()`**: Checks the integrity of the output directory.

### 2. State & Parameter Management

* **`AmalgkitParams`**: A standardized object for passing configuration to any tool in this module. Ensures type safety and consistent parameter handling.

### 3. Execution Engine

* **`run_amalgkit()`**: The core driver that builds CLI commands, manages subprocesses, and handles logging/monitoring (optional).
* **`ensure_cli_available()`**: Self-healing capability to check for and optionally install the underlying `amalgkit` tool via UV.
