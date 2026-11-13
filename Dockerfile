# ─────────────────────────────────────────────────────────────────────────
# MetaInformAnt — Amalgkit RNA-seq Pipeline Container
#
# Installs all bioinformatics tools and Python deps needed to run the
# streaming RNA-seq pipeline (download → quantify → merge → finalize).
#
# Python dependencies, including Amalgkit, come only from uv.lock. Amalgkit
# is pinned by commit (upstream tags and bioconda only publish x.y.0
# releases), so it is never taken from bioconda. Native tools are pinned to
# the versions that produced the existing quantifications: kallisto 0.52.0
# (index_version 13) and sra-tools 3.4.1.
#
# Build:  docker build -t metainformant-pipeline .
# Run:    docker run -v $(pwd)/output:/app/output metainformant-pipeline
# ─────────────────────────────────────────────────────────────────────────

FROM python:3.14-slim AS base

LABEL maintainer="docxology"
LABEL description="MetaInformAnt amalgkit RNA-seq pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    KALLISTO_VERSION=0.52.0 \
    SRA_TOOLKIT_VERSION=3.4.1

# ── System dependencies ─────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    curl \
    git \
    wget \
    bzip2 \
    pigz \
    zlib1g-dev \
    libhdf5-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ── R 4.5 ───────────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/* \
    && Rscript -e 'install.packages(c("BiocManager", "amap", "RColorBrewer", "colorspace", "dendextend", "NMF", "MASS", "pvclust", "Rtsne", "ggplot2", "patchwork", "optparse", "reshape2", "gridExtra"), repos="https://cloud.r-project.org", quiet=TRUE)'

# ── SRA Toolkit (fasterq-dump) ──────────────────────────────────────────
# The SRA toolkit uses a launcher pattern: fasterq-dump is a wrapper that
# execs fasterq-dump.VERSION. We install it to /opt and add to PATH to preserve its environment.
RUN cd /tmp && \
    wget -q "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_TOOLKIT_VERSION}/sratoolkit.${SRA_TOOLKIT_VERSION}-ubuntu64.tar.gz" && \
    tar xzf "sratoolkit.${SRA_TOOLKIT_VERSION}-ubuntu64.tar.gz" && \
    mv "sratoolkit.${SRA_TOOLKIT_VERSION}-ubuntu64" /opt/sratoolkit && \
    rm -rf /tmp/sratoolkit* && \
    mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n/libs/cloud/report_instance_identity = "false"\n/repository/user/main/public/root = "/tmp/sra-cache"\n' "$(uuidgen 2>/dev/null || echo 'docker-container')" > /root/.ncbi/user-settings.mkfg
ENV PATH="/opt/sratoolkit/bin:${PATH}"

# ── uv (pinned binary) for lock-exact Python environments ──────────────
COPY --from=ghcr.io/astral-sh/uv:0.13.7 /uv /uvx /usr/local/bin/

# ── Native tools from bioconda (never Amalgkit) ────────────────────────
# Ensure all conda-installed tools (kallisto, seqkit, fastp) are on PATH.
ENV PATH="/opt/conda/bin:${PATH}"
RUN cd /tmp && \
    wget -qO micromamba.tar.bz2 "https://micro.mamba.pm/api/micromamba/linux-64/latest" && \
    tar -xjf micromamba.tar.bz2 bin/micromamba && \
    export MAMBA_ROOT_PREFIX=/opt/conda && \
    bin/micromamba create -y -p /opt/conda -c conda-forge -c bioconda "kallisto=${KALLISTO_VERSION}" seqkit fastp && \
    rm -f /opt/conda/bin/kallisto_orig && \
    ln -s /opt/conda/bin/seqkit /usr/local/bin/seqkit && \
    ln -s /opt/conda/bin/fastp /usr/local/bin/fastp && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf /tmp/micromamba* bin && \
    kallisto version | grep -q "${KALLISTO_VERSION}"

# ── Application code and the locked Python environment ────────────────
WORKDIR /app
COPY pyproject.toml uv.lock README.md ./
COPY src/ ./src/
COPY scripts/ ./scripts/
COPY config/ ./config/

# Exactly the versions in uv.lock (Amalgkit at its pinned commit), using the
# image's Python rather than a downloaded interpreter. MetaInformAnt never
# installs packages at runtime, so this is the only install step.
ENV UV_PROJECT_ENVIRONMENT=/app/.venv \
    UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=/usr/local/bin/python3 \
    UV_CACHE_DIR=/tmp/uv-cache \
    UV_LINK_MODE=copy
RUN uv sync --frozen --no-dev --extra rna && \
    rm -rf /tmp/uv-cache && \
    /app/.venv/bin/python scripts/package/verify_amalgkit_runtime.py
ENV PATH="/app/.venv/bin:${PATH}"

# Pre-download NCBI taxdump for ete4 to prevent pipeline hanging during metadata initialization
RUN mkdir -p /root/.local/share/ete && \
    wget -qO /root/.local/share/ete/taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# ── Output volume ───────────────────────────────────────────────────────
RUN mkdir -p /app/output/amalgkit
VOLUME ["/app/output"]

# ── Healthcheck ─────────────────────────────────────────────────────────
HEALTHCHECK --interval=60s --timeout=10s --retries=3 \
    CMD python3 -c "import metainformant; print('ok')" || exit 1

# ── Entrypoint ──────────────────────────────────────────────────────────
# Default: run the full species pipeline
# Override workers/threads/max-gb via environment variables
ENV PIPELINE_MAX_GB=20.0 \
    PIPELINE_WORKERS=80 \
    PIPELINE_THREADS=96

CMD python3 scripts/rna/run_all_species.py \
    --max-gb ${PIPELINE_MAX_GB} \
    --workers ${PIPELINE_WORKERS} \
    --threads ${PIPELINE_THREADS}
