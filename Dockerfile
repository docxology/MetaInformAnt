# ─────────────────────────────────────────────────────────────────────────
# MetaInformAnt — Amalgkit RNA-seq Pipeline Container
#
# Installs all bioinformatics tools and Python deps needed to run the
# streaming RNA-seq pipeline (download → quantify → merge → curate).
#
# Build:  docker build -t metainformant-pipeline .
# Run:    docker run -v $(pwd)/output:/app/output metainformant-pipeline
# ─────────────────────────────────────────────────────────────────────────

FROM python:3.13-slim AS base

LABEL maintainer="docxology"
LABEL description="MetaInformAnt amalgkit RNA-seq pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

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
    && Rscript -e 'install.packages(c("ggplot2"), repos="https://cloud.r-project.org", quiet=TRUE)'

# ── kallisto 0.48.0 ─────────────────────────────────────────────────────
RUN cd /tmp && \
    wget -q https://github.com/pachterlab/kallisto/releases/download/v0.48.0/kallisto_linux-v0.48.0.tar.gz && \
    tar xzf kallisto_linux-v0.48.0.tar.gz && \
    cp kallisto/kallisto /usr/local/bin/ && \
    rm -rf /tmp/kallisto*

# ── fastp 0.24.0 ────────────────────────────────────────────────────────
RUN cd /tmp && \
    wget -q http://opengene.org/fastp/fastp && \
    chmod +x fastp && mv fastp /usr/local/bin/

# ── SRA Toolkit (fasterq-dump) ──────────────────────────────────────────
RUN cd /tmp && \
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz && \
    tar xzf sratoolkit.3.2.1-ubuntu64.tar.gz && \
    cp sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump /usr/local/bin/ && \
    cp sratoolkit.3.2.1-ubuntu64/bin/prefetch /usr/local/bin/ && \
    cp sratoolkit.3.2.1-ubuntu64/bin/vdb-config /usr/local/bin/ && \
    rm -rf /tmp/sratoolkit*

# ── UV for fast Python dependency resolution ────────────────────────────
RUN pip install uv

# ── Application code ────────────────────────────────────────────────────
WORKDIR /app
COPY pyproject.toml ./
COPY src/ ./src/
COPY scripts/ ./scripts/
COPY config/ ./config/

# Install Python deps
RUN uv pip install --system -e "." && \
    cd /tmp && \
    wget -qO micromamba.tar.bz2 "https://micro.mamba.pm/api/micromamba/linux-64/latest" && \
    tar -xjf micromamba.tar.bz2 bin/micromamba && \
    export MAMBA_ROOT_PREFIX=/opt/conda && \
    bin/micromamba create -y -p /opt/conda -c conda-forge -c bioconda amalgkit && \
    ln -s /opt/conda/bin/amalgkit /usr/local/bin/amalgkit && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf /tmp/micromamba* bin

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
