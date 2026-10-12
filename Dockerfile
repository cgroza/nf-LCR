FROM --platform=linux/amd64 ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MINIMAP2_VERSION=2.31
ARG SAMTOOLS_VERSION=1.21
ARG IMPG_VERSION=main

ENV PATH="/root/.cargo/bin:${PATH}"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        autoconf \
        automake \
        bash \
        build-essential \
        ca-certificates \
        cmake \
        coreutils \
        curl \
        gawk \
        git \
        libbz2-dev \
        libclang-dev \
        libcurl4-openssl-dev \
        libdeflate-dev \
        liblzma-dev \
        libncurses-dev \
        libssl-dev \
        libtool \
        openjdk-21-jre-headless \
        pkg-config \
        r-base-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 --branch "v${MINIMAP2_VERSION}" \
        https://github.com/lh3/minimap2.git /tmp/minimap2 \
    && make --directory /tmp/minimap2 --jobs "$(nproc)" \
    && install -m 0755 /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && rm -rf /tmp/minimap2

RUN git clone --depth 1 --branch "${SAMTOOLS_VERSION}" --recurse-submodules \
        https://github.com/samtools/samtools.git /tmp/samtools \
    && cd /tmp/samtools \
    && autoheader \
    && autoconf -Wno-syntax \
    && ./configure --prefix=/usr/local \
    && make --jobs "$(nproc)" \
    && make install \
    && rm -rf /tmp/samtools

RUN curl --fail --silent --show-error --location https://sh.rustup.rs | \
        sh -s -- -y --profile minimal

RUN git clone --depth 1 --branch "${IMPG_VERSION}" --recurse-submodules \
        https://github.com/pangenome/impg.git /tmp/impg \
    && cd /tmp/impg \
    && cargo install --path . --locked \
    && rm -rf /tmp/impg

RUN Rscript -e 'install.packages("ggplot2", repos = "https://cloud.r-project.org")'

RUN curl --fail --silent --show-error --location https://get.nextflow.io | bash \
    && install -m 0755 nextflow /usr/local/bin/nextflow \
    && rm nextflow

WORKDIR /opt/chr22q

COPY \
    main.nf \
    nextflow.config \
    run_chr22q_lcr_analysis.sh \
    prepare_chr22q_haplotypes.py \
    extract_chr22q_protein_coding_genes.py \
    plot_chr22q_full_lcr_layout.R \
    plot_chr22q_full_lcr_gene_layout.R \
    plot_chr22q_lcr_gene_configurations.R \
    /opt/chr22q/

WORKDIR /work

CMD ["nextflow", "run", "/opt/chr22q/main.nf", "-profile", "slurm"]
