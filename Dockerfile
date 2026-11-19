# syntax=docker/dockerfile:1.7

FROM cgroza/nf-chr22q@sha256:a2bfa5aa13a27746899722a60d862d2468acf38ff7e873509de1e8bd5a96e1cb AS impg-builder

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MINIMAP2_VERSION=2.31
ARG HTSLIB_VERSION=1.21
ARG SAMTOOLS_VERSION=1.21
ARG NEXTFLOW_VERSION=26.10.0
ARG WFMASH_VERSION=v0.24.2

ENV NXF_HOME=/opt/nextflow

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
        libgsl-dev \
        libjemalloc-dev \
        liblzma-dev \
        libncurses-dev \
        libssl-dev \
        libtool \
        openjdk-21-jre-headless \
        pkg-config \
        r-base-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then export GIT_SSL_CAINFO=/run/secrets/mac_ca; fi \
    && git clone --depth 1 --branch "v${MINIMAP2_VERSION}" \
        https://github.com/lh3/minimap2.git /tmp/minimap2 \
    && make --directory /tmp/minimap2 --jobs "$(nproc)" \
    && install -m 0755 /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && rm -rf /tmp/minimap2

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then export GIT_SSL_CAINFO=/run/secrets/mac_ca; fi \
    && git clone --depth 1 --branch "${HTSLIB_VERSION}" --recurse-submodules \
        https://github.com/samtools/htslib.git /tmp/htslib \
    && cd /tmp/htslib \
    && git submodule update --init --recursive \
    && autoreconf --install \
    && ./configure --prefix=/usr/local \
    && make --jobs "$(nproc)" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/htslib

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then export GIT_SSL_CAINFO=/run/secrets/mac_ca; fi \
    && git clone --depth 1 --branch "${SAMTOOLS_VERSION}" \
        https://github.com/samtools/samtools.git /tmp/samtools \
    && cd /tmp/samtools \
    && autoheader \
    && autoconf -Wno-syntax \
    && ./configure --prefix=/usr/local --with-htslib=/usr/local \
    && make --jobs "$(nproc)" \
    && make install \
    && rm -rf /tmp/samtools

COPY --from=impg-builder /usr/local/bin/impg /usr/local/bin/impg
COPY --from=impg-builder /usr/local/lib/libwfa2*.so* /usr/local/lib/

RUN ldconfig

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then export GIT_SSL_CAINFO=/run/secrets/mac_ca; fi \
    && git clone --depth 1 --branch "${WFMASH_VERSION}" --recurse-submodules \
        https://github.com/waveygang/wfmash.git /tmp/wfmash \
    && cmake -S /tmp/wfmash -B /tmp/wfmash/build -DCMAKE_BUILD_TYPE=Generic -DDISABLE_LTO=ON \
        -DCMAKE_INSTALL_RPATH=/opt/wfmash/lib -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    && cmake --build /tmp/wfmash/build --parallel "$(nproc)" \
    && install -m 0755 /tmp/wfmash/build/bin/wfmash /usr/local/bin/wfmash \
    && mkdir -p /opt/wfmash/lib \
    && cp -a /tmp/wfmash/build/lib/libwflign.so* /tmp/wfmash/build/lib/libwfa2*.so* /opt/wfmash/lib/ \
    && /usr/local/bin/wfmash --version \
    && rm -rf /tmp/wfmash

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then export CURL_CA_BUNDLE=/run/secrets/mac_ca; fi \
    && Rscript -e 'install.packages("ggplot2", repos = "https://cloud.r-project.org")'

RUN --mount=type=secret,id=mac_ca,required=false \
    if [ -s /run/secrets/mac_ca ]; then \
        curl --cacert /run/secrets/mac_ca --fail --silent --show-error --location \
            "https://www.nextflow.io/releases/v${NEXTFLOW_VERSION}/nextflow?a=install" \
            --output /usr/local/bin/nextflow; \
        chmod 0755 /usr/local/bin/nextflow; \
        cp /usr/local/bin/nextflow /tmp/nextflow-bootstrap; \
        sed -i 's/curl /curl --cacert \/run\/secrets\/mac_ca /g' /tmp/nextflow-bootstrap; \
        chmod 0755 /tmp/nextflow-bootstrap; \
        /tmp/nextflow-bootstrap -version; \
        rm /tmp/nextflow-bootstrap; \
    else \
        curl --fail --silent --show-error --location \
            "https://www.nextflow.io/releases/v${NEXTFLOW_VERSION}/nextflow?a=install" \
            --output /usr/local/bin/nextflow; \
        chmod 0755 /usr/local/bin/nextflow; \
        nextflow -version; \
    fi

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
