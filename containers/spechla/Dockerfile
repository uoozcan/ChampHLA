FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Umut Ozcan <u.onurozcan@gmail.com>"
LABEL description="SpecHLA v1.1 — HLA typing with variant calling support"
LABEL version="1.1"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        autoconf \
        automake \
        libdeflate-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN git clone --depth 1 https://github.com/deepomicslab/SpecHLA.git && \
    cd SpecHLA && \
    bash install.sh

ENV PATH="/app/SpecHLA/spechla_env/bin:/app/SpecHLA/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/SpecHLA/spechla_env/lib:${LD_LIBRARY_PATH}"

RUN cd /app/SpecHLA && \
    python -c "import sys; print('Python:', sys.version)" && \
    samtools --version | head -1

WORKDIR /data
