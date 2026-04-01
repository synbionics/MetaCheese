# ======================================
# BASE IMAGE E AMBIENTE
# ======================================

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# ======================================
# PACCHETTI DI SISTEMA
# ======================================

RUN apt-get update && \
    apt-get install -y \
        wget \
        curl \
        git \
        build-essential \
        python \
        python-dev \
        openjdk-11-jre-headless \
        unzip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        locales \
        perl \
        gnupg2 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ======================================
# LOCALE UTF-8
# ======================================

RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# ======================================
# R (DA CRAN UFFICIALE)
# ======================================

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /etc/apt/trusted.gpg.d/cran-r.gpg && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y r-base && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ======================================
# INSTALLA MINICONDA E CONFIGURA CANALI
# ======================================

ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -afy

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update -n base -c defaults conda -y

# ======================================
# TOOL YQ (YAML parser CLI)
# ======================================

RUN wget -qO /usr/local/bin/yq \
    https://github.com/mikefarah/yq/releases/download/v4.34.1/yq_linux_amd64 && \
    chmod +x /usr/local/bin/yq

# ======================================
# AMBIENTI CONDA PER TOOL BIOINFO
# ======================================

# Ambiente bioinfo generico con molti tool
RUN conda create -y -n bioenv python=3.8 && \
    conda run -n bioenv conda install -y \
        r-vegan \
        r-optparse \
        r-ape \
        r-rbiom \
        r-compositions \
        bioconductor-microbiome \
        adapterremoval=2.3.3 \
        metaphlan=4.1.0 \
        spades=3.13.0 \
        metabat2=2.17

# Link simbolici per Nextflow
RUN ln -s /opt/conda/envs/bioenv/bin/metaphlan /usr/local/bin/metaphlan
RUN ln -s /opt/conda/envs/bioenv/bin/merge_metaphlan_tables.py /usr/local/bin/merge_metaphlan_tables.py

# ======================================
# AMBIENTE HUMAnN + DATABASE (workaround bioconda/pip bug)
# ======================================
# Prepara le cartelle database HUMAnN
RUN mkdir -p /main/db/humann/chocophlan && \
    mkdir -p /main/db/humann/uniref

RUN conda create -y -n humannenv python=3.7 && \
    conda run -n humannenv bash -c "\
      conda install -y -c bioconda humann=3.9 && \
      pip install humann && \
      humann_databases --help && \
      humann_databases --download chocophlan full /main/db/humann/chocophlan && \
      humann_databases --download uniref uniref90_diamond /main/db/humann/uniref \
    "

# ======================================
# ALTRI AMBIENTI CONDA
# ======================================

# Bowtie2 + Samtools legacy
RUN conda create -y -n bowtieenv python=2.7 bowtie2=2.3.3.1 samtools=1.6

# CheckM legacy
RUN conda create -y -n checkmenv python=3.7 && \
    conda run -n checkmenv conda install -y -c conda-forge -c bioconda checkm-genome=1.2.2

# CheckM2
RUN conda create -y -n checkm2env python=3.8 && \
    conda run -n checkm2env conda install -y -c bioconda checkm2=1.0.2 && \
    conda run -n checkm2env pip install numpy==1.23.5 'pandas>=1.3,<2' --no-cache-dir && \
    conda run -n checkm2env pip install --upgrade checkm2 && \
    conda run -n checkm2env conda install -y -c bioconda prodigal diamond

# TORMES
RUN wget https://raw.githubusercontent.com/biobrad/Tormes-Meta-Create/main/tormes-1.3.0.yml && \
    CONDA_CHANNEL_PRIORITY=flexible conda env create -n tormes-1.3.0 --file tormes-1.3.0.yml && \
    rm tormes-1.3.0.yml

# ======================================
# INSTALLA PANDOC
# ======================================

RUN wget https://github.com/jgm/pandoc/releases/download/2.19.2/pandoc-2.19.2-1-amd64.deb && \
    dpkg -i pandoc-2.19.2-1-amd64.deb && \
    rm pandoc-2.19.2-1-amd64.deb

# ======================================
# SCARICA TUTTI I DATABASE BIOINFO RESTANTI
# ======================================

RUN mkdir -p /main/db/metaphlan && \
    mkdir -p /main/db/checkm2

# Metaphlan DB
RUN conda run -n bioenv metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /main/db/metaphlan

# CheckM2 database
RUN conda run -n checkm2env checkm2 database --download --path /main/db/checkm2

# TORMES setup (DB vengono gestiti dal tool interno)
RUN conda run -n tormes-1.3.0 tormes-setup

# ======================================
# COPIA LE CARTELLE DEL TUO PROGETTO (DOPO i DB, così non li sovrascrivi MAI!)
# ======================================

COPY scripts /main/scripts
COPY config /main/config
COPY data /main/data

# ======================================
# WORKDIR E DEFAULT SHELL
# ======================================

WORKDIR /main
CMD ["/bin/bash"]
