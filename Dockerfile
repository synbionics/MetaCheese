# ======================================
# MetaCheese image
#
# Every conda environment is created from docker/envs/<env>.txt: the exact
# package list (URL + MD5) of the image that passed the regression test
# (tests/regression, 2026-11), so a new build installs the same packages.
# Databases are downloaded at fixed versions. To change a tool, change its
# list (or the version below), rebuild and run the regression test.
#
# Tools and databases (diagram ids as in the README):
#   bioenv       AdapterRemoval 2.3.3, MetaPhlAn 4.1.1, MetaBAT2 2.17, R
#                (SPAdes 3.13.0 is also in it, not used)
#   humannenv    HUMAnN 3.9
#   bowtieenv    Bowtie2 2.5.5, samtools 1.25
#   spadesenv    metaSPAdes 4.3.0
#   checkm2env   CheckM2 1.1.0 (checkmenv: CheckM 1.2.2, installed, not run)
#   tormes-1.3.0 TORMES 1.3.0 + Kraken2 Standard-8 2026-10-22
#   qcenv        FastQC 0.12.1;  multiqcenv MultiQC 1.35;  quastenv QUAST 5.3.0
# ======================================

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# ======================================
# System packages
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

RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# R from CRAN (bionic-cran40)
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /etc/apt/trusted.gpg.d/cran-r.gpg && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y r-base && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ======================================
# Miniconda (conda 25.3.1, kept as installed: the environments come from
# explicit package lists, so no solver runs)
# ======================================

ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.3.1-1-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -afy

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# yq (YAML parser used by the pipeline scripts)
RUN wget -qO /usr/local/bin/yq \
    https://github.com/mikefarah/yq/releases/download/v4.34.1/yq_linux_amd64 && \
    chmod +x /usr/local/bin/yq

# pandoc (report rendering)
RUN wget https://github.com/jgm/pandoc/releases/download/2.19.2/pandoc-2.19.2-1-amd64.deb && \
    dpkg -i pandoc-2.19.2-1-amd64.deb && \
    rm pandoc-2.19.2-1-amd64.deb

# ======================================
# Conda environments, from the exact package lists
# ======================================

COPY docker/envs /tmp/envs
RUN set -e; for env in bioenv humannenv bowtieenv spadesenv checkmenv checkm2env \
                       tormes-1.3.0 qcenv multiqcenv quastenv; do \
        conda create -y -n "$env" --file "/tmp/envs/$env.txt"; \
    done && \
    rm -rf /tmp/envs && conda clean -afy

# HUMAnN 3.9 is used from pip, on Python 3.7 (the conda package of the list
# goes to a python3.12 folder that this env does not run)
RUN conda run -n humannenv pip install --no-deps --no-cache-dir humann==3.9 && \
    conda run -n humannenv humann --version 2>&1 | grep -qx "humann v3.9"

RUN ln -s /opt/conda/envs/bioenv/bin/metaphlan /usr/local/bin/metaphlan && \
    ln -s /opt/conda/envs/bioenv/bin/merge_metaphlan_tables.py /usr/local/bin/merge_metaphlan_tables.py

# Install checks: versions, and each tool's own test where it has one
RUN set -e; \
    conda run -n bioenv metaphlan --version 2>&1 | grep -q "version 4.1.1" && \
    conda run -n checkm2env checkm2 --version | grep -qx "1.1.0" && \
    conda run -n bowtieenv bowtie2 --version | head -1 | grep -q "version 2.5.5$" && \
    conda run -n bowtieenv samtools --version | head -1 | grep -qx "samtools 1.25" && \
    conda run -n spadesenv spades.py --version | grep -qx "SPAdes genome assembler v4.3.0" && \
    mkdir /tmp/bt && cd /tmp/bt && \
    python3 -c "import random; random.seed(1); s=''.join(random.choice('ACGT') for _ in range(5000)); open('r.fa','w').write('>r\n'+s+'\n'); open('q.fq','w').write(''.join('@q%d\n%s\n+\n%s\n'%(i,s[i*40:i*40+100],'I'*100) for i in range(100)))" && \
    conda run -n bowtieenv bowtie2-build -q r.fa r && \
    conda run -n bowtieenv bash -c "bowtie2 -x r -U q.fq 2>/dev/null | samtools view -c -F 4 -" | grep -qx 100 && \
    cd /tmp && rm -rf /tmp/bt && \
    conda run -n spadesenv spades.py --test > spades_test.log 2>&1 && \
    grep -q "TEST PASSED CORRECTLY" spades_test.log && \
    rm -rf /tmp/spades_test /tmp/spades_test.log

# ======================================
# Databases, at fixed versions
# ======================================

# HUMAnN 3.9 (step 04b / C.2-C.3): ChocoPhlAn, UniRef90 DIAMOND and the full
# utility mapping (UniRef90 -> KO, EC, GO, Pfam, eggNOG, MetaCyc), with the
# paths written to HUMAnN's config
RUN set -e; mkdir -p /main/db/humann/chocophlan /main/db/humann/uniref /main/db/humann/utility_mapping && \
    conda run -n humannenv humann_databases --download chocophlan full /main/db/humann/chocophlan --update-config yes && \
    conda run -n humannenv humann_databases --download uniref uniref90_diamond /main/db/humann/uniref --update-config yes && \
    conda run -n humannenv humann_databases --download utility_mapping full /main/db/humann/utility_mapping --update-config yes

# MetaPhlAn 4.1.1 (step 04 / C.1): vJun23_202403. Stays on vJun23: HUMAnN 3.9
# maps vJun23 SGBs to ChocoPhlAn and cannot use vJan25.
ARG MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403
RUN set -e; mkdir -p /main/db/metaphlan && \
    conda run -n bioenv metaphlan --install --index "$MPA_INDEX" --bowtie2db /main/db/metaphlan && \
    test -s "/main/db/metaphlan/$MPA_INDEX.pkl"

# CheckM2 1.1.0 (step 12-13 / A.5): DIAMOND database v3 (zenodo 14897628,
# the newest this CheckM2 accepts); `testrun` checks the install
RUN set -e; mkdir -p /main/db/checkm2 && \
    conda run -n checkm2env checkm2 database --download --path /main/db/checkm2 && \
    test -s /main/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd && \
    conda run -n checkm2env checkm2 testrun --threads 4

# TORMES 1.3.0 (step 15 / A.6-A.7): pinned replacement of tormes-setup (see
# docker/tormes/setup.sh)
COPY docker/tormes /tmp/tormes
RUN bash /tmp/tormes/setup.sh && rm -rf /tmp/tormes

# Kraken2 database of TORMES (A.7: MAG taxonomy): Standard-8 build of
# 2026-10-22 (RefSeq archaea/bacteria/viral/plasmid/human/UniVec_Core), from
# https://benlangmead.github.io/aws-indexes/k2, checked against the build's
# published MD5 file. TORMES's config_file.txt points to it.
ARG KRAKEN2_DB=k2_standard_08_GB_20261022
ARG KRAKEN2_MD5=standard_08_GB_20261022/standard_08_GB.md5
RUN set -e; db=/opt/conda/envs/tormes-1.3.0/db; mkdir -p "$db/$KRAKEN2_DB" && cd "$db/$KRAKEN2_DB" && \
    wget -q "https://genome-idx.s3.amazonaws.com/kraken/$KRAKEN2_MD5" -O checksums.md5 && \
    wget -q "https://genome-idx.s3.amazonaws.com/kraken/$KRAKEN2_DB.tar.gz" && \
    grep " $KRAKEN2_DB.tar.gz\$" checksums.md5 | md5sum -c --quiet && \
    tar -xzf "$KRAKEN2_DB.tar.gz" && rm "$KRAKEN2_DB.tar.gz" && \
    test -s hash.k2d -a -s opts.k2d -a -s taxo.k2d && \
    grep -E ' (hash|opts|taxo)\.k2d$' checksums.md5 | md5sum -c --quiet && \
    grep -q "^KRAKEN-DATABASE	$db/$KRAKEN2_DB$" /opt/conda/envs/tormes-1.3.0/files/config_file.txt

# ======================================
# Project files (after the databases, so changing them does not invalidate
# the database layers)
# ======================================

COPY scripts /main/scripts
COPY config /main/config
COPY data /main/data

WORKDIR /main
CMD ["/bin/bash"]
