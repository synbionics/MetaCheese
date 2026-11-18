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

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.3.1-1-Linux-x86_64.sh -O /tmp/miniconda.sh && \
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

RUN conda create -y -n humannenv python=3.7
RUN conda run -n humannenv conda install -y -c bioconda humann=3.9
RUN conda run -n humannenv pip install humann
RUN conda run -n humannenv humann --version
RUN conda run -n humannenv humann_databases --help
RUN conda run -n humannenv humann_databases --download chocophlan full /main/db/humann/chocophlan
RUN conda run -n humannenv humann_databases --download uniref uniref90_diamond /main/db/humann/uniref

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
# QUALITY CONTROL OF RAW READS (FastQC, step 00 / diagram Q.1)
# Own environment so it cannot change the versions pinned in bioenv.
# Placed after the databases so adding it does not invalidate their layers.
# ======================================

RUN conda create -y -n qcenv -c conda-forge -c bioconda fastqc=0.12.1 && \
    conda clean -afy

# ======================================
# HUMAnN utility mapping files (step 04b-last / diagram C.3)
# The default install only maps UniRef90 -> MetaCyc reactions. The full set
# adds UniRef90 -> KEGG Orthology, EC, GO, Pfam and eggNOG, needed to regroup
# gene families into KOs/ECs. --update-config makes HUMAnN find them.
# ======================================

RUN mkdir -p /main/db/humann/utility_mapping && \
    conda run -n humannenv humann_databases --download utility_mapping full \
        /main/db/humann/utility_mapping --update-config yes

# ======================================
# MultiQC (step 00 / diagram Q.1): merges the per-file FastQC reports into
# one page. Own environment, and its own layer placed after the databases so
# adding it does not invalidate them.
# ======================================

RUN conda create -y -n multiqcenv -c conda-forge -c bioconda multiqc=1.35 && \
    conda clean -afy

# ======================================
# QUAST (step 06b / diagram A.8): quality of the filtered assemblies. Own
# environment and layer, after the databases. Run without a reference, so it
# needs no download at run time.
# ======================================

RUN conda create -y -n quastenv -c conda-forge -c bioconda python=3.10 quast=5.3.0 && \
    conda clean -afy

# ======================================
# Kraken2 database of TORMES (step 15 / diagram A.7: MAG taxonomy).
# tormes-setup installs minikraken_8GB_20200312, which uses the Lactobacillus
# names from before the 2020 split. Replace it with the Standard-8 build of
# 2026-10-22 (same 8 GB size class, RefSeq archaea/bacteria/viral/plasmid/
# human/UniVec_Core), from https://benlangmead.github.io/aws-indexes/k2.
# The build's published checksum file (standard_08_GB.md5) checks the
# download and then the extracted files. Own layer after the databases, so the
# layers above stay as they are.
# ======================================

ARG KRAKEN2_DB=k2_standard_08_GB_20261022
ARG KRAKEN2_MD5=standard_08_GB_20261022/standard_08_GB.md5
RUN set -e; db=/opt/conda/envs/tormes-1.3.0/db; mkdir -p "$db/$KRAKEN2_DB" && cd "$db/$KRAKEN2_DB" && \
    wget -q "https://genome-idx.s3.amazonaws.com/kraken/$KRAKEN2_MD5" -O checksums.md5 && \
    wget -q "https://genome-idx.s3.amazonaws.com/kraken/$KRAKEN2_DB.tar.gz" && \
    grep " $KRAKEN2_DB.tar.gz\$" checksums.md5 | md5sum -c --quiet && \
    tar -xzf "$KRAKEN2_DB.tar.gz" && rm "$KRAKEN2_DB.tar.gz" && \
    test -s hash.k2d -a -s opts.k2d -a -s taxo.k2d && \
    grep -E ' (hash|opts|taxo)\.k2d$' checksums.md5 | md5sum -c --quiet && \
    rm -rf "$db/minikraken-DB" && \
    sed -i "s|^KRAKEN-DATABASE\t.*|KRAKEN-DATABASE\t$db/$KRAKEN2_DB|" /opt/conda/envs/tormes-1.3.0/files/config_file.txt && \
    grep -q "^KRAKEN-DATABASE	$db/$KRAKEN2_DB$" /opt/conda/envs/tormes-1.3.0/files/config_file.txt

# ======================================
# MetaPhlAn 4.1.1 + database mpa_vJun23_CHOCOPhlAnSGB_202403 (step 04 /
# diagram C.1). Replaces 4.1.0 + vJun23_202307. Stays on vJun23: HUMAnN 3.9
# maps vJun23 SGBs to ChocoPhlAn (humann3_sgb_map.Jun23.tsv) and cannot use
# vJan25. --freeze-installed keeps the other bioenv packages (SPAdes,
# AdapterRemoval, MetaBAT2, R) as they are; the diff of the package list is
# printed in the build log to check that. Own layer after the databases, so
# the layers above stay as they are.
# ======================================

ARG MPA_VERSION=4.1.1
ARG MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403
RUN set -e; conda list -n bioenv --export > /tmp/bioenv.before && \
    conda install -y -n bioenv --freeze-installed "metaphlan=$MPA_VERSION" && \
    conda list -n bioenv --export > /tmp/bioenv.after && \
    echo "bioenv package changes:" && (diff /tmp/bioenv.before /tmp/bioenv.after || true) && \
    conda run -n bioenv metaphlan --version 2>&1 | grep -q "version $MPA_VERSION" && \
    rm -rf /main/db/metaphlan/* && \
    conda run -n bioenv metaphlan --install --index "$MPA_INDEX" --bowtie2db /main/db/metaphlan && \
    test -s "/main/db/metaphlan/$MPA_INDEX.pkl" && \
    rm -f /tmp/bioenv.before /tmp/bioenv.after && conda clean -afy

# ======================================
# CheckM2 1.1.0 + DIAMOND database v3 (step 12-13 / diagram A.5). Replaces
# CheckM2 1.0.2 + database v2. checkm2env is rebuilt: the 1.1.0 build for
# Python 3.8 (pyh7e72e81_0) is broken (TensorFlow 2.4 cannot load its
# .keras model), so build pyh7e72e81_1 is used (Python 3.12, TensorFlow 2.17,
# DIAMOND 2.1.11). The base image has glibc 2.27 and conda does not report
# it, so CONDA_OVERRIDE_GLIBC states it. `database --download` takes the
# newest database this CheckM2 accepts (v3, zenodo 14897628) and records its
# path; `testrun` checks the install on CheckM2's own test genomes.
# ======================================

ARG CHECKM2_SPEC=checkm2=1.1.0=pyh7e72e81_1
RUN set -e; conda env remove -y -n checkm2env && \
    CONDA_OVERRIDE_GLIBC=2.27 conda create -y -n checkm2env -c conda-forge -c bioconda "$CHECKM2_SPEC" && \
    conda run -n checkm2env checkm2 --version | grep -qx "1.1.0" && \
    echo "checkm2env packages:" && conda list -n checkm2env | grep -E '^(python|tensorflow|keras|diamond|prodigal|lightgbm|scikit-learn|numpy|pandas) ' && \
    rm -rf /main/db/checkm2/* && \
    conda run -n checkm2env checkm2 database --download --path /main/db/checkm2 && \
    test -s /main/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd && \
    conda run -n checkm2env checkm2 testrun --threads 4 && \
    conda clean -afy

# ======================================
# Bowtie2 2.5.5 + samtools 1.25 (steps 03, 07, 08-09 and the host index).
# Replaces Bowtie2 2.3.3.1 + samtools 1.6 in bowtieenv (Python 2.7, no
# longer needed). MetaPhlAn, HUMAnN and TORMES use their own copies in their
# own envs and are not changed. Bowtie2 2.x indexes stay valid, so the host
# index in data/gene is not rebuilt. The last command builds a tiny index
# and maps reads to it, to check the install.
# ======================================

ARG BOWTIE2_VERSION=2.5.5
ARG SAMTOOLS_VERSION=1.25
RUN set -e; conda env remove -y -n bowtieenv && \
    CONDA_OVERRIDE_GLIBC=2.27 conda create -y -n bowtieenv -c conda-forge -c bioconda \
        "bowtie2=$BOWTIE2_VERSION" "samtools=$SAMTOOLS_VERSION" && \
    conda run -n bowtieenv bowtie2 --version | head -1 | grep -q "version $BOWTIE2_VERSION$" && \
    conda run -n bowtieenv samtools --version | head -1 | grep -qx "samtools $SAMTOOLS_VERSION" && \
    echo "bowtieenv packages:" && conda list -n bowtieenv | grep -E '^(python|perl|bowtie2|samtools|htslib|tbb) ' && \
    mkdir /tmp/bt && cd /tmp/bt && \
    python3 -c "import random; random.seed(1); s=''.join(random.choice('ACGT') for _ in range(5000)); open('r.fa','w').write('>r\n'+s+'\n'); open('q.fq','w').write(''.join('@q%d\n%s\n+\n%s\n'%(i,s[i*40:i*40+100],'I'*100) for i in range(100)))" && \
    conda run -n bowtieenv bowtie2-build -q r.fa r && \
    conda run -n bowtieenv bash -c "bowtie2 -x r -U q.fq 2>/dev/null | samtools view -c -F 4 -" | grep -qx 100 && \
    cd / && rm -rf /tmp/bt && conda clean -afy

# ======================================
# metaSPAdes 4.3.0 (step 05 / diagram A.1). Replaces SPAdes 3.13.0 of bioenv,
# which stays installed but is no longer used: SPAdes 4 gets its own env
# (spadesenv) so the other bioenv tools are not touched, and step 05
# activates spadesenv. `spades.py --test` assembles SPAdes's own toy data to
# check the install.
# ======================================

ARG SPADES_VERSION=4.3.0
RUN set -e; CONDA_OVERRIDE_GLIBC=2.27 conda create -y -n spadesenv -c conda-forge -c bioconda \
        "spades=$SPADES_VERSION" && \
    conda run -n spadesenv spades.py --version | grep -qx "SPAdes genome assembler v$SPADES_VERSION" && \
    echo "spadesenv packages:" && conda list -n spadesenv | grep -E '^(python|spades) ' && \
    cd /tmp && conda run -n spadesenv spades.py --test > spades_test.log 2>&1 && \
    grep -q "TEST PASSED CORRECTLY" spades_test.log && \
    rm -rf /tmp/spades_test /tmp/spades_test.log && cd / && conda clean -afy

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
