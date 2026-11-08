#!/usr/bin/env bash
# Print the tool and database versions found inside the image, one per line
# (name<TAB>version). Fed to the container on stdin by run.sh.
source /opt/conda/etc/profile.d/conda.sh

tool() {  # name, conda environment, command
    local out
    conda activate "$2" 2>/dev/null
    out="$(eval "$3" 2>&1 | grep -m1 -E '[0-9]+\.[0-9]+' | tr -s ' ' | cut -c1-200)"
    conda deactivate 2>/dev/null
    printf '%s\t%s\n' "$1" "${out:-unknown}"
}
db() {  # name, folder: its file names and total size identify the build
    printf '%s\t%s (%s bytes)\n' "$1" \
        "$(find "$2" -maxdepth 2 -type f -printf '%f\n' | sort | head -1)" \
        "$(find "$2" -type f -printf '%s\n' | awk '{n+=$1} END {printf "%.0f", n}')"
}

tool FastQC qcenv "fastqc --version"
tool MultiQC multiqcenv "multiqc --version"
tool AdapterRemoval bioenv "AdapterRemoval --version"
tool Bowtie2 bowtieenv "bowtie2 --version"
tool samtools bowtieenv "samtools --version"
tool MetaPhlAn bioenv "metaphlan --version"
tool HUMAnN humannenv "humann --version"
tool DIAMOND_humann humannenv "diamond --version"
tool metaSPAdes bioenv "metaspades.py --version"
tool QUAST quastenv "quast.py --version"
tool MetaBAT2 bioenv "metabat2 2>&1 | head -3"
tool CheckM2 checkm2env "checkm2 --version"
tool DIAMOND_checkm2 checkm2env "diamond --version"
tool TORMES tormes-1.3.0 "tormes --version"
tool Kraken2 tormes-1.3.0 "kraken2 --version"

db db_MetaPhlAn /main/db/metaphlan
db db_ChocoPhlAn /main/db/humann/chocophlan
db db_UniRef /main/db/humann/uniref
db db_CheckM2 /main/db/checkm2
db db_Kraken2 /opt/conda/envs/tormes-1.3.0/db/minikraken-DB
