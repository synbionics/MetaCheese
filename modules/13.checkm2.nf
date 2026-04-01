nextflow.enable.dsl=2

process checkm2_predict {

    tag "$sample"

    input:
    tuple val(sample), path(mag_files)


    output:
    tuple val(sample), path("*.tsv")


    publishDir "${params.output_base}/13_checkm2/", mode: 'copy'

    script:
    """
    set -e

    # Salva sample in variabile
    sample="${sample}"

    source /opt/conda/etc/profile.d/conda.sh
    conda activate checkm2env

    mkdir -p all_MAGs_\$sample

    # Trasforma lista Nextflow in array bash e copia file
    mag_files=(${mag_files})
    cp \${mag_files[@]} all_MAGs_\$sample/

    checkm2 predict \
        --threads ${task.cpus} \
        --input all_MAGs_$sample \
        --output-directory checkm2 \
        -x .fa

    mv checkm2/quality_report.tsv "${sample}.checkm2.tsv"
    """


}
