nextflow.enable.dsl=2

process bowtie_mag_database {
    tag "$sample"

    input:
    tuple val(sample), path(fasta_file)

    output:
    tuple val(sample), path("${sample}_index_base.1.bt2"), path("${sample}_index_base.2.bt2"), path("${sample}_index_base.3.bt2"), path("${sample}_index_base.4.bt2"), path("${sample}_index_base.rev.1.bt2"), path("${sample}_index_base.rev.2.bt2"), emit: db

    publishDir "${params.output_base}/07_Bowtie_Index", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bowtieenv

    bowtie2-build ${fasta_file} ${sample}_index_base -p ${task.cpus}

    conda deactivate
    """
}
