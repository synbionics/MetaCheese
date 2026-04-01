nextflow.enable.dsl=2

process bowtie2_map {
    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.unmapped.1.fq.gz"), path("${sample}.unmapped.2.fq.gz"), emit: unmapped

    publishDir "${params.output_base}/03_bowtie2_output", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bowtieenv

    bowtie2 -x ${params.bowtie2_db} \
      -1 ${read1} -2 ${read2} \
      --un-conc-gz ${sample}.unmapped.fq.gz \
      -p ${task.cpus} \
      -S ${sample}.sam

    mv ${sample}.unmapped.fq.1.gz ${sample}.unmapped.1.fq.gz
    mv ${sample}.unmapped.fq.2.gz ${sample}.unmapped.2.fq.gz

    rm -f ${sample}.sam

    conda deactivate
    """
}
