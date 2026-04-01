process metaphlan_sample {
    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.bowtie2.bz2"), path("${sample}.sam.bz2"), path("${sample}.txt"), emit: metaphlan_out

    publishDir "${params.output_base}/04_metaphlan_output", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    metaphlan \\
        "$read1","$read2" \\
        --input_type fastq \\
        --bowtie2db "/main/db/metaphlan" \\
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \\
        --nproc 4 \\
        --bowtie2out "${sample}.bowtie2.bz2" \\
        -s "${sample}.sam.bz2" \\
        -o "${sample}.txt"

    conda deactivate
    """
}
