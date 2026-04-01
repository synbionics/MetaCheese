process mapp_coverage {
    tag { sample }

    input:
    tuple val(sample), path(fq1), path(fq2), val(idx_prefix)

    output:
      tuple val(sample),
            path("${sample}.bam"),
            path("${sample}.sorted.bam"),
            path("${sample}.sorted.bam.bai"),
            path("${sample}.sorted.bam.idxstat"),
            path("${sample}_bowtie2.log")


    publishDir "${params.output_base}/08-09_mapping_coverage", mode: 'copy'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bowtieenv

    echo "Using index: ${idx_prefix}"
    
    bowtie2 -x "${idx_prefix}" \
      -1 "${fq1}" -2 "${fq2}" \
      -q --no-unal --very-sensitive-local -p ${task.cpus} \
      2> "${sample}_bowtie2.log" \
      | samtools view -bS -o "${sample}.bam" -

    samtools sort "${sample}.bam" -o "${sample}.sorted.bam"
    samtools index "${sample}.sorted.bam"
    samtools idxstats "${sample}.sorted.bam" > "${sample}.sorted.bam.idxstat"

    conda deactivate
    """
}
