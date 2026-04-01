process metabat_depth {
    tag "${sample}"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${sample}.depth.txt")

    publishDir "${params.output_base}/10_metabat_depth", mode: 'copy'
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    jgi_summarize_bam_contig_depths --outputDepth ${sample}.depth.txt $bam

    conda deactivate
    """
}
