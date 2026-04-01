process metabat_binning {
    tag "${sample}"

    input:
        tuple val(sample), path(contig), path(depth)

    output:
        tuple val(sample), path("${sample}.*.fa") // <-- emetti la cartella dei MAG

    publishDir "${params.output_base}/11b_metabat_MAG/", mode: 'copy'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv


    metabat2 -i "$contig" -a "$depth" -o "${sample}" -m ${params.metabat_binning.min_contig_length} -t "${task.cpus}"

    conda deactivate
    """
}