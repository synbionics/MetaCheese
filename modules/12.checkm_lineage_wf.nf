process checkm_lineage_wf {
    tag "${sample}"

    input:
        tuple val(sample), path(mag_dir)

    output:
        tuple val(sample), path("output")
        tuple val(sample), path("${sample}.checkm.txt")

    publishDir "${params.output_base}/12_checkm_lineage_wf", mode: 'copy'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate checkmenv

    mkdir input_MAGs
    for f in ${mag_dir}; do cp "\$f" input_MAGs/; done

    mkdir output

    checkm lineage_wf \
        -t "${task.cpus}" -x fa input_MAGs output \
        --pplacer_threads ${params.checkm.pplacer_threads} \
        -f "${sample}.checkm.txt" || { echo "Errore in CheckM"; exit 1; }

    conda deactivate
    """
}
