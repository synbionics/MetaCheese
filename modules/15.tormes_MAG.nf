process tormes_MAG {
    tag "tormes"

    input:
        path metadata

    output:
        path "tormes_output"

    publishDir "${params.output_base}/", mode: 'copy'

    script:
    """
    set +u

    source /opt/conda/etc/profile.d/conda.sh
    conda activate tormes-1.3.0

    tormes --metadata \$metadata --output tormes_output --threads "${task.cpus}"

    conda deactivate
    """
}
