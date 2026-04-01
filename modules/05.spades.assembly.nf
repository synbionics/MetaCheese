nextflow.enable.dsl=2

process spades_assembly {
    tag "$sample"

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    tuple val(sample), path("${sample}"), emit: assembly   // <--- Passa tutta la cartella di output!

    publishDir "${params.output_base}/05_spades_output", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    spades.py -1 ${fq1} -2 ${fq2} --meta -t ${task.cpus} --memory ${params.spades.memory} --only-assembler -o ${sample}

    conda deactivate
    """
}
