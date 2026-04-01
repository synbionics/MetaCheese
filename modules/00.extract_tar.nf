nextflow.enable.dsl=2

process extract_tar {
    tag "$tar_file"

    input:
    path tar_file

    output:
    path("*.fq.gz"), emit: fqgz

    publishDir "${params.output_base}/00_Extract", mode: 'copy'

    script:
    """
    set -euo pipefail
    outdir=\$(basename ${tar_file} .tar)
    mkdir \$outdir
    tar -xf ${tar_file} -C \$outdir
    cp \$outdir/*.fq.gz .
    """
}
