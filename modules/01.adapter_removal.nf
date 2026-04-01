nextflow.enable.dsl=2

process adapter_removal {
    tag "$base"

    input:
    tuple val(base), path(reads)  // reads è una lista: [fq1, fq2]

    output:
    tuple val(base), path("${base}_cleaned_L1.fq.gz"), path("${base}_cleaned_L2.fq.gz"), emit: cleaned

    publishDir "${params.output_base}/01_AdapterRemoval", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} \
      --output1 ${base}_cleaned_L1.fq.gz --output2 ${base}_cleaned_L2.fq.gz \
      --threads ${task.cpus} --gzip --minlength ${params.adapterRemoval.minlength} --trimqualities --minquality ${params.adapterRemoval.minquality} --trimns --maxns ${params.adapterRemoval.maxns} --trim5p ${params.adapterRemoval.trim5p} --trim3p ${params.adapterRemoval.trim3p}

    conda deactivate
    """
}
