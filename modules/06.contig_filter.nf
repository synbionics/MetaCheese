nextflow.enable.dsl=2

process contig_filter {
    tag "$sample"

    input:
    tuple val(sample), path(spades_dir)

    output:
    tuple val(sample), path("filtered/${sample}.fasta_sort.fasta"), emit: filtered

    publishDir "${params.output_base}/06_contig_filter", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    mkdir -p filtered

    # Copia il contig nel formato sample.fasta
    cp ${spades_dir}/contigs.fasta ${sample}.fasta

    # Filtra i contigs >500bp in filtered/
    awk 'BEGIN{RS=\">\"; ORS=\"\"} length(\$0) > ${params.contigFilter.Length} {print \">\"\$0}' ${sample}.fasta > filtered/${sample}.fasta_sort.fasta

    conda deactivate
    """
}
