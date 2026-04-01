process humann_sample {
    tag "$sample"

    input:
        tuple val(sample), path(read1), path(taxonomic_profile)

    output:
        path("${sample}_humann")                // puoi cambiare secondo le tue esigenze

    publishDir "${params.output_base}/04b_humann_output", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate humannenv

    NUC_DB="/main/db/humann/chocophlan/chocophlan"
    PROT_DB="/main/db/humann/uniref/uniref"

    # Check che i db ci siano
    if [ ! -d "\$NUC_DB" ] || [ -z "\$(ls -A "\$NUC_DB")" ]; then
        echo "ERRORE: DB ChocoPhlAn mancante" >&2; exit 1
    fi
    if [ ! -d "\$PROT_DB" ] || [ -z "\$(ls -A "\$PROT_DB")" ]; then
        echo "ERRORE: DB UniRef90 mancante" >&2; exit 1
    fi

    # Output in cartella per sample
    mkdir -p ${sample}_humann

    humann \\
      --input "$read1" \\
      --output ${sample}_humann \\
      --threads 4 \\
      --remove-temp-output \\
      --taxonomic-profile "$taxonomic_profile" \\
      --nucleotide-database "\$NUC_DB" \\
      --protein-database "\$PROT_DB"

    conda deactivate
    """
}
