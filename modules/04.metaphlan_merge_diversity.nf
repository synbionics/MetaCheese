process metaphlan_merge_diversity {
    tag "metaphlan_merge"

    input:
    path(txts)

    output:
    path("merged_abundance_table.txt")
    path("diversity/*")

    publishDir "${params.output_base}/04_metaphlan_output", mode: 'copy'

    script:
    """
    set -euo pipefail
    echo "DEBUG: Files in workdir:"; ls -lh
    echo "DEBUG: TXT files:"; ls -lh *.txt || true
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    N=\$(ls *.txt 2>/dev/null | wc -l)
    echo "DEBUG: numero file .txt = \$N"

    if [[ "\$N" -eq 1 ]]; then
        cp *.txt merged_abundance_table.txt
    elif [[ "\$N" -gt 1 ]]; then
        merge_metaphlan_tables.py *.txt > merged_abundance_table.txt
    else
        echo "No txt files found!" >&2
        exit 1
    fi

    # Analisi alpha diversity
    for metric in richness shannon simpson gini; do
        Rscript /main/data/calculate_diversity.R \\
            -f merged_abundance_table.txt \\
            -d alpha \\
            -m \$metric \\
            -o diversity
    done

    ## Analisi beta diversity
    #for metric in bray-curtis jaccard clr aitchison; do
    #    Rscript /main/data/calculate_diversity.R \\
    #        -f merged_abundance_table.txt \\
    #        -d beta \\
    #        -m \$metric \\
    #        -o diversity
    #done
    #
    ## UniFrac solo se tree.nwk esiste (opzionale, togli se vuoi semplificare)
    #if [ -f tree.nwk ]; then
    #    for metric in weighted-unifrac unweighted-unifrac; do
    #        Rscript /main/data/calculate_diversity.R \\
    #            -f merged_abundance_table.txt \\
    #            -d beta \\
    #            -m \$metric \\
    #            -t tree.nwk \\
    #            -o diversity
    #    done
    #fi

    conda deactivate
    """
}
