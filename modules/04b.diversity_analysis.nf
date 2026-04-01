process merge_diversity {
    publishDir "${params.output_base}/04b_Metaphlan", mode: 'copy'

    input:
    tuple val(dummy), path(profile_txts)
    // dummy è solo per trigger, puoi passare sempre 'all' o il nome di un sample

    output:
    path("merged_abundance_table.txt")
    path("diversity/")

    script:
    """
    merge_metaphlan_tables.py ${profile_txts.join(' ')} > merged_abundance_table.txt

    mkdir -p diversity

    for metric in richness shannon simpson gini; do
        echo "Calcolo alpha diversity (\$metric)..."
        Rscript ../../data/calculate_diversity.R \\
            -f merged_abundance_table.txt \\
            -d alpha \\
            -m "\$metric" \\
            -o diversity
    done

    for metric in bray-curtis jaccard clr aitchison; do
        echo "Calcolo beta diversity (\$metric)..."
        Rscript ../../data/calculate_diversity.R \\
            -f merged_abundance_table.txt \\
            -d beta \\
            -m "\$metric" \\
            -o diversity
    done

    if [[ -f tree.nwk ]]; then
        for metric in weighted-unifrac unweighted-unifrac; do
            echo "Calcolo beta diversity (\$metric)..."
            Rscript ../../data/calculate_diversity.R \\
                -f merged_abundance_table.txt \\
                -d beta \\
                -m "\$metric" \\
                -t tree.nwk \\
                -o diversity
        done
    else
        echo "Salto metriche UniFrac (weighted/unweighted) per mancanza di tree.nwk."
    fi
    """
}
