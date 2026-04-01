process humann_merge {
    tag "humann_merge"

    input:
        path(humann_out_dir)

    output:
        path("humann_joined_genefamilies.tsv")
        path("humann_joined_genefamilies_cpm.tsv")
        path("merged_table_genefamilies_rxn.txt")
        path("merged_table_ko.txt")
        path("merged_table_kegg_pathway.txt")

    publishDir "${params.output_base}/04b_humann_output/merge", mode: 'copy'

    script:
    """
    set -euo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate humannenv

    INPUT_DIR="\$PWD"

    # Controlla che ci siano file *_genefamilies.tsv
    if ! ls *_genefamilies.tsv 1>/dev/null 2>&1; then
        echo "Errore: Nessun file *_genefamilies.tsv in \$INPUT_DIR" >&2
        exit 2
    fi

    # 1. Join tables
    humann_join_tables \\
        --input "\$INPUT_DIR" \\
        --output "humann_joined_genefamilies.tsv"

    # 2. Normalize
    humann_renorm_table \\
        --input "humann_joined_genefamilies.tsv" \\
        --output "humann_joined_genefamilies_cpm.tsv" \\
        --units cpm

    # 3. Regroup: reazioni Metacyc (default HUMAnN 3.9)
    humann_regroup_table \\
        --input "humann_joined_genefamilies_cpm.tsv" \\
        --groups uniref90_rxn \\
        --output "merged_table_genefamilies_rxn.txt"

    # 4. Rename table: KO
    humann_rename_table \\
        --input "merged_table_genefamilies_rxn.txt" \\
        --names kegg-orthology \\
        --output "merged_table_ko.txt"

    # 5. Rename table: KEGG Pathway
    humann_rename_table \\
        --input "merged_table_genefamilies_rxn.txt" \\
        --names kegg-pathway \\
        --output "merged_table_kegg_pathway.txt"

    echo "HUMAnN merge e annotazione completate."
    conda deactivate
    """
}
