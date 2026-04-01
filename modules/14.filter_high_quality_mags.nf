process filter_high_quality_mags {
        tag "filter"

    input:
        path quality_report        // <--- il file quality_report.tsv unico (output di combine_quality_reports)
        path mag_files            // <--- tutti i MAG .fa come array/lista di path

    output:
        //path "MAGs_high_quality"  // <--- directory con i MAG filtrati
        path "MAGs_high_quality/my-metadata2.txt"
        path "tormes_output"

    publishDir "${params.output_base}/14_MAGs_high_quality/", pattern: "MAGs_high_quality/my-metadata2.txt", mode: 'copy'    
    publishDir "${params.output_base}/15_tormes_output/", pattern: "tormes_output", mode: 'copy'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bioenv

    mkdir -p MAGs_high_quality

    # Filtro: completeness > 50, contamination < 10 (modifica se vuoi)
    awk -F'\\t' 'NR==1 || (\$3 > ${params.filter.completeness} && \$4 < ${params.filter.contamination})' $quality_report > filtered_data.tsv

    # Ricava i nomi dei MAG filtrati
    awk -F'\\t' 'NR>1 {print \$2".fa"}' filtered_data.tsv > mags_to_keep.txt

    while read -r filename; do
        if [[ -e \$filename ]]; then
            cp \$filename MAGs_high_quality/
        else
            echo "ATTENZIONE: \$filename non trovato" >&2
        fi
    done < mags_to_keep.txt

    # Genera metadati filtrati con path relativo
    ls MAGs_high_quality/*.fa | sed "s|.*/||; s/.fa//" > uno.tmp
    while read -r filename; do
        echo "GENOME" >> dos.tmp
        echo "MAGs_high_quality/\${filename}.fa" >> tres.tmp
        echo "This is the genome for \${filename}" >> cuatro.tmp
    done < uno.tmp
    paste uno.tmp dos.tmp tres.tmp cuatro.tmp | sed "1iSamples\tRead1\tRead2\tDescription" > MAGs_high_quality/my-metadata2.txt
    rm uno.tmp dos.tmp tres.tmp cuatro.tmp

    # Metadata per tutti i MAG rimasti nella cartella filtrata (opzionale, come nel bash)
    #ls MAGs_high_quality/*.fa | sed "s|.*/||; s/.fa//" > uno.tmp
    #while read -r filename; do
    #    echo "GENOME" >> dos.tmp
    #    realpath "MAGs_high_quality/\${filename}.fa" >> tres.tmp
    #    echo "This is the genome for \${filename}" >> cuatro.tmp
    #done < uno.tmp
    #paste uno.tmp dos.tmp tres.tmp cuatro.tmp | sed "1iSamples\tRead1\tRead2\tDescription" > MAGs_high_quality/all_metadata.tsv
    #rm uno.tmp dos.tmp tres.tmp cuatro.tmp

    conda deactivate

    set +u

    source /opt/conda/etc/profile.d/conda.sh
    conda activate tormes-1.3.0

    tormes --metadata MAGs_high_quality/my-metadata2.txt --output tormes_output --threads ${task.cpus}

    conda deactivate
    """
}
