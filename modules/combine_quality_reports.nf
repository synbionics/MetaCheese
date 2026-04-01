process combine_quality_reports {

    tag "combine_reports"

    input:
    path(checkm2_files)

    output:
    path "quality_report.tsv"

    publishDir "${params.output_base}/13_checkm2/", mode: 'copy'

    script:
    """
    set -e

    echo -e "sample_name\tCompleteness\tContamination\tOtherColumns..." > quality_report.tsv

    for f in ${checkm2_files}; do
        sample=\$(basename "\$f" .checkm2.tsv)
        if [ -s "\$f" ]; then
            tail -n +2 "\$f" | awk -v s="\$sample" 'BEGIN{OFS="\\t"} {print s, \$0}' >> quality_report.tsv
        fi
    done
    """
}
