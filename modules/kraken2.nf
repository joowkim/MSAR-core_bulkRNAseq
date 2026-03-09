process kraken2 {
    tag "${sample_name}"
    label "process_high"

    publishDir "${projectDir}/analysis/kraken2/"

    module "kraken/2.1.2"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("${sample_name}_kraken2_report.txt"), emit: kraken2_report
    path("${sample_name}_kraken2_output.txt")

    script:
    def db = params.kraken2_db
    if (!is_SE) {
    """
    kraken2 --db ${db} \
        --threads ${task.cpus} \
        --paired ${reads[0]} ${reads[1]} \
        --report ${sample_name}_kraken2_report.txt \
        --output ${sample_name}_kraken2_output.txt
    """
    } else {
    """
    kraken2 --db ${db} \
        --threads ${task.cpus} \
        ${reads[0]} \
        --report ${sample_name}_kraken2_report.txt \
        --output ${sample_name}_kraken2_output.txt
    """
    }
}
