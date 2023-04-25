nextflow.enable.dsl=2

process fastqc {
    debug true
    tag "Fastqc on ${reads}"
    //label "universal" // getting cpu and memory usage from salmon.config - called universal
    cpus 8
    memory '8 GB'

    publishDir "${launchDir}/analysis/fastqc/", mode : "copy"

    module 'FastQC/0.11.9'

    input:
    path(reads)

    output:
    path ("*.zip"), emit: zips
    path ("*.html"), emit: htmls

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}


process multiqc {
    debug true
    tag "Multiqc on this project"

    cpus 1
    memory '1 GB'

    publishDir "${launchDir}/analysis/multiqc/", mode : "copy"

    //module 'python/3.6.2'

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    """
    multiqc ${files} --filename "multiqc_report.html"
    """
}


process fastq_screen {
    debug true
    tag "Fastq-screen on ${reads}", mode : "copy"

    cpus 4
    memory '8 GB'

    publishDir "${launchDir}/analysis/fastq_screen"

    module 'FastQScreen/0.14.1'
    module 'bowtie2/2.3.4.1'

    input:
    path(reads)

    output:
    path("*.html")
    path("*.txt"), emit: fastq_screen_out

    // threads option is already defined in fastq_screeN_conf
    script:
    conf = "/mnt/beegfs/kimj32/polymerase/polymeraseDependencies/FastQ_Screen_Genomes/fastq_screen.conf"
    """
    fastq_screen --aligner bowtie2 \
    --conf ${conf} \
    ${reads} \
    --outdir ./
    """
}


workflow {
    ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
    fastqc(ch_reads)
    fastq_screen(ch_reads)
    multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips).collect() )

}

workflow.onComplete {
	println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}