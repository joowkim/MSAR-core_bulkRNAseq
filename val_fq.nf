nextflow.enable.dsl=2

process validate_fastq {
    debug true
    tag "${meta.sample_name}"
    // label "universal"
    cpus 4
    memory '8 GB'
    time "2h"

    publishDir "${launchDir}/analysis/validate_fastq/", mode: "copy"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.sample_name}.{err,out}")

    script:
    def val_jar = "/mnt/beegfs/kimj32/tools/validatefastq-assembly-0.1.1.jar"
    if(!meta.single_end) {
        """
        java -jar ${val_jar} \
        -i ${reads[0]} \
        -j ${reads[1]} \
        -l warn 2>${meta.sample_name}.err 1>${meta.sample_name}.out

        """
    } else {
        """
        java -jar ${val_jar} \
        -i ${reads[0]} \
        -l warn 2>${meta.sample_name}.err 1>${meta.sample_name}.out
        """
    }
}


// See https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

// adapted from https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_reads = ch_samplesheet.splitCsv(header:true).map {

    // This is the read1 and read2 entry
    r1 = it['fq1']
    r2 = it['fq2']

    // Detect wiether single-end or paired-end
    is_singleEnd = r2.toString() =='' ? true : false

    // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
    meta = [sample_name: it['sample'], single_end: is_singleEnd]

    // We return a nested map, the first entry is the meta map, the second one is the read(s)
    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

}

workflow {

    validate_fastq(ch_reads)

}

workflow.onComplete {
	println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
