nextflow.enable.dsl=2

def getLibraryId( prefix ){
  // prefix = GC.NG.1429_1_S74_L001_R1_001.fastq.gz
  // Return the ID number, you can change for other file formats, here it just takes the first part before "_"
  prefix.split("_")[0] + "_" + prefix.split("_")[1]
}

//params.raw_data_dir = "rawdata/"

// Gather the pairs of R1/R2 according to sample ID
Channel
     .fromFilePairs(params.rawdata + '/*R{1,2}*.fastq.gz', flat: true, checkExists: true)
     .map { prefix, R1, R2 -> tuple(getLibraryId(prefix), R1, R2) }
     .groupTuple().set{ files_channel }


process merge_lane {
    debug true
    tag "merging ${sample}"
    cpus 1
    memory '1 GB'
    time '2h'

    publishDir "${launchDir}/analysis/merge_lane", mode : "copy"

    input:
        tuple val(sample), path(R1), path(R2)
    output:
        path("${sample}_R1.fastq.gz")
        path("${sample}_R2.fastq.gz")
        // .findAll{ it.contains("_R1_001")}
    script:
        """
        cat ${ R1.collect{ it }.join(" ") } > ${sample}_R1.fastq.gz
        cat ${ R2.collect{ it }.join(" ") } > ${sample}_R2.fastq.gz
        """
}

workflow {
    merge_lane(files_channel)
}

workflow.onComplete {
	println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}

// how to run
// nextflow run this.nf --raw_data_dir "some_path"
