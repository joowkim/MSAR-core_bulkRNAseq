nextflow.enable.dsl=2

process FASTQC {
    debug true
    //tag "Fastqc on ${meta.sample_name}"
    //label "universal" // getting cpu and memory usage from salmon.config - called universal
    cpus 8
    memory '4 GB'

    publishDir "${projectDir}/analysis/fastqc/"

    module 'FastQC/0.11.9'

    input:
    tuple val(meta), path(reads)

    output:
    path ("*.zip"), emit: zips
    path ("*.html"), emit: htmls

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}

// process FASTP {
//     debug true
//     tag "Fastp on ${sample_name}"
//     // label "universal"
//     cpus 8
//     memory '8 GB'
//
//     publishDir "${projectDir}/analysis/fastp/"
//
//     module 'fastp/0.21.0'
//
//     input:
//     tuple val(sample_name), path(reads)
//
//     output:
//     tuple val(sample_name), path("${sample_name}_trimmed.R{1,2}.fq.gz"), emit: trim_reads
//     path("${sample_name}.fastp.json"), emit: json
//
//     script:
//     """
//     fastp \
//     -i ${reads[0]} \
//     -I ${reads[1]} \
//     --thread ${task.cpus} \
//     --detect_adapter_for_pe \
//     --qualified_quality_phred 25 \
//     -o ${sample_name}_trimmed.R1.fq.gz \
//     -O ${sample_name}_trimmed.R2.fq.gz \
//     --json ${sample_name}.fastp.json
//     """
// }

process TRIM_GALORE {
    debug true
    //tag "trim_galore on ${meta.sample_name}"
    // label "universal"
    cpus 4
    memory '8 GB'

    publishDir "${projectDir}/analysis/trim_galore/"

    module 'fastp/0.21.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("*.gz"), emit: trim_reads
    path("*.html")
    path("*.zip")
    path("*.txt")

    script:
    if(meta.single_end) {
    """
    trim_galore \
    ${reads} \
    --cores ${task.cpus} \
    -q 20 \
    --fastqc \
    --output_dir ./
    """
    } else {
    """
    trim_galore \
    --paired \
    ${reads[0]} ${reads[1]} \
    --cores ${task.cpus} \
    -q 20 \
    --fastqc \
    --output_dir ./
    """
    }
}

process STAR {
    debug true
    tag "STAR on ${sample_name}"
    // label "universal"
    cpus 8
    memory '64 GB'

    publishDir "${projectDir}/analysis/star/", mode : "copy"

    module "STAR/2.7.10a"
    module 'samtools/1.16.1'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_name), path("${sample_name}.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(sample_name), path("${sample_name}.Log.final.out"), emit: log_final
    tuple val(sample_name), path("${sample_name}.Log.out"), emit: log_out
    tuple val(sample_name), path("${sample_name}.ReadsPerGene.out.tab"), emit: read_per_gene_out
    tuple val(sample_name), path("${sample_name}.SJ.out.tab"), emit: sj_out
    tuple val(sample_name), path("${sample_name}._STAR*"), emit: out_dir // STARgenome and STARpass1

    script:
    index = params.ref_fa.(params.genome)
    """
    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${index} \
    --readFilesIn ${reads} \
    --twopassMode Basic \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${sample_name}. \
    --quantMode GeneCounts \
    --outStd Log 2> ${sample_name}.log \

    samtools index "${sample_name}.Aligned.sortedByCoord.out.bam"
    """
}

process SEQTK {
    debug true
    tag "seqtk on ${sample_name}"

    cpus 2
    memory '8 GB'

    publishDir "${projectDir}/analysis/seqtk/"

    module 'seqtk/1.3'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}.subsample.100000.R1.fq.gz"), emit: subsample_reads

    script:
    """
    seqtk sample -s 100 ${reads[0]} 100000 | gzip -c > ${sample_name}.subsample.100000.R1.fq.gz
    """
}

 process SORTMERNA {
     debug true
     tag "SortMeRNA on ${sample_name}"

     cpus 8
     memory '8 GB'

     publishDir "${projectDir}/analysis/sortMeRNA/"

     module 'sortmerna/4.3.6'

     input:
     tuple val(sample_name), path(reads)

     output:
     path ("*"), emit: sortMeRNA_out

     script:
     idx = "/mnt/beegfs/kimj32/tools/sortmerna/idx"
     rfam5_8s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta"
     rfam5s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta"
     silva_arc_16s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta"
     silva_arc_23s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta"
     silva_euk_18s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta"
     silva_euk_28s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta"
     silva_bac_16s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta"
     silva_bac_23s = "/home/kimj32/beegfs/tools/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta"
     """
     sortmerna --threads ${task.cpus} \
     -reads ${reads[0]} \
     --workdir sortMeRNA_${sample_name}  \
     --idx-dir ${idx}  \
     --ref ${rfam5s}  \
     --ref ${rfam5_8s}  \
     --ref ${silva_arc_16s}  \
     --ref ${silva_arc_23s}  \
     --ref ${silva_bac_16s}  \
     --ref ${silva_bac_23s}  \
     --ref ${silva_euk_18s}  \
     --ref ${silva_euk_28s}
     """
}


process MULTIQC {
    debug true
    tag "Multiqc on this project"

    cpus 1
    memory '4 GB'

    publishDir "${projectDir}/analysis/multiqc/", mode : "copy"

    module 'python/3.6.2'

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    """
    multiqc ${projectDir}/analysis --filename "multiqc_report.html" --ignore '*STARpass1'
    """
}

process QUALIMAP {
    debug true
    tag "Qualimap on ${sample_name}"

    cpus 8
    memory '8 GB'

    publishDir "${projectDir}/analysis/qualimap/"

    module 'qualimap/2.2.1'

    input:
    tuple val(sample_name), path(bam)

    output:
    path("*"), emit: qualimap_out

    script:
    gtf = params.gtf.(params.genome)
    """
    qualimap rnaseq -bam ${bam} \
    -gtf ${gtf} \
    --paired \
    -outdir quailmap_${sample_name} \
    --java-mem-size=6G \
    --sequencing-protocol strand-specific-reverse
    """
}

process FASTQSCREEN {
    debug true
    tag "Fastq-screen on ${sample_name}"

    cpus 4
    memory '8 GB'

    publishDir "${projectDir}/analysis/fastq_screen"

    module 'FastQScreen/0.14.1'
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(reads)

    output:
    path("*.html")
    path("*.txt")

    // threads option is already defined in fastq_screeN_conf
    script:
    conf = params.fastq_screen_conf
    """
    fastq_screen --aligner bowtie2 \
    --conf ${conf} \
    ${reads[0]} \
    --outdir ./
    """
}

process TPMCALCULATOR {
    debug true
    tag "tpm_calculator on ${sample_name}"

    cpus 4
    memory '8 GB'

    publishDir "${projectDir}/analysis/tpm_calculator/"

    module "TPMCalculator/0.4.0"

    input:
    tuple val(sample_name), path(bam_file)

    output:
    path("*out"), emit : "tpm_calculator_out"

    script:
    // -p option is for paired end data
    gtf = params.gtf.(params.genome)
    """
    TPMCalculator -g ${gtf} \
    -b ${bam_file} \
    -p

    /home/kimj32/basic-tools/renamer 's/.Aligned.sortedByCoord.out_genes//' *out_genes*
    """
}


log.info """
bulkRNAseq Nextflow
=============================================
samplesheet                           : ${params.samplesheet}
reference                       : ${params.ref_fa.(params.genome)}
"""

reference = Channel.fromPath(params.ref_fa.(params.genome), checkIfExists: true)

// See https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

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

    FASTQC(ch_reads)
    TRIM_GALORE(ch_reads)
    SEQTK(TRIM_GALORE.out.trim_reads)
    FASTQSCREEN(SEQTK.out.subsample_reads)
    STAR(TRIM_GALORE.out.trim_reads)
    SORTMERNA(SEQTK.out.subsample_reads)
    QUALIMAP(STAR.out.bam)
    TPMCALCULATOR(STAR.out.bam)
    MULTIQC( QUALIMAP.out.qualimap_out.mix(SORTMERNA.out.sortMeRNA_out).collect() )

}

workflow.onComplete {
	println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}