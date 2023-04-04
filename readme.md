# bulkRNAseq in nextflow

This is a standard bulkRNAseq workflow in nextflow inspired by a bulkRNAseq pipeline from BBC at VAI and scwgbs_biscuit_nf from njspix at VAI

## How to run the pipeline
Adjust the slurm job setting and then, `sbatch run_bulk_rnaseq.slurm`

## To do
- [x] Support `single end`
- [x] Include `SortMeRNA`
- [ ] Utilize `singularity` rather than `env module`
- [ ] Capture `command line options` for each tools

## Miscellaneous

### General guidelines for bulkRNAseq analysis
1. General gene-level differential expression
   - ENCODE guidelines suggest 30 million SE reads per sample (stranded).
   - 15 million reads per sample is often sufficient, if there are a good number of replicates (>3).
   - Spend money on more **biological replicates**, if possible.
   - Generally recommended having read length >= 50 bp
2. Gene-level differential expressoin with detection of lowly-expressed genes
   - Similarly benefits from biological replicates more than sequencing depth.
   - Sequence deeper with at least 30-60 million reads depending on level of expression (start with 30 million with a good number of replicates).
   - Generally recommended to have read length >= 50 bp
3. Isoform-level differential expression:
   - Of known isoforms, suggested to have a depth of at least 30 million reads per sample and paired-end reads.
   - Of novel isoforms should have more depth (> 60 million reads per sample).
   - Choose biological replicates over paired/deeper sequencing.
   - Generally recommended having read length >= 50 bp, but longer is better as the reads will be more likely to cross exon junctions
   - Perform careful QC of RNA quality. Be careful to use high quality preparation methods and restrict analysis to high quality RIN # samples.

**This is from hbctraining.github.io** (https://hbctraining.github.io/Intro-to-rnaseq-fasrc-salmon-flipped/lessons/02_experimental_planning_considerations.html)

### Considerations for RNA Seq read length and coverage
Different RNA-Seq experiment types require different sequencing read lengths and depth (number of reads per sample). This bulletin reviews RNA sequencing considerations and offers resources for planning RNA-Seq experiments.

**What resources should I consult first?**

For RNA sequencing, read depth is typically used instead of coverage. Detecting low-expression genes can require an increase in read depth. The [ENCODE project](https://www.encodeproject.org/documents/cede0cbe-d324-4ce7-ace4-f0c3eddf5972/@@download/attachment/ENCODE%20Best%20Practices%20for%20RNA_v2.pdf) has data standards for and that are an excellent resource for many projects.
ENCODE project here RNA-Seq Small RNA sequencing Illumina recommends consulting the primary literature for your field and organism for the most up-to-date
guidance on experiment design.

**How many reads should I target per sample?**

Read depth varies depending on the goals of the RNA-Seq study. Most experiments require 5–200 million reads per sample, depending on organism complexity and size, along with project aims.

- Gene expression profiling experiments that are looking for a quick snapshot of highly expressed genes may only need 5–25 million reads per sample. In these cases, researchers can pool multiple RNA-Seq samples into one lane of a sequencing run, which allows for high multiplexing of samples.

- Experiments looking for a more global view of gene expression, and some information on alternative splicing, typically require 30–60 million reads per sample. This range encompasses most published RNA-Seq experiments for mRNA/whole transcriptome sequencing.

- Experiments looking to get an in-depth view of the transcriptome, or to assemble new transcripts, may require 100–200 million reads. In these cases, researchers may need to sequence multiple samples across several high output sequencing lanes.

- Targeted RNA expression requires fewer reads. For example, Illumina recommends 3 million reads per sample for TruSight RNA Pan Cancer TruSight RNA Fusion Panel, which are compatible with high plexity pooling of samples. 

- miRNA-Seq or small RNA Analysis experiments may require even fewer reads than whole transcriptome sequencing. This requirement varies significantly depending on the tissue type being sequenced. Illumina strongly recommends using the primary literature to determine how many reads are needed, with most applications ranging from 1–5 million reads per sample. 

To determine how many samples can be run at one time, divide the number of reads produced by the flow cell by the number of reads needed per sample: 

- number of reads per flow cell / number of reads per sample=number of samples per flow cell

**How long should my reads be?**

Read length depends on the application and final size of the library. The Library Prep Kit Selector provides read length guidance for each type of RNA-Seq library. Sequencing reads that are longer than the insert length do not provide additional useful data.

- Gene expression / RNA Profiling – Quantifying the coding transcriptome typically requires a short single read (often 50–75 bp) to minimize reading across splice junctions while counting all RNAs in the pool.

- Transcriptome Analysis – Novel transcriptome assembly and annotation projects tend to benefit from longer, paired-end reads (such as 2 x 75 bp or 2 x 100 bp) to enable more complete coverage of the transcripts and identification of novel variants or splice sites. Paired-end reads are required to get
information from both 5’ and 3’ ends of RNA species with stranded RNA-Seq library preparation kits.

- Small RNA Analysis – Due to the short length of small RNA, a single read (usually a 50 bp read)
typically covers the entire sequence. A read length of 50 bp sequences most small RNAs, plus enough of the adapter to be accurately identified and trimmed during data analysis.

Reference(https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000001243)

### Designing the right experiment - How many reads do we need?

The coverage is defined as:

Read length * Number of reads / Length of target sequence

The amount of sequencing needed for a given sample is determined by the goals of the experiment and the nature of the RNA sample.

- For a general view of differential expression: 5 ~ 25M reads per sample
- For alternative splicing and lowly expressed genes: 30 ~ 60M reads per sample
- In-depth view of the transcriptome/assemble new transcripts: 100 ~ 200M reads
- Targeted RNA expression requires fewer reads
- miRNA-Seq or Small RNA Analysis require even fewer reads

Reference(https://bioinformatics-core-shared-training.github.io/Bulk_RNASeq_Course_March23/Bulk_RNAseq_Course_Base/Markdowns/01_Introduction_to_RNAseq_Methods.html#6)