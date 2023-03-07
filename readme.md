# bulkRNAseq in nextflow

This is a standard bulkRNAseq workflow in nextflow inspired by a bulkRNAseq pipeline from BBC at VAI and scwgbs_biscuit_nf from njspix at VAI

## How to run the pipeline
Adjust the slurm job setting and then, `sbatch run_bulkRNAseq.sh`

## To do
- Support `single end`
- Include `SortMeRNA`
- Utilize `singularity` rather than `env module`
- Capture `command line options` for each tools

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