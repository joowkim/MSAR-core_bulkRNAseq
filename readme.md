# bulkRNAseq in nextflow

This is a standard bulkRNAseq workflow in nextflow inspired by a bulkRNAseq pipeline from BBC at VAI and scwgbs_biscuit_nf from njspix at VAI

## How to run the pipeline
Adjust the slurm job setting and then, `sbatch run_bulkRNAseq.sh`

## To do
- Support `single end`
- Include `SortMeRNA`
- Utilize `singularity` rather than `env module`
- Capture `command line options` for each tools