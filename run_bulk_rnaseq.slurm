#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=nf-bulk_rnaseq
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH -o mrnaseq.%A.o
#SBATCH -e mrnaseq.%A.e

module load nextflow/23.10.0
module load singularity/3.8.0
module load java/21.0.6

nextflow run bulk_rnaseq.nf -c ./bulk_rnaseq_conf/run.config -resume -profile slurm

module unload nextflow/23.10.0
module unload singularity/3.8.0
module unload java/21.0.6