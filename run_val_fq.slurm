#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=validate_fastq
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH -o val_fq.%A.o
#SBATCH -e val_fq.%A.e

module load nextflow/22.04.3

nextflow run val_fq.nf -c ./val_fq_conf/run.config -resume -profile slurm

module unload nextflow/22.04.3