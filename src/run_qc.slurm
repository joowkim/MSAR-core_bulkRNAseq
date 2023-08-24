#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=QCing
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH -o qc.%A.o
#SBATCH -e qc.%A.e

module load nextflow/22.04.3
module load singularity/3.8.0

nextflow run qc.nf -c ./conf/run.config -resume -profile slurm

module unload nextflow/22.04.3
module unload singularity/3.8.0