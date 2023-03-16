#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=bulk-rnaseq
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH -o mrnaseq.%A.o
#SBATCH -e mrnaseq.%A.e

module load nextflow/22.04.3
module load singularity/3.8.0

nextflow run bulk_rnaseq.nf -c ./conf/run.config -resume

module unload nextflow/22.04.3
module unload singularity/3.8.0