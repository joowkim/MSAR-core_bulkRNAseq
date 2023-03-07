#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=bulk-rnaseq
#SBATCH --ntasks=8
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30000
#SBATCH --time=12:00:00
#SBATCH -o bulk_rnaseq.%A.o
#SBATCH -e bulk_rnaseq.%A.e

module load nextflow/22.04.3

nextflow run bulk_rnaseq.nf

module unload nextflow/22.04.3