#! /usr/bin/env bash

set -euo pipefail

star_output=$1
gtf_path="~/beegfs/reference/mouse/gencode/GRCm38.p6/annotation/gencode.vM25.annotation.gtf"


module load R/4.2.2

Rscript ~/beegfs/src/star_to_mat-not-removal-gene-version.R $star_output/

sed -i '1 s/^gene_id//' star_read_cnt.csv

csvtk transpose star_read_cnt.csv > star_read_cnt.trans.csv

rnanorm tpm star_read_cnt.trans.csv --gtf $gtf_path > tpm.csv

sed -i '1 s/^,/gene_id,/' tpm.csv

csvtk transpose tpm.csv > tpm.trans.csv