#!/usr/bin/bash

if [ $# -ne 3 ]; then
  echo "bash this.sh [delimiter] [input_fq_dir] [merged_fq_dir]"
  exit 1
fi

set -eou pipefail

delimiter=$1
input_fq_dir=$2
merged_fq_dir=$3

mkdir -p $merged_fq_dir

ls $input_fq_dir/*R1* | cut -d $delimiter -f 1 | sort | uniq | sed "s/$input_fq_dir\///" \
    | while read id; do \
        echo $input_fq_dir/$id*R1*.fastq.gz --\> $merged_fq_dir/$id\_R1.fastq.gz;
        cat $input_fq_dir/$id*R1*.fastq.gz > $merged_fq_dir/$id\_R1.fastq.gz;
        echo $input_fq_dir/$id*R2*.fastq.gz --\>  $merged_fq_dir/$id\_R2.fastq.gz;
        cat $input_fq_dir/$id*R2*.fastq.gz > $merged_fq_dir/$id\_R2.fastq.gz;
      done
