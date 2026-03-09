# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bulk RNA-seq analysis pipeline for the MSAR core at Cleveland Clinic Foundation. Processes paired-end (or single-end) RNA-seq FASTQ files through QC, alignment, gene counting, and transcript quantification using Nextflow (DSL2) with Singularity containers on a SLURM cluster.

## Running the Pipeline

### Main workflow
```bash
sbatch run_bulk_rnaseq.slurm
```
Or directly:
```bash
module load nextflow/23.10.0 singularity/3.8.0 java/21.0.6
nextflow run bulk_rnaseq.nf -c ./bulk_rnaseq_conf/run.config -resume -profile slurm
```

### FASTQ validation (pre-check)
```bash
sbatch run_val_fq.slurm
```

### QC-only workflow
```bash
sbatch src/run_qc.slurm
```

### Generate samplesheet from a FASTQ directory
```bash
Rscript src/make_samplesheet_template.R -i /path/to/fastq_dir
```

## Configuration

All workflow config lives in `bulk_rnaseq_conf/`:

- **run.config** — primary entry point: samplesheet path, genome (`human`, `mouse`, `pig`, `zebra_fish`), and feature toggles (`run_salmon`, `run_tpm_calculator`, `run_ribodetector`, `run_kraken2`, `run_preseq`)
- **reference.config** — paths to STAR/Salmon indexes, GTF, FASTA, gene symbol lookup tables, sortMeRNA DB (`params.sortmerna_db`), MultiQC config (`params.multiqc_config`), and Kraken2 DB (`params.kraken2_db`) per genome
- **processes.config** — resource labels (CPUs/RAM) and Singularity image paths per process
- **cluster.config** — SLURM executor settings (queue `defq`, `queueSize = 35`, `submitRateLimit = "15 sec"`)

## Architecture

### Workflow (`bulk_rnaseq.nf`)
Main DSL2 Nextflow workflow. Key process order:

1. `fastqc` → raw read QC
2. `fastp` → adapter trimming (`--trim_front1 1 --trim_front2 1`)
3. `seqtk` → subsample 100k reads for screening
4. `fastq_screen` → contamination detection (R1 only)
5. `sortMeRNA` → rRNA detection on subsampled R1 reads (SE and PE use same command)
6. `ribo_detector` → rRNA removal from trimmed reads (optional, `run_ribodetector`)
7. `star` → alignment + gene counts (`--quantMode GeneCounts`, TwoPass mode); unmapped reads saved for Kraken2
8. `samtools_index` → BAM indexing
9. `qualimap` → RNA-seq QC metrics
10. `feature_count` → featureCounts (`-s 2`, reverse strand)
11. `gene_count_mat` → compiles STAR counts into matrix (calls `bin/star_to_mat.R` via Nextflow `bin/` PATH)
12. `gene_id_to_gene_symbol` → annotates with gene symbols (calls `bin/gene_symbol_annot.R`)
13. `multiqc` → aggregated QC report; inputs built incrementally based on which optional modules ran

Optional modules (all default `false` except `run_salmon`):
- `run_salmon` — transcript-level quantification; uses ribodetector-filtered reads if `run_ribodetector = true`
- `run_kraken2` — contamination check on STAR unmapped reads; report fed into MultiQC
- `run_tpm_calculator` — TPM calculation from BAM
- `run_preseq` — library complexity (currently throws errors, under investigation)

### Modules (`modules/`)
- **kraken2.nf** — runs Kraken2 on STAR unmapped reads; takes `tuple val(sample_name), path(reads), val(is_SE)`

### R scripts (`bin/`)
- **star_to_mat.R** — reads STAR `ReadsPerGene.out.tab` files, detects strandedness heuristically, outputs count matrix CSV
- **gene_symbol_annot.R** — joins count matrix with genome-specific gene symbol lookup table

### Supporting scripts (`src/`)
- **make_samplesheet_template.R** — generates `samplesheet.csv` from a FASTQ directory using regex-based sample naming
- **rename.py / concat_rename.py** — rename/concatenate demultiplexed FASTQ files
- **de_analysis_template.qmd** — Quarto template for downstream DESeq2 differential expression analysis

## Samplesheet Format

```csv
sample,fq1,fq2
sampleA,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz
```

For single-end data, omit the `fq2` column.

## Key Design Decisions

- **Strandedness**: Reverse strand (`-s 2` in featureCounts) for Illumina TruSeq Stranded libraries
- **Subsampling**: 100k reads used only for `fastq_screen` and `sortMeRNA` (not for alignment)
- **SE/PE output globs**: Use `R*.fastq.gz` pattern (not `R{1,2}`) so SE samples producing only R1 don't fail output collection
- **Error strategy**: `errorStrategy = "finish"` — all running jobs complete before pipeline exits on error
