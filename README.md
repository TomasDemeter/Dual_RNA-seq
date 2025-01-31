# Dual RNA-seq Analysis Pipeline

This Snakemake pipeline is designed for processing and analyzing dual RNA sequencing data, capable of handling both human and bacterial RNA-seq samples simultaneously.

## Overview

The Snakemake pipeline performs the following steps:
1. Quality control and trimming of raw reads (fastp)
2. Quality assessment of sequencing data (FastQC)
3. Alignment to human genome (HISAT2)
4. Alignment to bacterial genomes (Bowtie2)
5. Read counting (featureCounts)
6. Quality control report generation (MultiQC)

## Installation

1. Clone this git repository
2. Install all the conda environments in `workflow/envs` using `conda env create -f <environment>.yml` 

## Setup

1. Place your raw sequencing data in the directory specified by `raw_reads_dir` in the config file
2. Place your genomes in the directory specified by `genomes_dir` in the config file
   - make sure that the names of the **DIRECTORIES** are exactly the same as the genomes names in metadata.csv
   - the directory should contain **ONE** .fna file and **ONE** .gtf(or .gff) file
3. Prepare a metadata CSV file with the following columns:
   - sample_id
   - genome_id
   - **ONE** sample and **ONE** genome per row!!! If you want to map a sample to multiple genomes, create multiple rows with the same sample and different genomes. The order of the rows doesnt matter.
4. Configure the `config.yaml` file in the `config` directory with appropriate parameters for your experiment.
   - suffix of your raw reads
   - the pattern that distinguises forward and reverse reads
   - parameters for individual tools 


## Usage

### Running on the SIT Cluster using Slurm
```bash
snakemake -s workflow/Snakefile.py --workflow-profile ./profiles/dual_seq_pipeline/ -n
```

### Running Locally
```bash
snakemake -s workflow/Snakefile.py --profile profiles/default -n
```
Note: Remove the `-n` flag after verifying the dry run.

## Pipeline Output

The pipeline generates the following outputs:

1. **Quality Control**
   - Trimmed reads (fastp)
   - FastQC reports
   - MultiQC summary reports

2. **Alignment Results**
   - Human genome alignments (HISAT2)
   - Bacterial genome alignments (Bowtie2)

3. **Quantification**
   - Feature counts for human genes
   - Feature counts for bacterial genes
   - Combined feature counts
## Configuration

The pipeline requires a configuration file (`config.yaml`) with the following parameters:
- `raw_reads_dir`: Directory containing raw sequencing files
- `result_dir`: Directory for output files
- `genomes_dir`: Directory containing reference genomes
- `log_dir`: Directory for log files
- `metadata`: Path to the sample metadata CSV file

## Metadata Format

The metadata CSV file should contain information about samples and their corresponding genomes:

```csv
sample_id,genome_id
sample1,genome1
sample1,genome2
sample2,genome1
sample2,genome3
```

## Directory Structure

```plaintext
.
├── config/
│   └── config.yaml
├── inputs/
│   └── adapters/
│       └── adapters.fa
│   └── genomes/
│       └── genome1_folder/   # These must be the same as the genome names in metadata csv file
│           ├── genome1.fna
│           └── genomic.gtf
│       └── genome2_folder/   # These must be the same as the genome names in metadata csv file
│           ├── genome2.fna
│           └── genomic.gtf
│       └── genome3_folder/   # These must be the same as the genome names in metadata csv file
│           ├── genome3.fna
│           └── genomic.gtf
│   └── metadata/
│       └── metadata.csv
│   └── raw_reads/
│           ├── sample1.fastq
│           └── sample2.fastq
├── results/
│   └── bowtie2/
│       └── genome1_folder/
│           ├── sample1.bam
│           ├── sample1.bai
│           ├── sample2.bam
│           └── sample2.bami
│       └── genome2_folder/
│           ├── sample1.bam
│           └── sample1.bai
│       └── genome3_folder/
│           ├── sample2.bam
│           └── sample2.bami
│   └── fastp/
│       ├── logs/
│       ├── sample1.fastq
│       └── sample2.fastq
│   └── FastQC/
│   └── featureCounts/
│       └── genome1_folder/
│           ├── sample1.bam
│           ├── sample1.bai
│           ├── sample2.bam
│           └── sample2.bami
│       └── genome2_folder/
│           ├── sample1.bam
│           └── sample1.bai
│       └── genome3_folder/
│           ├── sample2.bam
│           └── sample2.bami
│   └── hisat2/
│       └── mapped/
│           ├── sample1.bam
│           ├── sample2.bam
│           └── sample3.bam
│       └── unmapped/
│           ├── sample1.1
│           ├── sample1.2
│           ├── sample2.1
│           ├── sample2.2
│           ├── sample3.1
│           └── sample3.2
│   └── logs/
│       ├── bowtie2/
│       ├── bowtie2_index/
│       ├── fastp/
│       ├── FastQC/
│       ├── hisat2/
│       ├── hisat2_index/
│       └── MultiQC/
│   └── MultiQC/
│       └── genome1_folder/
│           └── genome1_report.html
│       └── genome2_folder/
│           └── genome2_report.html
│       └── genome3_folder/
│           └── genome3_report.html
|   └── config_used_in_this_experiment.yaml
|   └── metadata_used_in_this_experiment.csv
├── workflow/
│   ├── Snakefile.py
│   └── rules/
│       ├── fastp.smk
│       ├── hisat2.smk
│       ├── bowtie2.smk
│       ├── fastqc.smk
│       ├── featuecounts.smk
│       └── multiqc.smk
|    └── envs/
│       ├── fastp_bowtie2.yml
│       ├── featurecounts_env.yml
│       ├── snakemake.yml
│       └── multiqc_env.yml 
└── profiles/
    ├── dual_seq_pipeline/
    └── default/
```