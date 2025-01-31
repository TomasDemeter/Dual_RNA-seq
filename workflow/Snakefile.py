
##########################################
# Snakemake pipeline for bacterial RNA-seq 
##########################################


#################################################
# to run the pipeline use one of these commands #
#################################################
# use this to run on the cluster
# snakemake -s workflow/Snakefile.py --workflow-profile ./profiles/dual_seq_pipeline/ -n

# use this to run locally on the laptop
# default resources are based on the withe lab notebook. If you run it on your own machine adjust the resources as needed
# snakemake -s workflow/Snakefile.py --profile profiles/default -n 

####################
# Python pacakages #
####################
import pandas as pd
import glob

#################
# Configuration #
#################
configfile: "../config/config.yaml"  
RAW_READS_DIR = config["raw_reads_dir"]
RESULT_DIR = config["result_dir"]
GENOMES_DIR = config["genomes_dir"]
LOG_DIR = config["log_dir"]

##########################
# Samples and conditions #
##########################
# Read sample information
samples_df = pd.read_csv(config["metadata"], header=0)
samples_df.columns = samples_df.columns.str.strip().str.lower()

# Get unique samples and genomes
SAMPLES = samples_df['sample_id'].unique().tolist()
GENOMES = samples_df['genome_id'].unique().tolist()

#########
# rules #
#########
include: "rules/fastp.smk"
include: "rules/hisat2.smk"
include: "rules/bowtie2.smk"
include: "rules/fastqc.smk"
include: "rules/featuecounts.smk"
include: "rules/multiqc.smk"
include: "rules/save_config.smk"

###################
# Desired outputs #
###################
fastp                   = expand(rules.fastp.output.trimmed_1, sample = SAMPLES)
hisat2                  = expand(rules.hisat2_H_sapiens.output.unmapped_1, sample = SAMPLES)
bowtie2_map_single      = expand(rules.bowtie2_map_single.output.flag, sample = SAMPLES, genome = GENOMES)
featureCounts_single    = expand(rules.featureCounts_single.output.counts, sample = SAMPLES, genome = GENOMES)
featureCounts_all       = expand(rules.featureCounts_all.output.flag, sample = SAMPLES)
featureCounts_H_sapiens = expand(rules.featureCounts_H_sapiens.output.counts, sample = SAMPLES)
fastqc                  = expand(rules.FastQC.output.html_1, sample = SAMPLES)
multiqc                 = expand(rules.multiqc.output.output, genome = GENOMES)
multiqc_human           = expand(rules.multiqc_H_sapiens.output.output, sample = SAMPLES)
save_config             = rules.save_config.output
rule all:
    input:
        fastp,
        bowtie2_map_single,
        featureCounts_single,
        featureCounts_all,
        featureCounts_H_sapiens,
        fastqc,
        multiqc,
        multiqc_human,
        save_config
    message:
        "RNA-seq pipeline finished successfully"