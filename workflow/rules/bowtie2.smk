################################
# create an index for a genome #
################################
rule bowtie2_index_single:
    input:
        genome      = lambda wildcards: glob.glob(os.path.join(config["genomes_dir"], wildcards.genome, "*.fna"))[0]
    output:
        flag        = os.path.join(config["genomes_dir"],"{genome}/bowtie_indices/.index_done")
    params:
        index_base  = os.path.join(config["genomes_dir"], "{genome}/bowtie_indices", "{genome}")
    log:
        logs        = os.path.join(LOG_DIR, "bowtie2_index/{genome}/bowtie_indices", "log.txt")
    conda:
        "fastp_bowtie2"
    shell:
        """
        mkdir -p $(dirname {output.flag})
        bowtie2-build {input.genome} {params.index_base} \
            --threads {resources.cpus_per_task} \
            2 > {log.logs}; \
        touch {output.flag}
        """

################################################
# execute genome indexing rule for all genomes #
################################################ 
rule bowtie2_index_all:
    input:
        flags   = expand(rules.bowtie2_index_single.output.flag, genome = GENOMES)
    output:
        flag    = os.path.join(GENOMES_DIR, ".all_indices_done")
    shell:
        """
        touch {output.flag}
        """


##################################################
# map one sample to one of the specified genomes #
##################################################
rule bowtie2_map_single:
    input:
        read_1      = rules.hisat2_H_sapiens.output.unmapped_1,
        read_2      = rules.hisat2_H_sapiens.output.unmapped_2,
        index_done  = rules.bowtie2_index_all.output.flag
    output:
        bam         = os.path.join(RESULT_DIR, "bowtie2/{genome}/{sample}.bam"),
        bai         = os.path.join(RESULT_DIR, "bowtie2/{genome}/{sample}.bam.bai"),
        flag        = os.path.join(RESULT_DIR, "bowtie2/{genome}/.{sample}.flag")
    params:
        idx         = os.path.join(config["genomes_dir"], "{genome}/bowtie_indices/{genome}")
    log:
        logs        = os.path.join(LOG_DIR, "bowtie2/{genome}/{sample}.log")
    conda:
        "fastp_bowtie2"
    shell:
        """
        (mkdir -p $(dirname {output.bam}) && \
        mkdir -p $(dirname {log.logs}) && \
        bowtie2 -x {params.idx} \
                -1 {input.read_1} \
                -2 {input.read_2} \
                --threads {resources.cpus_per_task} \
                --sensitive-local | \
        samtools view -bS - | \
        samtools sort -@ {resources.cpus_per_task} -m 2G -o {output.bam} - && \
        samtools index {output.bam} && \
        touch {output.flag}) 2> {log.logs}
        """  

#################################################################################
# execute the mapping rule for all required combinations of samples and genomes # 
#################################################################################
rule map_all_samples:
    input:
        flags   = expand(
            os.path.join(
                RESULT_DIR,
                "bowtie2/{genome}/.{sample}.flag"
            ), 
            zip, 
            sample = samples_df["sample_id"], 
            genome = samples_df["genome_id"]
        )
    output:
        flag    = os.path.join(RESULT_DIR, "bowtie2/.all_samples_mapped")
    shell:
        """
        touch {output.flag}
        """