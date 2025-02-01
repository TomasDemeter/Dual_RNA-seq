################################
# create an index for a genome #
################################
rule hisat2_index_single:
    input:
        genome      = lambda wildcards: glob.glob(os.path.join(config["genomes_dir"], wildcards.genome, "*.fna"))[0]
    output:
        flag        = os.path.join(config["genomes_dir"], "{genome}/hisat_indices/.index_done")
    params:
        index_base  = os.path.join(config["genomes_dir"], "{genome}/hisat_indices", "{genome}")
    log:
        logs        = os.path.join(LOG_DIR, "hisat2_index/{genome}/hisat_indices", "log.txt")
    conda:
        "hisat2_env"
    shell:
        """
        mkdir -p $(dirname {output.flag})
        hisat2-build {input.genome} {params.index_base} \
            --threads {resources.cpus_per_task} \
            2 > {log.logs}; \
        touch {output.flag}
        """

################################################
# execute genome indexing rule for all genomes #
################################################ 
rule hisat2_index_all:
    input:
        flags   = expand(rules.hisat2_index_single.output.flag, genome = GENOMES)
    output:
        flag    = os.path.join(GENOMES_DIR, ".all_hisat_indices_done")
    shell:
        """
        touch {output.flag}
        """


##################################################
# map one sample to one of the specified genomes #
##################################################
rule hisat2_map_single:
    input:
        read_1      = rules.hisat2_H_sapiens.output.unmapped_1,
        read_2      = rules.hisat2_H_sapiens.output.unmapped_2,
        index_done  = rules.hisat2_index_all.output.flag
    output:
        bam         = os.path.join(RESULT_DIR, "hisat2/bacteria/{genome}/{sample}.bam"),
        bai         = os.path.join(RESULT_DIR, "hisat2/bacteria/{genome}/{sample}.bam.bai"),
        flag        = os.path.join(RESULT_DIR, "hisat2/bacteria/{genome}/.{sample}.flag")
    params:
        idx         = os.path.join(config["genomes_dir"], "{genome}/hisat_indices/{genome}"),
        no_splice   = "--no-spliced-alignment" if config["hisat2"]["no_spliced_alignment"] else "",
        phred       = "--phred33" if config["hisat2"]["phred33"] else "",
        max_align   = "-k " + str(config["hisat2"]["max_alignments"]),
        mismatch    = "--mp " + config["hisat2"]["mismatch_penalty"],
        score_min   = "--score-min " + config["hisat2"]["score_min"],
        read_gap    = "--rdg " + config["hisat2"]["read_gap"],
        ref_gap     = "--rfg " + config["hisat2"]["ref_gap"],
        dta         = "--dta" if config["hisat2"]["dta"] else "",
        summary     = "--new-summary" if config["hisat2"]["new_summary"] else ""
    log:
        logs        = os.path.join(LOG_DIR, "hisat2/bacteria/{genome}/{sample}.log")
    conda:
        "hisat2_env"
    shell:
        """
        (hisat2 -x {params.idx} \
            -1 {input.read_1} \
            -2 {input.read_2} \
            --threads {resources.cpus_per_task} \
            {params.no_splice} \
            {params.phred} \
            {params.max_align} \
            {params.mismatch} \
            {params.score_min} \
            {params.read_gap} \
            {params.ref_gap} \
            {params.dta} \
            {params.summary} \
        | samtools view -bS - \
        | samtools sort -@ {resources.cpus_per_task} -m 2G -o {output.bam} - \
        && samtools index {output.bam} \
        && touch {output.flag}) 2> {log.logs}
        """

#################################################################################
# execute the mapping rule for all required combinations of samples and genomes # 
#################################################################################
rule map_all_samples:
    input:
        flags   = expand(
            os.path.join(
                RESULT_DIR,
                "hisat2/bacteria/{genome}/.{sample}.flag"
            ), 
            zip, 
            sample = samples_df["sample_id"], 
            genome = samples_df["genome_id"]
        )
    output:
        flag    = os.path.join(RESULT_DIR, "hisat2/.all_samples_mapped")
    shell:
        """
        touch {output.flag}
        """