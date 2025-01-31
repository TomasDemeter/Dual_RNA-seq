rule featureCounts_single:
    input:
        bowtie_flag = rules.map_all_samples.output.flag,
        bam         = os.path.join(RESULT_DIR, "bowtie2/{genome}/{sample}.bam")
    output:
        counts      = os.path.join(RESULT_DIR, "featureCounts/{genome}/{sample}_{genome}_counts.txt")
    log:
        log         = os.path.join(LOG_DIR, "featureCounts/{genome}/{sample}/{genome}.log")
    params:
        gtf_files   = lambda wildcards: glob.glob(os.path.join(config["genomes_dir"], wildcards.genome, "*.gtf"))[0],
        strand      = config["featureCounts"]["strand"],
        paired_end  = config["featureCounts"]["paired_end"],
        features    = config["featureCounts"]["features"],
        attribute   = config["featureCounts"]["attribute"]
    conda:
        "featurecounts_env"
    shell:
        """
        featureCounts \
            -a {params.gtf_files} \
            -o {output.counts} \
            {params.paired_end} \
            -s {params.strand} \
            {params.features} \
            {params.attribute} \
            -T {resources.cpus_per_task} \
            {input.bam} 2> {log.log}
        """

rule featureCounts_all:
    input:
        expand(
            os.path.join(RESULT_DIR, "featureCounts/{genome}/{sample}_{genome}_counts.txt"),
            zip, 
            sample = samples_df["sample_id"], 
            genome = samples_df["genome_id"]
        )
    output:
        flag = os.path.join(RESULT_DIR, "featureCounts/.all_counts_done")
    shell:
        """
        touch {output.flag}
        """


rule featureCounts_H_sapiens:
    input:
        hisat_bam   = expand(rules.hisat2_H_sapiens.output.bam, sample = SAMPLES)
    output:
        counts      = os.path.join(RESULT_DIR, "featureCounts/H_sapiens/{sample}_counts.txt")
    log:
        log         = os.path.join(LOG_DIR, "featureCounts/H_sapiens/{sample}.log")
    params:
        gtf_files   = config["H_sapiens_index_gtf"],
        strand      = config["featureCounts"]["strand"],
        paired_end  = config["featureCounts"]["paired_end"],
        features    = config["featureCounts"]["features"],
        attribute   = config["featureCounts"]["attribute"]
    conda:
        "featurecounts_env"
    shell:
        """
        featureCounts \
            -a {params.gtf_files} \
            -o {output.counts} \
            {params.paired_end} \
            -s {params.strand} \
            {params.features} \
            {params.attribute} \
            -T {resources.cpus_per_task} \
            {input.hisat_bam} 2> {log.log}
        """