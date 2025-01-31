#########################
# RNA-seq read trimming #
#########################
rule fastp:
    input:
        read_1              = os.path.join(RAW_READS_DIR, "{sample}" + config["pattern_FWD_read"] + config["suffix"]),
        read_2              = os.path.join(RAW_READS_DIR, "{sample}" + config["pattern_REV_read"] + config["suffix"])
    output:   
        trimmed_1           = os.path.join(RESULT_DIR, "fastp/{sample}" + config["pattern_FWD_read"] + config["suffix"]),
        trimmed_2           = os.path.join(RESULT_DIR, "fastp/{sample}" + config["pattern_REV_read"] + config["suffix"]),
        html                = os.path.join(RESULT_DIR, "fastp/logs/{sample}_fastp.html"),
        json                = os.path.join(RESULT_DIR, "fastp/logs/{sample}_fastp.json")
    log:
        logs                = os.path.join(LOG_DIR, "fastp/{sample}_fastp.log")
    message:
        "fastp trimming RNA-seq reads for {wildcards.sample}"
    conda: 
        "fastp_bowtie2"
    params:
        phread_quality      = config["fastp"]["phread_quality"],
        length_required     = config["fastp"]["length_required"],
        poly_x_min_len      = config["fastp"]["poly_x_min_len"],
        adapters            = config["adapters"]
    shell:
        """
        (mkdir -p $(dirname {output.json})
        fastp \
            -i {input.read_1} \
            -I {input.read_2} \
            -o {output.trimmed_1} \
            -O {output.trimmed_2} \
            --thread {resources.cpus_per_task} \
            --qualified_quality_phred {params.phread_quality} \
            --length_required {params.length_required} \
            --cut_tail \
            --trim_poly_x \
            --poly_x_min_len {params.poly_x_min_len} \
            --html {output.html} \
            --json {output.json} \
            --adapter_fasta {params.adapters}) 2> {log.logs}
        """