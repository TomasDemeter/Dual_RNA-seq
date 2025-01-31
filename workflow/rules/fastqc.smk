############################################
# Running FastQC on raw and filtered reads #
############################################
rule FastQC:
    input:
        trimmed_read_1  = rules.fastp.output.trimmed_1,
        trimmed_read_2  = rules.fastp.output.trimmed_2
    output:
        html_1          = os.path.join(RESULT_DIR, "FastQC/{sample}" + config["pattern_FWD_read"] + "_fastqc.html"),
        zip_file_1      = os.path.join(RESULT_DIR, "FastQC/{sample}" + config["pattern_FWD_read"] + "_fastqc.zip"),
        html_2          = os.path.join(RESULT_DIR, "FastQC/{sample}" + config["pattern_REV_read"] + "_fastqc.html"),
        zip_file_2      = os.path.join(RESULT_DIR, "FastQC/{sample}" + config["pattern_REV_read"] + "_fastqc.zip") 
    log:
        logs            = os.path.join(LOG_DIR, "FastQC/{sample}" + config["pattern_REV_read"] + "_fastqc.zip") 
    message:
        "Running FastQC on trimmed files"
    conda: 
        "multiqc_env"
    shell:
        """
        mkdir -p $(dirname {output.html_1})
        fastqc \
            {input.trimmed_read_1} {input.trimmed_read_2} \
            --outdir $(dirname {output.html_1}) \
            --threads {resources.cpus_per_task} \
            2> {log.logs}
        """