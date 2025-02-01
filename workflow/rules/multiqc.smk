rule multiqc:
    input:
        hisat               = rules.map_all_samples.output.flag,
        fastqc_input_1      = expand(rules.FastQC.output.zip_file_1, sample = SAMPLES),
        fastqc_input_2      = expand(rules.FastQC.output.zip_file_2, sample = SAMPLES),
        featurecounts       = expand(rules.featureCounts_single.output.counts, sample = SAMPLES, genome = GENOMES)
    output:
        output              = os.path.join(RESULT_DIR, "MultiQC/{genome}/{genome}_multiqc_report.html")
    params:
        hisat_logs         = os.path.join(LOG_DIR, "hisat2/bacteria/{genome}"),
        fastqc_zip          = os.path.join(RESULT_DIR, "FastQC/"),
        featurecounts_dir   = os.path.join(RESULT_DIR,"featureCounts/{genome}")
    log:
        logs                = os.path.join(LOG_DIR, "MultiQC/{genome}_log.logs")
    conda: 
        "multiqc_env"
    message:
        "Summarising reports with multiqc"
    shell:
        """
        (mkdir -p $(dirname {output.output}) && \
        multiqc \
            --force \
            --outdir $(dirname {output.output}) \
            {params.hisat_logs} \
            {params.fastqc_zip} \
            {params.featurecounts_dir} \
            2>&1 && \
        mv $(dirname {output.output})/multiqc_report.html {output.output}) > {log.logs} 2>&1
        """

rule multiqc_H_sapiens:
    input:
        hisat2              = expand(rules.hisat2_H_sapiens.output.bam, sample = SAMPLES),
        fastp_input_html    = expand(rules.fastp.output.html, sample = SAMPLES),
        fastp_input_json    = expand(rules.fastp.output.json, sample = SAMPLES),
        fastqc_input_1      = expand(rules.FastQC.output.zip_file_1, sample = SAMPLES),
        fastqc_input_2      = expand(rules.FastQC.output.zip_file_2, sample = SAMPLES),
        featurecounts       = expand(rules.featureCounts_H_sapiens.output.counts, sample = SAMPLES, genome = GENOMES)
    output:
        output              = os.path.join(RESULT_DIR, "MultiQC/H_sapiens/H_sapiens_multiqc_report.html")
    params:
        hisat2_logs         = os.path.join(LOG_DIR, "hisat2/H_sapiens/"),
        fastp_logs          = os.path.join(RESULT_DIR, "fastp/logs/"),
        fastqc_zip          = os.path.join(RESULT_DIR, "FastQC/"),
        featurecounts_dir   = os.path.join(RESULT_DIR,"featureCounts/H_sapiens")
    log:
        logs                = os.path.join(LOG_DIR, "MultiQC/H_sapiens_log.logs")
    conda: 
        "multiqc_env"
    message:
        "Summarising reports with multiqc for human genome"
    shell:
        """
        (mkdir -p $(dirname {output.output}) && \
        multiqc \
            --force \
            --outdir $(dirname {output.output}) \
            {params.hisat2_logs} \
            {params.fastp_logs} \
            {params.fastqc_zip} \
            {params.featurecounts_dir} \
            2>&1 && \
        mv $(dirname {output.output})/multiqc_report.html {output.output}) > {log.logs} 2>&1
        """