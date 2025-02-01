rule hisat2_index_H_sapiens:
    input:
        genome      = config["H_sapiens_index_fastq"],
        gtf         = config["H_sapiens_index_gtf"]
    output:
        flag        = os.path.join(config["H_sapiens_index"] + ".index_done"),
        ss          = temp(os.path.join(GENOMES_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly/splicesites.txt")),
        exon        = temp(os.path.join(GENOMES_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly/exons.txt"))
    params:
        index_base  = config["H_sapiens_index"] + "/H_sapiens"  
    log:
        logs        = os.path.join(LOG_DIR, "hisat2_index/hisat2_index.log")
    conda:
        "hisat2_env"
    shell:
        """
        # Decompress files first
        gunzip -c {input.gtf} > temp.gtf
        gunzip -c {input.genome} > temp.fa
        
        # Create index directory if it doesn't exist
        mkdir -p {params.index_base}
        
        # Extract splice sites
        hisat2_extract_splice_sites.py temp.gtf > {output.ss}
        
        # Extract exons
        hisat2_extract_exons.py temp.gtf > {output.exon}
        
        # Build the index
        hisat2-build \
            -p {resources.cpus_per_task} \
            --ss {output.ss} \
            --exon {output.exon} \
            temp.fa \
            {params.index_base} \
            2> {log.logs}
            
        # Clean up temporary files
        rm temp.gtf temp.fa
        
        # Create flag file
        touch {output.flag}
        """

rule hisat2_H_sapiens:
    input:
        indexes_done    = rules.hisat2_index_H_sapiens.output.flag,
        read_1          = rules.fastp.output.trimmed_1,
        read_2          = rules.fastp.output.trimmed_2
    output:
        bam             = os.path.join(RESULT_DIR, "hisat2/H_sapiens/mapped/{sample}.bam"),
        unmapped_1      = os.path.join(RESULT_DIR, "hisat2/H_sapiens/unmapped/{sample}_unmapped.1.fastq.gz"),
        unmapped_2      = os.path.join(RESULT_DIR, "hisat2/H_sapiens/unmapped/{sample}_unmapped.2.fastq.gz"),
        summary         = os.path.join(RESULT_DIR, "hisat2/H_sapiens/{sample}.summary.txt")
    params:
        unmapped_prefix = os.path.join(RESULT_DIR, "hisat2/H_sapiens/unmapped/{sample}_unmapped"),
        genome_index    = config["H_sapiens_index"] + "/H_sapiens"
    log:
        logs            = os.path.join(LOG_DIR, "hisat2/H_sapiens/{sample}_hisat2.log")
    conda:
        "hisat2_env"
    shell:
        """
        (
            # Run HISAT2
            hisat2 \
                -p {resources.cpus_per_task} \
                -x {params.genome_index} \
                -1 {input.read_1} \
                -2 {input.read_2} \
                --un-conc-gz {params.unmapped_prefix} \
                --summary-file {output.summary} \
                --no-temp-splicesite \
                --rna-strandness RF \
                | samtools view -bS - > {output.bam}
                
            # Rename the unmapped files to include .fastq.gz extension
            mv {params.unmapped_prefix}.1 {output.unmapped_1}
            mv {params.unmapped_prefix}.2 {output.unmapped_2}
        ) &> {log.logs}
        """