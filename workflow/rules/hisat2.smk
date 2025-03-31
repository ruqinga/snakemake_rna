rule align_reads_with_hisat2_se:
    input:
        read = get_trimmed_list
    output:
        sam = temp("Results/04_align_out/{sample}_se.sam"),
        sorted_bam = "Results/04_align_out/{sample}_se.Hisat_aln.sorted.bam"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["hisat2"]["params"],
        index= config["hisat2"]["index"][config["Spe"]]
    log:
        log = "Results/04_align_out/logs/{sample}.log"
    shell:
        """
        hisat2 {params.option} "{params.index}" -U "{input.read}" -S "{output.sam}" > {log.log} 2>&1
 
        samtools sort -@ 20 -o "{output.sorted_bam}" "{output.sam}" >> {log.log} 2>&1
        
        echo "alignment finished at $(date)" >> {log.log} 2>&1
        """

rule align_reads_with_hisat2_pe:
    input:
        read=get_trimmed_list
    output:
        sam=temp("Results/04_align_out/{sample}_pe.sam"),
        sorted_bam="Results/04_align_out/{sample}_pe.Hisat_aln.sorted.bam"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option=config["hisat2"]["params"],
        index=config["hisat2"]["index"][config["Spe"]]
    log:
        log="Results/04_align_out/logs/{sample}.log"
    shell:
        """
        hisat2 {params.option} "{params.index}" -1 {input.read[0]} -2 {input.read[1]} -S "{output.sam}" > {log.log} 2>&1

        samtools sort -@ 20 -o "{output.sorted_bam}" "{output.sam}" >> {log.log} 2>&1
        
        echo "alignment finished at $(date)" >> {log.log} 2>&1
        """

rule bam2bw:
    input: bam="Results/04_align_out/{sample}.Hisat_aln.sorted.bam"
    output: bw = "Results/04_align_out/bw/{sample}.Hisat_aln.sorted.bw"
    conda: config["conda_env"]
    group: "processing_group"
    params:
        option=config["bamCoverage"]["params"]
    log:
        log="Results/04_align_out/logs/bw_{sample}.log"
    shell:
        """
        bamCoverage --bam "{input.bam}" --outFileName "{output.bw}" --outFileFormat bigwig {params.option} >> {log.log} 2>&1
        """

rule bam_index:
    input: bam = "Results/04_align_out/{sample}.Hisat_aln.sorted.bam"
    output: bam_bai = "Results/04_align_out/{sample}.Hisat_aln.sorted.bam.bai"
    conda: config["conda_env"]
    group: "processing_group"
    log:
        log="Results/04_align_out/logs/index_{sample}.log"
    shell:
        """
        samtools index "{input.bam}" "{output.bam_bai}" >> {log.log} 2>&1
        """
