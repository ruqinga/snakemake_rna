rule featurecount_se:
    input:
        bam = "Results/04_align_out/{sample}_se.Hisat_aln.sorted.bam"
    output:
        counts = "Results/05_count/{sample}_se.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        gtf = config["gtf"][config["Spe"]]["genome"],
        thread = config["threads"]
    log:
        log = "Results/05_count/logs/{sample}.log"
    shell:
        """
        featureCounts -T {params.thread} -t exon -g gene_name \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1
        
        echo "featureCount finished at $(date)" >> {log.log} 2>&1
        """

rule featurecount_pe:
    input:
        bam="Results/04_align_out/{sample}_pe.Hisat_aln.sorted.bam"
    output:
        counts="Results/05_count/{sample}_pe.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        gtf=config["gtf"][config["Spe"]]["genome"],
        thread=config["threads"]
    log:
        log="Results/05_count/logs/{sample}.log"
    shell:
        """
        featureCounts -T {params.thread} -p --countReadPairs -t exon -g gene_name \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1

        echo "featureCount finished at $(date)" >> {log.log} 2>&1
        """

rule repeats_featurecount_se:
    input:
        bam="Results/04_align_out/{sample}_se.Hisat_aln.sorted.bam"
    output:
        counts="Results/07_repeats_count/{sample}_se.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        gtf=config["gtf"][config["Spe"]]["repeats"],
        thread=config["threads"]
    log:
        log="Results/07_repeats_count/logs/{sample}.log"
    shell:
        """
        featureCounts -T {params.thread} -t exon -g gene_id \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1

        echo "featureCount finished at $(date)" >> {log.log} 2>&1
        """

rule repeats_featurecount_pe:
    input:
        bam="Results/04_align_out/{sample}_pe.Hisat_aln.sorted.bam"
    output:
        counts="Results/07_repeats_count/{sample}_pe.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        gtf=config["gtf"][config["Spe"]]["repeats"],
        thread=config["threads"]
    log:
        log="Results/07_repeats_count/logs/{sample}.log"
    shell:
        """
        featureCounts -T {params.thread} -p --countReadPairs -t exon -g gene_id \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1
        echo "featureCount finished at $(date)" >> {log.log} 2>&1
        """
