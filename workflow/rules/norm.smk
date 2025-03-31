rule norm_cufflinks:
    input:
        bam = "Results/04_align_out/{sample}_{dt}.Hisat_aln.sorted.bam"
    output:
        fpkm = "Results/06_norm/{sample}_{dt}/genes.fpkm_tracking"
    conda:
        config["conda_cufflink"]
    group: "processing_group"
    params:
        norm_out = "Results/06_norm/{sample}_{dt}",
        gtf = config["gtf"][config["Spe"]]["genome"],
        thread = config["threads"]
    log:
        log="Results/06_norm/logs/{sample}_{dt}.log"
    shell:
        """
        cufflinks -p {params.thread} -G {params.gtf} -o {params.norm_out} {input.bam} > {log.log} 2>&1
        """