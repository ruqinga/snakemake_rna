rule norm_cufflinks:
    input:
        bam = get_alined_list
    output:
        fpkm = "Results/06_norm/{sample}/genes.fpkm_tracking"
    conda:
        config["conda_cufflink"]
    group: "processing_group"
    params:
        norm_out_sample = lambda wildcards: f"Results/06_norm/{wildcards.sample}",
        gtf = config["gtf"][config["Spe"]]["genome"],
        thread = config["threads"]
    log:
        log="Results/06_norm/logs/{sample}.log"
    shell:
        """
        cufflinks -p {params.thread} -G {params.gtf} -o {params.norm_out_sample} {input.bam} > {log.log} 2>&1
        """