rule count_feature_with_featurecount_se:
    input:
        bam = get_alined_list
    output:
        counts = "{count_out}/{sample}.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        count_out = directories["count_out"],
        gtf = config["gtf"][config["Spe"]],
        thread = config["threads"]
    log:
        log = "{count_out}/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
           featureCounts -T {params.thread} -t exon -g gene_name \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1
        else
           featureCounts -T {params.thread} -p --countReadPairs -t exon -g gene_name \
                      -a {params.gtf} -o {output.counts} {input.bam} > {log.log} 2>&1
        fi
        
        echo "featureCount finished at $(date)" > {log.log} 2>&1
        """