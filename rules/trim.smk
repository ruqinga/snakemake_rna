rule trim:
    input:
        read = get_fq_list
    output:
        trimmed_read = (
            "Results/01_cleandata/{sample}_trimmed.fq.gz"
            if config["dt"] == "SE"
            else [
                "Results/01_cleandata/{sample}_1_val_1.fq.gz",
                "Results/01_cleandata/{sample}_2_val_2.fq.gz"
            ]
        )
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["trim"]["params"],
        trim_out = "RNA_results/01_cleandata"
    log:
        log = "Results/01_cleandata/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
           trim_galore {params.option} {input.read} -o {params.trim_out} > {log.log} 2>&1
        else
           trim_galore {params.option} --paired {input.read[0]} {input.read[1]} -o {params.trim_out} > {log.log} 2>&1
        fi
        """
