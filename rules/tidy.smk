# 所有sample都跑完了再运行

# 合并fpkm
rule cufflinks_merge:
    input:
        files=expand("Results/06_norm/{sample}/genes.fpkm_tracking",
            sample=samples)
    output:
        merged_fpkm = "Results/06_norm/merged_fpkm.txt"
    conda:
        config["conda_env"]
    group: "processing_group"
    log:
        log="{norm_out}/logs/merge_fpkm.log"
    script:
        "../scripts/merge_fpkm.py"

# 合并count
rule count_merge:
    input:
        files=expand("{count_out}/{sample}.txt",
            sample=samples)
    output:
        merged_count = "Results/05_count/merged_count.txt"
    params:
        id_column="Geneid"
    log:
        log="Results/05_count/logs/merge_count.log"
    conda:
        config["conda_env"]
    group: "processing_group"
    script:
        "../scripts/merge_count.py"

rule rep_count_merge:
    input:
        files=expand("Results/07_rep_count/{sample}_rep.txt",
            sample=samples)
    output:
        merged_count = "Results/07_rep_count/merged_count_rep.txt"
    params:
        id_column="Geneid"
    log:
        log="Results/07_rep_count/logs/merge_count_rep.log"
    conda:
        config["conda_env"]
    group: "processing_group"
    script:
        "../scripts/merge_count.py"
