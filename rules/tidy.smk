# 所有sample都跑完了再运行

# 合并fpkm
rule cufflinks_merge:
    input:
        files=expand("{norm_out}/{sample}/genes.fpkm_tracking",
            norm_out=directories['norm_out'],
            sample=samples)
    output:
        merged_fpkm = "{norm_out}/merged_fpkm.txt"
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
            count_out=directories['count_out'],
            sample=samples)
    output:
        merged_count = "{count_out}/merged_count.txt"
    params:
        id_column="Geneid"
    log:
        log="{count_out}/logs/merge_count.log"
    conda:
        config["conda_env"]
    group: "processing_group"
    script:
        "../scripts/merge_count.py"

rule rep_count_merge:
    input:
        files=expand("{rep_out}/{sample}.txt",
            rep_out=directories['rep_out'],
            sample=samples)
    output:
        merged_count = "{rep_out}/merged_count.rep.txt"
    params:
        id_column="Geneid"
    log:
        log="{rep_out}/logs/merge_count_rep.log"
    conda:
        config["conda_env"]
    group: "processing_group"
    script:
        "../scripts/merge_count.py"
