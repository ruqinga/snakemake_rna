# 所有sample都跑完了再运行

# 合并fpkm
rule cufflinks_merge:
    input:
        files=lambda wildcards: [
            f"Results/06_norm/{sample}_{'pe' if sample_info[sample] == 'PE' else 'se'}/genes.fpkm_tracking"
            for sample in samples
        ]
    output:
        merged_fpkm = "Results/06_norm/merged_fpkm.txt"
    conda:
        config["conda_env"]
    group: "global_process"
    log:
        log="Results/06_norm/logs/merge_fpkm.log"
    script:
        "../scripts/merge_fpkm.py"

# 合并count
rule count_merge:
    input:
        files = lambda wildcards: [
            f"Results/05_count/{sample}_{'pe' if sample_info[sample] == 'PE' else 'se'}.txt"
            for sample in samples
        ]
    output:
        merged_count = "Results/05_count/merged_count.txt"
    params:
        id_column="Geneid"
    log:
        log="Results/05_count/logs/merge_count.log"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/merge_count.py"

rule rep_count_merge:
    input:
        files = lambda wildcards: [
            f"Results/07_repeats_count/{sample}_{'pe' if sample_info[sample] == 'PE' else 'se'}.txt"
            for sample in samples
        ]
    output:
        merged_count = "Results/07_repeats_count/merged_count.rep.txt"
    params:
        id_column="Geneid"
    log:
        log="Results/07_repeats_count/logs/merge_count_rep.log"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/merge_count.py"
