# 所有sample都跑完了再运行
rule cufflinks_merge:
    input:
        norm_out = directories['norm_out']
    output:
        fpkm = "{norm_out}/merged_fpkm.txt"
    conda:
        config["conda_cufflink"]
    group: "processing_group"
    params:
        norm_out=directories['norm_out']
    log:
        log="{norm_out}/logs/merge_fpkm.log"
    script:
        "../scripts/merge_FPKM.py"