# 分别跟踪trim和align是否完成会很复杂，直接用是否输出merged_fpkm.txt当作mian pipeline 的完成标记。通过params提供目录文件夹
rule summary:
    input:
        all_rule_finish = "Results/06_norm/merged_fpkm.txt"
    output:
        summary = "Results/summary.csv"
    params:
        trim_log_dir="Results/02_trim_out/logs/",
        align_log_dir="Results/04_align_out/logs/"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/summary.py"
