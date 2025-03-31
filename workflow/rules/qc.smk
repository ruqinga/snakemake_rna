rule fastqc_rawdata:
    input:
        rawdata = get_fq_list
    output:
        raw_qc_finish = temp("Results/03_qc/rawdata/{sample}_tmp.txt")
    group: "processing_group"
    conda:
        config["conda_env"]
    params:
        qc_out = "Results/03_qc/rawdata"
    log:
        "Results/03_qc/rawdata/logs/{sample}_fastqc.log"
    shell:
        """
        fastqc {input.rawdata} -o {params.qc_out} > {log} 2>&1
        touch {output.raw_qc_finish}
        """

rule fastqc_cleandata:
    input:
        trimmed_data = get_trimmed_list
    output:
        trim_qc_finish = temp("Results/03_qc/cleandata/{sample}_tmp.txt")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        qc_out_trim = "Results/03_qc/cleandata"
    log:
        "Results/03_qc/cleandata/logs/{sample}.log"
    shell:
        """
        fastqc {input.trimmed_data} -o {params.qc_out_trim} > {log} 2>&1
        touch {output.trim_qc_finish}
        """

rule multiqc_rawdata:
    input:
        qc_raw_out=expand("Results/03_qc/rawdata/{sample}_tmp.txt", sample = samples)
    output:
        report="Results/03_qc/rawdata/multiqc_report.html"
    conda:
        config["conda_env"]
    group: "global_process"
    log:
        "Results/03_qc/rawdata/logs/multiqc_raw.log"
    params:
        qc_out="Results/03_qc/rawdata"
    shell:
        """
        multiqc {params.qc_out} -o {params.qc_out} > {log} 2>&1
        """


rule multiqc_cleandata:
    input:
        qc_clean_out=expand("Results/03_qc/cleandata/{sample}_tmp.txt", sample = samples)
    output:
        report="Results/03_qc/cleandata/multiqc_report.html"
    conda:
        config["conda_env"]
    group: "global_process"
    log:
        "Results/03_qc/cleandata/logs/multiqc_clean.log"
    params:
        qc_out="Results/03_qc/cleandata"
    shell:
        """
        multiqc {params.qc_out} -o {params.qc_out} > {log} 2>&1
        """
